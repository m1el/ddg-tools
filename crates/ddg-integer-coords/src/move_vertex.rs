//! Moving an inserted vertex (paper Sec. 3.6; SPEC §4.7, I9).
//!
//! The paper's general route: the new place `x = exp_i(v)` is traced over the triangulation
//! `T¹` from the vertex `i`, a vertex is inserted at `x` (a face split, [`crate::insert`], or an
//! edge split when the path ends on an edge, [`crate::split`]), and `i` is removed
//! ([`crate::remove`]). The trace needs directions at `i`; this representation stores none, so
//! they are laid out from the corner angles (SPEC I9):
//! [`AngleFrames::of(&t).trace(...)`](ddg_intrinsic::AngleFrames::trace), the same tracing
//! kernel as the signposts'. **The frame of `i`**: direction `0` is that of
//! [`Complex::outgoing_halfedge`](ddg_intrinsic::Complex::outgoing_halfedge), the others the
//! running sum of the corner angles counter-clockwise, scaled by `2π / Θ` with `Θ` from the
//! lengths; it belongs to the current connectivity (a flip at `i` can change it), so a
//! displacement is to be computed right before the move.
//!
//! The removed vertex's id is taken by the inserted one (the removal renumbers the last
//! vertex into the freed id, [`VertexRemoval`](ddg_intrinsic::complex::VertexRemoval)), so the
//! vertex keeps its id; edges and faces are renumbered as the removal says.
//!
//! **A refused move changes nothing.** In particular a removal that gets stuck
//! ([`IntegerError::RemovalStuck`], e.g. the straight star of I27 when the path ends on an edge)
//! refuses the move and the state is restored (the reference ignores the failed removal and
//! leaves both vertices, `ic:1717-1738`, SPEC §11): the removal undoes its own flips, and the
//! inserted vertex, of degree three after a face split, is removed again, which restores every
//! id and value (I29: this replaced a copy of the whole state per move, `O(|E|)`, which made
//! the optimal Delaunay iteration quadratic); after an edge split (a path that ends on an edge)
//! the state is still copied first.
//! Inserted vertices on the boundary cannot be removed (I7) and so are not moved
//! ([`IntegerError::Unsupported`]).
//!
//! The optimal Delaunay iteration of IT-SPEC N21 runs on these moves ([`crate::odt`]).

use std::f64::consts::TAU;

use ddg_intrinsic::{
    AngleFrames, Insertion, IntrinsicError, MoveRefusal, SurfacePoint, TraceOptions, TraceStop,
};

use crate::coords::IntegerTriangulation;
use crate::error::IntegerError;
use crate::insert::{FaceSplitOptions, FaceSplitReport};
use crate::remove::Removal;

/// How the vertex was inserted at its new place.
#[derive(Debug, Clone, PartialEq)]
pub enum MoveInsertion {
    /// The path ended inside a face: a face split.
    Face(Box<FaceSplitReport>),
    /// The path ended on an edge: an edge split.
    Edge(Insertion),
}

/// What [`IntegerTriangulation::move_vertex`] did.
#[derive(Debug, Clone, PartialEq)]
pub struct VertexMoveReport {
    /// The vertex (the same id before and after the move).
    pub vertex: u32,
    /// Its place on the input mesh before the move.
    pub from: SurfacePoint,
    /// Its place on the input mesh after the move.
    pub to: SurfacePoint,
    /// The length of the displacement.
    pub distance: f64,
    /// Where the traced path ended on the triangulation before the move.
    pub end: SurfacePoint,
    /// The number of edges of the triangulation the path crossed.
    pub n_crossings: usize,
    /// The insertion at the new place (`None` for a zero displacement).
    pub insertion: Option<MoveInsertion>,
    /// The removal of the vertex at its old place (`None` for a zero displacement); its flips
    /// and renumbering refer to the state after the insertion.
    pub removal: Option<Removal>,
}

impl IntegerTriangulation {
    /// Moves the inserted vertex `v` by `displacement`, a vector in the frame of `v` laid out
    /// from the corner angles (see the [module documentation](crate::move_vertex)): traces the
    /// new place over the triangulation, inserts a vertex there (a face split with `options`,
    /// or an edge split), removes `v`. The moved vertex keeps the id `v`. A zero displacement
    /// changes nothing.
    ///
    /// # Errors
    /// [`IntegerError::IndexOutOfRange`]; [`IntegerError::Intrinsic`] with
    /// [`IntrinsicError::MoveRefused`] for an input vertex ([`MoveRefusal::InputVertex`]), a
    /// displacement that is not finite ([`MoveRefusal::NotFinite`]), a path that ends at a
    /// vertex ([`MoveRefusal::OnInputVertex`] for an input vertex), and
    /// [`IntrinsicError::MoveNotLocated`] for a path that leaves the surface or stops at a vertex
    /// before its length; [`IntegerError::Unsupported`] for a vertex on the boundary (I7); the
    /// errors of the insertion ([`split_face`](Self::split_face),
    /// [`split_edge`](Self::split_edge)) and of the removal
    /// ([`remove_vertex`](Self::remove_vertex), e.g. [`IntegerError::RemovalStuck`]). The state
    /// is unchanged in every case.
    pub fn move_vertex(
        &mut self,
        v: usize,
        displacement: [f64; 2],
        options: &FaceSplitOptions,
    ) -> Result<VertexMoveReport, IntegerError> {
        let nv = self.complex.n_vertices();
        if v >= nv {
            return Err(IntegerError::IndexOutOfRange {
                what: "vertex",
                index: v,
                bound: nv,
            });
        }
        let refused = |reason| {
            IntegerError::Intrinsic(IntrinsicError::MoveRefused {
                vertex: v as u32,
                reason,
            })
        };
        let not_located = |detail: String| {
            IntegerError::Intrinsic(IntrinsicError::MoveNotLocated {
                vertex: v as u32,
                detail,
            })
        };
        if self.is_input_vertex(v) {
            return Err(refused(MoveRefusal::InputVertex));
        }
        if !(displacement[0].is_finite() && displacement[1].is_finite()) {
            return Err(refused(MoveRefusal::NotFinite));
        }
        if self.complex.is_boundary_vertex(v) {
            return Err(IntegerError::Unsupported {
                operation: "move_vertex",
                reason: format!(
                    "vertex {v} lies on the boundary; inserted boundary vertices cannot be \
                     removed (SPEC I7), so they are not moved"
                ),
            });
        }
        let from = self.locations[v];
        let distance = displacement[0].hypot(displacement[1]);
        if distance == 0.0 {
            return Ok(VertexMoveReport {
                vertex: v as u32,
                from,
                to: from,
                distance,
                end: SurfacePoint::Vertex(v as u32),
                n_crossings: 0,
                insertion: None,
                removal: None,
            });
        }
        // The tracing kernel takes directions in [0, 2π).
        let direction = displacement[1].atan2(displacement[0]).rem_euclid(TAU);
        let direction = if direction >= TAU { 0.0 } else { direction };
        let trace = AngleFrames::of(&*self).trace(
            SurfacePoint::Vertex(v as u32),
            direction,
            distance,
            &TraceOptions::default(),
        )?;
        if trace.stop != TraceStop::Completed {
            return Err(not_located(format!(
                "the path from vertex {v} stopped ({:?}) at {} after {:?} of {distance:?}",
                trace.stop, trace.end, trace.length
            )));
        }
        // A refused insertion changes nothing. A refused removal (it restores what it flipped)
        // is undone by removing the inserted vertex again: after a face split it has degree
        // three, and the merge restores the face, its halfedges and every id exactly (the
        // split appended the vertex, its edges and faces, so nothing is renumbered). After an
        // edge split (a path that ends on an edge, rare) the whole state is copied first.
        let mut snapshot = None;
        let end = self.snap_move_end(trace.end);
        let insertion = match end {
            SurfacePoint::Face { face, bary } => {
                MoveInsertion::Face(Box::new(self.split_face(face as usize, bary, options)?))
            }
            SurfacePoint::Edge { edge, t } => {
                snapshot = Some(self.clone());
                match self.split_edge(edge as usize, t) {
                    Ok(insertion) => MoveInsertion::Edge(insertion),
                    Err(err) => {
                        if let Some(s) = snapshot {
                            *self = s;
                        }
                        return Err(err);
                    }
                }
            }
            SurfacePoint::Vertex(u) => {
                return Err(if self.is_input_vertex(u as usize) {
                    refused(MoveRefusal::OnInputVertex { vertex: u })
                } else {
                    not_located(format!("the path from vertex {v} ends at vertex {u}"))
                });
            }
        };
        match self.remove_vertex(v) {
            Ok(removal) => Ok(VertexMoveReport {
                vertex: v as u32,
                from,
                to: self.locations[v],
                distance,
                end,
                n_crossings: trace.n_crossings,
                insertion: Some(insertion),
                removal: Some(removal),
            }),
            Err(err) => {
                match (snapshot, &insertion) {
                    (Some(s), _) => *self = s,
                    (None, MoveInsertion::Face(split)) => {
                        let m = split.insertion.vertex as usize;
                        self.remove_vertex(m).map_err(|undo| {
                            IntegerError::Intrinsic(IntrinsicError::Inconsistent {
                                what: format!(
                                    "the move of vertex {v} could not be undone: removing the \
                                     inserted vertex {m} again failed ({undo}) after the \
                                     removal of {v} failed ({err})"
                                ),
                            })
                        })?;
                    }
                    (None, MoveInsertion::Edge(_)) => {}
                }
                Err(err)
            }
        }
    }

    /// The direction of halfedge `h` at its tail in the frame of [`move_vertex`]
    /// (`AngleFrames`, SPEC I9): for building displacements.
    ///
    /// [`move_vertex`]: Self::move_vertex
    pub fn frame_direction(&self, h: usize) -> f64 {
        AngleFrames::of(self).signpost(h)
    }

    /// The end of a move's path as the point to insert at: a point of a face with exactly one
    /// barycentric coordinate at most [`MOVE_END_SNAP`] is the point of the side opposite that
    /// corner, if that is a point strictly inside the side (a path along an edge ends a
    /// rounding error inside one of its faces, and a face split there would be refused as too
    /// thin). Other points as they are (near a corner the insertion is refused as too thin).
    fn snap_move_end(&self, end: SurfacePoint) -> SurfacePoint {
        let SurfacePoint::Face { face, bary } = end else {
            return end;
        };
        let small: Vec<usize> = (0..3).filter(|&k| bary[k] <= MOVE_END_SNAP).collect();
        let f = face as usize;
        match small[..] {
            [k] => {
                // The side opposite corner k: halfedge k + 1, from corner k + 1 to k + 2.
                let (i, j) = ((k + 1) % 3, (k + 2) % 3);
                let h = 3 * f + i;
                let th = bary[j] / (bary[i] + bary[j]);
                let e = self.complex.edge(h);
                let t = if self.complex.edge_halfedge(e) == h {
                    th
                } else {
                    1.0 - th
                };
                if t > 0.0 && t < 1.0 {
                    SurfacePoint::Edge { edge: e as u32, t }
                } else {
                    end
                }
            }
            _ => end,
        }
    }
}

/// The barycentric coordinate below which the end of a move's path is taken as a point of a
/// side of its face ([`IntegerTriangulation::move_vertex`]); the `LOCATE_SNAP` of the point
/// queries (SPEC §4.9).
pub const MOVE_END_SNAP: f64 = 1e-12;
