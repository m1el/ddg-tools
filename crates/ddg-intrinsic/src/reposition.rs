//! Vertex repositioning (paper Sec. 3.3.3, Alg. 10): an inserted vertex is moved inside its
//! flat star without changing the surface.
//!
//! An inserted vertex is flat (its angle sum is `2π`, or `π` on the boundary), so its star
//! unfolds into the plane, the vertex at the origin and the directions of its signposts as
//! they are. Moving it by a vector `v` inside that polygon changes only the lengths of the
//! edges at the vertex (the law of cosines, [`moved_length`]) and the directions at both ends
//! of those edges; the polygon, and so the surface, stays the same.
//!
//! # When the local update applies (SPEC N19, N20)
//!
//! Exactly when every triangle of the star keeps a positive area: the new place must lie in
//! the kernel of the star polygon. The area of the triangle on the link edge `d = k − j` is
//! affine in `v` ([`moved_double_area`]), so this is tested exactly, in addition to the shape
//! limit of N5 for each new triangle. The paper's test (`|v|` below the smallest height of a
//! triangle of the star over its link edge) is sufficient and stricter; whether it held is
//! reported ([`VertexMove::within_height_bound`]). A move that leaves the star is refused,
//! with the fraction of `v` that would still be inside ([`MoveRefusal::NotInStar`]); it is
//! never shortened here. (The paper's general route, removing the vertex and inserting it
//! again, is not implemented.)
//!
//! # The place on the input mesh
//!
//! Found, before anything is changed, by a tracing query over the input mesh from the
//! vertex's old place along `v` for `|v|` (the straight segment lies in the star, which
//! contains no other vertex). The paper's vertex update traces from a neighbour instead; in
//! exact arithmetic it is the same point. The direction back along the trace is the
//! direction of `−v` in the frame of the new place, which fixes the first signpost; the others
//! follow from the new corner angles (Sec. 3.2.3). A place within
//! [`LOCATION_SNAP`](crate::insert::LOCATION_SNAP) of an input edge is put on it, as for an
//! insertion.
//!
//! A vertex inserted on the boundary moves along the boundary only: its place is
//! interpolated on the input edge it lies on, as for an edge split.

use ddg_core::geometry::triangle_area;
use ddg_core::topology::NONE;

use crate::error::{IntrinsicError, MoveRefusal};
use crate::flip::triangle_is_valid;
use crate::formulas::{FLIP_SHAPE_LIMIT, moved_double_area, moved_length};
use crate::insert::wrap;
use crate::point::SurfacePoint;
use crate::trace::{TraceOptions, TraceStop, run};
use crate::triangulation::IntrinsicTriangulation;

/// What [`IntrinsicTriangulation::move_vertex`] did.
#[derive(Debug, Clone, PartialEq)]
pub struct VertexMove {
    /// The vertex.
    pub vertex: u32,
    /// Its place on the input mesh before the move.
    pub from: SurfacePoint,
    /// Its place on the input mesh after the move.
    pub to: SurfacePoint,
    /// The length of the displacement.
    pub distance: f64,
    /// The direction from the new place back to the old one, in the frame of the vertex after
    /// the move: moving by `distance` in this direction undoes the move.
    pub back_direction: f64,
    /// The barycentric coordinate that the traced place had opposite the input edge it was
    /// moved onto ([`LOCATION_SNAP`](crate::insert::LOCATION_SNAP)); `None` if it was not
    /// moved.
    pub snapped: Option<f64>,
    /// Whether the displacement was shorter than the smallest height of a triangle of the star
    /// over its side opposite the vertex: the paper's conservative test for the local update.
    /// The update is exact without it (see the [module documentation](self)).
    pub within_height_bound: bool,
    /// The smallest shape number (twice the area over the product of the longer new edge and
    /// the longest side, SPEC N5) of the triangles of the new star.
    pub min_shape: f64,
}

/// A triangle of the star of a vertex, unfolded around it.
struct StarFace {
    /// The halfedge of the face that leaves the vertex.
    h: usize,
    /// Its direction at the vertex (the signpost).
    phi: f64,
    /// The corner angle at the vertex.
    corner: f64,
    /// The lengths of the side leaving the vertex, of the side arriving at it, of the side
    /// opposite.
    lengths: [f64; 3],
    /// Twice the area.
    double_area: f64,
}

impl IntrinsicTriangulation {
    /// The faces of the star of `v` in counter-clockwise order from its first halfedge.
    fn star(&self, v: usize) -> Vec<StarFace> {
        self.outgoing_halfedges(v)
            .map(|h| {
                let lengths = [
                    self.edge_lengths[self.halfedge_edge[h] as usize],
                    self.edge_lengths[self.halfedge_edge[self.prev(h)] as usize],
                    self.edge_lengths[self.halfedge_edge[self.next(h)] as usize],
                ];
                let area = triangle_area(lengths[0], lengths[2], lengths[1]).unwrap_or(0.0);
                StarFace {
                    h,
                    phi: self.signposts[h],
                    corner: self.angle_at(h),
                    lengths,
                    double_area: 2.0 * area,
                }
            })
            .collect()
    }

    /// Moves the inserted vertex `v` by `displacement`, a vector in the frame of `v` (the
    /// frame its signposts refer to: angle zero is the direction
    /// [`signpost`](Self::signpost) `0`), inside its flat star (paper Sec. 3.3.3, Alg. 10).
    /// The lengths of the edges at `v`, the signposts at both ends of those edges and the
    /// place of `v` on the input mesh are updated; nothing else changes, and the surface
    /// stays the same. See the [module documentation](self) for when this applies.
    ///
    /// A vertex on the boundary moves along the boundary only: in its frame the boundary is
    /// the x-axis, so `displacement[1]` must be zero. A zero displacement changes nothing.
    ///
    /// # Errors
    /// [`IntrinsicError::IndexOutOfRange`]; [`IntrinsicError::MoveRefused`] for a vertex of
    /// the input mesh, a vertex with an edge to itself, a displacement that is not finite,
    /// leaves the boundary or leaves the star, a new triangle that would be too thin, a new
    /// place that is an input vertex; [`IntrinsicError::MoveNotLocated`] if the trace over the
    /// input mesh fails although the move is inside the star (an inconsistency, not repaired).
    /// The triangulation is unchanged in every case.
    pub fn move_vertex(
        &mut self,
        v: usize,
        displacement: [f64; 2],
    ) -> Result<VertexMove, IntrinsicError> {
        if v >= self.n_vertices() {
            return Err(IntrinsicError::IndexOutOfRange {
                what: "vertex",
                index: v,
                bound: self.n_vertices(),
            });
        }
        let refused = |reason| IntrinsicError::MoveRefused {
            vertex: v as u32,
            reason,
        };
        let not_located = |detail: String| IntrinsicError::MoveNotLocated {
            vertex: v as u32,
            detail,
        };
        if v < self.n_input_vertices() {
            return Err(refused(MoveRefusal::InputVertex));
        }
        if !(displacement[0].is_finite() && displacement[1].is_finite()) {
            return Err(refused(MoveRefusal::NotFinite));
        }
        let boundary = self.is_boundary_vertex(v);
        if boundary && displacement[1] != 0.0 {
            return Err(refused(MoveRefusal::LeavesBoundary {
                normal: displacement[1],
            }));
        }
        if let Some(h) = self.outgoing_halfedges(v).find(|&h| self.to_vertex(h) == v) {
            return Err(refused(MoveRefusal::SelfEdge {
                edge: self.halfedge_edge[h],
            }));
        }
        let from = self.locations[v];
        let r = displacement[0].hypot(displacement[1]);
        let star = self.star(v);
        if r == 0.0 {
            let min_shape = star
                .iter()
                .map(|s| shape_number(s.double_area, s.lengths))
                .fold(f64::INFINITY, f64::min);
            return Ok(VertexMove {
                vertex: v as u32,
                from,
                to: from,
                distance: 0.0,
                back_direction: 0.0,
                snapped: None,
                within_height_bound: true,
                min_shape,
            });
        }
        let phi_v = displacement[1].atan2(displacement[0]);

        // The new star: areas, shapes, lengths.
        let mut fraction = f64::INFINITY;
        let mut turned: Option<u32> = None;
        let mut thinnest = (f64::INFINITY, 0u32);
        let mut min_height = f64::INFINITY;
        // New length of the side leaving v in every face, and of the side arriving in the last.
        let mut leaving = Vec::with_capacity(star.len());
        let mut arriving_last = 0.0;
        for s in &star {
            let [la, lb, lc] = s.lengths;
            let (pa, pb) = (s.phi, s.phi + s.corner);
            let a = [la * pa.cos(), la * pa.sin()];
            let b = [lb * pb.cos(), lb * pb.sin()];
            let d = [b[0] - a[0], b[1] - a[1]];
            let new_double_area = moved_double_area(s.double_area, displacement, d);
            let push = displacement[0] * d[1] - displacement[1] * d[0];
            if push > 0.0 {
                fraction = fraction.min(s.double_area / push);
            }
            min_height = min_height.min(s.double_area / lc);
            let face = (s.h / 3) as u32;
            if new_double_area <= 0.0 && turned.is_none() {
                turned = Some(face);
            }
            let la_new = moved_length(la, r, pa - phi_v);
            let lb_new = moved_length(lb, r, pb - phi_v);
            let shape = if triangle_is_valid(la_new, lc, lb_new) {
                shape_number(new_double_area, [la_new, lb_new, lc])
            } else {
                0.0
            };
            if shape < thinnest.0 || shape.is_nan() {
                thinnest = (shape, face);
            }
            leaving.push(la_new);
            arriving_last = lb_new;
        }
        if let Some(face) = turned {
            return Err(refused(MoveRefusal::NotInStar { face, fraction }));
        }
        if !(thinnest.0 >= FLIP_SHAPE_LIMIT) {
            return Err(refused(MoveRefusal::TooThin {
                face: thinnest.1,
                shape: thinnest.0,
            }));
        }

        // The new place on the input mesh, and the direction back from it.
        let (to, back, snapped) = if boundary {
            let SurfacePoint::Edge { edge, t } = from else {
                return Err(not_located(format!(
                    "a vertex on the boundary is located at {from}, not on an input edge"
                )));
            };
            let first = star[0].h;
            let j = self.to_vertex(first);
            let Some(tj) = self.parameter_on_input_edge(j, edge as usize) else {
                return Err(not_located(format!(
                    "its boundary neighbour {j} does not lie on input edge {edge}"
                )));
            };
            let t_new = t + (displacement[0] / star[0].lengths[0]) * (tj - t);
            if !(t_new > 0.0 && t_new < 1.0) {
                return Err(not_located(format!(
                    "the new parameter {t_new:?} on input edge {edge} is not inside the edge"
                )));
            }
            // In the boundary frame the old place is in the direction opposite the motion.
            let back = if displacement[0] > 0.0 { std::f64::consts::PI } else { 0.0 };
            (SurfacePoint::Edge { edge, t: t_new }, back, None)
        } else {
            let direction = wrap(phi_v);
            let outcome = run(
                &*self.input,
                from,
                self.direction_on_input(v, direction),
                r,
                &TraceOptions::default(),
                None,
            )
            .map_err(|e| not_located(e.to_string()))?;
            let trace = outcome.trace;
            if trace.stop != TraceStop::Completed {
                if let (TraceStop::Vertex, SurfacePoint::Vertex(u)) = (trace.stop, trace.end) {
                    return Err(refused(MoveRefusal::OnInputVertex { vertex: u }));
                }
                return Err(not_located(format!(
                    "the path from {from} stopped at {} ({:?}) after a length of {:?} of {r:?}",
                    trace.end, trace.stop, trace.length
                )));
            }
            if let SurfacePoint::Vertex(u) = trace.end {
                return Err(refused(MoveRefusal::OnInputVertex { vertex: u }));
            }
            let placed = self
                .place_at_end(&trace)
                .map_err(|e| not_located(e.to_string()))?;
            (placed.location, placed.back, placed.snapped)
        };

        // Apply: lengths, place, signposts.
        for (s, &l) in star.iter().zip(&leaving) {
            let e = self.halfedge_edge[s.h] as usize;
            self.edge_lengths[e] = l;
            if !boundary {
                self.edge_on_input[e] = NONE;
            }
        }
        if boundary {
            let last = star[star.len() - 1].h;
            let e = self.halfedge_edge[self.prev(last)] as usize;
            self.edge_lengths[e] = arriving_last;
        }
        self.locations[v] = to;
        // Outgoing: the first from the trace (interior) or the boundary frame, then corner by
        // corner.
        if !boundary {
            let s = &star[0];
            let a = [s.lengths[0] * s.phi.cos(), s.lengths[0] * s.phi.sin()];
            let to_first = (a[1] - displacement[1]).atan2(a[0] - displacement[0]);
            let to_old = (-displacement[1]).atan2(-displacement[0]);
            self.signposts[s.h] = wrap(back + (to_first - to_old));
        }
        for n in 0..star.len().saturating_sub(1) {
            self.signposts[star[n + 1].h] = self.signpost_after(star[n].h);
        }
        // Incoming: each from the side of the link that leaves the neighbour.
        for s in &star {
            let t = self.twin[s.h];
            if t != NONE {
                self.signposts[t as usize] = self.signpost_after(self.next(s.h));
            }
        }
        Ok(VertexMove {
            vertex: v as u32,
            from,
            to,
            distance: r,
            back_direction: back,
            snapped,
            within_height_bound: r < min_height,
            min_shape: thinnest.0,
        })
    }

    /// The vector from `v` to the point the optimal Delaunay update (Chen and Holst 2011,
    /// Eq. 4.13, step 1) moves it to, in the frame of `v`: the mean of the vectors from `v` to
    /// the circumcentres of the triangles of its star, weighted by their areas (paper Sec.
    /// 4.3). With `barycenters_at_boundary`, a triangle that has a vertex on the boundary
    /// contributes its barycentre instead (Chen and Holst, Sec. 4). The vectors are those of
    /// the planar layout of each triangle around `v`, which for an inserted vertex (angle sum
    /// `2π`, or `π` on the boundary) is the unfolding of its star.
    ///
    /// # Errors
    /// [`IntrinsicError::IndexOutOfRange`]; [`IntrinsicError::MoveRefused`] for a vertex of the
    /// input mesh (its frame is scaled and it cannot be moved) and for a vertex with an edge
    /// to itself.
    pub fn odt_target(
        &self,
        v: usize,
        barycenters_at_boundary: bool,
    ) -> Result<[f64; 2], IntrinsicError> {
        use crate::formulas::{
            odt_displacement, weighted_barycenter_vector, weighted_circumcenter_vector,
        };
        if v >= self.n_vertices() {
            return Err(IntrinsicError::IndexOutOfRange {
                what: "vertex",
                index: v,
                bound: self.n_vertices(),
            });
        }
        let refused = |reason| IntrinsicError::MoveRefused {
            vertex: v as u32,
            reason,
        };
        if v < self.n_input_vertices() {
            return Err(refused(MoveRefusal::InputVertex));
        }
        if let Some(h) = self.outgoing_halfedges(v).find(|&h| self.to_vertex(h) == v) {
            return Err(refused(MoveRefusal::SelfEdge {
                edge: self.halfedge_edge[h],
            }));
        }
        let mut weighted = [0.0; 2];
        let mut total = 0.0;
        for s in self.star(v) {
            let [la, lb, lc] = s.lengths;
            let area = 0.5 * s.double_area;
            let f = s.h / 3;
            let touches_boundary = barycenters_at_boundary
                && self.faces[f]
                    .iter()
                    .any(|&u| self.is_boundary_vertex(u as usize));
            let w = if touches_boundary {
                weighted_barycenter_vector(la, lb, lc, area)
            } else {
                weighted_circumcenter_vector(la, lb, lc, area)
            };
            let (sin, cos) = s.phi.sin_cos();
            weighted[0] += cos * w[0] - sin * w[1];
            weighted[1] += sin * w[0] + cos * w[1];
            total += area;
        }
        Ok(odt_displacement(weighted, total, 1.0))
    }
}

/// Twice the area over the product of the longer of the first two sides (the edges at the
/// moved vertex) and the longest side.
fn shape_number(double_area: f64, l: [f64; 3]) -> f64 {
    let at_vertex = l[0].max(l[1]);
    double_area / (at_vertex * at_vertex.max(l[2]))
}

