//! Removal of an inserted vertex (paper Sec. 3.5, Alg. 8, Thm. D.1; SPEC §4.6, I7, I8, I24,
//! I27).
//!
//! Edges at the vertex are flipped ([`IntegerTriangulation::flip_edge`], Eq. 3, 4) until its
//! degree is three, then its star is merged into one face
//! ([`Complex::remove_vertex`](ddg_intrinsic::Complex::remove_vertex)); no integer changes
//! beyond the flips: the three outer edges keep their lengths and normal coordinates and the
//! outer halfedges their roundabouts, and the renumbering of
//! [`VertexRemoval`] moves `n`, `r`, the lengths and
//! `q⁰`.
//!
//! **The choice of flip** (I8, I24, as the reference's `removeInsertedVertex`): among the edges
//! at the vertex that can be flipped ([`ddg_intrinsic::plan_flip`]) and whose flip lowers the
//! vertex's degree, an edge from the vertex to itself is preferred; then the largest score
//! `min(A₁, A₂) / (A₁ + A₂)` of the two triangles the flip creates (twice their signed areas,
//! [`Diamond::flipped_areas`](ddg_intrinsic::formulas::Diamond::flipped_areas)); ties by the
//! smaller edge id. Every flip lowers the degree, so at most `deg − 3` flips happen (the bound
//! `10 deg` of I24 is kept as a guard).
//!
//! **Straight stars (I27, resolved).** An interior vertex made by an edge split has degree four
//! with two angles of `π`. The two edges to the opposite corners then have a quadrilateral with
//! an angle of `π` at the vertex and are never flipped (the flip would create a triangle
//! without area, which [`plan_flip`](ddg_intrinsic::plan_flip) refuses as too thin); the two
//! halves of the split edge flip only where their quadrilateral is convex at its far end. When
//! Alg. 8 removes the vertex, nothing changes. When it is stuck and the vertex has two edges on
//! one straight line ([`ddg_intrinsic::straight_pair`]), its flips are undone, the other edges
//! at the vertex are flipped by the same choice (never the two halves) until its degree is
//! four, and the halves are merged back into one edge ([`IntegerTriangulation::merge_edge`],
//! exact integers, `crate::merge`). An inserted vertex on the boundary (its two boundary edges
//! are straight; I7) is removed the same way, flipping down to degree three. The decisions are
//! `ddg_intrinsic`'s ([`merge_removal_halves`], [`removal_flip_choice_except`]), so the
//! signposts remove the same vertices by the same flips and merges.
//!
//! **Refusals.** An input vertex ([`RemovalRefusal::InputVertex`]), a vertex with a constrained
//! edge ([`RemovalRefusal::ConstrainedEdge`]): refused before anything changes. A boundary
//! vertex whose corner angles do not sum to `π` ([`RemovalRefusal::BoundaryVertex`]). A star
//! where no edge can be flipped to lower the degree and no merge applies:
//! [`IntegerError::RemovalStuck`]; a star of degree three that cannot be merged (an edge from
//! the vertex to itself, a merged triangle without area): the refusal of
//! `Complex::remove_vertex`; in all cases every flip done is undone and the state is
//! **exactly** what it was (from a journal of what the flips overwrote: the connectivity by
//! [`Complex::undo_flip`](ddg_intrinsic::Complex::undo_flip), the integers and lengths here; a
//! copy of the whole complex, as before, cost `O(n)` per removal).

use ddg_core::geometry::triangle_area;
use ddg_intrinsic::complex::{ComplexFlipUndo, VertexRemoval};
use ddg_intrinsic::remove::{merge_removal_halves, removal_flip_choice_except};
use ddg_intrinsic::{IntrinsicError, RemovalRefusal, RemovalStep};

use crate::coords::IntegerTriangulation;
use crate::error::IntegerError;

/// What [`IntegerTriangulation::remove_vertex`] did.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Removal {
    /// The removed vertex (its id before the removal).
    pub vertex: u32,
    /// Its degree before the flips.
    pub degree: usize,
    /// The edges flipped, in order (ids before the removal's renumbering). For a removal by a
    /// merge after Alg. 8 was stuck, the flips of the merge only (those of Alg. 8 were undone).
    pub flips: Vec<u32>,
    /// The last step and its renumbering: the star merged into one face (Alg. 8), or the two
    /// halves merged into one edge (I27).
    pub removal: RemovalStep,
}

/// What a flip overwrote, to undo it.
struct Undo {
    /// What the flip overwrote in the connectivity.
    complex: ComplexFlipUndo,
    edge: usize,
    normal: i64,
    length: f64,
    /// The roundabouts of the six halfedges of the two faces of the flip.
    roundabouts: [(usize, u32); 6],
}

impl IntegerTriangulation {
    /// Removes the inserted vertex `v` (paper Alg. 8): flips edges at it until its degree is
    /// three, then merges its three faces; where Alg. 8 is stuck at a straight star, and for an
    /// inserted boundary vertex, by merging two halves instead (I27). See the
    /// [module documentation](crate::remove) for the choice of flips, the refusals and I27. Ids
    /// are renumbered as
    /// [`Complex::remove_vertex`](ddg_intrinsic::Complex::remove_vertex) or
    /// [`Complex::merge_edge`](ddg_intrinsic::Complex::merge_edge) say (the last vertex, edges
    /// and faces take the freed ids).
    ///
    /// # Errors
    /// [`IntegerError::IndexOutOfRange`]; [`IntegerError::Intrinsic`] with
    /// [`IntrinsicError::RemovalRefused`] for an input vertex, a constrained edge at `v`, a
    /// degree-three star that cannot be merged, a boundary vertex that cannot be merged;
    /// [`IntegerError::RemovalStuck`] when no edge can be flipped and no merge applies; the
    /// errors of [`flip_edge`](Self::flip_edge). The state is unchanged in every case.
    pub fn remove_vertex(&mut self, v: usize) -> Result<Removal, IntegerError> {
        let nv = self.complex.n_vertices();
        if v >= nv {
            return Err(IntegerError::IndexOutOfRange {
                what: "vertex",
                index: v,
                bound: nv,
            });
        }
        let refused = |reason| {
            IntegerError::Intrinsic(IntrinsicError::RemovalRefused {
                vertex: v as u32,
                reason,
            })
        };
        if self.is_input_vertex(v) {
            return Err(refused(RemovalRefusal::InputVertex));
        }
        if let Some(h) = self
            .complex
            .outgoing_halfedges(v)
            .find(|&h| self.complex.is_constrained(self.complex.edge(h)))
        {
            return Err(refused(RemovalRefusal::ConstrainedEdge {
                edge: self.complex.edge(h) as u32,
            }));
        }
        let degree = self.complex.vertex_degree(v);
        if self.complex.is_boundary_vertex(v) {
            return match merge_removal_halves(&self.complex, &self.lengths, v) {
                Some((halves, target)) => self.remove_by_merge(v, degree, halves, target),
                None => Err(refused(RemovalRefusal::BoundaryVertex)),
            };
        }
        let mut journal: Vec<Undo> = Vec::new();
        match self.flip_down_and_remove(v, degree, &mut journal) {
            Ok(removal) => Ok(Removal {
                vertex: v as u32,
                degree,
                flips: journal.iter().map(|u| u.edge as u32).collect(),
                removal: RemovalStep::Star(removal),
            }),
            Err(err) => {
                self.undo(&journal);
                if matches!(err, IntegerError::RemovalStuck { .. })
                    && let Some((halves, target)) =
                        merge_removal_halves(&self.complex, &self.lengths, v)
                    && let Ok(removal) = self.remove_by_merge(v, degree, halves, target)
                {
                    return Ok(removal);
                }
                Err(err)
            }
        }
    }

    /// Undoes the journaled flips, last first.
    fn undo(&mut self, journal: &[Undo]) {
        for undo in journal.iter().rev() {
            self.complex.undo_flip(&undo.complex);
            self.normal[undo.edge] = undo.normal;
            self.lengths[undo.edge] = undo.length;
            for &(h, r) in &undo.roundabouts {
                self.roundabouts[h] = r;
            }
        }
    }

    /// Flips edge `e` (interior), journaling what the flip overwrites first; on an error
    /// nothing was changed by this flip and its journal entry restores the same values.
    fn journaled_flip(&mut self, e: usize, journal: &mut Vec<Undo>) -> Result<(), IntegerError> {
        let h0 = self.complex.edge_halfedge(e);
        let h1 = self.complex.twin(h0).expect("a flippable edge is interior");
        let slots = [h0 / 3, h1 / 3].map(|f| [3 * f, 3 * f + 1, 3 * f + 2]);
        let roundabouts = [
            slots[0][0],
            slots[0][1],
            slots[0][2],
            slots[1][0],
            slots[1][1],
            slots[1][2],
        ]
        .map(|h| (h, self.roundabouts[h]));
        journal.push(Undo {
            complex: self.complex.flip_undo(e),
            edge: e,
            normal: self.normal[e],
            length: self.lengths[e],
            roundabouts,
        });
        self.flip_edge(e)?;
        Ok(())
    }

    /// The flips and the merge of [`remove_vertex`](Self::remove_vertex), journaling what the
    /// flips overwrite; on an error the caller restores.
    fn flip_down_and_remove(
        &mut self,
        v: usize,
        degree: usize,
        journal: &mut Vec<Undo>,
    ) -> Result<VertexRemoval, IntegerError> {
        let bound = 10 * degree;
        loop {
            let d = self.complex.vertex_degree(v);
            if d == 3 {
                break;
            }
            let Some(e) = (journal.len() < bound)
                .then(|| self.removal_flip_choice(v))
                .flatten()
            else {
                return Err(IntegerError::RemovalStuck {
                    vertex: v as u32,
                    degree: d,
                });
            };
            self.journaled_flip(e, journal)?;
        }
        let refused = |reason| {
            IntegerError::Intrinsic(IntrinsicError::RemovalRefused {
                vertex: v as u32,
                reason,
            })
        };
        let out: Vec<usize> = self.complex.outgoing_halfedges(v).collect();
        if let Some(&h) = out.iter().find(|&&h| self.complex.to_vertex(h) == v) {
            return Err(refused(RemovalRefusal::SelfEdge {
                edge: self.complex.edge(h) as u32,
            }));
        }
        let lengths =
            [out[0], out[1], out[2]].map(|h| self.lengths[self.complex.edge(self.complex.next(h))]);
        if triangle_area(lengths[0], lengths[1], lengths[2]).is_none() {
            return Err(refused(RemovalRefusal::Degenerate { lengths }));
        }
        let removal = self.complex.remove_vertex(v)?;
        removal.apply_halfedge_data(&mut self.roundabouts);
        removal.apply_edge_data(&mut self.lengths);
        removal.apply_edge_data(&mut self.normal);
        removal.apply_vertex_data(&mut self.locations);
        Ok(removal)
    }

    /// Removal by a merge (I27): flips the edges at `v` other than `halves` by
    /// `ddg_intrinsic`'s [`removal_flip_choice_except`] until the degree of `v` is `target`
    /// (four, three on the boundary), then [`merge_edge`](Self::merge_edge) of the halves;
    /// everything undone on failure.
    fn remove_by_merge(
        &mut self,
        v: usize,
        degree: usize,
        halves: [u32; 2],
        target: usize,
    ) -> Result<Removal, IntegerError> {
        let mut journal: Vec<Undo> = Vec::new();
        let result = (|| {
            let bound = 10 * degree;
            loop {
                let d = self.complex.vertex_degree(v);
                if d <= target {
                    break;
                }
                let Some(e) = (journal.len() < bound)
                    .then(|| removal_flip_choice_except(&self.complex, &self.lengths, v, &halves))
                    .flatten()
                else {
                    return Err(IntegerError::RemovalStuck {
                        vertex: v as u32,
                        degree: d,
                    });
                };
                self.journaled_flip(e, &mut journal)?;
            }
            self.merge_edge(v, halves[0] as usize)
        })();
        match result {
            Ok(merge) => Ok(Removal {
                vertex: v as u32,
                degree,
                flips: journal.iter().map(|u| u.edge as u32).collect(),
                removal: RemovalStep::Merge(merge),
            }),
            Err(err) => {
                self.undo(&journal);
                Err(err)
            }
        }
    }

    /// The edge at `v` that removal flips next (I24), or `None` if no edge at `v` can be
    /// flipped so that the degree of `v` drops: `ddg_intrinsic`'s
    /// [`removal_flip_choice`](ddg_intrinsic::removal_flip_choice), which the signpost
    /// representation uses too, so both remove a vertex by the same flips (I11, WP6).
    fn removal_flip_choice(&self, v: usize) -> Option<usize> {
        ddg_intrinsic::removal_flip_choice(&self.complex, &self.lengths, v)
    }
}
