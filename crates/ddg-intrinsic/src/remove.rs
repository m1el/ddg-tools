//! Removal of an inserted vertex of any degree (Gillespie, Sharp, Crane 2021, Alg. 8; SPEC
//! N51, N53), shared by the representations so that they flip the same edges.
//!
//! Edges at the vertex are flipped until its degree is three, then its star is merged into one
//! face ([`Complex::remove_vertex`]). **The choice of flip** is that of the reference's
//! `removeInsertedVertex` (IC-SPEC I8, I24): among the edges at the vertex that can be flipped
//! ([`plan_flip`]) and whose flip lowers the vertex's degree, an edge from the vertex to itself
//! is preferred; then the largest score `min(A₁, A₂) / (A₁ + A₂)` of the two triangles the flip
//! creates; ties by the smaller edge id ([`removal_flip_choice`]). Every such flip lowers the
//! degree, so at most `deg − 3` flips happen (the bound `10 deg` is kept as a guard).
//!
//! **Straight stars** (N53, IC-SPEC I27). A star in which no flip lowers the degree used to be
//! refused with [`RemovalRefusal::Stuck`]. Now, if the vertex has two edges on one straight
//! line ([`straight_pair`], as a vertex made by an edge split has), the flips are undone and the
//! vertex is removed by merging those two edges back into one ([`Complex::merge_edge`]): the
//! other edges at it are flipped by the same choice, never the two halves, until its degree is
//! four, then the halves are merged ([`plan_edge_merge`](crate::merge::plan_edge_merge)). An inserted vertex on the boundary
//! (its two boundary edges are always straight) is removed the same way, flipping down to
//! degree three. Where this fails too, the refusal is the one before (`Stuck`), the
//! triangulation restored. Removals that Alg. 8 completes are unchanged.

use crate::complex::{Complex, ComplexFlipUndo, EdgeMerge, VertexRemoval};
use crate::error::{IntrinsicError, RemovalRefusal};
use crate::flip::plan_flip;
use crate::formulas::unfold_diamond;
use crate::merge::straight_pair;
use crate::triangulation::IntrinsicTriangulation;

/// The last step of a removal: the star of degree three merged into one face (Alg. 8), or the
/// two halves of a split edge merged into one edge (N53). Both renumber ids the same way (the
/// last elements take the freed ids); the accessors below give what both have.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum RemovalStep {
    /// [`Complex::remove_vertex`].
    Star(VertexRemoval),
    /// [`Complex::merge_edge`].
    Merge(EdgeMerge),
}

impl RemovalStep {
    /// The removed edges (ids before the removal).
    pub fn removed_edges(&self) -> Vec<u32> {
        match self {
            RemovalStep::Star(r) => r.removed_edges.to_vec(),
            RemovalStep::Merge(m) => m
                .removed_edges
                .iter()
                .copied()
                .filter(|&e| e != ddg_core::topology::NONE)
                .collect(),
        }
    }

    /// The faces the removal made (ids after renumbering): the merged face, or the merged faces
    /// on the sides of the merged edge.
    pub fn faces(&self) -> Vec<u32> {
        match self {
            RemovalStep::Star(r) => vec![r.face],
            RemovalStep::Merge(m) => m
                .faces
                .iter()
                .copied()
                .filter(|&f| f != ddg_core::topology::NONE)
                .collect(),
        }
    }

    /// Every surviving halfedge whose id changed, `(old, new)`.
    pub fn halfedge_moves(&self) -> &[(u32, u32)] {
        match self {
            RemovalStep::Star(r) => &r.halfedge_moves,
            RemovalStep::Merge(m) => &m.halfedge_moves,
        }
    }

    /// Every surviving edge whose id changed, `(old, new)`.
    pub fn edge_moves(&self) -> &[(u32, u32)] {
        match self {
            RemovalStep::Star(r) => &r.edge_moves,
            RemovalStep::Merge(m) => &m.edge_moves,
        }
    }

    /// Every surviving face whose id changed, `(old, new)`.
    pub fn face_moves(&self) -> &[(u32, u32)] {
        match self {
            RemovalStep::Star(r) => &r.face_moves,
            RemovalStep::Merge(m) => &m.face_moves,
        }
    }

    /// The vertex that took the freed id, `(old, new)`.
    pub fn vertex_move(&self) -> Option<(u32, u32)> {
        match self {
            RemovalStep::Star(r) => r.vertex_move,
            RemovalStep::Merge(m) => m.vertex_move,
        }
    }

    /// The numbers of halfedges, edges, faces and vertices after the removal.
    pub fn counts(&self) -> [u32; 4] {
        match self {
            RemovalStep::Star(r) => r.counts,
            RemovalStep::Merge(m) => m.counts,
        }
    }

    /// Brings per-halfedge data along (see [`VertexRemoval::apply_halfedge_data`],
    /// [`EdgeMerge::apply_halfedge_data`]).
    pub fn apply_halfedge_data<T: Copy>(&self, data: &mut Vec<T>) {
        match self {
            RemovalStep::Star(r) => r.apply_halfedge_data(data),
            RemovalStep::Merge(m) => m.apply_halfedge_data(data),
        }
    }

    /// Brings per-edge data along.
    pub fn apply_edge_data<T: Copy>(&self, data: &mut Vec<T>) {
        match self {
            RemovalStep::Star(r) => r.apply_edge_data(data),
            RemovalStep::Merge(m) => m.apply_edge_data(data),
        }
    }

    /// Brings per-face data along.
    pub fn apply_face_data<T: Copy>(&self, data: &mut Vec<T>) {
        match self {
            RemovalStep::Star(r) => r.apply_face_data(data),
            RemovalStep::Merge(m) => m.apply_face_data(data),
        }
    }

    /// Brings per-vertex data along.
    pub fn apply_vertex_data<T: Copy>(&self, data: &mut Vec<T>) {
        match self {
            RemovalStep::Star(r) => r.apply_vertex_data(data),
            RemovalStep::Merge(m) => m.apply_vertex_data(data),
        }
    }
}

/// What a removal of an inserted vertex of any degree did
/// ([`Triangulation::remove_inserted_vertex`](crate::Triangulation::remove_inserted_vertex)).
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct InsertedVertexRemoval {
    /// The removed vertex (its id before the removal).
    pub vertex: u32,
    /// Its degree before the flips.
    pub degree: usize,
    /// The edges flipped, in order (ids before the final renumbering). For a removal by a
    /// merge after Alg. 8 was stuck, the flips of the merge only (those of Alg. 8 were undone).
    pub flips: Vec<u32>,
    /// The last step and its renumbering.
    pub removal: RemovalStep,
}

/// The edge at vertex `v` that removal flips next (IC-SPEC I24; see the
/// [module documentation](self)), or `None` if no edge at `v` can be flipped so that the
/// degree of `v` drops. Decided from the connectivity and the lengths alone, so every
/// representation makes the same choice.
pub fn removal_flip_choice(complex: &Complex, lengths: &[f64], v: usize) -> Option<usize> {
    removal_flip_choice_except(complex, lengths, v, &[])
}

/// [`removal_flip_choice`] among the edges at `v` other than `except` (the halves of a merge,
/// N53).
pub fn removal_flip_choice_except(
    complex: &Complex,
    lengths: &[f64],
    v: usize,
    except: &[u32],
) -> Option<usize> {
    let c = complex;
    // (is a loop, score, edge): the best is the largest loop flag, then score, then the
    // smallest edge id.
    let mut best: Option<(bool, f64, usize)> = None;
    for h in c.outgoing_halfedges(v) {
        let e = c.edge(h);
        if except.contains(&(e as u32)) {
            continue;
        }
        let Ok(plan) = plan_flip(c, lengths, e) else {
            continue;
        };
        let h0 = c.edge_halfedge(e);
        let Some(h1) = c.twin(h0) else { continue };
        let (i, j) = (c.from_vertex(h0), c.to_vertex(h0));
        let (k, l) = (c.from_vertex(c.prev(h0)), c.from_vertex(c.prev(h1)));
        let ends = |a: usize, b: usize| usize::from(a == v) + usize::from(b == v);
        if ends(k, l) >= ends(i, j) {
            continue;
        }
        let [l_ij, l_jk, l_ki, l_il, l_lj] = plan.lengths;
        let Some(diamond) = unfold_diamond(l_ij, l_jk, l_ki, l_il, l_lj) else {
            continue;
        };
        let [a1, a2] = diamond.flipped_areas();
        let score = a1.min(a2) / (a1 + a2);
        let candidate = (i == j, score, e);
        let better = match best {
            None => true,
            Some((loop_b, score_b, e_b)) => {
                (candidate.0, candidate.1) > (loop_b, score_b)
                    || (candidate.0 == loop_b && candidate.1 == score_b && e < e_b)
            }
        };
        if better {
            best = Some(candidate);
        }
    }
    best.map(|(_, _, e)| e)
}

/// The plan of a removal by merging (N53) for vertex `v` of a triangulation with this
/// connectivity and these lengths: the two halves to keep, and the degree to flip down to
/// (four, or three on the boundary); `None` if `v` has no two edges on a straight line.
pub fn merge_removal_halves(
    complex: &Complex,
    lengths: &[f64],
    v: usize,
) -> Option<([u32; 2], usize)> {
    let pair = straight_pair(complex, lengths, v)?;
    let target = if complex.is_boundary_vertex(v) { 3 } else { 4 };
    Some((pair, target))
}

impl IntrinsicTriangulation {
    /// Removes the inserted vertex `v` of any degree (Gillespie, Sharp, Crane 2021, Alg. 8):
    /// flips edges at it by [`removal_flip_choice`] until its degree is three, then
    /// [`remove_vertex`](Self::remove_vertex). Where no flip lowers the degree and `v` has two
    /// edges on a straight line, and for an inserted boundary vertex, by merging those two
    /// edges back into one instead (SPEC N53, see the [module documentation](self)). Ids are
    /// renumbered as [`Complex::remove_vertex`] or [`Complex::merge_edge`] say. None of this
    /// crate's own operations calls it, except refinement by
    /// [`RefineRules::Chew`](crate::RefineRules::Chew) (SPEC N51).
    ///
    /// # Errors
    /// [`IntrinsicError::IndexOutOfRange`]; [`IntrinsicError::RemovalRefused`] for an input
    /// vertex, a constrained edge at `v`, a star where no flip lowers the degree and no merge
    /// applies ([`RemovalRefusal::Stuck`]), the refusals of
    /// [`remove_vertex`](Self::remove_vertex), and for a boundary vertex
    /// [`RemovalRefusal::BoundaryVertex`] (its corner angles do not sum to `π`), `Stuck` or the
    /// refusals of [`merge_edge`](Self::merge_edge). The triangulation is unchanged in every case (the flips are undone exactly from a
    /// journal of what they overwrote).
    pub fn remove_inserted_vertex(
        &mut self,
        v: usize,
    ) -> Result<InsertedVertexRemoval, IntrinsicError> {
        let nv = self.n_vertices();
        if v >= nv {
            return Err(IntrinsicError::IndexOutOfRange {
                what: "vertex",
                index: v,
                bound: nv,
            });
        }
        let refused = |reason| IntrinsicError::RemovalRefused {
            vertex: v as u32,
            reason,
        };
        if v < self.n_input_vertices() {
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
            return match merge_removal_halves(&self.complex, &self.edge_lengths, v) {
                Some((halves, target)) => self.remove_by_merge(v, degree, halves, target),
                // The corner angles do not sum to π (or a self-edge): as before N53.
                None => Err(refused(RemovalRefusal::BoundaryVertex)),
            };
        }
        // What the flips overwrite, to undo them if the removal fails (a copy of the whole
        // triangulation would cost O(n) per removal).
        let mut journal: Vec<FlipUndo> = Vec::new();
        let mut flips = Vec::new();
        let result = (|| {
            let bound = 10 * degree;
            loop {
                let d = self.complex.vertex_degree(v);
                if d == 3 {
                    break;
                }
                let choice = (flips.len() < bound)
                    .then(|| removal_flip_choice(&self.complex, &self.edge_lengths, v))
                    .flatten();
                let Some(e) = choice else {
                    return Err(refused(RemovalRefusal::Stuck { degree: d }));
                };
                let undo = self.flip_undo(e);
                self.flip_edge(e)?;
                journal.push(undo);
                flips.push(e as u32);
            }
            self.remove_vertex(v)
        })();
        match result {
            Ok(removal) => Ok(InsertedVertexRemoval {
                vertex: v as u32,
                degree,
                flips,
                removal: RemovalStep::Star(removal),
            }),
            Err(err) => {
                for undo in journal.iter().rev() {
                    self.undo_flip(undo);
                }
                let stuck = matches!(
                    err,
                    IntrinsicError::RemovalRefused {
                        reason: RemovalRefusal::Stuck { .. },
                        ..
                    }
                );
                if stuck
                    && let Some((halves, target)) =
                        merge_removal_halves(&self.complex, &self.edge_lengths, v)
                    && let Ok(removal) = self.remove_by_merge(v, degree, halves, target)
                {
                    return Ok(removal);
                }
                Err(err)
            }
        }
    }

    /// Flips the edges at `v` other than `halves` by [`removal_flip_choice_except`] until the
    /// degree of `v` is `target`, then merges the halves ([`merge_edge`](Self::merge_edge));
    /// everything undone on failure.
    fn remove_by_merge(
        &mut self,
        v: usize,
        degree: usize,
        halves: [u32; 2],
        target: usize,
    ) -> Result<InsertedVertexRemoval, IntrinsicError> {
        let mut journal: Vec<FlipUndo> = Vec::new();
        let mut flips = Vec::new();
        let result = (|| {
            let bound = 10 * degree;
            loop {
                let d = self.complex.vertex_degree(v);
                if d <= target {
                    break;
                }
                let choice = (flips.len() < bound)
                    .then(|| {
                        removal_flip_choice_except(&self.complex, &self.edge_lengths, v, &halves)
                    })
                    .flatten();
                let Some(e) = choice else {
                    return Err(IntrinsicError::RemovalRefused {
                        vertex: v as u32,
                        reason: RemovalRefusal::Stuck { degree: d },
                    });
                };
                let undo = self.flip_undo(e);
                self.flip_edge(e)?;
                journal.push(undo);
                flips.push(e as u32);
            }
            self.merge_edge(v, halves[0] as usize)
        })();
        match result {
            Ok(merge) => Ok(InsertedVertexRemoval {
                vertex: v as u32,
                degree,
                flips,
                removal: RemovalStep::Merge(merge),
            }),
            Err(err) => {
                for undo in journal.iter().rev() {
                    self.undo_flip(undo);
                }
                Err(err)
            }
        }
    }
}

/// What a signpost flip overwrites: the connectivity, the signposts of the six halfedges of
/// its two faces, the edge's length and its input edge.
struct FlipUndo {
    complex: ComplexFlipUndo,
    edge: usize,
    length: f64,
    on_input: u32,
    signposts: [(usize, f64); 6],
}

impl IntrinsicTriangulation {
    /// What flipping edge `e` (interior) overwrites.
    fn flip_undo(&self, e: usize) -> FlipUndo {
        let h0 = self.complex.edge_halfedge[e] as usize;
        let h1 = self.complex.twin[h0] as usize;
        let [f0, f1] = [h0 / 3, h1 / 3];
        let slots = [
            3 * f0,
            3 * f0 + 1,
            3 * f0 + 2,
            3 * f1,
            3 * f1 + 1,
            3 * f1 + 2,
        ];
        FlipUndo {
            complex: self.complex.flip_undo(e),
            edge: e,
            length: self.edge_lengths[e],
            on_input: self.edge_on_input[e],
            signposts: slots.map(|h| (h, self.signposts[h])),
        }
    }

    /// Restores what a flip overwrote.
    fn undo_flip(&mut self, undo: &FlipUndo) {
        self.complex.undo_flip(&undo.complex);
        self.edge_lengths[undo.edge] = undo.length;
        self.edge_on_input[undo.edge] = undo.on_input;
        for &(h, x) in &undo.signposts {
            self.signposts[h] = x;
        }
    }
}
