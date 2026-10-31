//! Removal of an inserted vertex of any degree (Gillespie, Sharp, Crane 2021, Alg. 8; SPEC
//! N51), shared by the representations so that they flip the same edges.
//!
//! Edges at the vertex are flipped until its degree is three, then its star is merged into one
//! face ([`Complex::remove_vertex`]). **The choice of flip** is that of the reference's
//! `removeInsertedVertex` (IC-SPEC I8, I24): among the edges at the vertex that can be flipped
//! ([`plan_flip`]) and whose flip lowers the vertex's degree, an edge from the vertex to itself
//! is preferred; then the largest score `min(A₁, A₂) / (A₁ + A₂)` of the two triangles the flip
//! creates; ties by the smaller edge id ([`removal_flip_choice`]). Every such flip lowers the
//! degree, so at most `deg − 3` flips happen (the bound `10 deg` is kept as a guard). A star in
//! which no flip lowers the degree is refused with [`RemovalRefusal::Stuck`], everything
//! restored.

use crate::complex::{Complex, ComplexFlipUndo, VertexRemoval};
use crate::error::{IntrinsicError, RemovalRefusal};
use crate::flip::plan_flip;
use crate::formulas::unfold_diamond;
use crate::triangulation::IntrinsicTriangulation;

/// What a removal of an inserted vertex of any degree did
/// ([`Triangulation::remove_inserted_vertex`](crate::Triangulation::remove_inserted_vertex)).
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct InsertedVertexRemoval {
    /// The removed vertex (its id before the removal).
    pub vertex: u32,
    /// Its degree before the flips.
    pub degree: usize,
    /// The edges flipped, in order (ids before the final renumbering).
    pub flips: Vec<u32>,
    /// The renumbering of the final step (see [`VertexRemoval`]).
    pub removal: VertexRemoval,
}

/// The edge at vertex `v` that removal flips next (IC-SPEC I24; see the
/// [module documentation](self)), or `None` if no edge at `v` can be flipped so that the
/// degree of `v` drops. Decided from the connectivity and the lengths alone, so every
/// representation makes the same choice.
pub fn removal_flip_choice(complex: &Complex, lengths: &[f64], v: usize) -> Option<usize> {
    let c = complex;
    // (is a loop, score, edge): the best is the largest loop flag, then score, then the
    // smallest edge id.
    let mut best: Option<(bool, f64, usize)> = None;
    for h in c.outgoing_halfedges(v) {
        let e = c.edge(h);
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

impl IntrinsicTriangulation {
    /// Removes the inserted interior vertex `v` of any degree (Gillespie, Sharp, Crane 2021,
    /// Alg. 8): flips edges at it by [`removal_flip_choice`] until its degree is three, then
    /// [`remove_vertex`](Self::remove_vertex). Ids are renumbered as
    /// [`Complex::remove_vertex`] says. None of this crate's own operations calls it, except
    /// refinement by [`RefineRules::Chew`](crate::RefineRules::Chew) (SPEC N51).
    ///
    /// # Errors
    /// [`IntrinsicError::IndexOutOfRange`]; [`IntrinsicError::RemovalRefused`] for an input
    /// vertex, a boundary vertex, a constrained edge at `v`, a star where no flip lowers the
    /// degree ([`RemovalRefusal::Stuck`]) and the refusals of
    /// [`remove_vertex`](Self::remove_vertex). The triangulation is unchanged in every case
    /// (the flips are undone exactly from a journal of what they overwrote).
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
        if self.complex.is_boundary_vertex(v) {
            return Err(refused(RemovalRefusal::BoundaryVertex));
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
                removal,
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
