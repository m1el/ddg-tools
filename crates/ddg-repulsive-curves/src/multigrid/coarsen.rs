//! Coarsening and prolongation of curve networks (`SPEC.md` §8.4, App. B.3.2; D17).
//!
//! # Coarsening [formula: rc-coarsen]
//!
//! Per component, breadth-first from its smallest vertex index, a vertex is **black** (kept) iff
//! the search proposes it (the proposals alternate along the search), or `deg ≠ 2`, or it is
//! *forced* (it carries a `Point` or `Tangent(Vertex)` constraint), or an already explored
//! neighbour is white. Closed loops with fewer than 6 vertices are not coarsened. So no two
//! adjacent vertices are white and every white vertex has two black neighbours (proved in
//! `lean/ddg/DDG/Repulsive/Curves/Multigrid.lean`).
//!
//! The coarse network has the black vertices (numbered in increasing fine index) at their fine
//! positions (the paper's "removing white vertices"; the reference uses `J⁺γ`, D17). Its edges
//! are the fine edges with two black ends and one edge `(a, b)` per white vertex with neighbours
//! `a, b`, in the order in which they first appear along the fine edge list. A white vertex whose
//! coarse edge would join two vertices that are already joined (a fine triangle, or a cycle of
//! four vertices through junctions; the frozen small loops are the pure-loop case of this) is
//! turned black, so that no level gets two edges between the same vertices.
//!
//! The breadth-first rule is the reference's (`ref: src/poly_curve_network.cpp`, `Coarsen`,
//! MIT, Copyright (c) 2019 Christopher Yu); the frozen small loops and the duplicate-edge check
//! are ours.
//!
//! # Prolongation [formula: rc-prolongation]
//!
//! `J ∈ R^{|V_f| × |V_c|}`: 1 at a black vertex, ½ and ½ from the two black neighbours of a white
//! one (paper: "linear interpolation by edge length did not improve"), so `J 𝟙 = 𝟙`. Vectors in
//! the interleaved layout are prolonged with `J ⊗ I_d` and restricted with `Jᵀ ⊗ I_d`.

use std::collections::HashSet;

use ddg_core::sparse::SparseMatrix;

use crate::constraints::{Constraint, ConstraintSet, ResolvedConstraint, TangentSite};
use crate::error::RcError;
use crate::network::CurveNetwork;
use crate::solve::sparse_error;

/// Closed loops with fewer vertices than this are not coarsened.
pub const MIN_COARSENED_LOOP: usize = 6;

/// One row of the prolongation `J`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ProlongationRow {
    /// A black vertex: its coarse index (weight 1).
    Black(u32),
    /// A white vertex: the coarse indices of its two black neighbours (weights ½, ½).
    White(u32, u32),
}

/// The vertices that must stay black on the next level: those with a `Point` or a
/// `Tangent(Vertex)` constraint (`SPEC.md` §8.4). Vertices of degree other than 2 are black by
/// the rule itself.
pub fn forced_black(n_vertices: usize, constraints: Option<&ConstraintSet>) -> Vec<bool> {
    let mut forced = vec![false; n_vertices];
    if let Some(set) = constraints {
        for c in set.constraints() {
            match c {
                ResolvedConstraint::Point { vertex, .. }
                | ResolvedConstraint::Tangent {
                    site: TangentSite::Vertex(vertex),
                    ..
                } => {
                    if let Some(f) = forced.get_mut(*vertex as usize) {
                        *f = true;
                    }
                }
                _ => {}
            }
        }
    }
    forced
}

/// The black/white marking of `SPEC.md` §8.4 before the duplicate-edge check: per component,
/// breadth-first from the smallest vertex; black iff proposed, `deg ≠ 2`, forced, or an explored
/// neighbour is white; closed loops of fewer than [`MIN_COARSENED_LOOP`] vertices all black.
///
/// [formula: rc-coarsen]
///
/// # Errors
/// [`RcError::InvalidOption`] if `forced` does not have one entry per vertex.
pub fn mark_black(network: &CurveNetwork, forced: &[bool]) -> Result<Vec<bool>, RcError> {
    let n = network.n_vertices();
    if forced.len() != n {
        return Err(RcError::InvalidOption {
            name: "forced",
            value: forced.len().to_string(),
            requirement: "one entry per vertex",
        });
    }
    let edges = network.edges();
    let mut black = vec![true; n];
    let mut explored = vec![false; n];
    let mut queue = std::collections::VecDeque::new();
    for comp in network.components() {
        if comp.closed && comp.vertices.len() < MIN_COARSENED_LOOP {
            continue; // all black
        }
        let start = comp.vertices[0];
        queue.push_back((start, true));
        while let Some((v, proposed)) = queue.pop_front() {
            let vi = v as usize;
            if explored[vi] {
                continue;
            }
            explored[vi] = true;
            let inc = network.incident_edges(v);
            let other = |e: u32| {
                let [a, b] = edges[e as usize];
                if a == v { b } else { a }
            };
            let keep = proposed
                || inc.len() != 2
                || forced[vi]
                || inc.iter().any(|&e| {
                    let w = other(e) as usize;
                    explored[w] && !black[w]
                });
            black[vi] = keep;
            for &e in inc {
                let w = other(e);
                if !explored[w as usize] {
                    queue.push_back((w, !keep));
                }
            }
        }
    }
    Ok(black)
}

/// One coarsening step: the coarse network, the prolongation and the maps between the levels.
#[derive(Debug, Clone)]
pub struct Coarsening {
    coarse: CurveNetwork,
    /// Per fine vertex: its coarse index if black.
    fine_to_coarse: Vec<Option<u32>>,
    /// Per coarse vertex: its fine index.
    coarse_to_fine: Vec<u32>,
    /// The rows of `J`, one per fine vertex.
    rows: Vec<ProlongationRow>,
    /// Per fine edge: the coarse edge that contains it.
    edge_to_coarse: Vec<u32>,
    /// Per coarse edge: the fine edges it is made of (one or two).
    coarse_edge_parts: Vec<Vec<u32>>,
}

/// Coarsens `network` once (`SPEC.md` §8.4) with the vertices in `forced` kept black. Returns
/// `None` when no vertex can be removed.
///
/// [formula: rc-coarsen]
///
/// # Errors
/// [`RcError::InvalidOption`] for a `forced` of the wrong length; the errors of
/// [`CurveNetwork::new`] for the coarse network (e.g. [`RcError::NoInteractingPairs`] when the
/// coarse network is too small).
pub fn coarsen(network: &CurveNetwork, forced: &[bool]) -> Result<Option<Coarsening>, RcError> {
    let n = network.n_vertices();
    let edges = network.edges();
    let mut black = mark_black(network, forced)?;
    let neighbours = |v: u32| -> [u32; 2] {
        let inc = network.incident_edges(v);
        let other = |e: u32| {
            let [a, b] = edges[e as usize];
            if a == v { b } else { a }
        };
        [other(inc[0]), other(inc[1])]
    };
    // Duplicate-edge check: a white vertex whose (a, b) is already a coarse edge turns black.
    let key = |a: u32, b: u32| if a < b { (a, b) } else { (b, a) };
    let mut pairs: HashSet<(u32, u32)> = edges
        .iter()
        .filter(|&&[a, b]| black[a as usize] && black[b as usize])
        .map(|&[a, b]| key(a, b))
        .collect();
    for v in 0..n as u32 {
        if black[v as usize] {
            continue;
        }
        let [a, b] = neighbours(v);
        if !pairs.insert(key(a, b)) {
            black[v as usize] = true;
            pairs.insert(key(a, v));
            pairs.insert(key(v, b));
        }
    }
    if black.iter().all(|&b| b) {
        return Ok(None);
    }
    let mut fine_to_coarse = vec![None; n];
    let mut coarse_to_fine = Vec::new();
    for v in 0..n {
        if black[v] {
            fine_to_coarse[v] = Some(coarse_to_fine.len() as u32);
            coarse_to_fine.push(v as u32);
        }
    }
    let c = |v: u32| fine_to_coarse[v as usize].expect("black vertex");
    let rows: Vec<ProlongationRow> = (0..n as u32)
        .map(|v| match fine_to_coarse[v as usize] {
            Some(cv) => ProlongationRow::Black(cv),
            None => {
                let [a, b] = neighbours(v);
                ProlongationRow::White(c(a), c(b))
            }
        })
        .collect();
    let mut coarse_edges: Vec<[u32; 2]> = Vec::new();
    let mut coarse_edge_parts: Vec<Vec<u32>> = Vec::new();
    let mut edge_to_coarse = vec![u32::MAX; edges.len()];
    for (e, &[a, b]) in edges.iter().enumerate() {
        if edge_to_coarse[e] != u32::MAX {
            continue;
        }
        let ce = coarse_edges.len() as u32;
        if black[a as usize] && black[b as usize] {
            coarse_edges.push([c(a), c(b)]);
            coarse_edge_parts.push(vec![e as u32]);
            edge_to_coarse[e] = ce;
        } else {
            let w = if black[a as usize] { b } else { a };
            let inc = network.incident_edges(w);
            let [na, nb] = neighbours(w);
            coarse_edges.push([c(na), c(nb)]);
            coarse_edge_parts.push(inc.to_vec());
            for &f in inc {
                edge_to_coarse[f as usize] = ce;
            }
        }
    }
    let coarse = CurveNetwork::new(coarse_to_fine.len(), coarse_edges)?;
    Ok(Some(Coarsening {
        coarse,
        fine_to_coarse,
        coarse_to_fine,
        rows,
        edge_to_coarse,
        coarse_edge_parts,
    }))
}

impl Coarsening {
    /// The coarse network.
    pub fn coarse(&self) -> &CurveNetwork {
        &self.coarse
    }

    /// Number of fine vertices.
    pub fn n_fine(&self) -> usize {
        self.rows.len()
    }

    /// Number of coarse vertices.
    pub fn n_coarse(&self) -> usize {
        self.coarse_to_fine.len()
    }

    /// Per fine vertex: its coarse index if it is black.
    pub fn fine_to_coarse(&self) -> &[Option<u32>] {
        &self.fine_to_coarse
    }

    /// Per coarse vertex: the fine vertex it is.
    pub fn coarse_to_fine(&self) -> &[u32] {
        &self.coarse_to_fine
    }

    /// The rows of `J`, one per fine vertex.
    pub fn rows(&self) -> &[ProlongationRow] {
        &self.rows
    }

    /// Per fine edge: the coarse edge that contains it.
    pub fn edge_to_coarse(&self) -> &[u32] {
        &self.edge_to_coarse
    }

    /// Per coarse edge: the fine edges it replaces (one, or the two edges of a white vertex).
    pub fn coarse_edge_parts(&self) -> &[Vec<u32>] {
        &self.coarse_edge_parts
    }

    /// The coarse positions: the fine positions of the black vertices (D17).
    ///
    /// # Errors
    /// [`RcError::PositionCountMismatch`].
    pub fn coarse_positions(&self, fine: &[[f64; 3]]) -> Result<Vec<[f64; 3]>, RcError> {
        if fine.len() != self.n_fine() {
            return Err(RcError::PositionCountMismatch {
                expected: self.n_fine(),
                found: fine.len(),
            });
        }
        Ok(self
            .coarse_to_fine
            .iter()
            .map(|&v| fine[v as usize])
            .collect())
    }

    /// `(J ⊗ I_d) x` for a coarse vector `x` with `d` entries per vertex (interleaved).
    ///
    /// [formula: rc-prolongation]
    ///
    /// # Errors
    /// [`RcError::InvalidOption`] for `d = 0` or a vector of the wrong length.
    pub fn prolong(&self, coarse: &[f64], d: usize) -> Result<Vec<f64>, RcError> {
        self.check_len(coarse.len(), self.n_coarse(), d)?;
        Ok(self.prolong_unchecked(coarse, d))
    }

    pub(crate) fn prolong_unchecked(&self, coarse: &[f64], d: usize) -> Vec<f64> {
        let mut out = vec![0.0; self.n_fine() * d];
        for (v, row) in self.rows.iter().enumerate() {
            let o = &mut out[v * d..(v + 1) * d];
            match *row {
                ProlongationRow::Black(c) => {
                    o.copy_from_slice(&coarse[c as usize * d..(c as usize + 1) * d]);
                }
                ProlongationRow::White(a, b) => {
                    for s in 0..d {
                        o[s] = 0.5 * (coarse[a as usize * d + s] + coarse[b as usize * d + s]);
                    }
                }
            }
        }
        out
    }

    /// `(Jᵀ ⊗ I_d) y` for a fine vector `y` (summed in increasing fine vertex order).
    ///
    /// [formula: rc-prolongation]
    ///
    /// # Errors
    /// [`RcError::InvalidOption`] for `d = 0` or a vector of the wrong length.
    pub fn restrict(&self, fine: &[f64], d: usize) -> Result<Vec<f64>, RcError> {
        self.check_len(fine.len(), self.n_fine(), d)?;
        Ok(self.restrict_unchecked(fine, d))
    }

    pub(crate) fn restrict_unchecked(&self, fine: &[f64], d: usize) -> Vec<f64> {
        let mut out = vec![0.0; self.n_coarse() * d];
        for (v, row) in self.rows.iter().enumerate() {
            let y = &fine[v * d..(v + 1) * d];
            match *row {
                ProlongationRow::Black(c) => {
                    for s in 0..d {
                        out[c as usize * d + s] += y[s];
                    }
                }
                ProlongationRow::White(a, b) => {
                    for s in 0..d {
                        out[a as usize * d + s] += 0.5 * y[s];
                        out[b as usize * d + s] += 0.5 * y[s];
                    }
                }
            }
        }
        out
    }

    fn check_len(&self, len: usize, n: usize, d: usize) -> Result<(), RcError> {
        if d == 0 || len != n * d {
            return Err(RcError::InvalidOption {
                name: "multigrid vector",
                value: format!("{len} entries with d = {d}"),
                requirement: "d ≥ 1 entries per vertex of the level",
            });
        }
        Ok(())
    }

    /// The scalar prolongation `J` as a sparse `|V_f| × |V_c|` matrix.
    ///
    /// [formula: rc-prolongation]
    ///
    /// # Errors
    /// Resource errors of the sparse matrix.
    pub fn prolongation_matrix(&self) -> Result<SparseMatrix, RcError> {
        let mut t = Vec::with_capacity(2 * self.n_fine());
        for (v, row) in self.rows.iter().enumerate() {
            match *row {
                ProlongationRow::Black(c) => t.push((v as u32, c, 1.0)),
                ProlongationRow::White(a, b) => {
                    t.push((v as u32, a, 0.5));
                    t.push((v as u32, b, 0.5));
                }
            }
        }
        SparseMatrix::from_triplets(self.n_fine(), self.n_coarse(), &t)
            .map_err(|e| sparse_error(e, "prolongation"))
    }

    /// The constraints of the coarse level (`SPEC.md` §8.4), as a list for
    /// [`ConstraintSet::new`] on the coarse network: `Barycenter`, `TotalLength`; `EdgeLengths`
    /// on all coarse edges when the fine level constrains all its edges, else on the coarse edges
    /// whose fine edges were all constrained; `Point` and `Tangent(Vertex)` on the surviving
    /// vertex; `Tangent(Edge)` on the coarse edge that contains the fine edge; `Surface` on the
    /// surviving constrained vertices (the reference drops surface rows on coarse levels).
    /// Targets are left `None` (the coarse levels only use the Jacobian). Constraints with no
    /// row left are omitted.
    pub fn coarse_constraints(&self, fine: &ConstraintSet) -> Vec<Constraint> {
        let m_fine = self.edge_to_coarse.len();
        let mut out = Vec::new();
        for c in fine.constraints() {
            match c {
                ResolvedConstraint::Barycenter { .. } => {
                    out.push(Constraint::Barycenter { target: None });
                }
                ResolvedConstraint::TotalLength { .. } => {
                    out.push(Constraint::TotalLength { target: None });
                }
                ResolvedConstraint::EdgeLengths { edges, .. } => {
                    let mut listed = vec![false; m_fine];
                    for &e in edges {
                        if let Some(l) = listed.get_mut(e as usize) {
                            *l = true;
                        }
                    }
                    if listed.iter().all(|&l| l) {
                        out.push(Constraint::EdgeLengths {
                            edges: None,
                            targets: None,
                        });
                    } else {
                        let list: Vec<u32> = self
                            .coarse_edge_parts
                            .iter()
                            .enumerate()
                            .filter(|(_, parts)| parts.iter().all(|&f| listed[f as usize]))
                            .map(|(ce, _)| ce as u32)
                            .collect();
                        if !list.is_empty() {
                            out.push(Constraint::EdgeLengths {
                                edges: Some(list),
                                targets: None,
                            });
                        }
                    }
                }
                ResolvedConstraint::Point { vertex, .. } => {
                    if let Some(Some(cv)) = self.fine_to_coarse.get(*vertex as usize) {
                        out.push(Constraint::Point {
                            vertex: *cv,
                            target: None,
                        });
                    }
                }
                ResolvedConstraint::Surface { vertices, surface } => {
                    let list: Vec<u32> = vertices
                        .iter()
                        .filter_map(|&v| self.fine_to_coarse.get(v as usize).copied().flatten())
                        .collect();
                    if !list.is_empty() {
                        out.push(Constraint::Surface {
                            vertices: Some(list),
                            surface: surface.clone(),
                        });
                    }
                }
                ResolvedConstraint::Tangent {
                    site, direction, ..
                } => {
                    let site = match *site {
                        TangentSite::Vertex(v) => self
                            .fine_to_coarse
                            .get(v as usize)
                            .copied()
                            .flatten()
                            .map(TangentSite::Vertex),
                        TangentSite::Edge(e) => self
                            .edge_to_coarse
                            .get(e as usize)
                            .map(|&ce| TangentSite::Edge(ce)),
                    };
                    if let Some(site) = site {
                        out.push(Constraint::Tangent {
                            site,
                            direction: Some(*direction),
                        });
                    }
                }
            }
        }
        out
    }
}
