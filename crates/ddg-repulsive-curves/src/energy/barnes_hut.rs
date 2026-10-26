//! Barnes–Hut energy and differential (`SPEC.md` §8.2, App. B.1.2, Eq. 26–27, D15, D25, D26).
//!
//! The edges are the elements of a [`Bvh<6>`] over their tangent-points `(x_I, T_I) ∈ R⁶` with
//! mass `ℓ_I` (`SPEC.md` §8.1). For every edge `I` the tree is traversed from the root:
//!
//! * an **internal node** `N` that is admissible for `I` [formula: rc-bh-admissibility]
//!
//!   ```text
//!   (guard) N's endpoint box contains neither γ_{i₁} nor γ_{i₂},
//!   r_x^N ≤ θ_x |x_I − x̄_N|,    r_T^N ≤ θ_T
//!   ```
//!
//!   is replaced by one lumped term `ℓ_I L_N k(x_I, x̄_N, T_I)` (Eq. 26);
//! * otherwise its children are visited; **leaves are always exact**:
//!   `Σ_{J ∈ leaf, (I,J) ∈ 𝒫} k̂_IJ ℓ_I ℓ_J` (Eq. 18).
//!
//! The endpoint-box guard makes the exclusion of neighbouring edges exact: an edge sharing a
//! vertex with `I` has that vertex in the endpoint box of every node that contains it, so no
//! such node is ever lumped (the paper relies on `r_x/d` being small instead).
//!
//! The lumped data are `L_N = Σ ℓ_I`, `x̄_N = Σ ℓ_I x_I / L_N` and the **normalized** averaged
//! tangent `T̂_N = Σ ℓ_I T_I / |Σ ℓ_I T_I|` (the reference normalizes; the paper's `T̄_N` does not,
//! D15) [formula: rc-bvh-lumped]. The structure (a tree over tangent-points, depth-first
//! traversal with a per-node opening test, exact leaves) follows the authors' MIT code
//! (`src/spatial/tpe_bvh.cpp`, Copyright (c) 2019 Christopher Yu); the energy is the edge-pair
//! form of `SPEC.md` §3 (D1), not the reference's vertex form, and the criteria are those of
//! Eq. 27 plus the guard (D15).
//!
//! **Differential** [formula: rc-bh-differential]. Lumping the pairs `(I, J)` and `(J, I)` for
//! `J ∈ N` gives `ℓ_I L_N (k(x_I, x̄_N, T_I) + k(x̄_N, x_I, T̂_N))` with `L_N, x̄_N, T̂_N` held
//! fixed; its derivative with respect to the two endpoints of `I` is added for an admissible
//! node, and the exact `G_{I,a}` of `SPEC.md` §3.2 restricted to the leaf's edges for a leaf.
//! With `s₁ = −1, s₂ = +1`, `g_T = ∇_T k(x_I, x̄_N, T_I)`:
//!
//! ```text
//! ∂/∂γ_{i_a} = s_a L_N (k₁ T_I + g_T + k₂ T_I) + ½ ℓ_I L_N (∇_p k(x_I, x̄_N, T_I) − ∇_p k(x̄_N, x_I, T̂_N))
//! ```
//!
//! (the second kernel's argument `x_I` is its `q`, hence `∇_q = −∇_p`.) This is the exact
//! derivative of the lumped energy, not the derivative of the *evaluated* Barnes–Hut energy: the
//! latter is the energy of Eq. 26 only (lumping the pairs `(I, J)`), whose derivative would miss
//! the `(J, I)` terms. The two are different approximations of the same differential; the test
//! suite documents the tolerance between finite differences of one and the other.
//!
//! **Determinism (D25).** The tree is built and refitted sequentially; each edge's traversal is
//! a sequential depth-first walk; rows are computed in parallel into per-edge slots and reduced
//! in increasing edge index, exactly as in [`super::exact`]. The result is bitwise independent of
//! the thread count.
//!
//! **Line searches (D26).** [`BhTree::refit`] keeps the topology and recomputes the node data
//! only, so every trial energy of one line search uses the same approximation.

use rayon::prelude::*;

use crate::error::RcError;
use crate::kernel::v3::{axpy, dot, norm, scale, sub};
use crate::kernel::{Exponents, kernel, kernel_with_gradient};
use crate::network::CurveNetwork;
use crate::spatial::aabb::Aabb;
use crate::spatial::bvh::{Bvh, DEFAULT_MIN_SPLIT_FRACTION, Element};

use super::BhOptions;

impl BhOptions {
    /// Checks the options.
    ///
    /// # Errors
    /// [`RcError::InvalidOption`] for a non-finite or non-positive `theta_x` or `theta_t`
    /// (`0` would make the evaluation exact and slow; use [`super::EnergyEval::Exact`]), and for
    /// `leaf_size = 0`.
    pub fn validate(&self) -> Result<(), RcError> {
        let bad = |name: &'static str, value: String, requirement: &'static str| {
            Err(RcError::InvalidOption {
                name,
                value,
                requirement,
            })
        };
        if !(self.theta_x.is_finite() && self.theta_x > 0.0) {
            return bad(
                "theta_x",
                self.theta_x.to_string(),
                "a finite number greater than 0",
            );
        }
        if !(self.theta_t.is_finite() && self.theta_t > 0.0) {
            return bad(
                "theta_t",
                self.theta_t.to_string(),
                "a finite number greater than 0",
            );
        }
        if self.leaf_size == 0 {
            return bad("leaf_size", "0".into(), "a leaf holds at least 1 edge");
        }
        Ok(())
    }
}

/// Work counters of one evaluation.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub struct BhStats {
    /// Ordered edge pairs `(I, J) ∈ 𝒫` evaluated exactly (in leaves).
    pub pairs_exact: u64,
    /// Lumped (admissible) nodes used, summed over all edges `I`.
    pub nodes_lumped: u64,
}

/// An edge as an element of the tree: the tangent-point `(x, T)` with mass `ℓ`.
struct TangentPoint {
    p: [f64; 6],
    mass: f64,
}

impl Element<6> for TangentPoint {
    fn point(&self) -> [f64; 6] {
        self.p
    }
    fn mass(&self) -> f64 {
        self.mass
    }
}

/// Derived data of a node, recomputed by every build and refit.
#[derive(Debug, Clone, Copy)]
struct NodeData {
    /// `x̄_N`.
    xbar: [f64; 3],
    /// `T̂_N`.
    that: [f64; 3],
    /// `L_N`.
    mass: f64,
    /// `r_x²`.
    rx2: f64,
    /// `r_T`.
    rt: f64,
    /// Box of the endpoints of the node's edges (the guard).
    ends: Aabb<3>,
}

/// Per-edge geometry snapshot.
#[derive(Debug, Clone, Default)]
struct Snapshot {
    ends: Vec<[[f64; 3]; 2]>,
    length: Vec<f64>,
    tangent: Vec<[f64; 3]>,
}

/// The Barnes–Hut tree of a curve network: topology of the BVH plus the node data of one
/// geometry. Build it once per step, [`refit`](BhTree::refit) it for the trial points of the
/// line search (D26).
#[derive(Debug, Clone)]
pub struct BhTree {
    opts: BhOptions,
    edges: Vec<[u32; 2]>,
    n_vertices: usize,
    bvh: Bvh<6>,
    data: Vec<NodeData>,
    snap: Snapshot,
    lump_leaves: bool,
}

impl BhTree {
    /// Builds the tree for `positions`. Cheap geometric checks only (count, finite coordinates,
    /// zero-length edges); the full validation (intersections) is the caller's
    /// ([`CurveNetwork::validate_positions`], done by [`super::energy`]).
    ///
    /// # Errors
    /// [`RcError::InvalidOption`] for bad options (`theta_x, theta_t ≤ 0`, `leaf_size = 0`), and
    /// [`RcError::PositionCountMismatch`], [`RcError::NonFinitePosition`],
    /// [`RcError::ZeroLengthEdge`].
    pub fn build(
        network: &CurveNetwork,
        positions: &[[f64; 3]],
        opts: BhOptions,
    ) -> Result<BhTree, RcError> {
        Self::build_with(network, positions, opts, DEFAULT_MIN_SPLIT_FRACTION)
    }

    /// As [`BhTree::build`] with an explicit balance bound of the BVH splits
    /// ([`Bvh::build_with`]; `0` is the paper's rule as printed).
    ///
    /// # Errors
    /// As [`BhTree::build`]; [`RcError::InvalidOption`] for `min_split_fraction ∉ [0, ½]`.
    pub fn build_with(
        network: &CurveNetwork,
        positions: &[[f64; 3]],
        opts: BhOptions,
        min_split_fraction: f64,
    ) -> Result<BhTree, RcError> {
        opts.validate()?;
        let edges = network.edges().to_vec();
        let snap = snapshot(&edges, network.n_vertices(), positions)?;
        let elems = elements(&snap);
        let bvh = Bvh::build_with(&elems, opts.leaf_size, min_split_fraction)?;
        let mut tree = BhTree {
            opts,
            edges,
            n_vertices: network.n_vertices(),
            bvh,
            data: Vec::new(),
            snap,
            lump_leaves: false,
        };
        tree.compute_node_data();
        Ok(tree)
    }

    /// Moves the tree to `positions` keeping its topology (D26).
    ///
    /// # Errors
    /// As [`BhTree::build`] for the positions.
    pub fn refit(&mut self, positions: &[[f64; 3]]) -> Result<(), RcError> {
        self.snap = snapshot(&self.edges, self.n_vertices, positions)?;
        let elems = elements(&self.snap);
        self.bvh.refit(&elems)?;
        self.compute_node_data();
        Ok(())
    }

    /// Experimental, off by default (the contract is "leaves are always exact", §8.2): also lump
    /// leaves that are admissible. Their node data are those of any node.
    pub fn set_lump_leaves(&mut self, on: bool) {
        self.lump_leaves = on;
    }

    /// The underlying hierarchy.
    pub fn bvh(&self) -> &Bvh<6> {
        &self.bvh
    }

    /// The options the tree was built with.
    pub fn options(&self) -> &BhOptions {
        &self.opts
    }

    /// `[formula: rc-bvh-lumped]`: `L_N, x̄_N, T̂_N`, the radii and the endpoint boxes.
    fn compute_node_data(&mut self) {
        let snap = &self.snap;
        let ends = self.bvh.fold_up(
            |idx| {
                let mut b = Aabb::<3>::empty();
                for &e in idx {
                    for p in snap.ends[e as usize] {
                        b.grow(p);
                    }
                }
                b
            },
            |a, b| a.union(b),
        );
        self.data = self
            .bvh
            .nodes()
            .iter()
            .zip(ends)
            .map(|(nd, ends)| {
                let c = nd.centroid();
                let xbar = [c[0], c[1], c[2]];
                let tsum = [nd.moment[3], nd.moment[4], nd.moment[5]];
                let tn = norm(tsum);
                // Opposite tangents may cancel; r_T is then large, so the node is never used
                // unless theta_t is huge, and a zero tangent is a harmless kernel argument.
                let that = if tn > 0.0 {
                    scale(tsum, 1.0 / tn)
                } else {
                    [0.0; 3]
                };
                NodeData {
                    xbar,
                    that,
                    mass: nd.mass,
                    rx2: {
                        let r = nd.bbox.half_diagonal_of(0..3);
                        r * r
                    },
                    rt: nd.bbox.half_diagonal_of(3..6),
                    ends,
                }
            })
            .collect();
    }

    /// `[formula: rc-bh-admissibility]`
    #[inline]
    fn admissible(&self, nd: &NodeData, xi: [f64; 3], pi: &[[f64; 3]; 2]) -> bool {
        if nd.rt > self.opts.theta_t {
            return false;
        }
        let d = sub(xi, nd.xbar);
        let d2 = dot(d, d);
        if d2 <= 0.0 || nd.rx2 > self.opts.theta_x * self.opts.theta_x * d2 {
            return false;
        }
        !(nd.ends.contains_point(pi[0]) || nd.ends.contains_point(pi[1]))
    }

    /// Energy of the row of edge `i` and (optionally) the derivative of the row with respect to
    /// the two endpoints of `i`.
    fn row(&self, i: usize, e: Exponents, with_gradient: bool) -> Row {
        let s = &self.snap;
        let (pi, ti, li) = (&s.ends[i], s.tangent[i], s.length[i]);
        let xi = scale(
            [
                pi[0][0] + pi[1][0],
                pi[0][1] + pi[1][1],
                pi[0][2] + pi[1][2],
            ],
            0.5,
        );
        let mut row = Row::default();
        let nodes = self.bvh.nodes();
        let mut stack: Vec<usize> = Vec::with_capacity(64);
        stack.push(0);
        let adjacent = |j: usize| {
            let [a, b] = self.edges[i];
            let [c, d] = self.edges[j];
            a == c || a == d || b == c || b == d
        };
        while let Some(n) = stack.pop() {
            let node = &nodes[n];
            let lump = !node.is_leaf() || self.lump_leaves;
            if lump && self.admissible(&self.data[n], xi, pi) {
                let nd = &self.data[n];
                row.lumped += 1;
                let k1 = kernel(xi, nd.xbar, ti, e);
                row.energy += li * nd.mass * k1;
                if with_gradient {
                    lumped_gradient(&mut row, ti, li, xi, nd, e);
                }
                continue;
            }
            if node.is_leaf() {
                for &j in self.bvh.elements_of(n) {
                    let j = j as usize;
                    if adjacent(j) {
                        continue;
                    }
                    row.pairs += 1;
                    if with_gradient {
                        pair_with_gradient(s, i, j, e, &mut row);
                    } else {
                        row.energy += pair_energy(s, i, j, e);
                    }
                }
                continue;
            }
            // Left child first: push the right child, then the left one.
            stack.push(node.right as usize);
            stack.push(n + 1);
        }
        row
    }

    /// Diagnostic: how edge `i` is evaluated. Returns the indices of the tree nodes that are
    /// lumped and the edges `J` evaluated exactly (leaf edges not adjacent to `i`), in
    /// traversal order. The two sets cover every `J` with `(I, J) ∈ 𝒫` exactly once and never
    /// contain an edge adjacent to `i` (the endpoint-box guard, D15).
    pub fn classify(&self, i: usize) -> (Vec<usize>, Vec<usize>) {
        let s = &self.snap;
        let pi = &s.ends[i];
        let xi = scale(
            [
                pi[0][0] + pi[1][0],
                pi[0][1] + pi[1][1],
                pi[0][2] + pi[1][2],
            ],
            0.5,
        );
        let (mut lumped, mut exact) = (Vec::new(), Vec::new());
        let mut stack = vec![0usize];
        while let Some(n) = stack.pop() {
            let node = &self.bvh.nodes()[n];
            if node.is_leaf() {
                for &j in self.bvh.elements_of(n) {
                    let [a, b] = self.edges[i];
                    let [c, d] = self.edges[j as usize];
                    if !(a == c || a == d || b == c || b == d) {
                        exact.push(j as usize);
                    }
                }
            } else if self.admissible(&self.data[n], xi, pi) {
                lumped.push(n);
            } else {
                stack.push(node.right as usize);
                stack.push(n + 1);
            }
        }
        (lumped, exact)
    }

    /// Barnes–Hut energy only.
    ///
    /// # Errors
    /// [`RcError::InvalidOption`] never for a tree that was built; the signature mirrors
    /// [`BhTree::energy_and_differential`].
    pub fn energy(&self, pool: &rayon::ThreadPool, exponents: Exponents) -> (f64, BhStats) {
        let rows = self.rows(pool, exponents, false);
        let mut energy = 0.0;
        let mut stats = BhStats::default();
        for r in &rows {
            energy += r.energy;
            stats.pairs_exact += r.pairs;
            stats.nodes_lumped += r.lumped;
        }
        (energy, stats)
    }

    /// Barnes–Hut energy and differential `∂Ê/∂γ_v` (one entry per vertex).
    ///
    /// [formula: rc-bh-differential]
    pub fn energy_and_differential(
        &self,
        pool: &rayon::ThreadPool,
        exponents: Exponents,
    ) -> (f64, Vec<[f64; 3]>, BhStats) {
        let rows = self.rows(pool, exponents, true);
        let mut energy = 0.0;
        let mut stats = BhStats::default();
        let mut diff = vec![[0.0; 3]; self.n_vertices];
        for (i, r) in rows.iter().enumerate() {
            energy += r.energy;
            stats.pairs_exact += r.pairs;
            stats.nodes_lumped += r.lumped;
            let [a, b] = self.edges[i];
            diff[a as usize] = axpy(diff[a as usize], 1.0, r.grad[0]);
            diff[b as usize] = axpy(diff[b as usize], 1.0, r.grad[1]);
        }
        (energy, diff, stats)
    }

    /// Rows for all edges, indexed by edge. Computed in tree order (neighbouring traversals share
    /// cache lines) into per-edge slots; the reductions are done by the caller in edge order.
    fn rows(&self, pool: &rayon::ThreadPool, e: Exponents, with_gradient: bool) -> Vec<Row> {
        let order = self.bvh.order();
        let computed: Vec<Row> = pool.install(|| {
            order
                .par_iter()
                .with_min_len(16)
                .map(|&i| self.row(i as usize, e, with_gradient))
                .collect()
        });
        let mut rows = vec![Row::default(); order.len()];
        for (r, &i) in computed.into_iter().zip(order) {
            rows[i as usize] = r;
        }
        rows
    }
}

/// Free-function form: build the tree and evaluate once.
///
/// # Errors
/// As [`BhTree::build`].
pub fn energy_and_differential(
    pool: &rayon::ThreadPool,
    network: &CurveNetwork,
    positions: &[[f64; 3]],
    exponents: Exponents,
    opts: BhOptions,
) -> Result<(f64, Vec<[f64; 3]>, BhStats), RcError> {
    let tree = BhTree::build(network, positions, opts)?;
    Ok(tree.energy_and_differential(pool, exponents))
}

#[derive(Debug, Clone, Copy, Default)]
struct Row {
    energy: f64,
    grad: [[f64; 3]; 2],
    pairs: u64,
    lumped: u64,
}

fn snapshot(
    edges: &[[u32; 2]],
    n_vertices: usize,
    positions: &[[f64; 3]],
) -> Result<Snapshot, RcError> {
    if positions.len() != n_vertices {
        return Err(RcError::PositionCountMismatch {
            expected: n_vertices,
            found: positions.len(),
        });
    }
    for (v, p) in positions.iter().enumerate() {
        if !p.iter().all(|c| c.is_finite()) {
            return Err(RcError::NonFinitePosition { vertex: v as u32 });
        }
    }
    let mut s = Snapshot::default();
    for (i, &[a, b]) in edges.iter().enumerate() {
        let (pa, pb) = (positions[a as usize], positions[b as usize]);
        let d = sub(pb, pa);
        let l = norm(d);
        if l == 0.0 {
            return Err(RcError::ZeroLengthEdge { edge: i, a, b });
        }
        s.ends.push([pa, pb]);
        s.length.push(l);
        s.tangent.push(scale(d, 1.0 / l));
    }
    Ok(s)
}

fn elements(s: &Snapshot) -> Vec<TangentPoint> {
    (0..s.length.len())
        .map(|i| {
            let [p, q] = s.ends[i];
            let t = s.tangent[i];
            TangentPoint {
                p: [
                    0.5 * (p[0] + q[0]),
                    0.5 * (p[1] + q[1]),
                    0.5 * (p[2] + q[2]),
                    t[0],
                    t[1],
                    t[2],
                ],
                mass: s.length[i],
            }
        })
        .collect()
}

/// `e_IJ = k̂_IJ ℓ_I ℓ_J`.
///
/// [formula: rc-discrete-kernel]
#[inline]
fn pair_energy(s: &Snapshot, i: usize, j: usize, e: Exponents) -> f64 {
    let (pi, pj, ti) = (&s.ends[i], &s.ends[j], s.tangent[i]);
    let mut k = 0.0;
    for a in pi {
        for b in pj {
            k += kernel(*a, *b, ti, e);
        }
    }
    0.25 * k * s.length[i] * s.length[j]
}

/// Adds `e_IJ` and the derivative of `e_IJ + e_JI` with respect to the endpoints of `I` to `row`
/// (the same formulas as `exact::row_with_gradient`, `SPEC.md` §3.2).
///
/// [formula: rc-energy-differential]
#[inline]
fn pair_with_gradient(s: &Snapshot, i: usize, j: usize, e: Exponents, row: &mut Row) {
    let (pi, ti, li) = (&s.ends[i], s.tangent[i], s.length[i]);
    let (pj, tj, lj) = (&s.ends[j], s.tangent[j], s.length[j]);
    let signs = [-1.0, 1.0];
    let mut khat_ij = 0.0;
    let mut gp = [[0.0; 3]; 2];
    let mut gt = [0.0; 3];
    for a in 0..2 {
        for b in 0..2 {
            let ke = kernel_with_gradient(pi[a], pj[b], ti, e);
            khat_ij += ke.value;
            gp[a] = axpy(gp[a], 1.0, ke.grad_p);
            gt = axpy(gt, 1.0, ke.grad_t);
        }
    }
    khat_ij *= 0.25;
    let mut khat_ji = 0.0;
    let mut gq = [[0.0; 3]; 2];
    for a in 0..2 {
        for b in 0..2 {
            let ke = kernel_with_gradient(pj[b], pi[a], tj, e);
            khat_ji += ke.value;
            gq[a] = axpy(gq[a], -1.0, ke.grad_p);
        }
    }
    khat_ji *= 0.25;
    row.energy += khat_ij * li * lj;
    let ll = 0.25 * li * lj;
    for a in 0..2 {
        let first = axpy(scale(ti, khat_ij), 0.25, gt);
        let mut v = scale(first, signs[a] * lj);
        v = axpy(v, ll, gp[a]);
        v = axpy(v, signs[a] * lj * khat_ji, ti);
        v = axpy(v, ll, gq[a]);
        row.grad[a] = axpy(row.grad[a], 1.0, v);
    }
}

/// Derivative of `ℓ_I L_N (k(x_I, x̄_N, T_I) + k(x̄_N, x_I, T̂_N))` with respect to the endpoints
/// of `I`, `L_N, x̄_N, T̂_N` fixed.
///
/// [formula: rc-bh-differential]
#[inline]
fn lumped_gradient(
    row: &mut Row,
    ti: [f64; 3],
    li: f64,
    xi: [f64; 3],
    nd: &NodeData,
    e: Exponents,
) {
    let k1 = kernel_with_gradient(xi, nd.xbar, ti, e);
    let k2 = kernel_with_gradient(nd.xbar, xi, nd.that, e);
    let signs = [-1.0, 1.0];
    // ℓ_I-derivative and T_I-derivative (through δT_I = P(δγ₂ − δγ₁)/ℓ_I) of the first term;
    // ℓ_I-derivative of the second term (its kernel does not depend on T_I).
    let along = axpy(scale(ti, k1.value + k2.value), 1.0, k1.grad_t);
    // x_I = (γ₁ + γ₂)/2: ½ ∇_p k₁ and ½ ∇_q k₂ = −½ ∇_p k₂.
    let pos = scale(axpy(k1.grad_p, -1.0, k2.grad_p), 0.5 * li);
    for a in 0..2 {
        let v = axpy(pos, signs[a], along);
        row.grad[a] = axpy(row.grad[a], nd.mass, v);
    }
}
