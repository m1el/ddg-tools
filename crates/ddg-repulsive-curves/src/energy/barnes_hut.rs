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
//! **Differential: the discrete adjoint** (`SPEC.md` §17i) [formula: rc-bh-adjoint]. On fixed
//! interaction lists ([`BhTree::interaction_lists`]: per edge `I` the lumped nodes `𝓛(I)` and the
//! exactly evaluated edges `𝓧(I)`) the evaluated energy is
//!
//! ```text
//! Ê_bh = Σ_I [ Σ_{N ∈ 𝓛(I)} ℓ_I L_N k(x_I, x̄_N, T_I) + Σ_{J ∈ 𝓧(I)} e_IJ ],
//! L_N = Σ_{J ∈ N} ℓ_J,   x̄_N = Σ_{J ∈ N} ℓ_J x_J / L_N
//! ```
//!
//! (`T̂_N` does not enter it). [`BhTree::energy_and_differential`] returns its exact derivative:
//!
//! * the **own part** of every lumped term, with `s₁ = −1, s₂ = +1`, `k, ∇_p k, g_T = ∇_T k` at
//!   `(x_I, x̄_N, T_I)`: `∂/∂γ_{i_a} = L_N (s_a (k T_I + g_T) + ½ ℓ_I ∇_p k)`;
//! * the **node part**: each lumped term adds `b = ℓ_I k` to `∂Ê/∂L_N` (at fixed `x̄_N`) and
//!   `u = −ℓ_I ∇_p k` to `∂Ê/∂x̄_N · L_N`. Through `L_N` and `x̄_N` an edge `J ∈ N` sees
//!   `∂/∂ℓ_J = b_N + u_N·(x_J − x̄_N)` and `∂/∂x_J = ℓ_J u_N` (sums `b_N, u_N` over all `I`
//!   lumping `N`). The nodes containing `J` are the ancestors of its leaf; one downward pass
//!   accumulates `β_C = β_P + U_P·(x̄_C − x̄_P) + b_C`, `U_C = U_P + u_C` (re-centred at each node,
//!   so no large coordinates cancel) [formula: rc-bh-adjoint-shift], and at the leaf `Λ` of `J`
//!   `λ_J = β_Λ + U_Λ·(x_J − x̄_Λ)`, `μ_J = U_Λ`, `∂/∂γ_{j_b} = s_b λ_J T_J + ½ ℓ_J μ_J`;
//! * the **exact pairs**: `∂e_IJ/∂γ_I` in row `I` and `∂e_IJ/∂γ_J` in row `J` (the transposed lists;
//!   the lists are not symmetric) [formula: rc-energy-differential].
//!
//! So the differential is the derivative of exactly the function that the line search evaluates
//! with [`BhTree::energy_with_lists`] (the base curve's lists, D26); its cost is the energy's plus
//! one sequential pass over the lumped entries and one over the nodes (`O(n log n)` in all).
//! The previous differential of §8.2 (WP4: the derivative of the lumped `(I, J) + (J, I)` pair with
//! `L_N, x̄_N, T̂_N` held fixed) [formula: rc-bh-differential] remains as
//! [`BhTree::energy_and_differential_lumped_pairs`]; it is not the derivative of `Ê_bh` (it differs
//! from finite differences by the approximation error, `SPEC.md` §17d, §17f).
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
        // The exact pair energies are added after the lumped terms, so that the sum is the one
        // of `energy_with_lists` on this tree's own lists, bit for bit.
        let mut exact_energies: Vec<f64> = Vec::new();
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
                    exact_energies.push(if with_gradient {
                        pair_with_gradient(s, i, j, e, &mut row)
                    } else {
                        pair_energy(s, i, j, e)
                    });
                }
                continue;
            }
            // Left child first: push the right child, then the left one.
            stack.push(node.right as usize);
            stack.push(n + 1);
        }
        for v in exact_energies {
            row.energy += v;
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
            // The same decisions, in the same order, as `row` (leaves only lumped when
            // `set_lump_leaves` is on).
            let lump = !node.is_leaf() || self.lump_leaves;
            if lump && self.admissible(&self.data[n], xi, pi) {
                lumped.push(n);
            } else if node.is_leaf() {
                for &j in self.bvh.elements_of(n) {
                    let [a, b] = self.edges[i];
                    let [c, d] = self.edges[j as usize];
                    if !(a == c || a == d || b == c || b == d) {
                        exact.push(j as usize);
                    }
                }
            } else {
                stack.push(node.right as usize);
                stack.push(n + 1);
            }
        }
        (lumped, exact)
    }

    /// The interaction lists of the current tree and positions: for every edge `I`, the lumped
    /// nodes and the exactly evaluated edges (what [`BhTree::classify`] returns).
    ///
    /// Barnes–Hut energies at different positions are not comparable at the level of a line
    /// search's decrease: the admissibility tests flip as the curve moves, and each flip changes the
    /// energy by about the approximation error (1e-6 relative on Freedman 600), which swamps
    /// the decrease near a minimum. Evaluating the trial energies with the lists of the base
    /// curve ([`BhTree::energy_with_lists`], after a [`BhTree::refit`]) removes that noise: the
    /// energy is then a smooth function of the positions.
    pub fn interaction_lists(&self, pool: &rayon::ThreadPool) -> InteractionLists {
        let (lumped, exact): (Vec<_>, Vec<_>) = pool.install(|| {
            (0..self.edges.len())
                .into_par_iter()
                .with_min_len(16)
                .map(|i| {
                    let (l, x) = self.classify(i);
                    (
                        l.into_iter().map(|n| n as u32).collect::<Vec<_>>(),
                        x.into_iter().map(|j| j as u32).collect::<Vec<_>>(),
                    )
                })
                .unzip()
        });
        InteractionLists {
            lumped,
            exact,
            n_nodes: self.bvh.nodes().len(),
        }
    }

    /// The Barnes–Hut energy with the given (frozen) interaction lists at the current positions
    /// (set by [`BhTree::refit`]); summed in edge order, so independent of the thread count.
    ///
    /// # Panics
    /// If `lists` was made for a tree with a different number of edges.
    pub fn energy_with_lists(
        &self,
        lists: &InteractionLists,
        pool: &rayon::ThreadPool,
        e: Exponents,
    ) -> f64 {
        assert_eq!(lists.lumped.len(), self.edges.len());
        let s = &self.snap;
        let rows: Vec<f64> = pool.install(|| {
            (0..self.edges.len())
                .into_par_iter()
                .with_min_len(16)
                .map(|i| {
                    let (pi, ti, li) = (&s.ends[i], s.tangent[i], s.length[i]);
                    let xi = scale(
                        [
                            pi[0][0] + pi[1][0],
                            pi[0][1] + pi[1][1],
                            pi[0][2] + pi[1][2],
                        ],
                        0.5,
                    );
                    let mut sum = 0.0;
                    for &n in &lists.lumped[i] {
                        let nd = &self.data[n as usize];
                        sum += li * nd.mass * kernel(xi, nd.xbar, ti, e);
                    }
                    for &j in &lists.exact[i] {
                        sum += pair_energy(s, i, j as usize, e);
                    }
                    sum
                })
                .collect()
        });
        rows.iter().sum()
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

    /// Barnes–Hut energy and its exact differential `∂Ê_bh/∂γ_v` (one entry per vertex): the
    /// discrete adjoint on the interaction lists of the current tree and positions (module docs,
    /// `SPEC.md` §17i). The energy is [`BhTree::energy_with_lists`] on those lists, bit for bit,
    /// and the differential is its derivative with the lists held fixed, so a line search that
    /// evaluates its trials with these lists sees a descent direction whose directional
    /// derivative is exact. Computes the lists itself; a caller that has them already uses
    /// [`BhTree::energy_and_differential_with_lists`].
    ///
    /// [formula: rc-bh-adjoint]
    pub fn energy_and_differential(
        &self,
        pool: &rayon::ThreadPool,
        exponents: Exponents,
    ) -> (f64, Vec<[f64; 3]>, BhStats) {
        let lists = self.interaction_lists(pool);
        self.adjoint(&lists, pool, exponents)
    }

    /// As [`BhTree::energy_and_differential`] with given interaction lists (normally the ones
    /// [`BhTree::interaction_lists`] made for this tree at the current positions; any lists of
    /// this tree's topology give the energy on those lists and its exact derivative).
    ///
    /// # Errors
    /// [`RcError::InvalidOption`] (`lists`) if the lists were made for a tree with a different
    /// number of edges or nodes.
    ///
    /// [formula: rc-bh-adjoint]
    pub fn energy_and_differential_with_lists(
        &self,
        lists: &InteractionLists,
        pool: &rayon::ThreadPool,
        exponents: Exponents,
    ) -> Result<(f64, Vec<[f64; 3]>, BhStats), RcError> {
        let (m, nn) = (self.edges.len(), self.bvh.nodes().len());
        if lists.lumped.len() != m || lists.exact.len() != m || lists.n_nodes != nn {
            return Err(RcError::InvalidOption {
                name: "lists",
                value: format!(
                    "interaction lists of a tree with {} edges and {} nodes",
                    lists.lumped.len(),
                    lists.n_nodes
                ),
                requirement: "the lists of this tree (BhTree::interaction_lists): same edge and \
                              node count",
            });
        }
        Ok(self.adjoint(lists, pool, exponents))
    }

    /// The discrete adjoint of [`BhTree::energy_with_lists`] (module docs, `SPEC.md` §17i).
    ///
    /// 1. rows in parallel (tree order, per-edge slots): the energy, the own part of every term
    ///    of the row, the node terms `(b, u)` of its lumped entries, and `∂e_KI/∂γ_I` for the
    ///    edges `K` whose exact list contains `I`;
    /// 2. node sums `b_N, u_N` in increasing edge order (sequential);
    /// 3. the downward pass over the nodes in index order (parents precede their children);
    /// 4. per edge, the node part at its leaf; vertex sums in increasing edge order.
    ///
    /// Every reduction has a fixed order: bitwise independent of the thread count (D25).
    ///
    /// [formula: rc-bh-adjoint]
    fn adjoint(
        &self,
        lists: &InteractionLists,
        pool: &rayon::ThreadPool,
        e: Exponents,
    ) -> (f64, Vec<[f64; 3]>, BhStats) {
        let m = self.edges.len();
        let s = &self.snap;
        let nodes = self.bvh.nodes();
        let nn = nodes.len();
        let (toff, tidx) = transpose(&lists.exact, m);
        let order = self.bvh.order();
        let computed: Vec<AdjRow> = pool.install(|| {
            order
                .par_iter()
                .with_min_len(16)
                .map(|&i| {
                    let i = i as usize;
                    self.adjoint_row(i, lists, &tidx[toff[i]..toff[i + 1]], e)
                })
                .collect()
        });
        let mut rows: Vec<AdjRow> = (0..m).map(|_| AdjRow::default()).collect();
        for (r, &i) in computed.into_iter().zip(order) {
            rows[i as usize] = r;
        }
        // 2. Node sums: b_N = ∂Ê/∂L_N at fixed x̄_N, u_N = ∂Ê/∂x̄_N / L_N.
        let mut b = vec![0.0; nn];
        let mut u = vec![[0.0; 3]; nn];
        for (i, r) in rows.iter().enumerate() {
            for (&n, &(bk, uk)) in lists.lumped[i].iter().zip(&r.node_terms) {
                b[n as usize] += bk;
                u[n as usize] = axpy(u[n as usize], 1.0, uk);
            }
        }
        // 3. Downward pass [formula: rc-bh-adjoint-shift]: Σ_{A ⊇ N} (b_A + u_A·(x − x̄_A)) =
        //    beta_N + ubig_N·(x − x̄_N).
        let mut beta = vec![0.0; nn];
        let mut ubig = vec![[0.0; 3]; nn];
        if nn > 0 {
            beta[0] = b[0];
            ubig[0] = u[0];
        }
        for (n, node) in nodes.iter().enumerate() {
            if node.is_leaf() {
                continue;
            }
            for c in [n + 1, node.right as usize] {
                let shift = sub(self.data[c].xbar, self.data[n].xbar);
                beta[c] = beta[n] + dot(ubig[n], shift) + b[c];
                ubig[c] = axpy(ubig[n], 1.0, u[c]);
            }
        }
        // 4. The node part of every edge, at its leaf: s_b λ_J T_J + ½ ℓ_J μ_J.
        let mut through_nodes = vec![[[0.0; 3]; 2]; m];
        for (n, node) in nodes.iter().enumerate() {
            if !node.is_leaf() {
                continue;
            }
            for &j in self.bvh.elements_of(n) {
                let j = j as usize;
                let [p, q] = s.ends[j];
                let xj = scale([p[0] + q[0], p[1] + q[1], p[2] + q[2]], 0.5);
                let lambda = beta[n] + dot(ubig[n], sub(xj, self.data[n].xbar));
                let half = scale(ubig[n], 0.5 * s.length[j]);
                through_nodes[j] = [
                    axpy(half, -lambda, s.tangent[j]),
                    axpy(half, lambda, s.tangent[j]),
                ];
            }
        }
        let energy: f64 = rows.iter().map(|r| r.energy).sum();
        let mut stats = BhStats::default();
        let mut diff = vec![[0.0; 3]; self.n_vertices];
        for (i, r) in rows.iter().enumerate() {
            stats.pairs_exact += lists.exact[i].len() as u64;
            stats.nodes_lumped += lists.lumped[i].len() as u64;
            let [a, bb] = self.edges[i];
            let g0 = axpy(r.grad[0], 1.0, through_nodes[i][0]);
            let g1 = axpy(r.grad[1], 1.0, through_nodes[i][1]);
            diff[a as usize] = axpy(diff[a as usize], 1.0, g0);
            diff[bb as usize] = axpy(diff[bb as usize], 1.0, g1);
        }
        (energy, diff, stats)
    }

    /// Row `i` of the adjoint: the energy of its lists (the same operations, in the same order,
    /// as [`BhTree::energy_with_lists`]), the own derivative of its lumped terms and exact pairs,
    /// the partner derivatives `∂e_KI/∂γ_I` for `K ∈ partners`, and the node terms.
    ///
    /// [formula: rc-bh-adjoint]
    fn adjoint_row(
        &self,
        i: usize,
        lists: &InteractionLists,
        partners: &[u32],
        e: Exponents,
    ) -> AdjRow {
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
        let mut row = AdjRow {
            node_terms: Vec::with_capacity(lists.lumped[i].len()),
            ..AdjRow::default()
        };
        let signs = [-1.0, 1.0];
        for &n in &lists.lumped[i] {
            let nd = &self.data[n as usize];
            let k = kernel_with_gradient(xi, nd.xbar, ti, e);
            row.energy += li * nd.mass * k.value;
            // Own part: L_N (s_a (k T_I + g_T) + ½ ℓ_I ∇_p k).
            let along = axpy(scale(ti, k.value), 1.0, k.grad_t);
            let pos = scale(k.grad_p, 0.5 * li);
            for a in 0..2 {
                row.grad[a] = axpy(row.grad[a], nd.mass, axpy(pos, signs[a], along));
            }
            // Node part: b = ℓ_I k, u = ℓ_I ∇_q k = −ℓ_I ∇_p k.
            row.node_terms.push((li * k.value, scale(k.grad_p, -li)));
        }
        for &j in &lists.exact[i] {
            pair_own_gradient(s, i, j as usize, e, &mut row);
        }
        for &k in partners {
            pair_partner_gradient(s, k as usize, i, e, &mut row.grad);
        }
        row
    }

    /// The differential of `SPEC.md` §8.2 as implemented by WP4: the derivative of the lumped
    /// `(I, J) + (J, I)` pair with `L_N, x̄_N, T̂_N` held fixed. An approximation of the exact
    /// differential, but **not** the derivative of the evaluated Barnes–Hut energy (it differs
    /// from finite differences by about the approximation error, §17d, §17f); kept for
    /// comparison. [`BhTree::energy_and_differential`] is the exact derivative.
    ///
    /// [formula: rc-bh-differential]
    pub fn energy_and_differential_lumped_pairs(
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

/// Per-edge interaction lists of a [`BhTree`] ([`BhTree::interaction_lists`]).
#[derive(Debug, Clone, Default)]
pub struct InteractionLists {
    lumped: Vec<Vec<u32>>,
    exact: Vec<Vec<u32>>,
    /// Node count of the tree the lists were made for.
    n_nodes: usize,
}

impl InteractionLists {
    /// Number of lumped (edge, node) entries.
    pub fn n_lumped(&self) -> usize {
        self.lumped.iter().map(Vec::len).sum()
    }

    /// Number of exactly evaluated ordered edge pairs.
    pub fn n_exact(&self) -> usize {
        self.exact.iter().map(Vec::len).sum()
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

/// One row of the adjoint.
#[derive(Debug, Clone, Default)]
struct AdjRow {
    energy: f64,
    /// Derivative with respect to the two endpoints of the row's edge (own and partner parts).
    grad: [[f64; 3]; 2],
    /// `(b, u)` per lumped entry, in the order of the row's lumped list.
    node_terms: Vec<(f64, [f64; 3])>,
}

/// The transpose of per-edge lists: `idx[off[j]..off[j + 1]]` are the `i` with `j ∈ lists[i]`, in
/// increasing order.
fn transpose(lists: &[Vec<u32>], m: usize) -> (Vec<usize>, Vec<u32>) {
    let mut off = vec![0usize; m + 1];
    for l in lists {
        for &j in l {
            off[j as usize + 1] += 1;
        }
    }
    for k in 0..m {
        off[k + 1] += off[k];
    }
    let mut fill = off.clone();
    let mut idx = vec![0u32; off[m]];
    for (i, l) in lists.iter().enumerate() {
        for &j in l {
            let f = &mut fill[j as usize];
            idx[*f] = i as u32;
            *f += 1;
        }
    }
    (off, idx)
}

/// Adds `e_IJ` (the operations of [`pair_energy`]) and `∂e_IJ/∂γ_{i_a}` to `row`:
/// `s_a ℓ_J (k̂_IJ T_I + g_T) + ℓ_I ℓ_J ¼ Σ_b ∇_p k(γ_{i_a}, γ_{j_b}, T_I)`.
///
/// [formula: rc-energy-differential]
#[inline]
fn pair_own_gradient(s: &Snapshot, i: usize, j: usize, e: Exponents, row: &mut AdjRow) {
    let (pi, ti, li) = (&s.ends[i], s.tangent[i], s.length[i]);
    let (pj, lj) = (&s.ends[j], s.length[j]);
    let mut k = 0.0;
    let mut gp = [[0.0; 3]; 2];
    let mut gt = [0.0; 3];
    for a in 0..2 {
        for b in 0..2 {
            let ke = kernel_with_gradient(pi[a], pj[b], ti, e);
            k += ke.value;
            gp[a] = axpy(gp[a], 1.0, ke.grad_p);
            gt = axpy(gt, 1.0, ke.grad_t);
        }
    }
    row.energy += 0.25 * k * li * lj;
    let khat = 0.25 * k;
    let first = axpy(scale(ti, khat), 0.25, gt);
    let ll = 0.25 * li * lj;
    let signs = [-1.0, 1.0];
    for a in 0..2 {
        let v = axpy(scale(first, signs[a] * lj), ll, gp[a]);
        row.grad[a] = axpy(row.grad[a], 1.0, v);
    }
}

/// Adds `∂e_KI/∂γ_{i_b} = s_b ℓ_K k̂_KI T_I − ℓ_K ℓ_I ¼ Σ_a ∇_p k(γ_{k_a}, γ_{i_b}, T_K)` to `grad`
/// (the derivative of a pair term of row `K` with respect to the endpoints of its partner `I`).
///
/// [formula: rc-energy-differential]
#[inline]
fn pair_partner_gradient(s: &Snapshot, k: usize, i: usize, e: Exponents, grad: &mut [[f64; 3]; 2]) {
    let (pk, tk, lk) = (&s.ends[k], s.tangent[k], s.length[k]);
    let (pi, ti, li) = (&s.ends[i], s.tangent[i], s.length[i]);
    let mut ksum = 0.0;
    let mut gq = [[0.0; 3]; 2];
    for a in 0..2 {
        for b in 0..2 {
            let ke = kernel_with_gradient(pk[a], pi[b], tk, e);
            ksum += ke.value;
            gq[b] = axpy(gq[b], -1.0, ke.grad_p);
        }
    }
    let khat = 0.25 * ksum;
    let ll = 0.25 * lk * li;
    let signs = [-1.0, 1.0];
    for b in 0..2 {
        let v = axpy(scale(ti, signs[b] * lk * khat), ll, gq[b]);
        grad[b] = axpy(grad[b], 1.0, v);
    }
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

/// Adds the derivative of `e_IJ + e_JI` with respect to the endpoints of `I` to `row` and returns
/// `e_IJ` (the same formulas as `exact::row_with_gradient`, `SPEC.md` §3.2).
///
/// [formula: rc-energy-differential]
#[inline]
fn pair_with_gradient(s: &Snapshot, i: usize, j: usize, e: Exponents, row: &mut Row) -> f64 {
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
    let e_ij = khat_ij * li * lj;
    let ll = 0.25 * li * lj;
    for a in 0..2 {
        let first = axpy(scale(ti, khat_ij), 0.25, gt);
        let mut v = scale(first, signs[a] * lj);
        v = axpy(v, ll, gp[a]);
        v = axpy(v, signs[a] * lj * khat_ji, ti);
        v = axpy(v, ll, gq[a]);
        row.grad[a] = axpy(row.grad[a], 1.0, v);
    }
    e_ij
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
