//! Block cluster tree and its 0th-order fast-multipole matrix–vector product (`SPEC.md` §8.3,
//! App. B.2, Eq. 28; D16, D25; §16).
//!
//! The tree is generic: it is built on a [`Bvh`] of any dimension and evaluates a family of `N`
//! symmetric kernel matrices `K_k ∈ R^{n×n}` (`k < N`) given by a [`ClusterKernel`]. There is no
//! curve type in any signature (`SPEC.md` §16), so the surface papers can reuse it; the curve
//! kernels (`K^B`, `K⁰` of Eq. 31–32) live in [`crate::metric::hierarchical`].
//!
//! # Blocks
//!
//! Blocks are pairs `(𝒜, ℬ)` of BVH nodes, starting from `(root, root)`. A block is
//!
//! * an **admissible leaf** if `𝒜 ≠ ℬ` and [`ClusterKernel::admissible`] holds (for curves:
//!   Eq. 28 plus the endpoint-box guard, D16) [formula: rc-bct-admissibility];
//! * otherwise an **inadmissible leaf** if `|𝒜|·|ℬ| ≤ leaf_pairs` or both nodes are BVH leaves;
//! * otherwise it is split into the pairs of children (a BVH leaf is kept as it is, so only the
//!   non-leaf side splits).
//!
//! Starting from `(root, root)`, the two nodes of a block are either equal or disjoint, the blocks
//! partition the `n × n` index set, and the partition is symmetric (`(𝒜, ℬ)` is a block of a kind
//! iff `(ℬ, 𝒜)` is) when `admissible` is symmetric.
//!
//! # Product
//!
//! An admissible block uses the rank-1 form `K̂_𝒜ℬ = w[𝒜] κ(𝒜, ℬ) w[ℬ]ᵀ` with the element
//! weights `w` ([`ClusterKernel::weight`]; for curves `ℓ_I`) and one coefficient per kernel
//! ([`ClusterKernel::far`]) [formula: rc-bct-lowrank]; an inadmissible block uses the exact
//! entries ([`ClusterKernel::near`]). `φ = K ψ` for `C` channels, channel `c` multiplied by kernel
//! `map[c]`, is computed as a 0th-order fast multipole pass:
//!
//! 1. upward: `U_ℬ = Σ_{J ∈ ℬ} w_J ψ_J` per node (children summed);
//! 2. far field: `F_𝒜 = Σ_{(𝒜, ℬ) admissible} κ(𝒜, ℬ) U_ℬ` per node;
//! 3. downward: `G_𝒜 = Σ_{ancestors 𝒩 ⊇ 𝒜} F_𝒩`, and `φ_I = w_I G_{leaf(I)}`;
//! 4. near field: `φ_I += Σ_J K_IJ ψ_J` over the inadmissible blocks, stored as a sparse row
//!    matrix computed once at build time.
//!
//! The admissible coefficients and the near-field entries are computed once per tree, so a
//! product costs no kernel evaluation. The structure (cluster pairs split level by level, the
//! upward/downward passes of `SetVIs`/`PropagateBIs`, `K𝟙` computed once per tree) follows the
//! authors' MIT code (`src/product/block_cluster_tree.cpp`, `include/product/block_cluster_tree.h`,
//! Copyright (c) 2019 Christopher Yu); the admissibility differs (D16).
//!
//! # Determinism (D25)
//!
//! The build classifies the pairs of one level in parallel and collects them in order; the
//! product runs every reduction in a fixed order (per node for the far field, per element row for
//! the near field and the downward pass). Results are bitwise independent of the thread count.

use rayon::prelude::*;

use crate::error::RcError;
use crate::spatial::bvh::Bvh;

/// Parameters of the block cluster tree of the curve metric (`SPEC.md` §8.3, §9.1, §10).
///
/// The generic tree only uses `leaf_pairs`; `theta_x`, `theta_t` are the admissibility
/// parameters of Eq. 28 used by the curve kernel, `leaf_size` the leaf size of the BVH it is
/// built on.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct BctOptions {
    /// Spatial parameter: `max(r_x^𝒜, r_x^ℬ) ≤ θ_x |x̄_𝒜 − x̄_ℬ|`. `0` makes every block
    /// inadmissible (the exact product).
    pub theta_x: f64,
    /// Tangent parameter: `max(r_T^𝒜, r_T^ℬ) ≤ θ_T`. `0` makes every block inadmissible.
    pub theta_t: f64,
    /// A block with `|𝒜|·|ℬ| ≤ leaf_pairs` is not split further.
    pub leaf_pairs: usize,
    /// Maximal number of elements in a BVH leaf.
    pub leaf_size: usize,
}

impl Default for BctOptions {
    /// `theta_x 0.25, theta_t 0.25, leaf_pairs 64, leaf_size 8` (`SPEC.md` §10, §17e: 0.5
    /// missed the 1e-2 operator error target of §13.4).
    fn default() -> Self {
        BctOptions {
            theta_x: 0.25,
            theta_t: 0.25,
            leaf_pairs: 64,
            leaf_size: 8,
        }
    }
}

impl BctOptions {
    /// Checks the options.
    ///
    /// # Errors
    /// [`RcError::InvalidOption`] for a negative or non-finite `theta_x` or `theta_t`, and for
    /// `leaf_size = 0`.
    pub fn validate(&self) -> Result<(), RcError> {
        for (name, v) in [("theta_x", self.theta_x), ("theta_t", self.theta_t)] {
            if !(v.is_finite() && v >= 0.0) {
                return Err(RcError::InvalidOption {
                    name,
                    value: v.to_string(),
                    requirement: "a finite number ≥ 0 (0: exact product)",
                });
            }
        }
        if self.leaf_size == 0 {
            return Err(RcError::InvalidOption {
                name: "leaf_size",
                value: "0".into(),
                requirement: "a leaf holds at least 1 element",
            });
        }
        Ok(())
    }
}

/// A family of `N` symmetric kernel matrices over the elements of a [`Bvh`], with the cluster
/// data needed by the block cluster tree. Node indices are those of [`Bvh::nodes`], element
/// indices those of the elements the BVH was built from.
pub trait ClusterKernel<const N: usize>: Sync {
    /// Whether the block `(a, b)` of two distinct, disjoint nodes may use the rank-1 form. Must
    /// be symmetric in `(a, b)` for the hierarchical matrix to be symmetric.
    ///
    /// [formula: rc-bct-admissibility]
    fn admissible(&self, a: usize, b: usize) -> bool;
    /// The coefficients `κ_k(a, b)` of the rank-1 blocks `w[a] κ_k w[b]ᵀ`. Must be bitwise
    /// symmetric, `far(a, b) == far(b, a)`, for the hierarchical matrix to be symmetric.
    ///
    /// [formula: rc-bct-lowrank]
    fn far(&self, a: usize, b: usize) -> [f64; N];
    /// The exact entries `K_k(i, j)` for two distinct elements, `None` for a structurally zero
    /// entry (for curves: adjacent edges). Must be bitwise symmetric.
    fn near(&self, i: usize, j: usize) -> Option<[f64; N]>;
    /// The weight `w_i` of the rank-1 factors (for curves: `ℓ_I`).
    fn weight(&self, i: usize) -> f64;
}

/// Counts of a block cluster tree.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub struct BctStats {
    /// Number of admissible (rank-1) blocks.
    pub admissible_blocks: usize,
    /// Number of inadmissible (exact) blocks.
    pub inadmissible_blocks: usize,
    /// Element pairs covered by admissible blocks.
    pub admissible_pairs: u64,
    /// Stored exact entries (pairs of inadmissible blocks that are not structurally zero).
    pub near_entries: usize,
    /// Number of BVH nodes.
    pub nodes: usize,
    /// Levels of the block tree.
    pub depth: usize,
}

/// A block cluster tree with the precomputed far-field coefficients and near-field entries of
/// `N` kernels (`SPEC.md` §8.3). See the module documentation.
#[derive(Debug, Clone)]
pub struct BlockClusterTree<const N: usize> {
    /// Per node: `start, end` in `order`, and the right child (0 for a leaf).
    start: Vec<u32>,
    end: Vec<u32>,
    right: Vec<u32>,
    order: Vec<u32>,
    /// The leaf containing each element.
    leaf_of: Vec<u32>,
    weight: Vec<f64>,
    /// Admissible blocks grouped by row node: `far_ptr[a]..far_ptr[a + 1]` index `far_col`,
    /// `far_coef`.
    far_ptr: Vec<usize>,
    far_col: Vec<u32>,
    far_coef: Vec<[f64; N]>,
    /// Inadmissible blocks `(a, b)` in creation order.
    near_blocks: Vec<(u32, u32)>,
    /// Near field as sparse rows over the elements.
    near_ptr: Vec<usize>,
    near_col: Vec<u32>,
    near_val: Vec<[f64; N]>,
    /// `K_k 𝟙` of the hierarchical matrix.
    row_sums: Vec<[f64; N]>,
    stats: BctStats,
}

enum Class {
    Admissible,
    Inadmissible,
    Split,
}

impl<const N: usize> BlockClusterTree<N> {
    /// Builds the block tree of `bvh` for `kernel`, precomputes the far-field coefficients, the
    /// near-field entries and `K𝟙`, on `pool`. Bitwise independent of the pool size (D25).
    ///
    /// # Errors
    /// [`RcError::InvalidOption`] if the BVH is empty.
    pub fn build<const D: usize, K: ClusterKernel<N>>(
        pool: &rayon::ThreadPool,
        bvh: &Bvh<D>,
        kernel: &K,
        leaf_pairs: usize,
    ) -> Result<Self, RcError> {
        let nodes = bvh.nodes();
        if nodes.is_empty() || bvh.n_elements() == 0 {
            return Err(RcError::InvalidOption {
                name: "bvh",
                value: "empty".into(),
                requirement: "at least one element",
            });
        }
        let n_nodes = nodes.len();
        let n = bvh.n_elements();
        let start: Vec<u32> = nodes.iter().map(|nd| nd.start).collect();
        let end: Vec<u32> = nodes.iter().map(|nd| nd.end).collect();
        let right: Vec<u32> = nodes.iter().map(|nd| nd.right).collect();
        let order = bvh.order().to_vec();
        let mut leaf_of = vec![0u32; n];
        for (i, nd) in nodes.iter().enumerate() {
            if nd.is_leaf() {
                for &e in bvh.elements_of(i) {
                    leaf_of[e as usize] = i as u32;
                }
            }
        }
        let weight: Vec<f64> = (0..n).map(|i| kernel.weight(i)).collect();

        // Level-by-level splitting; each level is classified in parallel and collected in order.
        let is_leaf = |a: usize| right[a] == 0;
        let len = |a: usize| (end[a] - start[a]) as usize;
        let classify = |&(a, b): &(u32, u32)| {
            let (a, b) = (a as usize, b as usize);
            if a != b && kernel.admissible(a, b) {
                Class::Admissible
            } else if (is_leaf(a) && is_leaf(b)) || len(a).saturating_mul(len(b)) <= leaf_pairs {
                Class::Inadmissible
            } else {
                Class::Split
            }
        };
        let mut adm: Vec<(u32, u32)> = Vec::new();
        let mut near_blocks: Vec<(u32, u32)> = Vec::new();
        let mut level: Vec<(u32, u32)> = vec![(0, 0)];
        let mut depth = 0;
        while !level.is_empty() {
            depth += 1;
            let classes: Vec<Class> =
                pool.install(|| level.par_iter().with_min_len(64).map(classify).collect());
            let mut next = Vec::with_capacity(4 * level.len());
            for (&(a, b), c) in level.iter().zip(classes) {
                match c {
                    Class::Admissible => adm.push((a, b)),
                    Class::Inadmissible => near_blocks.push((a, b)),
                    Class::Split => {
                        let kids = |x: u32| -> ([u32; 2], usize) {
                            if is_leaf(x as usize) {
                                ([x, x], 1)
                            } else {
                                ([x + 1, right[x as usize]], 2)
                            }
                        };
                        let (ka, na) = kids(a);
                        let (kb, nb) = kids(b);
                        for &ca in &ka[..na] {
                            for &cb in &kb[..nb] {
                                next.push((ca, cb));
                            }
                        }
                    }
                }
            }
            level = next;
        }

        // Far field grouped by row node (stable: creation order within a node).
        let mut far_ptr = vec![0usize; n_nodes + 1];
        for &(a, _) in &adm {
            far_ptr[a as usize + 1] += 1;
        }
        for i in 0..n_nodes {
            far_ptr[i + 1] += far_ptr[i];
        }
        let mut fill = far_ptr.clone();
        let mut far_col = vec![0u32; adm.len()];
        for &(a, b) in &adm {
            far_col[fill[a as usize]] = b;
            fill[a as usize] += 1;
        }
        let far_coef: Vec<[f64; N]> = pool.install(|| {
            (0..n_nodes)
                .into_par_iter()
                .flat_map_iter(|a| {
                    far_col[far_ptr[a]..far_ptr[a + 1]]
                        .iter()
                        .map(move |&b| kernel.far(a, b as usize))
                })
                .collect()
        });
        let admissible_pairs = adm
            .iter()
            .map(|&(a, b)| (len(a as usize) * len(b as usize)) as u64)
            .sum();

        // Near field: inadmissible blocks grouped by row node, then rows per element over the
        // chain root → leaf(I).
        let mut parent = vec![u32::MAX; n_nodes];
        for i in 0..n_nodes {
            if !is_leaf(i) {
                parent[i + 1] = i as u32;
                parent[right[i] as usize] = i as u32;
            }
        }
        let mut nb_ptr = vec![0usize; n_nodes + 1];
        for &(a, _) in &near_blocks {
            nb_ptr[a as usize + 1] += 1;
        }
        for i in 0..n_nodes {
            nb_ptr[i + 1] += nb_ptr[i];
        }
        let mut fill = nb_ptr.clone();
        let mut nb_col = vec![0u32; near_blocks.len()];
        for &(a, b) in &near_blocks {
            nb_col[fill[a as usize]] = b;
            fill[a as usize] += 1;
        }
        let rows: Vec<Vec<(u32, [f64; N])>> = pool.install(|| {
            (0..n)
                .into_par_iter()
                .with_min_len(16)
                .map(|i| {
                    let mut chain = Vec::with_capacity(32);
                    let mut a = leaf_of[i];
                    loop {
                        chain.push(a);
                        if a == 0 {
                            break;
                        }
                        a = parent[a as usize];
                    }
                    let mut row = Vec::new();
                    for &a in chain.iter().rev() {
                        for &b in &nb_col[nb_ptr[a as usize]..nb_ptr[a as usize + 1]] {
                            let b = b as usize;
                            for &j in &order[start[b] as usize..end[b] as usize] {
                                if j as usize == i {
                                    continue;
                                }
                                if let Some(k) = kernel.near(i, j as usize) {
                                    row.push((j, k));
                                }
                            }
                        }
                    }
                    row
                })
                .collect()
        });
        let mut near_ptr = Vec::with_capacity(n + 1);
        near_ptr.push(0);
        let total: usize = rows.iter().map(Vec::len).sum();
        let mut near_col = Vec::with_capacity(total);
        let mut near_val = Vec::with_capacity(total);
        for r in rows {
            for (j, k) in r {
                near_col.push(j);
                near_val.push(k);
            }
            near_ptr.push(near_col.len());
        }

        let stats = BctStats {
            admissible_blocks: adm.len(),
            inadmissible_blocks: near_blocks.len(),
            admissible_pairs,
            near_entries: near_col.len(),
            nodes: n_nodes,
            depth,
        };
        let mut tree = BlockClusterTree {
            start,
            end,
            right,
            order,
            leaf_of,
            weight,
            far_ptr,
            far_col,
            far_coef,
            near_blocks,
            near_ptr,
            near_col,
            near_val,
            row_sums: Vec::new(),
            stats,
        };
        let mut map = [0usize; N];
        for (k, m) in map.iter_mut().enumerate() {
            *m = k;
        }
        tree.row_sums = tree.multiply(pool, &map, &vec![[1.0; N]; n])?;
        Ok(tree)
    }

    /// Number of elements `n`.
    pub fn n_elements(&self) -> usize {
        self.order.len()
    }

    /// Counts of the tree.
    pub fn stats(&self) -> BctStats {
        self.stats
    }

    /// `K_k 𝟙` of the hierarchical matrix, per element (computed once by the build).
    pub fn row_sums(&self) -> &[[f64; N]] {
        &self.row_sums
    }

    /// The elements of node `a`.
    pub fn elements_of(&self, a: usize) -> &[u32] {
        &self.order[self.start[a] as usize..self.end[a] as usize]
    }

    /// The admissible blocks `(a, b, κ(a, b))`, grouped by `a`.
    pub fn admissible_blocks(&self) -> impl Iterator<Item = (usize, usize, [f64; N])> + '_ {
        (0..self.far_ptr.len() - 1).flat_map(move |a| {
            (self.far_ptr[a]..self.far_ptr[a + 1])
                .map(move |p| (a, self.far_col[p] as usize, self.far_coef[p]))
        })
    }

    /// The inadmissible blocks `(a, b)` in creation order.
    pub fn inadmissible_blocks(&self) -> impl Iterator<Item = (usize, usize)> + '_ {
        self.near_blocks
            .iter()
            .map(|&(a, b)| (a as usize, b as usize))
    }

    /// `φ = K ψ` for `C` channels: channel `c` of `φ` is `K_{map[c]}` applied to channel `c` of
    /// `ψ` (one `C`-vector per element). Parallel on `pool`, bitwise independent of its size.
    ///
    /// [formula: rc-bct-lowrank]
    ///
    /// # Errors
    /// [`RcError::InvalidOption`] if `psi` has the wrong length or `map` names a kernel `≥ N`.
    pub fn multiply<const C: usize>(
        &self,
        pool: &rayon::ThreadPool,
        map: &[usize; C],
        psi: &[[f64; C]],
    ) -> Result<Vec<[f64; C]>, RcError> {
        let n = self.n_elements();
        if psi.len() != n {
            return Err(RcError::InvalidOption {
                name: "block cluster tree argument",
                value: psi.len().to_string(),
                requirement: "one value per element",
            });
        }
        if let Some(&m) = map.iter().find(|&&m| m >= N) {
            return Err(RcError::InvalidOption {
                name: "channel map",
                value: m.to_string(),
                requirement: "kernel indices below the number of kernels",
            });
        }
        Ok(self.multiply_unchecked(pool, map, psi))
    }

    /// [`BlockClusterTree::multiply`] without the argument checks.
    pub(crate) fn multiply_unchecked<const C: usize>(
        &self,
        pool: &rayon::ThreadPool,
        map: &[usize; C],
        psi: &[[f64; C]],
    ) -> Vec<[f64; C]> {
        let n = self.n_elements();
        let n_nodes = self.start.len();
        // 1. Upward sums U_B = Σ w_J ψ_J (sequential, bottom up).
        let mut up = vec![[0.0f64; C]; n_nodes];
        for a in (0..n_nodes).rev() {
            let r = self.right[a] as usize;
            let mut u = [0.0; C];
            if r == 0 {
                for &j in self.elements_of(a) {
                    let (w, p) = (self.weight[j as usize], &psi[j as usize]);
                    for c in 0..C {
                        u[c] += w * p[c];
                    }
                }
            } else {
                let (l, rr) = (&up[a + 1], &up[r]);
                for c in 0..C {
                    u[c] = l[c] + rr[c];
                }
            }
            up[a] = u;
        }
        // 2. Far field per node, in block order.
        let mut acc: Vec<[f64; C]> = pool.install(|| {
            (0..n_nodes)
                .into_par_iter()
                .with_min_len(32)
                .map(|a| {
                    let mut f = [0.0; C];
                    for p in self.far_ptr[a]..self.far_ptr[a + 1] {
                        let (k, u) = (&self.far_coef[p], &up[self.far_col[p] as usize]);
                        for c in 0..C {
                            f[c] += k[map[c]] * u[c];
                        }
                    }
                    f
                })
                .collect()
        });
        // 3. Downward accumulation in pre-order (a parent precedes its children).
        for a in 0..n_nodes {
            let r = self.right[a] as usize;
            if r != 0 {
                let f = acc[a];
                for child in [a + 1, r] {
                    for c in 0..C {
                        acc[child][c] += f[c];
                    }
                }
            }
        }
        // 4. Rows: far part w_I G_leaf(I), then the near field.
        pool.install(|| {
            (0..n)
                .into_par_iter()
                .with_min_len(64)
                .map(|i| {
                    let g = &acc[self.leaf_of[i] as usize];
                    let w = self.weight[i];
                    let mut out = [0.0; C];
                    for c in 0..C {
                        out[c] = w * g[c];
                    }
                    for p in self.near_ptr[i]..self.near_ptr[i + 1] {
                        let (k, x) = (&self.near_val[p], &psi[self.near_col[p] as usize]);
                        for c in 0..C {
                            out[c] += k[map[c]] * x[c];
                        }
                    }
                    out
                })
                .collect()
        })
    }
}
