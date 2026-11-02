//! A block cluster tree with a cache-friendly near field (`SPEC.md` §8.3; WP17).
//!
//! The same blocks, far-field coefficients, near-field entries and product as the curves
//! crate's [`ddg_repulsive_curves::spatial::bct::BlockClusterTree`] (its module documentation
//! describes them; the block structure follows the authors' MIT code,
//! `src/product/block_cluster_tree.cpp`, Copyright (c) 2019 Christopher Yu), stored for speed:
//!
//! * the elements are numbered in **BVH order** (the order of [`Bvh::order`]), so every node is a
//!   contiguous range and the arguments of a near block are read contiguously;
//! * every inadmissible block `(𝒜, ℬ)` is a **dense** `|𝒜| × |ℬ|` array of the `N` kernels (no
//!   column indices); the diagonal of a block `(𝒜, 𝒜)` is not used;
//! * the kernel of an element pair is evaluated **once**: the entries of `(ℬ, 𝒜)` are the
//!   transpose of those of `(𝒜, ℬ)` (the kernels are bitwise symmetric, [`ClusterKernel::near`]),
//!   and the build is parallel over blocks (the curves crate's flattening is sequential);
//! * the product takes any number of channels `C` (several right-hand sides through one pass
//!   over the stored entries) and runs in parallel over the BVH leaves.
//!
//! **The arithmetic is the curves crate's, operation for operation**: per element row the far
//! part `w_I G_leaf(I)` first, then the near entries over the inadmissible blocks of the chain
//! root → leaf, each block's columns in BVH order — the order in which that crate's flattened rows
//! hold them. Products are therefore bitwise equal to the curves crate's (a structurally zero
//! entry, `near = None`, is stored as `0`, which differs at most in the sign of a zero), and
//! bitwise independent of the pool size; a channel's result does not depend on the other
//! channels, so a batch of right-hand sides gives the same bits as one at a time.

use std::time::Instant;

use ddg_repulsive_curves::spatial::Bvh;
use ddg_repulsive_curves::spatial::bct::{BctStats, ClusterKernel};
use rayon::prelude::*;

use crate::error::RsError;

enum Class {
    Admissible,
    Inadmissible,
    Split,
}

/// A block cluster tree of `N` symmetric kernels (module docs).
#[derive(Debug, Clone)]
pub struct BlockTree<const N: usize> {
    start: Vec<u32>,
    end: Vec<u32>,
    right: Vec<u32>,
    /// `order[p]` = the element at BVH position `p`.
    order: Vec<u32>,
    /// `position[i]` = the BVH position of element `i`.
    position: Vec<u32>,
    /// The leaves in pre-order (ascending ranges), and the chain root → leaf of each, flattened.
    leaves: Vec<u32>,
    chain_ptr: Vec<usize>,
    chain: Vec<u32>,
    /// Weights in BVH order.
    weight: Vec<f64>,
    far_ptr: Vec<usize>,
    far_col: Vec<u32>,
    far_coef: Vec<[f64; N]>,
    /// Inadmissible blocks grouped by row node (creation order within a node): column node and
    /// offset of the dense row-major block in `near_val`.
    near_ptr: Vec<usize>,
    near_col: Vec<u32>,
    near_off: Vec<usize>,
    near_val: Vec<[f64; N]>,
    /// `K_k 𝟙` per element (original order).
    row_sums: Vec<[f64; N]>,
    stats: BctStats,
    build_seconds: f64,
}

impl<const N: usize> BlockTree<N> {
    /// Builds the blocks of `bvh` for `kernel` (split until admissible, both BVH leaves, or
    /// `|𝒜|·|ℬ| ≤ leaf_pairs`), the far coefficients, the near entries and `K𝟙`, on `pool`.
    /// Bitwise independent of the pool size.
    ///
    /// # Errors
    /// [`RsError::InvalidOption`] if the BVH is empty.
    pub fn build<K: ClusterKernel<N>>(
        pool: &rayon::ThreadPool,
        bvh: &Bvh<3>,
        kernel: &K,
        leaf_pairs: usize,
    ) -> Result<Self, RsError> {
        let t0 = Instant::now();
        let nodes = bvh.nodes();
        let n = bvh.n_elements();
        if nodes.is_empty() || n == 0 {
            return Err(RsError::InvalidOption {
                name: "bvh",
                value: "empty".into(),
                requirement: "at least one element",
            });
        }
        let n_nodes = nodes.len();
        let start: Vec<u32> = nodes.iter().map(|nd| nd.start).collect();
        let end: Vec<u32> = nodes.iter().map(|nd| nd.end).collect();
        let right: Vec<u32> = nodes.iter().map(|nd| nd.right).collect();
        let order = bvh.order().to_vec();
        let mut position = vec![0u32; n];
        for (p, &e) in order.iter().enumerate() {
            position[e as usize] = p as u32;
        }
        let weight: Vec<f64> = order.iter().map(|&e| kernel.weight(e as usize)).collect();
        let is_leaf = |a: usize| right[a] == 0;
        let len = |a: usize| (end[a] - start[a]) as usize;

        // Blocks, level by level (the curves crate's classification and order).
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

        // Far field grouped by row node (stable).
        let (far_ptr, far_col) = group_by_row(&adm, n_nodes);
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

        // Near blocks grouped by row node; storage: canonical blocks (a ≤ b) first, then the
        // mirrored ones (a > b).
        let (near_ptr, near_col) = group_by_row(&near_blocks, n_nodes);
        let n_blocks = near_col.len();
        let mut near_off = vec![0usize; n_blocks];
        let mut canonical_total = 0usize;
        let mut block_row = vec![0u32; n_blocks];
        for a in 0..n_nodes {
            for q in near_ptr[a]..near_ptr[a + 1] {
                block_row[q] = a as u32;
            }
        }
        for q in 0..n_blocks {
            let (a, b) = (block_row[q] as usize, near_col[q] as usize);
            if a <= b {
                near_off[q] = canonical_total;
                canonical_total += len(a) * len(b);
            }
        }
        let mut total = canonical_total;
        for q in 0..n_blocks {
            let (a, b) = (block_row[q] as usize, near_col[q] as usize);
            if a > b {
                near_off[q] = total;
                total += len(a) * len(b);
            }
        }
        // (a, b) → block index, for the mirrors.
        let mut keys: Vec<(u64, u32)> = (0..n_blocks)
            .map(|q| (((block_row[q] as u64) << 32) | near_col[q] as u64, q as u32))
            .collect();
        keys.sort_unstable();
        let find = |a: usize, b: usize| -> Option<usize> {
            let k = ((a as u64) << 32) | b as u64;
            keys.binary_search_by_key(&k, |e| e.0).ok().map(|i| keys[i].1 as usize)
        };

        let mut near_val: Vec<[f64; N]> = vec![[0.0; N]; total];
        let (canon_buf, mirror_buf) = near_val.split_at_mut(canonical_total);
        // Canonical blocks (a ≤ b): one evaluation per unordered pair.
        {
            let mut tasks: Vec<(usize, &mut [[f64; N]])> = Vec::new();
            let mut rest: &mut [[f64; N]] = canon_buf;
            for q in 0..n_blocks {
                let (a, b) = (block_row[q] as usize, near_col[q] as usize);
                if a <= b {
                    let (head, tail) = rest.split_at_mut(len(a) * len(b));
                    tasks.push((q, head));
                    rest = tail;
                }
            }
            pool.install(|| {
                tasks.into_par_iter().with_min_len(16).for_each(|(q, blk)| {
                    let (a, b) = (block_row[q] as usize, near_col[q] as usize);
                    let (sa, sb) = (start[a] as usize, start[b] as usize);
                    let (ra, cb) = (len(a), len(b));
                    for r in 0..ra {
                        let i = order[sa + r] as usize;
                        let c0 = if a == b { r + 1 } else { 0 };
                        for c in c0..cb {
                            let j = order[sb + c] as usize;
                            blk[r * cb + c] = kernel.near(i, j).unwrap_or([0.0; N]);
                        }
                    }
                    if a == b {
                        for r in 0..ra {
                            for c in 0..r {
                                blk[r * cb + c] = blk[c * cb + r];
                            }
                        }
                    }
                });
            });
        }
        // Mirrored blocks (a > b): the transpose of (b, a) (computed directly if the partition
        // is not symmetric).
        {
            let canon: &[[f64; N]] = canon_buf;
            let mut tasks: Vec<(usize, &mut [[f64; N]])> = Vec::new();
            let mut rest: &mut [[f64; N]] = mirror_buf;
            for q in 0..n_blocks {
                let (a, b) = (block_row[q] as usize, near_col[q] as usize);
                if a > b {
                    let (head, tail) = rest.split_at_mut(len(a) * len(b));
                    tasks.push((q, head));
                    rest = tail;
                }
            }
            pool.install(|| {
                tasks.into_par_iter().with_min_len(16).for_each(|(q, blk)| {
                    let (a, b) = (block_row[q] as usize, near_col[q] as usize);
                    let (ra, cb) = (len(a), len(b));
                    if let Some(m) = find(b, a) {
                        let src = &canon[near_off[m]..near_off[m] + ra * cb];
                        for r in 0..ra {
                            for c in 0..cb {
                                blk[r * cb + c] = src[c * ra + r];
                            }
                        }
                    } else {
                        let (sa, sb) = (start[a] as usize, start[b] as usize);
                        for r in 0..ra {
                            let i = order[sa + r] as usize;
                            for c in 0..cb {
                                let j = order[sb + c] as usize;
                                blk[r * cb + c] = kernel.near(i, j).unwrap_or([0.0; N]);
                            }
                        }
                    }
                });
            });
        }
        // Entries that are not structurally zero (all of them for the surface kernels).
        let near_entries: usize = near_blocks
            .iter()
            .map(|&(a, b)| {
                let (a, b) = (a as usize, b as usize);
                if a == b { len(a) * (len(a) - 1) } else { len(a) * len(b) }
            })
            .sum();

        // Leaves and their chains root → leaf.
        let mut parent = vec![u32::MAX; n_nodes];
        for i in 0..n_nodes {
            if !is_leaf(i) {
                parent[i + 1] = i as u32;
                parent[right[i] as usize] = i as u32;
            }
        }
        let leaves: Vec<u32> = (0..n_nodes).filter(|&a| is_leaf(a)).map(|a| a as u32).collect();
        let mut chain_ptr = Vec::with_capacity(leaves.len() + 1);
        chain_ptr.push(0);
        let mut chain = Vec::new();
        for &l in &leaves {
            let s = chain.len();
            let mut a = l;
            loop {
                chain.push(a);
                if a == 0 {
                    break;
                }
                a = parent[a as usize];
            }
            chain[s..].reverse();
            chain_ptr.push(chain.len());
        }

        let stats = BctStats {
            admissible_blocks: adm.len(),
            inadmissible_blocks: near_blocks.len(),
            admissible_pairs,
            near_entries,
            nodes: n_nodes,
            depth,
        };
        let mut tree = BlockTree {
            start,
            end,
            right,
            order,
            position,
            leaves,
            chain_ptr,
            chain,
            weight,
            far_ptr,
            far_col,
            far_coef,
            near_ptr,
            near_col,
            near_off,
            near_val,
            row_sums: Vec::new(),
            stats,
            build_seconds: 0.0,
        };
        let mut map = [0usize; N];
        for (k, m) in map.iter_mut().enumerate() {
            *m = k;
        }
        tree.row_sums = tree.multiply(pool, &map, &vec![[1.0; N]; n]);
        tree.build_seconds = t0.elapsed().as_secs_f64();
        Ok(tree)
    }

    /// Number of elements.
    pub fn n_elements(&self) -> usize {
        self.order.len()
    }

    /// Counts of the tree (as the curves crate's).
    pub fn stats(&self) -> BctStats {
        self.stats
    }

    /// Seconds of the build.
    pub fn build_seconds(&self) -> f64 {
        self.build_seconds
    }

    /// `order()[p]` = the element at BVH position `p`.
    pub fn order(&self) -> &[u32] {
        &self.order
    }

    /// `position()[i]` = the BVH position of element `i`.
    pub fn position(&self) -> &[u32] {
        &self.position
    }

    /// `K_k 𝟙` per element (original order).
    pub fn row_sums(&self) -> &[[f64; N]] {
        &self.row_sums
    }

    /// `φ = K ψ` (original element order), channel `c` with kernel `map[c]` (all `< N`).
    ///
    /// # Panics
    /// If `psi` does not have one entry per element or `map` names a kernel `≥ N`.
    pub fn multiply<const C: usize>(
        &self,
        pool: &rayon::ThreadPool,
        map: &[usize; C],
        psi: &[[f64; C]],
    ) -> Vec<[f64; C]> {
        assert_eq!(psi.len(), self.n_elements(), "one value per element");
        let pp: Vec<[f64; C]> = self.order.iter().map(|&e| psi[e as usize]).collect();
        let out = self.multiply_permuted(pool, map, &pp);
        self.position.iter().map(|&p| out[p as usize]).collect()
    }

    /// `φ = K ψ` with `ψ` and `φ` in BVH order (`psi[p]` belongs to element `order()[p]`).
    ///
    /// [formula: rc-bct-lowrank]
    ///
    /// # Panics
    /// If `psi` does not have one entry per element or `map` names a kernel `≥ N`.
    pub fn multiply_permuted<const C: usize>(
        &self,
        pool: &rayon::ThreadPool,
        map: &[usize; C],
        psi: &[[f64; C]],
    ) -> Vec<[f64; C]> {
        assert!(map.iter().all(|&m| m < N), "kernel index below N");
        self.multiply_core::<C, 0>(pool, Some(map), psi)
    }

    /// [`BlockTree::multiply_permuted`] with the channels `0..H` on kernel 0 and `H..C` on
    /// kernel `N − 1` (`N ≤ 2`): the same operations with the channel split known at compile
    /// time (vectorized).
    ///
    /// # Panics
    /// If `psi` does not have one entry per element, `N > 2` or `H > C`.
    pub fn multiply_split<const C: usize, const H: usize>(
        &self,
        pool: &rayon::ThreadPool,
        psi: &[[f64; C]],
    ) -> Vec<[f64; C]> {
        assert!(N <= 2 && H <= C, "a split of at most two kernels");
        self.multiply_core::<C, H>(pool, None, psi)
    }

    fn multiply_core<const C: usize, const H: usize>(
        &self,
        pool: &rayon::ThreadPool,
        map: Option<&[usize; C]>,
        psi: &[[f64; C]],
    ) -> Vec<[f64; C]> {
        let n = self.n_elements();
        assert_eq!(psi.len(), n, "one value per element");
        let n_nodes = self.start.len();
        // 1. Upward sums.
        let mut up = vec![[0.0f64; C]; n_nodes];
        for a in (0..n_nodes).rev() {
            let r = self.right[a] as usize;
            let mut u = [0.0; C];
            if r == 0 {
                for p in self.start[a] as usize..self.end[a] as usize {
                    let (w, x) = (self.weight[p], &psi[p]);
                    for c in 0..C {
                        u[c] += w * x[c];
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
        // 2. Far field per node.
        let mut acc: Vec<[f64; C]> = pool.install(|| {
            (0..n_nodes)
                .into_par_iter()
                .with_min_len(32)
                .map(|a| {
                    let mut f = [0.0; C];
                    for p in self.far_ptr[a]..self.far_ptr[a + 1] {
                        let (k, u) = (&self.far_coef[p], &up[self.far_col[p] as usize]);
                        accumulate::<N, C, H>(&mut f, std::slice::from_ref(k), std::slice::from_ref(u), map);
                    }
                    f
                })
                .collect()
        });
        // 3. Downward accumulation (pre-order).
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
        // 4. Rows, per leaf: far part, then the near blocks of the chain root → leaf.
        let mut out = vec![[0.0f64; C]; n];
        let mut chunks: Vec<(usize, &mut [[f64; C]])> = Vec::with_capacity(self.leaves.len());
        let mut rest: &mut [[f64; C]] = &mut out;
        for (li, &l) in self.leaves.iter().enumerate() {
            let k = (self.end[l as usize] - self.start[l as usize]) as usize;
            let (head, tail) = rest.split_at_mut(k);
            chunks.push((li, head));
            rest = tail;
        }
        pool.install(|| {
            chunks.into_par_iter().with_min_len(8).for_each(|(li, rows)| {
                let l = self.leaves[li] as usize;
                let s = self.start[l] as usize;
                let g = &acc[l];
                for (r, o) in rows.iter_mut().enumerate() {
                    let w = self.weight[s + r];
                    for c in 0..C {
                        o[c] = w * g[c];
                    }
                }
                for &a in &self.chain[self.chain_ptr[li]..self.chain_ptr[li + 1]] {
                    let a = a as usize;
                    let ra = (self.end[a] - self.start[a]) as usize;
                    let r0 = s - self.start[a] as usize;
                    for q in self.near_ptr[a]..self.near_ptr[a + 1] {
                        let b = self.near_col[q] as usize;
                        let sb = self.start[b] as usize;
                        let cb = (self.end[b] - self.start[b]) as usize;
                        let blk = &self.near_val[self.near_off[q]..self.near_off[q] + ra * cb];
                        let xs = &psi[sb..sb + cb];
                        for (r, o) in rows.iter_mut().enumerate() {
                            let row = &blk[(r0 + r) * cb..(r0 + r + 1) * cb];
                            let mut acc = *o;
                            if a == b {
                                // Skip the diagonal entry (the element itself).
                                let d = r0 + r;
                                accumulate::<N, C, H>(&mut acc, &row[..d], &xs[..d], map);
                                accumulate::<N, C, H>(&mut acc, &row[d + 1..], &xs[d + 1..], map);
                            } else {
                                accumulate::<N, C, H>(&mut acc, row, xs, map);
                            }
                            *o = acc;
                        }
                    }
                }
            });
        });
        out
    }
}

/// `acc[c] += k[map[c]] · x[c]` over the entries, in order; without `map`, kernel 0 on the
/// channels `0..H` and kernel `N − 1` on `H..C` (the same operations, vectorizable).
#[inline(always)]
fn accumulate<const N: usize, const C: usize, const H: usize>(
    acc: &mut [f64; C],
    ks: &[[f64; N]],
    xs: &[[f64; C]],
    map: Option<&[usize; C]>,
) {
    match map {
        None => {
            for (k, x) in ks.iter().zip(xs) {
                let (k0, k1) = (k[0], k[N - 1]);
                for c in 0..H {
                    acc[c] += k0 * x[c];
                }
                for c in H..C {
                    acc[c] += k1 * x[c];
                }
            }
        }
        Some(map) => {
            for (k, x) in ks.iter().zip(xs) {
                for c in 0..C {
                    acc[c] += k[map[c]] * x[c];
                }
            }
        }
    }
}

/// Groups `(a, b)` pairs by `a` (stable): CSR pointers over `n` rows and the `b`s.
fn group_by_row(pairs: &[(u32, u32)], n: usize) -> (Vec<usize>, Vec<u32>) {
    let mut ptr = vec![0usize; n + 1];
    for &(a, _) in pairs {
        ptr[a as usize + 1] += 1;
    }
    for i in 0..n {
        ptr[i + 1] += ptr[i];
    }
    let mut fill = ptr.clone();
    let mut col = vec![0u32; pairs.len()];
    for &(a, b) in pairs {
        col[fill[a as usize]] = b;
        fill[a as usize] += 1;
    }
    (ptr, col)
}
