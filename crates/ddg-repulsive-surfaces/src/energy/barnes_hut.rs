//! Barnes–Hut energy `Ẽ` and its exact gradient for frozen interaction lists (`SPEC.md` §3.4,
//! §8.2, D4, D5; paper Eq. 17, Sec. 4).
//!
//! The faces are the elements of a [`FaceTree`] (a [`Bvh<3>`](ddg_repulsive_curves::spatial::Bvh)
//! with cluster data). For every **leaf** `ℓ` the tree is traversed from the root, depth first,
//! left child first (`ref: tpe_barnes_hut_0.cpp:51–186`, MIT, Copyright (c) 2019 Nick Sharp):
//! a node `C` with `max(r_ℓ, r_C) < θ d_box(ℓ, C)` ([`bh_admissible`]) is **lumped**, an
//! inadmissible leaf is evaluated **exactly**, an inadmissible internal node is opened. This
//! gives the interaction lists `𝓛(ℓ)` (lumped nodes) and `𝓧(ℓ)` (exact leaves), built once
//! per step ([`BhTree::build`]) and frozen for the line search ([`BhTree::refit`]). With them
//!
//! ```text
//! Ẽ = Σ_ℓ [ Σ_{C ∈ 𝓛(ℓ)} Σ_{S ∈ ℓ} K(X_S, N_S; X_C) a_S A_C  +  Σ_{C ∈ 𝓧(ℓ)} Σ_{S ∈ ℓ} Σ_{T ∈ C, T ≠ S} K(X_S, N_S; X_T) a_S a_T ],
//! A_C = Σ_{T ∈ C} a_T,   X_C = M_C / A_C,   M_C = Σ_{T ∈ C} a_T X_T.
//! ```
//!
//! Every ordered pair `(S, T)`, `S ≠ T`, is covered exactly once (the leaves below the lumped
//! and exact nodes of `ℓ` partition the faces). `θ = 0` admits nothing, so `Ẽ = ℰ̂` up to the
//! order of summation. [formula: rs-bh-energy]
//!
//! **Gradient (D4).** [`BhTree::energy_and_differential`] returns the exact derivative of `Ẽ`
//! with the lists held fixed, not the paper's forward formula of Sec. 4.2:
//!
//! * own part (the face `S` of a term): `∂/∂a_S`, `∂/∂X_S`, `∂/∂N_S` of every term of its leaf;
//! * exact pairs, partner part: `∂/∂a_T = K a_S`, `∂/∂X_T = −a_S a_T ∂K/∂v`, gathered per leaf
//!   `C` over the transposed exact lists (the lists are not symmetric);
//! * lumped terms, cluster part: with `b = Σ_S a_S K` (`= ∂/∂A_C` at fixed `X_C`) and
//!   `u = −Σ_S a_S ∂K/∂v` (`= ∂/∂M_C`, i.e. `∂/∂X_C / A_C`) summed per node over the
//!   transposed lumped lists, every member face `T ∈ C` receives
//!   `∂/∂a_T = b_C + u_C·(X_T − X_C)` and `∂/∂X_T = a_T u_C`, i.e. `∂Ẽ/∂a_T = Σ_{C∋T} ∂Ẽ/∂A_C`
//!   and `∂Ẽ/∂(a_T X_T) = Σ_{C∋T} ∂Ẽ/∂M_C` (`ref: optimized_cluster_tree.cpp:1246–1311`,
//!   `PercolateDown`). The nodes containing `T` are the ancestors of its leaf; one downward pass
//!   accumulates `β_c = β_p + U_p·(X_c − X_p) + b_c`, `U_c = U_p + u_c` (re-centred at every
//!   node, as the curves crate's `rc-bh-adjoint-shift`), and at the leaf `Λ` of `T`:
//!   `∂/∂a_T = β_Λ + U_Λ·(X_T − X_Λ)`, `∂/∂X_T = a_T U_Λ`. [formula: rs-bh-aggregate-adjoint]
//!
//! The per-face partials `(∂/∂a, ∂/∂X, ∂/∂N)` go to the vertices by WP1's chain rule
//! ([`crate::geometry::pull_back`]).
//!
//! **Determinism (`SPEC.md` §8.2).** The tree and the lists are built sequentially or per leaf
//! into per-leaf slots; the evaluation runs in parallel over leaves into per-leaf results; every
//! reduction (energy over leaves, node sums over the transposed lists, partner sums, the
//! downward pass, the vertex sums) has a fixed order. The result is bitwise independent of the
//! thread count. [`BhTree::energy`] performs the same floating-point operations as
//! [`BhTree::energy_and_differential`] for the value, so a line search compares like with like.

use rayon::prelude::*;

use crate::energy::{BhOptions, EnergyEvaluation};
use crate::error::RsError;
use crate::geometry::FaceGeometry;
use crate::geometry::v3::{axpy, dot, sub};
use crate::kernel::{Exponents, kernel_of_v, kernel_with_gradient_of_v};
use crate::spatial::{Cluster, FaceTree, bh_admissible};
use crate::surface::Surface;

/// The frozen Barnes–Hut interaction lists of a [`FaceTree`] (`SPEC.md` §3.4): per leaf (in
/// the order of [`FaceTree::leaves`]) the lumped nodes `𝓛(ℓ)` and the exactly evaluated leaves
/// `𝓧(ℓ)`, in traversal order, plus their transposes.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct InteractionLists {
    /// `lumped[lumped_off[k]..lumped_off[k + 1]]`: node indices lumped by leaf `k`.
    lumped_off: Vec<usize>,
    lumped: Vec<u32>,
    /// `exact[exact_off[k]..exact_off[k + 1]]`: node indices of the leaves evaluated exactly by
    /// leaf `k`.
    exact_off: Vec<usize>,
    exact: Vec<u32>,
    /// Per exact entry: offset of its partner block (`|C|` faces) in the flat partner array.
    partner_off: Vec<usize>,
    /// Transpose of `lumped`: for every node, the entries (indices into `lumped`) that lump it,
    /// increasing.
    lumped_t_off: Vec<usize>,
    lumped_t: Vec<u32>,
    /// Transpose of `exact`: for every leaf (ordinal), the entries (indices into `exact`) that
    /// open it, increasing.
    exact_t_off: Vec<usize>,
    exact_t: Vec<u32>,
    pairs_exact: u64,
    n_nodes: usize,
}

impl InteractionLists {
    /// The lumped nodes of the leaf with ordinal `k` (position in [`FaceTree::leaves`]).
    pub fn lumped_of(&self, k: usize) -> &[u32] {
        &self.lumped[self.lumped_off[k]..self.lumped_off[k + 1]]
    }

    /// The leaves (node indices) evaluated exactly by the leaf with ordinal `k`.
    pub fn exact_of(&self, k: usize) -> &[u32] {
        &self.exact[self.exact_off[k]..self.exact_off[k + 1]]
    }

    /// Number of (leaf, node) lumped interactions.
    pub fn n_lumped(&self) -> u64 {
        self.lumped.len() as u64
    }

    /// Number of (leaf, leaf) exact blocks.
    pub fn n_exact_blocks(&self) -> u64 {
        self.exact.len() as u64
    }

    /// Number of ordered face pairs `(S, T)`, `S ≠ T`, evaluated exactly.
    pub fn pairs_exact(&self) -> u64 {
        self.pairs_exact
    }
}

/// The Barnes–Hut evaluator of one step: the face tree, the frozen interaction lists and the
/// face data of the current positions (`SPEC.md` §3.4, §7.2). Build it once per step,
/// [`refit`](BhTree::refit) it for the line-search trials.
#[derive(Debug, Clone)]
pub struct BhTree {
    opts: BhOptions,
    faces: Vec<[u32; 3]>,
    n_vertices: usize,
    positions: Vec<[f64; 3]>,
    geom: FaceGeometry,
    tree: FaceTree,
    lists: InteractionLists,
    /// The face data in BVH order (`tree.bvh().order()`), so that a leaf's faces are contiguous
    /// in memory (the same values as `geom`).
    packed: Vec<Packed>,
    /// The lists on the GPU ([`BhTree::attach_gpu`]); with the determinism level.
    gpu: Option<(std::sync::Arc<crate::gpu::GpuBh>, crate::gpu::Determinism)>,
}

/// Barycenter, normal and area of a face, stored in BVH order.
#[derive(Debug, Clone, Copy, Default)]
struct Packed {
    x: [f64; 3],
    n: [f64; 3],
    a: f64,
}

fn pack(tree: &FaceTree, geom: &FaceGeometry, out: &mut Vec<Packed>) {
    out.clear();
    out.extend(tree.bvh().order().iter().map(|&f| {
        let f = f as usize;
        Packed {
            x: geom.barycenters[f],
            n: geom.normals[f],
            a: geom.areas[f],
        }
    }));
}

impl BhTree {
    /// Builds the face tree and the interaction lists at `positions` (validated with
    /// [`Surface::validate_positions`]). The lists are computed in parallel on `pool`.
    ///
    /// # Errors
    /// [`RsError::InvalidOption`] for bad options; the errors of
    /// [`Surface::validate_positions`].
    pub fn build(
        pool: &rayon::ThreadPool,
        surface: &Surface,
        positions: &[[f64; 3]],
        opts: BhOptions,
    ) -> Result<BhTree, RsError> {
        let geom = FaceGeometry::new(surface, positions)?;
        Self::build_with_geometry(pool, surface, positions, geom, opts)
    }

    /// As [`BhTree::build`] with the face geometry of `positions` already computed (positions
    /// must be valid).
    ///
    /// # Errors
    /// [`RsError::InvalidOption`] for bad options or a `geom` of another face count;
    /// [`RsError::PositionCountMismatch`].
    pub fn build_with_geometry(
        pool: &rayon::ThreadPool,
        surface: &Surface,
        positions: &[[f64; 3]],
        geom: FaceGeometry,
        opts: BhOptions,
    ) -> Result<BhTree, RsError> {
        opts.validate()?;
        let tree = FaceTree::build(surface, positions, &geom, opts.leaf_size)?;
        let lists = build_lists(pool, &tree, opts.theta);
        let mut packed = Vec::new();
        pack(&tree, &geom, &mut packed);
        Ok(BhTree {
            opts,
            faces: surface.faces().to_vec(),
            n_vertices: surface.n_vertices(),
            positions: positions.to_vec(),
            geom,
            tree,
            lists,
            packed,
            gpu: None,
        })
    }

    /// Moves the evaluator to `positions` keeping the tree topology, its boxes and radii and the
    /// interaction lists (`SPEC.md` §7.2): refreshes the face data and the node aggregates
    /// `A_C`, `X_C`, `P̄_C`. A degenerate face is not an error here: it makes the energy
    /// non-finite (a rejected line-search trial).
    ///
    /// # Errors
    /// [`RsError::PositionCountMismatch`], [`RsError::NonFinitePosition`].
    pub fn refit(&mut self, positions: &[[f64; 3]]) -> Result<(), RsError> {
        if positions.len() != self.n_vertices {
            return Err(RsError::PositionCountMismatch {
                expected: self.n_vertices,
                found: positions.len(),
            });
        }
        if let Some(v) = positions
            .iter()
            .position(|p| !p.iter().all(|x| x.is_finite()))
        {
            return Err(RsError::NonFinitePosition { vertex: v as u32 });
        }
        self.positions.clear();
        self.positions.extend_from_slice(positions);
        let g = &mut self.geom;
        for (f, t) in self.faces.iter().enumerate() {
            let (a, n, x) =
                crate::geometry::face_quantities(crate::geometry::corners(positions, t));
            g.areas[f] = a;
            g.normals[f] = n;
            g.barycenters[f] = x;
        }
        pack(&self.tree, &self.geom, &mut self.packed);
        self.tree.refit_faces(&self.faces, positions, &self.geom)
    }

    /// The options of the build.
    pub fn options(&self) -> &BhOptions {
        &self.opts
    }

    /// The face tree (topology, cluster data; WP6 builds the block cluster tree on it).
    pub fn face_tree(&self) -> &FaceTree {
        &self.tree
    }

    /// The frozen interaction lists.
    pub fn lists(&self) -> &InteractionLists {
        &self.lists
    }

    /// The face geometry of the current positions.
    pub fn geometry(&self) -> &FaceGeometry {
        &self.geom
    }

    /// The current positions.
    pub fn positions(&self) -> &[[f64; 3]] {
        &self.positions
    }

    /// Evaluates [`BhTree::try_energy`] and [`BhTree::energy_and_differential`] on the GPU of
    /// `gpu` from now on (uploads the frozen lists; the face data and the node aggregates go up
    /// per evaluation, so [`BhTree::refit`] needs nothing more; `NOTES-GPU.md`). The kernels
    /// repeat the CPU's operations in its order: bitwise equal at
    /// [`crate::gpu::Determinism::MatchCpu`] for integer exponents; at that level non-integer
    /// exponents (whose `powf` the device cannot reproduce) keep the CPU path.
    ///
    /// # Errors
    /// [`RsError::Gpu`].
    pub fn attach_gpu(&mut self, gpu: &std::sync::Arc<crate::gpu::GpuContext>) -> Result<(), RsError> {
        #[cfg(feature = "cuda")]
        {
            let tree = &self.tree;
            let nodes = tree.bvh().nodes();
            let nstart: Vec<u32> = nodes.iter().map(|n| n.start).collect();
            let nend: Vec<u32> = nodes.iter().map(|n| n.end).collect();
            let mut leaf_of = vec![0u32; self.faces.len()];
            for (k, &l) in tree.leaves().iter().enumerate() {
                for p in nstart[l as usize]..nend[l as usize] {
                    leaf_of[p as usize] = k as u32;
                }
            }
            let l = &self.lists;
            let host = crate::gpu::cuda::BhHost {
                leaf_of: &leaf_of,
                leaves: tree.leaves(),
                nstart: &nstart,
                nend: &nend,
                lumped_off: &l.lumped_off,
                lumped: &l.lumped,
                exact_off: &l.exact_off,
                exact: &l.exact,
                partner_off: &l.partner_off,
                lumped_t_off: &l.lumped_t_off,
                lumped_t: &l.lumped_t,
                exact_t_off: &l.exact_t_off,
                exact_t: &l.exact_t,
                order: tree.bvh().order(),
            };
            let dev = crate::gpu::cuda::DeviceBh::new(gpu, &host)?;
            self.gpu = Some((std::sync::Arc::new(crate::gpu::GpuBh { inner: dev }), gpu.determinism()));
            Ok(())
        }
        #[cfg(not(feature = "cuda"))]
        {
            let _ = gpu;
            Err(RsError::Gpu {
                operation: "upload",
                detail: "built without the `cuda` feature".into(),
            })
        }
    }

    /// Whether the evaluations run on the GPU (for exponents `e`).
    pub fn on_gpu(&self, e: &Exponents) -> bool {
        match &self.gpu {
            None => false,
            Some((_, d)) => {
                let (a, b) = e.int_exponents();
                *d == crate::gpu::Determinism::RunToRun || (a.is_some() && b.is_some())
            }
        }
    }

    /// The face data and cluster aggregates as flat arrays (BVH positions, nodes).
    #[cfg(feature = "cuda")]
    fn device_inputs(&self) -> [Vec<f64>; 5] {
        let pk = &self.packed;
        let px = pk.iter().flat_map(|p| p.x).collect();
        let pn = pk.iter().flat_map(|p| p.n).collect();
        let pa = pk.iter().map(|p| p.a).collect();
        let cl = self.tree.clusters();
        let cc = cl.iter().flat_map(|c| c.centroid).collect();
        let ca = cl.iter().map(|c| c.area).collect();
        [px, pn, pa, cc, ca]
    }

    #[cfg(feature = "cuda")]
    fn exps(e: &Exponents) -> crate::gpu::cuda::Exps {
        let (a, b) = e.int_exponents();
        crate::gpu::cuda::Exps {
            alpha: e.alpha(),
            beta: e.beta(),
            alpha_int: a.map_or(-1, |x| x as i32),
            beta_int: b.map_or(-1, |x| x as i32),
        }
    }

    /// [`BhTree::energy`], on the GPU if one is attached ([`BhTree::attach_gpu`]).
    ///
    /// [formula: rs-bh-energy]
    ///
    /// # Errors
    /// [`RsError::Gpu`] on a device failure.
    pub fn try_energy(&self, pool: &rayon::ThreadPool, e: &Exponents) -> Result<f64, RsError> {
        #[cfg(feature = "cuda")]
        if self.on_gpu(e)
            && let Some((g, _)) = &self.gpu
        {
            let [px, pn, pa, cc, ca] = self.device_inputs();
            let per_leaf = g.inner.energy(Self::exps(e), &px, &pn, &pa, &cc, &ca)?;
            return Ok(per_leaf.iter().sum());
        }
        Ok(self.energy(pool, e))
    }

    /// `Ẽ` at the current positions on the frozen lists. May be non-finite (degenerate faces,
    /// coincident barycenters): the caller decides (a line search rejects the trial).
    ///
    /// [formula: rs-bh-energy]
    pub fn energy(&self, pool: &rayon::ThreadPool, e: &Exponents) -> f64 {
        let nl = self.tree.leaves().len();
        let per_leaf: Vec<f64> = pool.install(|| {
            (0..nl)
                .into_par_iter()
                .map(|k| self.leaf_energy(k, e))
                .collect()
        });
        per_leaf.iter().sum()
    }

    /// `Ẽ` and its exact differential `∂Ẽ/∂f_v` on the frozen lists (module docs, D4), with
    /// the pair and node counts.
    ///
    /// [formula: rs-bh-energy]
    /// [formula: rs-bh-aggregate-adjoint]
    ///
    /// # Errors
    /// [`RsError::NonFiniteEnergy`] naming a pair of faces when the energy or a differential
    /// entry is not finite: the first exact pair (in leaf and list order) with a non-finite
    /// term, else the exact pair with the largest term.
    pub fn energy_and_differential(
        &self,
        pool: &rayon::ThreadPool,
        e: &Exponents,
    ) -> Result<EnergyEvaluation, RsError> {
        let (energy, d_area, d_bary, d_normal) = self.try_face_partials(pool, e)?;
        let differential = pull_back_faces(
            &self.faces,
            self.n_vertices,
            &self.positions,
            &d_area,
            &d_bary,
            &d_normal,
        );
        if !energy.is_finite() || !differential.iter().flatten().all(|x| x.is_finite()) {
            return Err(RsError::NonFiniteEnergy {
                faces: self.culprit_pair(e),
            });
        }
        Ok(EnergyEvaluation {
            energy,
            differential,
            pairs_exact: self.lists.pairs_exact,
            nodes_lumped: self.lists.n_lumped(),
        })
    }

    /// [`BhTree::face_partials`], on the GPU if one is attached ([`BhTree::attach_gpu`]).
    ///
    /// [formula: rs-bh-aggregate-adjoint]
    ///
    /// # Errors
    /// [`RsError::Gpu`] on a device failure.
    #[allow(clippy::type_complexity)]
    pub fn try_face_partials(
        &self,
        pool: &rayon::ThreadPool,
        e: &Exponents,
    ) -> Result<(f64, Vec<f64>, Vec<[f64; 3]>, Vec<[f64; 3]>), RsError> {
        #[cfg(feature = "cuda")]
        if self.on_gpu(e)
            && let Some((g, _)) = &self.gpu
        {
            let tree = &self.tree;
            let clusters = tree.clusters();
            let nodes = tree.bvh().nodes();
            let nn = clusters.len();
            // The CPU's downward pass (step 3 of face_partials) on the node sums.
            let downward = |sums: &[f64]| -> (Vec<f64>, Vec<f64>) {
                let su = |n: usize| [sums[4 * n + 1], sums[4 * n + 2], sums[4 * n + 3]];
                let mut beta = vec![0.0; nn];
                let mut ubig = vec![[0.0; 3]; nn];
                if nn > 0 {
                    beta[0] = sums[0];
                    ubig[0] = su(0);
                }
                for (n, node) in nodes.iter().enumerate() {
                    if node.is_leaf() {
                        continue;
                    }
                    for c in [n + 1, node.right as usize] {
                        let shift = sub(clusters[c].centroid, clusters[n].centroid);
                        beta[c] = beta[n] + dot(ubig[n], shift) + sums[4 * c];
                        ubig[c] = axpy(ubig[n], 1.0, su(c));
                    }
                }
                (beta, ubig.into_iter().flatten().collect())
            };
            let [px, pn, pa, cc, ca] = self.device_inputs();
            let (per_leaf, d) = g.inner.partials(Self::exps(e), &px, &pn, &pa, &cc, &ca, downward)?;
            let energy: f64 = per_leaf.iter().sum();
            let nf = self.faces.len();
            let mut d_area = vec![0.0; nf];
            let mut d_bary = vec![[0.0; 3]; nf];
            let mut d_normal = vec![[0.0; 3]; nf];
            for f in 0..nf {
                let r = &d[7 * f..7 * f + 7];
                d_area[f] = r[0];
                d_bary[f] = [r[1], r[2], r[3]];
                d_normal[f] = [r[4], r[5], r[6]];
            }
            return Ok((energy, d_area, d_bary, d_normal));
        }
        Ok(self.face_partials(pool, e))
    }

    /// `Ẽ` and its per-face partials `(∂Ẽ/∂a_S, ∂Ẽ/∂X_S, ∂Ẽ/∂N_S)` (all paths: own, partner,
    /// cluster), indexed by face. Non-finite values are returned as they are.
    ///
    /// [formula: rs-bh-aggregate-adjoint]
    #[allow(clippy::type_complexity)]
    pub fn face_partials(
        &self,
        pool: &rayon::ThreadPool,
        e: &Exponents,
    ) -> (f64, Vec<f64>, Vec<[f64; 3]>, Vec<[f64; 3]>) {
        let tree = &self.tree;
        let lists = &self.lists;
        let leaves = tree.leaves();
        let nl = leaves.len();
        let nn = tree.clusters().len();
        let nf = self.faces.len();
        // 1. Per leaf, in parallel: energy, own partials, node terms, partner blocks.
        let outs: Vec<LeafOut> = pool.install(|| {
            (0..nl)
                .into_par_iter()
                .map(|k| self.leaf_with_gradient(k, e))
                .collect()
        });
        let energy: f64 = outs.iter().map(|o| o.energy).sum();
        // Flatten the node terms and partner blocks in entry order.
        let mut node_terms: Vec<(f64, [f64; 3])> = Vec::with_capacity(lists.lumped.len());
        let mut partner: Vec<(f64, [f64; 3])> =
            Vec::with_capacity(lists.partner_off.last().copied().unwrap_or(0));
        for o in &outs {
            node_terms.extend_from_slice(&o.node);
            partner.extend_from_slice(&o.partner);
        }
        // 2. Node sums b_C, u_C over the transposed lumped lists (increasing entry).
        let sums: Vec<(f64, [f64; 3])> = pool.install(|| {
            (0..nn)
                .into_par_iter()
                .with_min_len(64)
                .map(|n| {
                    let (mut b, mut u) = (0.0, [0.0; 3]);
                    for &k in &lists.lumped_t[lists.lumped_t_off[n]..lists.lumped_t_off[n + 1]] {
                        let (bk, uk) = node_terms[k as usize];
                        b += bk;
                        u = axpy(u, 1.0, uk);
                    }
                    (b, u)
                })
                .collect()
        });
        // 3. Downward pass (pre-order: parents precede children).
        let clusters = tree.clusters();
        let mut beta = vec![0.0; nn];
        let mut ubig = vec![[0.0; 3]; nn];
        if nn > 0 {
            beta[0] = sums[0].0;
            ubig[0] = sums[0].1;
        }
        for (n, node) in tree.bvh().nodes().iter().enumerate() {
            if node.is_leaf() {
                continue;
            }
            for c in [n + 1, node.right as usize] {
                let shift = sub(clusters[c].centroid, clusters[n].centroid);
                beta[c] = beta[n] + dot(ubig[n], shift) + sums[c].0;
                ubig[c] = axpy(ubig[n], 1.0, sums[c].1);
            }
        }
        // 4. Per leaf C, in parallel: own + partner (transposed exact lists) + cluster part.
        let per_leaf: Vec<Vec<(f64, [f64; 3], [f64; 3])>> = pool.install(|| {
            (0..nl)
                .into_par_iter()
                .map(|kc| {
                    let c = leaves[kc] as usize;
                    let faces = tree.faces_of(c);
                    let o = &outs[kc];
                    let mut res: Vec<(f64, [f64; 3], [f64; 3])> = o.own.clone();
                    for &ent in &lists.exact_t[lists.exact_t_off[kc]..lists.exact_t_off[kc + 1]] {
                        let off = lists.partner_off[ent as usize];
                        for (t, r) in res.iter_mut().enumerate() {
                            let (pa, px) = partner[off + t];
                            r.0 += pa;
                            r.1 = axpy(r.1, 1.0, px);
                        }
                    }
                    let xc = clusters[c].centroid;
                    for (t, r) in res.iter_mut().enumerate() {
                        let f = faces[t] as usize;
                        let a = self.geom.areas[f];
                        let x = self.geom.barycenters[f];
                        r.0 += beta[c] + dot(ubig[c], sub(x, xc));
                        r.1 = axpy(r.1, a, ubig[c]);
                    }
                    res
                })
                .collect()
        });
        let mut d_area = vec![0.0; nf];
        let mut d_bary = vec![[0.0; 3]; nf];
        let mut d_normal = vec![[0.0; 3]; nf];
        for (kc, res) in per_leaf.into_iter().enumerate() {
            for (&f, (da, dx, dn)) in tree.faces_of(leaves[kc] as usize).iter().zip(res) {
                d_area[f as usize] = da;
                d_bary[f as usize] = dx;
                d_normal[f as usize] = dn;
            }
        }
        (energy, d_area, d_bary, d_normal)
    }

    /// The energy of leaf `k`'s terms: the lumped nodes, then the exact leaves, in list order.
    /// The same operations as [`BhTree::leaf_with_gradient`] for the value.
    fn leaf_energy(&self, k: usize, e: &Exponents) -> f64 {
        let tree = &self.tree;
        let pk = &self.packed;
        let srange = self.range(tree.leaves()[k] as usize);
        let mut energy = 0.0;
        for &c in self.lists.lumped_of(k) {
            let cl = &tree.clusters()[c as usize];
            let mut b = 0.0;
            for fs in &pk[srange.clone()] {
                let v = sub(fs.x, cl.centroid);
                b += fs.a * kernel_of_v(e, v, fs.n);
            }
            energy += b * cl.area;
        }
        for &c in self.lists.exact_of(k) {
            let trange = self.range(c as usize);
            for ps in srange.clone() {
                let fs = &pk[ps];
                let (xs, ns) = (fs.x, fs.n);
                let mut row = 0.0;
                for pt in trange.clone() {
                    if pt == ps {
                        continue;
                    }
                    let ft = &pk[pt];
                    row += ft.a * kernel_of_v(e, sub(xs, ft.x), ns);
                }
                energy += fs.a * row;
            }
        }
        energy
    }

    /// The BVH positions of node `n`'s faces.
    #[inline]
    fn range(&self, n: usize) -> std::ops::Range<usize> {
        let nd = &self.tree.bvh().nodes()[n];
        nd.start as usize..nd.end as usize
    }

    /// Leaf `k`: its energy (as [`BhTree::leaf_energy`]), the own partials of its faces, the
    /// node terms `(b, u)` of its lumped entries and the partner partials of its exact blocks.
    ///
    /// [formula: rs-bh-aggregate-adjoint]
    fn leaf_with_gradient(&self, k: usize, e: &Exponents) -> LeafOut {
        let tree = &self.tree;
        let pk = &self.packed;
        let srange = self.range(tree.leaves()[k] as usize);
        let ns_len = srange.len();
        let mut out = LeafOut {
            energy: 0.0,
            own: vec![(0.0, [0.0; 3], [0.0; 3]); ns_len],
            node: Vec::with_capacity(self.lists.lumped_of(k).len()),
            partner: Vec::new(),
        };
        for &c in self.lists.lumped_of(k) {
            let cl: &Cluster = &tree.clusters()[c as usize];
            let ac = cl.area;
            let mut b = 0.0;
            let mut u = [0.0; 3];
            for (i, fs) in pk[srange.clone()].iter().enumerate() {
                let a_s = fs.a;
                let v = sub(fs.x, cl.centroid);
                let kv = kernel_with_gradient_of_v(e, v, fs.n);
                b += a_s * kv.value;
                let o = &mut out.own[i];
                o.0 += kv.value * ac;
                o.1 = axpy(o.1, a_s * ac, kv.d_v);
                o.2 = axpy(o.2, a_s * ac, kv.d_n);
                u = axpy(u, -a_s, kv.d_v);
            }
            out.energy += b * ac;
            out.node.push((b, u));
        }
        for &c in self.lists.exact_of(k) {
            let trange = self.range(c as usize);
            let base = out.partner.len();
            out.partner.resize(base + trange.len(), (0.0, [0.0; 3]));
            for (i, ps) in srange.clone().enumerate() {
                let fs = &pk[ps];
                let (xs, ns, a_s) = (fs.x, fs.n, fs.a);
                let mut row = 0.0;
                let mut dx = [0.0; 3];
                let mut dn = [0.0; 3];
                for (j, pt) in trange.clone().enumerate() {
                    if pt == ps {
                        continue;
                    }
                    let ft = &pk[pt];
                    let a_t = ft.a;
                    let kv = kernel_with_gradient_of_v(e, sub(xs, ft.x), ns);
                    row += a_t * kv.value;
                    dx = axpy(dx, a_t, kv.d_v);
                    dn = axpy(dn, a_t, kv.d_n);
                    let p = &mut out.partner[base + j];
                    p.0 += a_s * kv.value;
                    p.1 = axpy(p.1, -a_s * a_t, kv.d_v);
                }
                out.energy += a_s * row;
                let o = &mut out.own[i];
                o.0 += row;
                o.1 = axpy(o.1, a_s, dx);
                o.2 = axpy(o.2, a_s, dn);
            }
        }
        out
    }

    /// The exact pair `(S, T)` of the lists with the first non-finite term `K a_S a_T` or
    /// kernel gradient (leaf and list order), else the one with the largest term.
    fn culprit_pair(&self, e: &Exponents) -> [u32; 2] {
        let tree = &self.tree;
        let g = &self.geom;
        let mut best = ([0u32, 0u32], f64::NEG_INFINITY);
        for (k, &l) in tree.leaves().iter().enumerate() {
            for &c in self.lists.exact_of(k) {
                for &s in tree.faces_of(l as usize) {
                    for &t in tree.faces_of(c as usize) {
                        if s == t {
                            continue;
                        }
                        let (su, tu) = (s as usize, t as usize);
                        let kv = kernel_with_gradient_of_v(
                            e,
                            sub(g.barycenters[su], g.barycenters[tu]),
                            g.normals[su],
                        );
                        let term = kv.value * g.areas[su] * g.areas[tu];
                        let finite =
                            term.is_finite() && kv.d_v.iter().chain(&kv.d_n).all(|x| x.is_finite());
                        if !finite {
                            return [s, t];
                        }
                        if term > best.1 {
                            best = ([s, t], term);
                        }
                    }
                }
            }
        }
        best.0
    }
}

/// Per-leaf result of the gradient evaluation.
struct LeafOut {
    energy: f64,
    /// Own partials `(∂/∂a, ∂/∂X, ∂/∂N)` of the leaf's faces, in [`FaceTree::faces_of`] order.
    own: Vec<(f64, [f64; 3], [f64; 3])>,
    /// `(b, u)` per lumped entry.
    node: Vec<(f64, [f64; 3])>,
    /// `(∂/∂a_T, ∂/∂X_T)` of the partner faces, one block of `|C|` per exact entry.
    partner: Vec<(f64, [f64; 3])>,
}

/// The interaction lists of `tree` for `θ` (`SPEC.md` §3.4): per leaf a depth-first traversal
/// from the root, left child first.
///
/// [formula: rs-bh-energy]
/// [formula: rs-box-distance]
pub fn build_lists(pool: &rayon::ThreadPool, tree: &FaceTree, theta: f64) -> InteractionLists {
    let leaves = tree.leaves();
    let nl = leaves.len();
    let nodes = tree.bvh().nodes();
    let nn = nodes.len();
    let clusters = tree.clusters();
    let per_leaf: Vec<(Vec<u32>, Vec<u32>)> = pool.install(|| {
        leaves
            .par_iter()
            .map(|&l| {
                let cl = &clusters[l as usize];
                let (mut lumped, mut exact) = (Vec::new(), Vec::new());
                let mut stack: Vec<usize> = Vec::with_capacity(64);
                stack.push(0);
                while let Some(c) = stack.pop() {
                    if bh_admissible(cl, &clusters[c], theta) {
                        lumped.push(c as u32);
                    } else if nodes[c].is_leaf() {
                        exact.push(c as u32);
                    } else {
                        stack.push(nodes[c].right as usize);
                        stack.push(c + 1);
                    }
                }
                (lumped, exact)
            })
            .collect()
    });
    let mut lists = InteractionLists {
        lumped_off: Vec::with_capacity(nl + 1),
        exact_off: Vec::with_capacity(nl + 1),
        n_nodes: nn,
        ..Default::default()
    };
    lists.lumped_off.push(0);
    lists.exact_off.push(0);
    let mut poff = 0usize;
    for (k, (lumped, exact)) in per_leaf.into_iter().enumerate() {
        let ns = tree.faces_of(leaves[k] as usize).len() as u64;
        for &c in &exact {
            let nc = tree.faces_of(c as usize).len();
            lists.partner_off.push(poff);
            poff += nc;
            lists.pairs_exact += if c == leaves[k] {
                ns * (ns - 1)
            } else {
                ns * nc as u64
            };
        }
        lists.lumped.extend(lumped);
        lists.exact.extend(exact);
        lists.lumped_off.push(lists.lumped.len());
        lists.exact_off.push(lists.exact.len());
    }
    lists.partner_off.push(poff);
    let (o, t) = transpose(&lists.lumped, nn, |n| n as usize);
    lists.lumped_t_off = o;
    lists.lumped_t = t;
    let (o, t) = transpose(&lists.exact, nl, |n| {
        tree.leaf_ordinal(n as usize)
            .expect("exact entries are leaves")
    });
    lists.exact_t_off = o;
    lists.exact_t = t;
    lists
}

/// For a flat list of targets, the entries grouped by target (`key`), increasing within a group.
fn transpose(entries: &[u32], m: usize, key: impl Fn(u32) -> usize) -> (Vec<usize>, Vec<u32>) {
    let mut off = vec![0usize; m + 1];
    for &x in entries {
        off[key(x) + 1] += 1;
    }
    for i in 0..m {
        off[i + 1] += off[i];
    }
    let mut fill = off.clone();
    let mut idx = vec![0u32; entries.len()];
    for (k, &x) in entries.iter().enumerate() {
        let f = &mut fill[key(x)];
        idx[*f] = k as u32;
        *f += 1;
    }
    (off, idx)
}

/// [`crate::geometry::pull_back`] on a face list (the evaluator keeps faces, not the `Surface`).
fn pull_back_faces(
    faces: &[[u32; 3]],
    n_vertices: usize,
    positions: &[[f64; 3]],
    d_area: &[f64],
    d_bary: &[[f64; 3]],
    d_normal: &[[f64; 3]],
) -> Vec<[f64; 3]> {
    let mut out = vec![[0.0; 3]; n_vertices];
    for (f, t) in faces.iter().enumerate() {
        let gr = crate::geometry::face_pull_back(
            crate::geometry::corners(positions, t),
            d_area[f],
            d_bary[f],
            d_normal[f],
        );
        for c in 0..3 {
            let v = t[c] as usize;
            out[v] = axpy(out[v], 1.0, gr[c]);
        }
    }
    out
}

/// Hook for [`crate::energy::energy`]: builds the evaluator and returns `Ẽ`, `∂Ẽ/∂f_v` and the
/// pair/node counts. `geom` holds the validated face quantities of `positions`; `options` are
/// validated.
///
/// # Errors
/// As [`BhTree::build_with_geometry`] and [`BhTree::energy_and_differential`].
pub(crate) fn energy_and_differential(
    pool: &rayon::ThreadPool,
    surface: &Surface,
    positions: &[[f64; 3]],
    geom: &FaceGeometry,
    exponents: &Exponents,
    options: &BhOptions,
) -> Result<EnergyEvaluation, RsError> {
    let tree = BhTree::build_with_geometry(pool, surface, positions, geom.clone(), *options)?;
    tree.energy_and_differential(pool, exponents)
}
