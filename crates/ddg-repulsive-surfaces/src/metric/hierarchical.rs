//! Hierarchical (block cluster tree) products with `B`, `B⁰` and the fractional operators
//! (`SPEC.md` §4.3, §6.3, §8.3, App. A, B).
//!
//! Every operator here is a kernel-matrix form (App. A, `rs-kernel-matrix-form`)
//!
//! ```text
//! L u = 2 Xᵀ (diag(a)⁻¹ diag(H a) − H) X u = 2 Mᵀ diag(a) (diag(H a) M u − H diag(a) M u),   X = diag(a) M
//! ```
//!
//! with a symmetric kernel matrix `H` (zero diagonal) replaced by its block-cluster-tree
//! approximation `H̃` (the curves crate's block cluster tree, stored as [`BlockTree`]: 0th-order
//! multipole, rank one per
//! admissible block, exact near field), and `H̃ a` computed **once per tree with the same `H̃`**
//! (`ref: optimized_bct.cpp:936–1019`), so that `L 1 = 0` holds exactly for the compressed
//! operator when `M 1 = 1` (`rs-hierarchical-null-space`): the product evaluates
//! `diag(H̃a) M u − H̃ (a ⊙ M u)`, which for `u = 1` is `H̃a − H̃a` with bitwise equal arguments.
//!
//! * [`FaceMetricTree`]: the face pairs of §4.2 on the face BVH of WP3 ([`FaceTree`]). Two
//!   kernels on one tree (`N = 2`): `h_hi = ρ^{−2s}` (for `B`, through `M = D_f`, three channels
//!   per input column) and `h_lo = ½(K₂(S,T) + K₂(T,S)) ρ^{−2s}` (for `B⁰`, through the face
//!   average). Far blocks use the cluster data: `ρ = |X_I − X_J|` and, for `h_lo`, the
//!   **area-averaged projectors** `½(vᵀP̄_I v + vᵀP̄_J v) ρ^{−4}` (`rs-projector-average`: `vᵀP̄v`
//!   is the area average of `(N_S·v)²`); near blocks use the face data, with the same expressions
//!   as the dense kernels of [`super::dense`] (so `χ = 0` reproduces them to round-off).
//! * [`FaceFractionalTree`]: the face form of `L_{2−s}` (`h_fr = ρ^{2s−6}`, face averages; the
//!   paper's and the reference's preconditioner factor, `FractionalForm::FacePairs`).
//! * [`VertexFractionalTree`]: the vertex-pair form `L^V_{2−s}` of D32 (`SPEC.md` §6.3): a BVH over
//!   the vertices (point `f_v`, mass the dual area `a_v`), kernel `|f_v − f_w|^{2s−6}`,
//!   `M = I`, weights `a_v`.
//!
//! **Separation** (§8.3): blocks `(I, J)` of distinct nodes with
//! `max(r_I, r_J)² ≤ χ² d_box(I, J)²` ([`separated`]; `≤`, unlike Barnes–Hut's `<`, as the
//! reference), where `r_I` is the distance from the cluster centroid to the farthest corner of the
//! cluster's box (for faces: the box of the triangles' vertices; for vertices: of the points).
//! Touching or overlapping boxes are never separated; `χ = 0` makes every block exact. The
//! blocks are those of the curves crate's tree (both non-leaf sides split; the reference splits
//! the larger cluster), which changes the approximation only within the separation bound.
//!
//! Vector inputs are the interleaved `3|V|` position vectors (`3 v + axis`, as `A′`): the three
//! columns go through one tree pass (App. B; 12 channels for `B + B⁰`). All parallel loops write
//! disjoint outputs, the tree's product is deterministic: results are bitwise independent of the
//! size of the pool.

// `!(x > tol)` is deliberate throughout: a NaN must fail the test.
#![allow(clippy::neg_cmp_op_on_partial_ord)]

use std::time::Instant;

use ddg_repulsive_curves::spatial::bct::{BctStats, ClusterKernel};
use ddg_repulsive_curves::spatial::{Aabb, Bvh, Element};
use rayon::prelude::*;

use super::{FaceDiagonal, FaceOperators};
use crate::error::RsError;
use crate::geometry::v3::{dot, sub};
use crate::geometry::{FaceGeometry, vertex_areas};
use crate::kernel::Exponents;
use crate::spatial::block_tree::BlockTree;
use crate::spatial::{FaceTree, box_distance_sq, farthest_corner_distance};
use crate::surface::Surface;

/// Parameters of the block cluster trees (`SPEC.md` §8.3).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ClusterTreeOptions {
    /// The separation parameter `χ` of [`separated`] (default 0.5, as the reference's `θ`;
    /// 0 = every block exact).
    pub chi: f64,
    /// Maximal number of elements (faces, vertices) in a BVH leaf (default 8).
    pub leaf_size: usize,
    /// A block of `|I|·|J| ≤ leaf_pairs` element pairs that is not separated is not split further
    /// (default 64, i.e. split down to the leaves for leaf size 8).
    pub leaf_pairs: usize,
}

impl Default for ClusterTreeOptions {
    fn default() -> Self {
        ClusterTreeOptions {
            chi: 0.5,
            leaf_size: 8,
            leaf_pairs: 64,
        }
    }
}

impl ClusterTreeOptions {
    /// Checks the options.
    ///
    /// # Errors
    /// [`RsError::InvalidOption`] for a negative or non-finite `chi` or `leaf_size = 0`.
    pub fn validate(&self) -> Result<(), RsError> {
        if !(self.chi.is_finite() && self.chi >= 0.0) {
            return Err(RsError::InvalidOption {
                name: "chi",
                value: self.chi.to_string(),
                requirement: "a finite number ≥ 0 (0: exact products)",
            });
        }
        if self.leaf_size == 0 {
            return Err(RsError::InvalidOption {
                name: "leaf_size",
                value: "0".into(),
                requirement: "at least 1",
            });
        }
        Ok(())
    }
}

/// The separation of two clusters (`SPEC.md` §8.3): `max(r_I, r_J)² ≤ χ² d²` with `d` the box
/// distance ([`box_distance_sq`]), and `d > 0` (touching or overlapping boxes are never
/// separated); `χ = 0` separates nothing.
///
/// [formula: rs-box-distance]
#[inline]
pub fn separated(r_a: f64, box_a: &Aabb<3>, r_b: f64, box_b: &Aabb<3>, chi: f64) -> bool {
    if !(chi > 0.0) {
        return false;
    }
    let d2 = box_distance_sq(box_a, box_b);
    let h = r_a.max(r_b);
    d2 > 0.0 && h * h <= chi * chi * d2
}

/// Vertex → (face, corner) incidence, faces ascending (for deterministic gathers).
#[derive(Debug, Clone)]
struct Incidence {
    ptr: Vec<usize>,
    entries: Vec<(u32, u8)>,
}

impl Incidence {
    fn new(faces: &[[u32; 3]], n_vertices: usize) -> Incidence {
        let mut ptr = vec![0usize; n_vertices + 1];
        for t in faces {
            for &v in t {
                ptr[v as usize + 1] += 1;
            }
        }
        for v in 0..n_vertices {
            ptr[v + 1] += ptr[v];
        }
        let mut fill = ptr.clone();
        let mut entries = vec![(0u32, 0u8); ptr[n_vertices]];
        for (s, t) in faces.iter().enumerate() {
            for (c, &v) in t.iter().enumerate() {
                entries[fill[v as usize]] = (s as u32, c as u8);
                fill[v as usize] += 1;
            }
        }
        Incidence { ptr, entries }
    }

    fn of(&self, v: usize) -> &[(u32, u8)] {
        &self.entries[self.ptr[v]..self.ptr[v + 1]]
    }
}

fn check_len(name: &'static str, got: usize, want: usize) -> Result<(), RsError> {
    if got != want {
        return Err(RsError::InvalidOption {
            name,
            value: format!("{got} entries"),
            requirement: "3|V| interleaved entries (3 v + axis)",
        });
    }
    Ok(())
}

/// The first face whose `H̃ a` is not finite, with a face at the same barycenter (or the face
/// itself), as [`RsError::NonFiniteEnergy`].
fn coincident_faces(ops: &FaceOperators, bad: usize) -> RsError {
    let x = ops.barycenters[bad];
    let t = (0..ops.n_faces())
        .find(|&t| t != bad && ops.barycenters[t] == x)
        .unwrap_or(bad);
    let (a, b) = if bad < t { (bad, t) } else { (t, bad) };
    RsError::NonFiniteEnergy {
        faces: [a as u32, b as u32],
    }
}

/// `ρ²`-based face kernels, shared by the near and far fields.
#[inline]
fn hi_lo(d: [f64; 3], qa: f64, qb: f64, neg_s: f64) -> [f64; 2] {
    let r2 = dot(d, d);
    let hi = r2.powf(neg_s);
    [hi, 0.5 * (qa + qb) / (r2 * r2) * hi]
}

/// `h_hi`, `h_lo` on the face tree (build time only).
struct FaceKernels<'a> {
    tree: &'a FaceTree,
    ops: &'a FaceOperators,
    chi: f64,
    /// `−s`, applied to `ρ²`.
    neg_s: f64,
}

impl ClusterKernel<2> for FaceKernels<'_> {
    fn admissible(&self, a: usize, b: usize) -> bool {
        let (ca, cb) = (&self.tree.clusters()[a], &self.tree.clusters()[b]);
        separated(ca.radius, &ca.bbox, cb.radius, &cb.bbox, self.chi)
    }

    /// `[ρ^{−2s}, ½(vᵀP̄_I v + vᵀP̄_J v) ρ^{−4} ρ^{−2s}]` at the cluster centroids, evaluated for
    /// `(min, max)` so that it is bitwise symmetric.
    ///
    /// [formula: rs-projector-average]
    fn far(&self, a: usize, b: usize) -> [f64; 2] {
        let (a, b) = if a <= b { (a, b) } else { (b, a) };
        let (ca, cb) = (&self.tree.clusters()[a], &self.tree.clusters()[b]);
        let d = sub(ca.centroid, cb.centroid);
        hi_lo(
            d,
            ca.projector_quadratic_form(d),
            cb.projector_quadratic_form(d),
            self.neg_s,
        )
    }

    /// The dense kernels of [`super::dense::high_order_kernel`] and
    /// [`super::dense::low_order_kernel`] (same expressions, `(min, max)` order).
    fn near(&self, i: usize, j: usize) -> Option<[f64; 2]> {
        let (a, b) = if i <= j { (i, j) } else { (j, i) };
        let ops = self.ops;
        let d = sub(ops.barycenters[a], ops.barycenters[b]);
        let ns = dot(ops.normals[a], d);
        let nt = dot(ops.normals[b], d);
        Some(hi_lo(d, ns * ns, nt * nt, self.neg_s))
    }

    fn weight(&self, _i: usize) -> f64 {
        1.0
    }
}

/// `h_fr = ρ^{2s−6}` on the face tree.
struct FaceFractionalKernel<'a> {
    tree: &'a FaceTree,
    ops: &'a FaceOperators,
    chi: f64,
    /// `−(τ + 1)` for `L_τ`, applied to `ρ²`.
    exponent: f64,
}

impl ClusterKernel<1> for FaceFractionalKernel<'_> {
    fn admissible(&self, a: usize, b: usize) -> bool {
        let (ca, cb) = (&self.tree.clusters()[a], &self.tree.clusters()[b]);
        separated(ca.radius, &ca.bbox, cb.radius, &cb.bbox, self.chi)
    }

    fn far(&self, a: usize, b: usize) -> [f64; 1] {
        let (a, b) = if a <= b { (a, b) } else { (b, a) };
        let (ca, cb) = (&self.tree.clusters()[a], &self.tree.clusters()[b]);
        let d = sub(ca.centroid, cb.centroid);
        [dot(d, d).powf(self.exponent)]
    }

    fn near(&self, i: usize, j: usize) -> Option<[f64; 1]> {
        let (a, b) = if i <= j { (i, j) } else { (j, i) };
        let d = sub(self.ops.barycenters[a], self.ops.barycenters[b]);
        Some([dot(d, d).powf(self.exponent)])
    }

    fn weight(&self, _i: usize) -> f64 {
        1.0
    }
}

/// Which kernels a face product uses, with weights and the obstacle diagonal (§4.4).
struct FaceTerms<'a> {
    w_hi: f64,
    w_lo: f64,
    d_hi: &'a [f64],
    d_lo: &'a [f64],
}

/// The hierarchical `B` and `B⁰` of one configuration (`SPEC.md` §8.3): the face tree, its
/// block cluster tree with the kernels `h_hi`, `h_lo`, and `H̃ a` for both. Built once per step
/// (before the line search); a product costs no kernel evaluation.
#[derive(Debug, Clone)]
pub struct FaceMetricTree {
    face_tree: FaceTree,
    ops: FaceOperators,
    bct: BlockTree<2>,
    /// `H̃ a` per face for `h_hi`, `h_lo`.
    ha: Vec<[f64; 2]>,
    /// `H̃ a` in BVH order.
    ha_perm: Vec<[f64; 2]>,
    incidence: Incidence,
    options: ClusterTreeOptions,
    build_seconds: f64,
}

impl FaceMetricTree {
    /// Builds the face tree ([`FaceTree::build`]) and the block cluster tree at `positions`.
    ///
    /// # Errors
    /// The errors of [`Surface::validate_positions`]; [`RsError::InvalidOption`] for invalid
    /// options; [`RsError::NonFiniteEnergy`] for coincident barycenters.
    pub fn build(
        pool: &rayon::ThreadPool,
        surface: &Surface,
        positions: &[[f64; 3]],
        exponents: &Exponents,
        options: &ClusterTreeOptions,
    ) -> Result<FaceMetricTree, RsError> {
        options.validate()?;
        let geom = FaceGeometry::new(surface, positions)?;
        let tree = FaceTree::build(surface, positions, &geom, options.leaf_size)?;
        let ops = FaceOperators::from_geometry(surface, positions, geom);
        Self::from_parts(pool, tree, ops, exponents, options)
    }

    /// As [`FaceMetricTree::build`] from an existing face tree (e.g. the Barnes–Hut tree's
    /// [`crate::energy::barnes_hut::BhTree::face_tree`], cloned) and the face operators of the
    /// **same** positions. `options.leaf_size` is not used (the tree's own leaves are).
    ///
    /// # Errors
    /// [`RsError::InvalidOption`] for invalid options or a tree of another face count;
    /// [`RsError::NonFiniteEnergy`] for coincident barycenters.
    pub fn from_parts(
        pool: &rayon::ThreadPool,
        face_tree: FaceTree,
        ops: FaceOperators,
        exponents: &Exponents,
        options: &ClusterTreeOptions,
    ) -> Result<FaceMetricTree, RsError> {
        let t0 = Instant::now();
        options.validate()?;
        if face_tree.n_faces() != ops.n_faces() {
            return Err(RsError::InvalidOption {
                name: "face tree",
                value: format!("{} faces", face_tree.n_faces()),
                requirement: "the face count of the face operators",
            });
        }
        let kernel = FaceKernels {
            tree: &face_tree,
            ops: &ops,
            chi: options.chi,
            neg_s: -exponents.s(),
        };
        let bct = BlockTree::<2>::build(pool, face_tree.bvh(), &kernel, options.leaf_pairs)?;
        let psi: Vec<[f64; 2]> = ops.areas.iter().map(|&a| [a, a]).collect();
        let ha = bct.multiply(pool, &[0, 1], &psi);
        if let Some(s) = ha.iter().position(|h| !(h[0].is_finite() && h[1].is_finite())) {
            return Err(coincident_faces(&ops, s));
        }
        let ha_perm: Vec<[f64; 2]> = bct.order().iter().map(|&s| ha[s as usize]).collect();
        let incidence = Incidence::new(&ops.faces, ops.n_vertices);
        Ok(FaceMetricTree {
            face_tree,
            ops,
            bct,
            ha,
            ha_perm,
            incidence,
            options: *options,
            build_seconds: t0.elapsed().as_secs_f64(),
        })
    }

    /// The face tree.
    pub fn face_tree(&self) -> &FaceTree {
        &self.face_tree
    }

    /// The face operators.
    pub fn operators(&self) -> &FaceOperators {
        &self.ops
    }

    /// The counts of the block cluster tree.
    pub fn stats(&self) -> BctStats {
        self.bct.stats()
    }

    /// The options used.
    pub fn options(&self) -> &ClusterTreeOptions {
        &self.options
    }

    /// Seconds of the block-tree build (classification, coefficients, near field, `H̃ a`).
    pub fn build_seconds(&self) -> f64 {
        self.build_seconds
    }

    /// `|V|`.
    pub fn n_vertices(&self) -> usize {
        self.ops.n_vertices
    }

    /// `H̃ a` per face for `[h_hi, h_lo]`.
    pub fn kernel_times_areas(&self) -> &[[f64; 2]] {
        &self.ha
    }

    /// `y = (w_hi B̃ + w_lo B̃⁰ + O_diag) x` for an interleaved `3|V|` vector (`y` overwritten),
    /// where `O_diag` is the obstacle term `2 Vᵀ diag(hi/a) V + 2 Uᵀ diag(lo/a) U` of §4.4
    /// ([`FaceDiagonal`], unweighted; empty parts are skipped). Terms with zero weight and no
    /// diagonal are not evaluated.
    ///
    /// [formula: rs-kernel-matrix-form]
    /// [formula: rs-hierarchical-null-space]
    ///
    /// # Errors
    /// [`RsError::InvalidOption`] for wrong lengths (also of a non-empty diagonal) or
    /// non-finite weights.
    pub fn apply(
        &self,
        pool: &rayon::ThreadPool,
        w_hi: f64,
        w_lo: f64,
        diagonal: Option<&FaceDiagonal>,
        x: &[f64],
        y: &mut [f64],
    ) -> Result<(), RsError> {
        let n = 3 * self.ops.n_vertices;
        check_len("hierarchical product input", x.len(), n)?;
        check_len("hierarchical product output", y.len(), n)?;
        if !(w_hi.is_finite() && w_lo.is_finite()) {
            return Err(RsError::InvalidOption {
                name: "metric weights",
                value: format!("({w_hi}, {w_lo})"),
                requirement: "finite weights",
            });
        }
        let empty: &[f64] = &[];
        let (d_hi, d_lo) = match diagonal {
            Some(d) => (d.hi.as_slice(), d.lo.as_slice()),
            None => (empty, empty),
        };
        let nf = self.ops.n_faces();
        for (name, d) in [("obstacle diagonal hi", d_hi), ("obstacle diagonal lo", d_lo)] {
            if !d.is_empty() && d.len() != nf {
                return Err(RsError::InvalidOption {
                    name,
                    value: format!("{} entries", d.len()),
                    requirement: "one entry per face",
                });
            }
        }
        let terms = FaceTerms {
            w_hi,
            w_lo,
            d_hi,
            d_lo,
        };
        let hi = w_hi != 0.0 || !d_hi.is_empty();
        let lo = w_lo != 0.0 || !d_lo.is_empty();
        match (hi, lo) {
            (true, true) => self.product::<12>(pool, true, true, &terms, x, y),
            (true, false) => self.product::<9>(pool, true, false, &terms, x, y),
            (false, true) => self.product::<3>(pool, false, true, &terms, x, y),
            (false, false) => y.iter_mut().for_each(|v| *v = 0.0),
        }
        Ok(())
    }

    /// `B̃ x` (`w_hi = 1`, `w_lo = 0`).
    ///
    /// # Errors
    /// As [`FaceMetricTree::apply`].
    pub fn apply_high(&self, pool: &rayon::ThreadPool, x: &[f64], y: &mut [f64]) -> Result<(), RsError> {
        self.apply(pool, 1.0, 0.0, None, x, y)
    }

    /// `B̃⁰ x` (`w_hi = 0`, `w_lo = 1`).
    ///
    /// # Errors
    /// As [`FaceMetricTree::apply`].
    pub fn apply_low(&self, pool: &rayon::ThreadPool, x: &[f64], y: &mut [f64]) -> Result<(), RsError> {
        self.apply(pool, 0.0, 1.0, None, x, y)
    }

    /// The face product with `C` channels: `9` gradient channels (`3k + component` for axis `k`)
    /// if `hi`, then `3` average channels if `lo`.
    fn product<const C: usize>(
        &self,
        pool: &rayon::ThreadPool,
        hi: bool,
        lo: bool,
        terms: &FaceTerms<'_>,
        x: &[f64],
        y: &mut [f64],
    ) {
        debug_assert_eq!(C, 9 * hi as usize + 3 * lo as usize);
        let ops = &self.ops;
        let nf = ops.n_faces();
        let off = if hi { 9 } else { 0 };
        let mut map = [0usize; C];
        for (c, m) in map.iter_mut().enumerate() {
            *m = if c < off { 0 } else { 1 };
        }
        // m = M x per face (D_f per axis, face averages per axis), psi = a m; in BVH order.
        let order = self.bct.order();
        let m: Vec<[f64; C]> = pool.install(|| {
            order
                .par_iter()
                .with_min_len(256)
                .map(|&s| {
                    let s = s as usize;
                    let t = ops.faces[s];
                    let mut out = [0.0; C];
                    if hi {
                        let g = &ops.df[s];
                        for k in 0..3 {
                            for c in 0..3 {
                                let xc = x[3 * t[c] as usize + k];
                                for r in 0..3 {
                                    out[3 * k + r] += xc * g[c][r];
                                }
                            }
                        }
                    }
                    if lo {
                        for k in 0..3 {
                            out[off + k] = (x[3 * t[0] as usize + k]
                                + x[3 * t[1] as usize + k]
                                + x[3 * t[2] as usize + k])
                                / 3.0;
                        }
                    }
                    out
                })
                .collect()
        });
        let psi: Vec<[f64; C]> = m
            .iter()
            .zip(order)
            .map(|(mm, &s)| {
                let a = ops.areas[s as usize];
                let mut p = *mm;
                p.iter_mut().for_each(|v| *v *= a);
                p
            })
            .collect();
        let hx = self.bct.multiply_permuted(pool, &map, &psi);
        // z = w (diag(H̃a) m − H̃ (a m)) + d m, in BVH order.
        let z: Vec<[f64; C]> = pool.install(|| {
            (0..nf)
                .into_par_iter()
                .with_min_len(256)
                .map(|p| {
                    let s = order[p] as usize;
                    let mut out = [0.0; C];
                    for c in 0..C {
                        let (w, kern, d) = if c < off {
                            (terms.w_hi, 0, terms.d_hi.get(s).copied().unwrap_or(0.0))
                        } else {
                            (terms.w_lo, 1, terms.d_lo.get(s).copied().unwrap_or(0.0))
                        };
                        out[c] = w * (self.ha_perm[p][kern] * m[p][c] - hx[p][c]) + d * m[p][c];
                    }
                    out
                })
                .collect()
        });
        let position = self.bct.position();
        // y_v = 2 Σ_{S∋v} a_S (Mᵀ row of S at v) · z_S, faces ascending.
        let inc = &self.incidence;
        pool.install(|| {
            y.par_chunks_mut(3)
                .with_min_len(256)
                .enumerate()
                .for_each(|(v, yv)| {
                    let mut acc = [0.0; 3];
                    for &(s, c) in inc.of(v) {
                        let (s, c) = (s as usize, c as usize);
                        let a = ops.areas[s];
                        let zs = &z[position[s] as usize];
                        if hi {
                            let g = ops.df[s][c];
                            for k in 0..3 {
                                acc[k] += a
                                    * (g[0] * zs[3 * k] + g[1] * zs[3 * k + 1] + g[2] * zs[3 * k + 2]);
                            }
                        }
                        if lo {
                            for k in 0..3 {
                                acc[k] += a / 3.0 * zs[off + k];
                            }
                        }
                    }
                    for k in 0..3 {
                        yv[k] = 2.0 * acc[k];
                    }
                });
        });
    }
}

/// The face form of the fractional operator `L_τ` (`SPEC.md` §4.2, §6.3; `τ = 2 − s` in the
/// preconditioner, `FractionalForm::FacePairs`): kernel `ρ^{−(2τ+2)}` on face pairs, through the face
/// averages, hierarchical on a face tree.
#[derive(Debug, Clone)]
pub struct FaceFractionalTree {
    bct: BlockTree<1>,
    faces: Vec<[u32; 3]>,
    areas: Vec<f64>,
    ha: Vec<f64>,
    incidence: Incidence,
    n_vertices: usize,
    build_seconds: f64,
}

impl FaceFractionalTree {
    /// Builds the block cluster tree of `L_τ` on `face_tree` (built at the positions of `ops`).
    ///
    /// # Errors
    /// [`RsError::InvalidOption`] for invalid options, a non-finite `τ` or a tree of another
    /// face count; [`RsError::NonFiniteEnergy`] for coincident barycenters.
    pub fn build(
        pool: &rayon::ThreadPool,
        face_tree: &FaceTree,
        ops: &FaceOperators,
        tau: f64,
        options: &ClusterTreeOptions,
    ) -> Result<FaceFractionalTree, RsError> {
        let t0 = Instant::now();
        options.validate()?;
        if !tau.is_finite() {
            return Err(RsError::InvalidOption {
                name: "tau",
                value: tau.to_string(),
                requirement: "a finite order",
            });
        }
        if face_tree.n_faces() != ops.n_faces() {
            return Err(RsError::InvalidOption {
                name: "face tree",
                value: format!("{} faces", face_tree.n_faces()),
                requirement: "the face count of the face operators",
            });
        }
        let kernel = FaceFractionalKernel {
            tree: face_tree,
            ops,
            chi: options.chi,
            exponent: -(tau + 1.0),
        };
        let bct = BlockTree::<1>::build(pool, face_tree.bvh(), &kernel, options.leaf_pairs)?;
        let psi: Vec<[f64; 1]> = ops.areas.iter().map(|&a| [a]).collect();
        let ha: Vec<f64> = bct.multiply(pool, &[0], &psi).into_iter().map(|h| h[0]).collect();
        if let Some(s) = ha.iter().position(|h| !h.is_finite()) {
            return Err(coincident_faces(ops, s));
        }
        Ok(FaceFractionalTree {
            bct,
            faces: ops.faces.clone(),
            areas: ops.areas.clone(),
            ha,
            incidence: Incidence::new(&ops.faces, ops.n_vertices),
            n_vertices: ops.n_vertices,
            build_seconds: t0.elapsed().as_secs_f64(),
        })
    }

    /// The counts of the block cluster tree.
    pub fn stats(&self) -> BctStats {
        self.bct.stats()
    }

    /// Seconds of the build.
    pub fn build_seconds(&self) -> f64 {
        self.build_seconds
    }

    /// `y = L̃_τ x` for an interleaved `3|V|` vector (`y` overwritten).
    ///
    /// [formula: rs-kernel-matrix-form]
    /// [formula: rs-hierarchical-null-space]
    ///
    /// # Errors
    /// [`RsError::InvalidOption`] for wrong lengths.
    pub fn apply(&self, pool: &rayon::ThreadPool, x: &[f64], y: &mut [f64]) -> Result<(), RsError> {
        let n = 3 * self.n_vertices;
        check_len("fractional product input", x.len(), n)?;
        check_len("fractional product output", y.len(), n)?;
        let nf = self.faces.len();
        let m: Vec<[f64; 3]> = (0..nf)
            .map(|s| {
                let t = self.faces[s];
                std::array::from_fn(|k| {
                    (x[3 * t[0] as usize + k] + x[3 * t[1] as usize + k] + x[3 * t[2] as usize + k])
                        / 3.0
                })
            })
            .collect();
        let psi: Vec<[f64; 3]> = m
            .iter()
            .zip(&self.areas)
            .map(|(mm, &a)| [a * mm[0], a * mm[1], a * mm[2]])
            .collect();
        let hx = self.bct.multiply(pool, &[0, 0, 0], &psi);
        let z: Vec<[f64; 3]> = (0..nf)
            .map(|s| std::array::from_fn(|k| self.ha[s] * m[s][k] - hx[s][k]))
            .collect();
        let inc = &self.incidence;
        pool.install(|| {
            y.par_chunks_mut(3)
                .with_min_len(256)
                .enumerate()
                .for_each(|(v, yv)| {
                    let mut acc = [0.0; 3];
                    for &(s, _) in inc.of(v) {
                        let s = s as usize;
                        let a = self.areas[s];
                        for k in 0..3 {
                            acc[k] += a / 3.0 * z[s][k];
                        }
                    }
                    for k in 0..3 {
                        yv[k] = 2.0 * acc[k];
                    }
                });
        });
        Ok(())
    }
}

/// A vertex as an element of the vertex BVH: position, dual area.
struct VertexElement {
    p: [f64; 3],
    mass: f64,
}

impl Element<3> for VertexElement {
    fn point(&self) -> [f64; 3] {
        self.p
    }
    fn mass(&self) -> f64 {
        self.mass
    }
}

/// Cluster data of the vertex tree: centroid (dual-area weighted), box of the points, radius.
#[derive(Debug, Clone, Copy)]
struct VertexCluster {
    centroid: [f64; 3],
    bbox: Aabb<3>,
    radius: f64,
}

struct VertexKernel<'a> {
    clusters: &'a [VertexCluster],
    positions: &'a [[f64; 3]],
    chi: f64,
    exponent: f64,
}

impl ClusterKernel<1> for VertexKernel<'_> {
    fn admissible(&self, a: usize, b: usize) -> bool {
        let (ca, cb) = (&self.clusters[a], &self.clusters[b]);
        separated(ca.radius, &ca.bbox, cb.radius, &cb.bbox, self.chi)
    }

    fn far(&self, a: usize, b: usize) -> [f64; 1] {
        let (a, b) = if a <= b { (a, b) } else { (b, a) };
        let d = sub(self.clusters[a].centroid, self.clusters[b].centroid);
        [dot(d, d).powf(self.exponent)]
    }

    /// The entries of [`super::dense::vertex_fractional_operator`] (same expression).
    fn near(&self, i: usize, j: usize) -> Option<[f64; 1]> {
        let (a, b) = if i <= j { (i, j) } else { (j, i) };
        let d = sub(self.positions[a], self.positions[b]);
        Some([dot(d, d).powf(self.exponent)])
    }

    fn weight(&self, _i: usize) -> f64 {
        1.0
    }
}

/// The vertex-pair fractional operator `L^V_τ` of D32 (`SPEC.md` §6.3):
/// `uᵀ L^V_τ u = Σ_{v≠w} (u_v − u_w)² a_v a_w |f_v − f_w|^{−(2τ+2)}` (dual areas, ordered
/// pairs), i.e. `2 (diag(a ⊙ H a) − diag(a) H diag(a))`, hierarchical on a BVH over the vertices
/// (point `f_v`, mass `a_v`, leaf size and `χ` from [`ClusterTreeOptions`]).
#[derive(Debug, Clone)]
pub struct VertexFractionalTree {
    bct: BlockTree<1>,
    mass: Vec<f64>,
    ha: Vec<f64>,
    build_seconds: f64,
}

impl VertexFractionalTree {
    /// Builds the vertex BVH and the block cluster tree of `L^V_τ` for the masses `mass`
    /// (the dual areas; [`VertexFractionalTree::for_surface`] computes them).
    ///
    /// # Errors
    /// [`RsError::InvalidOption`] for invalid options, a non-finite `τ`, mismatched lengths,
    /// non-positive masses or two vertices at the same position.
    pub fn build(
        pool: &rayon::ThreadPool,
        positions: &[[f64; 3]],
        mass: &[f64],
        tau: f64,
        options: &ClusterTreeOptions,
    ) -> Result<VertexFractionalTree, RsError> {
        let t0 = Instant::now();
        options.validate()?;
        let n = positions.len();
        if !tau.is_finite() || mass.len() != n || n == 0 {
            return Err(RsError::InvalidOption {
                name: "vertex fractional operator",
                value: format!("tau {tau}, {} masses for {n} vertices", mass.len()),
                requirement: "a finite order and one mass per vertex (at least one vertex)",
            });
        }
        if let Some(v) = mass.iter().position(|m| !(m.is_finite() && *m > 0.0)) {
            return Err(RsError::InvalidOption {
                name: "mass",
                value: format!("vertex {v}: {:e}", mass[v]),
                requirement: "positive finite vertex areas",
            });
        }
        let elems: Vec<VertexElement> = positions
            .iter()
            .zip(mass)
            .map(|(&p, &m)| VertexElement { p, mass: m })
            .collect();
        let bvh = Bvh::build(&elems, options.leaf_size)?;
        let clusters: Vec<VertexCluster> = bvh
            .nodes()
            .iter()
            .map(|nd| {
                let c = nd.centroid();
                VertexCluster {
                    centroid: c,
                    bbox: nd.bbox,
                    radius: farthest_corner_distance(&nd.bbox, c),
                }
            })
            .collect();
        let kernel = VertexKernel {
            clusters: &clusters,
            positions,
            chi: options.chi,
            exponent: -(tau + 1.0),
        };
        let bct = BlockTree::<1>::build(pool, &bvh, &kernel, options.leaf_pairs)?;
        let psi: Vec<[f64; 1]> = mass.iter().map(|&a| [a]).collect();
        let ha: Vec<f64> = bct.multiply(pool, &[0], &psi).into_iter().map(|h| h[0]).collect();
        if let Some(v) = ha.iter().position(|h| !h.is_finite()) {
            let w = (0..n)
                .find(|&w| w != v && positions[w] == positions[v])
                .unwrap_or(v);
            return Err(RsError::InvalidOption {
                name: "positions",
                value: format!("vertices {} and {} coincide", v.min(w), v.max(w)),
                requirement: "distinct vertex positions for the vertex-pair operator",
            });
        }
        Ok(VertexFractionalTree {
            bct,
            mass: mass.to_vec(),
            ha,
            build_seconds: t0.elapsed().as_secs_f64(),
        })
    }

    /// [`VertexFractionalTree::build`] with the dual areas of `surface` at `positions`.
    ///
    /// # Errors
    /// The errors of [`Surface::validate_positions`] and of [`VertexFractionalTree::build`].
    pub fn for_surface(
        pool: &rayon::ThreadPool,
        surface: &Surface,
        positions: &[[f64; 3]],
        tau: f64,
        options: &ClusterTreeOptions,
    ) -> Result<VertexFractionalTree, RsError> {
        let geom = FaceGeometry::new(surface, positions)?;
        let mass = vertex_areas(surface, &geom.areas);
        Self::build(pool, positions, &mass, tau, options)
    }

    /// The counts of the block cluster tree.
    pub fn stats(&self) -> BctStats {
        self.bct.stats()
    }

    /// Seconds of the build (BVH and block tree).
    pub fn build_seconds(&self) -> f64 {
        self.build_seconds
    }

    /// `y = L̃^V_τ x = 2 (a ⊙ H̃a ⊙ x − a ⊙ H̃(a ⊙ x))` per axis, for an interleaved `3|V|` vector
    /// (`y` overwritten).
    ///
    /// [formula: rs-kernel-matrix-form]
    /// [formula: rs-hierarchical-null-space]
    ///
    /// # Errors
    /// [`RsError::InvalidOption`] for wrong lengths.
    pub fn apply(&self, pool: &rayon::ThreadPool, x: &[f64], y: &mut [f64]) -> Result<(), RsError> {
        let n = 3 * self.mass.len();
        check_len("fractional product input", x.len(), n)?;
        check_len("fractional product output", y.len(), n)?;
        let psi: Vec<[f64; 3]> = x
            .chunks(3)
            .zip(&self.mass)
            .map(|(xv, &a)| [a * xv[0], a * xv[1], a * xv[2]])
            .collect();
        let hx = self.bct.multiply(pool, &[0, 0, 0], &psi);
        pool.install(|| {
            y.par_chunks_mut(3)
                .with_min_len(512)
                .enumerate()
                .for_each(|(v, yv)| {
                    let a = self.mass[v];
                    for k in 0..3 {
                        yv[k] = 2.0 * a * (self.ha[v] * x[3 * v + k] - hx[v][k]);
                    }
                });
        });
        Ok(())
    }
}
