//! Hierarchical (block cluster tree) product with the metric `Ā` (`SPEC.md` §8.3, App. B.2,
//! Eq. 28, 31, 32; D2, D16, D25).
//!
//! [`HierarchicalMetric`] applies
//!
//! ```text
//! Ā u = Dᵀ((diag(K^B 𝟙) − K^B) ⊗ I₃) D u + Eᵀ(diag(K⁰ 𝟙) − K⁰) E u      (per coordinate, Eq. 31)
//! ```
//!
//! with `K^B` and `K⁰` replaced by their hierarchical approximations from one
//! [`BlockClusterTree`] over the edges. The edges are the elements of a [`Bvh<6>`] over their
//! tangent-points `(x_I, T_I)` with mass `ℓ_I` (§8.1). Per node: `L_N`, `x̄_N`, the normalized
//! `T̂_N`, `r_x`, `r_T` and the box of the edge endpoints [formula: rc-bvh-lumped] (the same data
//! as the Barnes–Hut tree of [`crate::energy::barnes_hut`]).
//!
//! * **Admissible** blocks (Eq. 28 plus the guard, D16) [formula: rc-bct-admissibility]:
//!   `𝒜 ≠ ℬ`, the endpoint boxes are disjoint (so no pair of adjacent edges is ever in a
//!   rank-1 block), `max(r_x^𝒜, r_x^ℬ) ≤ θ_x |x̄_𝒜 − x̄_ℬ|`, `max(r_T^𝒜, r_T^ℬ) ≤ θ_T`. `θ = 0`
//!   makes every block inadmissible. The reference uses the spatial criterion only (separation
//!   1.0) and lets admissible blocks contain neighbouring pairs.
//! * **Rank-1 blocks** [formula: rc-bct-lowrank]: `K̂_𝒜ℬ = ℓ[𝒜] κ ℓ[ℬ]ᵀ` with
//!   `κ^B = 2|x̄_𝒜 − x̄_ℬ|^{−(2σ+1)}` and
//!   `κ⁰ = k²_{2σ+5}(x̄_𝒜, x̄_ℬ, T̂_𝒜) + k²_{2σ+5}(x̄_ℬ, x̄_𝒜, T̂_ℬ)` (midpoint rule in the far field
//!   for either quadrature; `|T̂ × d|²` by the cross product, so `κ⁰ ≥ 0` exactly).
//! * **Exact blocks** use [`GramGeometry::pair_kernels`] with the session's quadrature, the same
//!   (bitwise symmetric) entries as the dense matrices of [`super::sobolev_dense`], so the
//!   hierarchical and the dense products agree as `θ → 0` and to rounding at `θ = 0`.
//!
//! One tree pass carries 12 channels: the three components of `D u` for each of the three
//! coordinates (multiplied by `K^B`) and `E u` for each coordinate (multiplied by `K⁰`). In planar
//! mode ([`DofLayout`]) the held coordinate's channels are zero. `K𝟙` is computed once per tree,
//! with the hierarchical `K`, so `Ā 𝟙 = 0` holds to rounding as for the dense matrix.
//!
//! Ported in structure from the authors' MIT code (`src/product/block_cluster_tree.cpp`,
//! `MultiplyVector`, `MultiplyAdmissibleFast`; Copyright (c) 2019 Christopher Yu); the
//! admissibility (D16), the leaf rule and the 12-channel single pass are ours.
//!
//! [`HierarchicalCgSaddle`] is a [`SaddleSolver`] on top of the operator: conjugate gradients on
//! the projected system of Eq. 30 through [`Projector`], with Eq. 33 for constraint right-hand
//! sides. It is a correctness check and a fallback; the multigrid of §8.4 (WP6) is the fast
//! solver.

use std::sync::Arc;
use std::time::Instant;

use ddg_core::sparse::SparseMatrix;

use super::sobolev_dense::GramGeometry;
use super::{DofLayout, GramQuadrature};
use crate::error::RcError;
use crate::kernel::Exponents;
use crate::kernel::v3::{cross, dot, norm, scale, sub};
use crate::network::CurveNetwork;
use crate::solve::{LinearReport, Projector, SaddleSolution, SaddleSolver, check_rhs};
use crate::spatial::aabb::Aabb;
use crate::spatial::bct::{BctOptions, BctStats, BlockClusterTree, ClusterKernel};
use crate::spatial::bvh::{Bvh, Element};

/// Channels of one tree pass: `3 × 3` components of `D u`, then `3` of `E u`.
const CHANNELS: usize = 12;
/// Kernel of each channel: 0 = `K^B`, 1 = `K⁰`.
const CHANNEL_KERNEL: [usize; CHANNELS] = [0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1];

/// An edge as a BVH element: the tangent-point `(x_I, T_I) ∈ R⁶` with mass `ℓ_I`.
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

/// Lumped data of a cluster [formula: rc-bvh-lumped].
#[derive(Debug, Clone, Copy)]
struct ClusterData {
    xbar: [f64; 3],
    that: [f64; 3],
    rx2: f64,
    rt: f64,
    ends: Aabb<3>,
}

/// `K^B`, `K⁰` of Eq. 31–32 on the edges, as a [`ClusterKernel`] (build time only).
struct CurveKernels<'a> {
    geom: &'a GramGeometry,
    edges: &'a [[u32; 2]],
    quadrature: GramQuadrature,
    clusters: Vec<ClusterData>,
    theta_x2: f64,
    theta_t: f64,
    /// `−(2σ+1)/2`, applied to squared distances.
    half_exponent: f64,
}

fn boxes_disjoint(a: &Aabb<3>, b: &Aabb<3>) -> bool {
    (0..3).any(|k| a.hi[k] < b.lo[k] || b.hi[k] < a.lo[k])
}

impl ClusterKernel<2> for CurveKernels<'_> {
    /// Eq. 28 with the endpoint-box guard (D16).
    ///
    /// [formula: rc-bct-admissibility]
    fn admissible(&self, a: usize, b: usize) -> bool {
        if self.theta_x2 <= 0.0 || self.theta_t <= 0.0 {
            return false;
        }
        let (ca, cb) = (&self.clusters[a], &self.clusters[b]);
        if ca.rt.max(cb.rt) > self.theta_t {
            return false;
        }
        let d = sub(ca.xbar, cb.xbar);
        let d2 = dot(d, d);
        if d2.is_nan() || d2 <= 0.0 || ca.rx2.max(cb.rx2) > self.theta_x2 * d2 {
            return false;
        }
        boxes_disjoint(&ca.ends, &cb.ends)
    }

    /// `[κ^B, κ⁰]`, computed for `(min(a, b), max(a, b))` so that it is bitwise symmetric.
    ///
    /// [formula: rc-bct-lowrank]
    fn far(&self, a: usize, b: usize) -> [f64; 2] {
        let (a, b) = if a <= b { (a, b) } else { (b, a) };
        let (ca, cb) = (&self.clusters[a], &self.clusters[b]);
        let d = sub(ca.xbar, cb.xbar);
        let r2 = dot(d, d);
        let p1 = r2.powf(self.half_exponent);
        let p5 = p1 / (r2 * r2);
        // |T̂ × d|² with the cross product: never negative (no cancellation when T̂ ∥ d).
        let (ca2, cb2) = (cross(ca.that, d), cross(cb.that, d));
        [2.0 * p1, (dot(ca2, ca2) + dot(cb2, cb2)) * p5]
    }

    fn near(&self, i: usize, j: usize) -> Option<[f64; 2]> {
        let [a, b] = self.edges[i];
        let [c, d] = self.edges[j];
        if a == c || a == d || b == c || b == d {
            return None;
        }
        Some(self.geom.pair_kernels(i, j, self.quadrature))
    }

    fn weight(&self, i: usize) -> f64 {
        self.geom.length[i]
    }
}

/// Cheap geometric checks: count, finite coordinates, zero-length edges.
fn check_positions(network: &CurveNetwork, positions: &[[f64; 3]]) -> Result<(), RcError> {
    if positions.len() != network.n_vertices() {
        return Err(RcError::PositionCountMismatch {
            expected: network.n_vertices(),
            found: positions.len(),
        });
    }
    for (v, p) in positions.iter().enumerate() {
        if !p.iter().all(|c| c.is_finite()) {
            return Err(RcError::NonFinitePosition { vertex: v as u32 });
        }
    }
    for (i, &[a, b]) in network.edges().iter().enumerate() {
        if positions[a as usize] == positions[b as usize] {
            return Err(RcError::ZeroLengthEdge { edge: i, a, b });
        }
    }
    Ok(())
}

/// The metric `Ā` of one curve, applied through a block cluster tree (`SPEC.md` §8.3).
#[derive(Debug, Clone)]
pub struct HierarchicalMetric {
    pool: Arc<rayon::ThreadPool>,
    layout: DofLayout,
    n_vertices: usize,
    n_components: usize,
    edges: Vec<[u32; 2]>,
    /// `T_I / ℓ_I`.
    t_over_l: Vec<[f64; 3]>,
    bvh: Bvh<6>,
    bct: BlockClusterTree<2>,
    options: BctOptions,
    quadrature: GramQuadrature,
    build_seconds: f64,
}

impl HierarchicalMetric {
    /// Builds the BVH and the block cluster tree of `positions` on `pool` and precomputes the
    /// far-field coefficients, the exact near-field entries and `K𝟙`. Bitwise independent of
    /// the pool size (D25). Cheap geometric checks only; the full validation (intersections) is
    /// the caller's ([`CurveNetwork::validate_positions`]).
    ///
    /// # Errors
    /// [`RcError::InvalidOption`] for bad options ([`BctOptions::validate`]);
    /// [`RcError::PositionCountMismatch`], [`RcError::NonFinitePosition`],
    /// [`RcError::ZeroLengthEdge`].
    pub fn build(
        pool: Arc<rayon::ThreadPool>,
        network: &CurveNetwork,
        positions: &[[f64; 3]],
        exponents: Exponents,
        quadrature: GramQuadrature,
        layout: DofLayout,
        options: BctOptions,
    ) -> Result<HierarchicalMetric, RcError> {
        let t0 = Instant::now();
        options.validate()?;
        check_positions(network, positions)?;
        let geom = GramGeometry::new(network, positions, exponents);
        let m = geom.n_edges();
        let elems: Vec<TangentPoint> = (0..m)
            .map(|i| {
                let (x, t) = (geom.midpoint[i], geom.tangent[i]);
                TangentPoint {
                    p: [x[0], x[1], x[2], t[0], t[1], t[2]],
                    mass: geom.length[i],
                }
            })
            .collect();
        let bvh = Bvh::build(&elems, options.leaf_size)?;
        let ends = bvh.fold_up(
            |idx| {
                let mut b = Aabb::<3>::empty();
                for &e in idx {
                    for p in geom.ends[e as usize] {
                        b.grow(p);
                    }
                }
                b
            },
            |a, b| a.union(b),
        );
        let clusters = bvh
            .nodes()
            .iter()
            .zip(ends)
            .map(|(nd, ends)| {
                let c = nd.centroid();
                let tsum = [nd.moment[3], nd.moment[4], nd.moment[5]];
                let tn = norm(tsum);
                let rx = nd.bbox.half_diagonal_of(0..3);
                ClusterData {
                    xbar: [c[0], c[1], c[2]],
                    // Cancelling tangents: r_T is then large and the node is never admissible.
                    that: if tn > 0.0 {
                        scale(tsum, 1.0 / tn)
                    } else {
                        [0.0; 3]
                    },
                    rx2: rx * rx,
                    rt: nd.bbox.half_diagonal_of(3..6),
                    ends,
                }
            })
            .collect();
        let kernels = CurveKernels {
            geom: &geom,
            edges: network.edges(),
            quadrature,
            clusters,
            theta_x2: options.theta_x * options.theta_x,
            theta_t: options.theta_t,
            half_exponent: -(2.0 * exponents.sigma() + 1.0) / 2.0,
        };
        let bct = BlockClusterTree::build(&pool, &bvh, &kernels, options.leaf_pairs)?;
        let t_over_l = (0..m)
            .map(|i| scale(geom.tangent[i], 1.0 / geom.length[i]))
            .collect();
        Ok(HierarchicalMetric {
            pool,
            layout,
            n_vertices: network.n_vertices(),
            n_components: network.components().len(),
            edges: network.edges().to_vec(),
            t_over_l,
            bvh,
            bct,
            options,
            quadrature,
            build_seconds: t0.elapsed().as_secs_f64(),
        })
    }

    /// `d|V|`.
    pub fn n_dofs(&self) -> usize {
        self.layout.n_dofs(self.n_vertices)
    }

    /// The layout of the vectors `apply` takes and returns.
    pub fn layout(&self) -> DofLayout {
        self.layout
    }

    /// The options of the tree.
    pub fn options(&self) -> &BctOptions {
        &self.options
    }

    /// The quadrature of the exact blocks.
    pub fn quadrature(&self) -> GramQuadrature {
        self.quadrature
    }

    /// Counts of the block cluster tree.
    pub fn stats(&self) -> BctStats {
        self.bct.stats()
    }

    /// Wall time of [`HierarchicalMetric::build`].
    pub fn build_seconds(&self) -> f64 {
        self.build_seconds
    }

    /// The edge BVH.
    pub fn bvh(&self) -> &Bvh<6> {
        &self.bvh
    }

    /// The block cluster tree of `[K^B, K⁰]` (kernel 0 is `K^B`, kernel 1 is `K⁰`).
    pub fn tree(&self) -> &BlockClusterTree<2> {
        &self.bct
    }

    /// The thread pool the products run on.
    pub fn pool(&self) -> &Arc<rayon::ThreadPool> {
        &self.pool
    }

    /// `Ā x` for `x` in the interleaved layout (D29), parallel and bitwise independent of the
    /// thread count (D25).
    ///
    /// [formula: rc-kernel-matrix-form]
    ///
    /// # Errors
    /// [`RcError::InvalidOption`] if `x.len() != d|V|`.
    pub fn apply(&self, x: &[f64]) -> Result<Vec<f64>, RcError> {
        if x.len() != self.n_dofs() {
            return Err(RcError::InvalidOption {
                name: "hierarchical metric argument",
                value: x.len().to_string(),
                requirement: "d|V| entries",
            });
        }
        Ok(self.apply_unchecked(x))
    }

    /// [`HierarchicalMetric::apply`] without the length check.
    pub(crate) fn apply_unchecked(&self, x: &[f64]) -> Vec<f64> {
        let d = self.layout.dim();
        // ψ_I: (D x_s)_I for each slot s (3 channels each), then (E x_s)_I.
        let psi: Vec<[f64; CHANNELS]> = self
            .edges
            .iter()
            .zip(&self.t_over_l)
            .map(|(&[a, b], t)| {
                let (a, b) = (a as usize * d, b as usize * d);
                let mut p = [0.0; CHANNELS];
                for s in 0..d {
                    let du = x[b + s] - x[a + s];
                    for k in 0..3 {
                        p[3 * s + k] = du * t[k];
                    }
                    p[9 + s] = 0.5 * (x[a + s] + x[b + s]);
                }
                p
            })
            .collect();
        let phi = self
            .bct
            .multiply_unchecked(&self.pool, &CHANNEL_KERNEL, &psi);
        let rs = self.bct.row_sums();
        let mut out = vec![0.0; self.n_dofs()];
        for (i, &[a, b]) in self.edges.iter().enumerate() {
            let (a, b) = (a as usize * d, b as usize * d);
            let t = self.t_over_l[i];
            let (p, f, r) = (&psi[i], &phi[i], &rs[i]);
            for s in 0..d {
                let mut v = 0.0;
                for k in 0..3 {
                    let c = 3 * s + k;
                    v += t[k] * (r[0] * p[c] - f[c]);
                }
                let e = 0.5 * (r[1] * p[9 + s] - f[9 + s]);
                out[b + s] += v + e;
                out[a + s] += e - v;
            }
        }
        out
    }
}

/// Options of [`HierarchicalCgSaddle`].
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CgOptions {
    /// Stop when `‖r_k‖ ≤ tolerance · ‖r_0‖` on the projected system (Eq. 30).
    pub tolerance: f64,
    /// Maximal number of CG iterations per solve.
    pub max_iterations: usize,
}

impl Default for CgOptions {
    /// `tolerance 1e-10`, `max_iterations 20000`.
    fn default() -> Self {
        CgOptions {
            tolerance: 1e-10,
            max_iterations: 20_000,
        }
    }
}

/// A [`SaddleSolver`] for `[Ā Cᵀ; C 0][x; λ] = [a; b]` with the hierarchical `Ā`: conjugate
/// gradients on `P Ā P y = P r` (Eq. 30, `x = P y`) through [`Projector`], and Eq. 33 for
/// `b ≠ 0`; then `λ = (C Cᵀ)⁻¹ C (a − Ā x)`. Unpreconditioned, so the iteration count grows
/// with the size; the report carries the iterations and the residual of the full saddle
/// system. A solve that reaches `max_iterations` is not an error: check
/// `report.relative_residual`.
///
/// [formula: rc-projected-system]
/// [formula: rc-projection-transform]
#[derive(Debug)]
pub struct HierarchicalCgSaddle {
    metric: HierarchicalMetric,
    projector: Projector,
    c: SparseMatrix,
    ct: SparseMatrix,
    options: CgOptions,
}

impl HierarchicalCgSaddle {
    /// The solver for `metric` and the constraint Jacobian `c` (`k × d|V|`).
    ///
    /// # Errors
    /// [`RcError::InvalidOption`] for a Jacobian of the wrong width or a non-positive or
    /// non-finite tolerance; the errors of [`Projector::new`].
    pub fn new(
        metric: HierarchicalMetric,
        c: &SparseMatrix,
        options: CgOptions,
    ) -> Result<HierarchicalCgSaddle, RcError> {
        if c.ncols() != metric.n_dofs() {
            return Err(RcError::InvalidOption {
                name: "constraint Jacobian",
                value: format!("{} columns", c.ncols()),
                requirement: "d|V| columns, matching the metric",
            });
        }
        if !(options.tolerance.is_finite() && options.tolerance > 0.0) {
            return Err(RcError::InvalidOption {
                name: "cg.tolerance",
                value: options.tolerance.to_string(),
                requirement: "a finite positive tolerance",
            });
        }
        let projector = Projector::new(c, None)?;
        Ok(HierarchicalCgSaddle {
            metric,
            projector,
            c: c.clone(),
            ct: c.transpose(),
            options,
        })
    }

    /// The operator.
    pub fn metric(&self) -> &HierarchicalMetric {
        &self.metric
    }

    /// The projector of `C`.
    pub fn projector(&self) -> &Projector {
        &self.projector
    }

    /// CG on `P Ā P y = P r` from `y = 0`; returns `P y` and the iteration count.
    fn solve_kernel(&self, r: &[f64]) -> Result<(Vec<f64>, usize), RcError> {
        let p_op = |v: &[f64]| {
            self.projector
                .projected_operator(|w| self.metric.apply_unchecked(w), v)
        };
        let dotv = |a: &[f64], b: &[f64]| a.iter().zip(b).map(|(x, y)| x * y).sum::<f64>();
        let mut res = self.projector.project(r)?;
        let r0 = dotv(&res, &res).sqrt();
        let mut y = vec![0.0; r.len()];
        if r0 == 0.0 {
            return Ok((y, 0));
        }
        let mut p = res.clone();
        let mut rr = r0 * r0;
        let mut it = 0;
        while it < self.options.max_iterations && rr.sqrt() > self.options.tolerance * r0 {
            let q = p_op(&p)?;
            let pq = dotv(&p, &q);
            if pq.is_nan() || pq <= 0.0 {
                return Err(RcError::SingularMetric {
                    components: self.metric.n_components,
                });
            }
            let alpha = rr / pq;
            for k in 0..y.len() {
                y[k] += alpha * p[k];
                res[k] -= alpha * q[k];
            }
            let rr_new = dotv(&res, &res);
            let beta = rr_new / rr;
            rr = rr_new;
            for k in 0..p.len() {
                p[k] = res[k] + beta * p[k];
            }
            it += 1;
        }
        Ok((self.projector.project(&y)?, it))
    }
}

impl SaddleSolver for HierarchicalCgSaddle {
    fn n_dofs(&self) -> usize {
        self.metric.n_dofs()
    }

    fn n_constraints(&self) -> usize {
        self.c.nrows()
    }

    fn solve(&mut self, a: &[f64], b: &[f64]) -> Result<SaddleSolution, RcError> {
        check_rhs(&*self, a, b)?;
        let t0 = Instant::now();
        let mut iterations = 0;
        let x = self.projector.solve_with_kernel_solver(
            |v| self.metric.apply_unchecked(v),
            |r| {
                let (y, it) = self.solve_kernel(r)?;
                iterations += it;
                Ok(y)
            },
            a,
            b,
        )?;
        let ax = self.metric.apply_unchecked(&x);
        let resid: Vec<f64> = a.iter().zip(&ax).map(|(p, q)| p - q).collect();
        let mut cr = vec![0.0; self.c.nrows()];
        self.c
            .mul_vec(&resid, &mut cr)
            .map_err(|e| crate::solve::sparse_error(e, "C"))?;
        let lambda = self.projector.solve_gram(&cr)?;
        // Residual of the full saddle system.
        let mut ctl = vec![0.0; a.len()];
        self.ct
            .mul_vec(&lambda, &mut ctl)
            .map_err(|e| crate::solve::sparse_error(e, "C^T"))?;
        let mut cx = vec![0.0; b.len()];
        self.c
            .mul_vec(&x, &mut cx)
            .map_err(|e| crate::solve::sparse_error(e, "C"))?;
        let top: f64 = resid.iter().zip(&ctl).map(|(r, c)| (c - r) * (c - r)).sum();
        let bot: f64 = cx.iter().zip(b).map(|(p, q)| (p - q) * (p - q)).sum();
        let rhs: f64 = a.iter().chain(b).map(|v| v * v).sum();
        let res = (top + bot).sqrt();
        let relative_residual = if rhs > 0.0 { res / rhs.sqrt() } else { res };
        Ok(SaddleSolution {
            x,
            lambda,
            report: LinearReport {
                solver: "projected CG (hierarchical)",
                dimension: a.len() + b.len(),
                iterations,
                relative_residual,
                fallback_used: false,
                levels: 0,
                pivot_ratio: None,
                seconds: t0.elapsed().as_secs_f64(),
            },
        })
    }
}
