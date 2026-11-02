//! The saddle solvers of one step (`SPEC.md` §6.1, §6.3, §4.5).
//!
//! * The Sobolev metric: [`Backend::Dense`] (dense `A`, LU) or the hierarchical back end (GMRES on
//!   the block-cluster-tree product with the sparse preconditioner), chosen by
//!   [`Backend::resolve`]; the hierarchical back end reuses the Barnes–Hut face tree when the leaf
//!   sizes agree.
//! * The comparison metrics `L2`, `H1`, `H2` (§4.5): [`SparseSaddle`], a sparse `LDLᵀ` of
//!   `[[K ⊗ I₃, C_Qᵀ], [C_Q, 0]]` with iterative refinement against the exact matrix.

// `!(x > tol)` is deliberate throughout: a NaN must fail the test.
#![allow(clippy::neg_cmp_op_on_partial_ord)]

use std::sync::Arc;
use std::time::Instant;

use ddg_core::SparseMatrix;

use crate::error::RsError;
use crate::geometry::{FaceGeometry, vertex_areas};
use crate::kernel::Exponents;
use crate::metric::dense::comparison_metric;
use crate::metric::{Metric, SobolevWeights};
use crate::solve::dense::DenseSaddle;
use crate::solve::iterative::{
    Backend, HierarchicalOptions, HierarchicalSaddle, MetricExtras, dense_metric_with_extras, hierarchical_saddle,
};
use crate::solve::precond::LaplacianSaddle;
use crate::solve::{LinearReport, SaddleSolution, SaddleSolver, SolveRecord};
use crate::spatial::FaceTree;
use crate::surface::Surface;

/// The relative `ε` of the regularization `K + ε (tr K / tr M) M` of the singular comparison
/// metrics `H1`, `H2` inside the factorization (the refinement solves the exact system).
pub const COMPARISON_EPS_REL: f64 = 1e-10;

/// Iterative-refinement steps of [`SparseSaddle`] at most.
pub const REFINEMENT_STEPS: usize = 20;

/// The refinement stops at this relative residual of the exact saddle system.
pub const REFINEMENT_TOLERANCE: f64 = 1e-13;

/// A direct solver for `A′ = [[K ⊗ I₃, C_Qᵀ], [C_Q, 0]]` with a sparse `K` (the comparison
/// metrics of `SPEC.md` §4.5): `K + εM` (ε = 0 for the positive definite L² metrics,
/// [`COMPARISON_EPS_REL`]` · tr K / tr M` for `H1`, `H2`, whose kernel holds the constants) is
/// factored once by sparse `LDLᵀ` (fill-reducing order, multipliers last), and every solve is
/// refined against the exact `A′` until [`REFINEMENT_TOLERANCE`].
#[derive(Debug)]
pub struct SparseSaddle {
    k: SparseMatrix,
    c: SparseMatrix,
    ct: SparseMatrix,
    factor: LaplacianSaddle,
    n_vertices: usize,
    factor_seconds: f64,
    epsilon: f64,
}

impl SparseSaddle {
    /// Assembles and factors the saddle of a comparison metric (`metric` must not be
    /// `Sobolev`).
    ///
    /// # Errors
    /// The errors of [`comparison_metric`] and [`LaplacianSaddle::new`].
    pub fn new(
        surface: &Surface,
        positions: &[[f64; 3]],
        metric: Metric,
        c: &SparseMatrix,
    ) -> Result<SparseSaddle, RsError> {
        let t0 = Instant::now();
        let k = comparison_metric(surface, positions, metric)?;
        let epsilon = match metric {
            Metric::H1 | Metric::H2 => {
                let geom = FaceGeometry::new(surface, positions)?;
                let mass = vertex_areas(surface, &geom.areas);
                let trk: f64 = k.diagonal().iter().sum();
                let trm: f64 = mass.iter().sum();
                COMPARISON_EPS_REL * trk / trm
            }
            _ => 0.0,
        };
        let kreg = if epsilon > 0.0 {
            let geom = FaceGeometry::new(surface, positions)?;
            let mass = vertex_areas(surface, &geom.areas);
            let n = mass.len();
            let mut trip: Vec<(u32, u32, f64)> = k.triplets().collect();
            trip.extend((0..n).map(|v| (v as u32, v as u32, epsilon * mass[v])));
            SparseMatrix::from_triplets(n, n, &trip).map_err(|e| RsError::InvalidOption {
                name: "comparison metric",
                value: e.to_string(),
                requirement: "a representable sparse matrix",
            })?
        } else {
            k.clone()
        };
        let factor = LaplacianSaddle::new(&kreg, c)?;
        Ok(SparseSaddle {
            n_vertices: k.nrows(),
            ct: c.transpose(),
            c: c.clone(),
            k,
            factor,
            factor_seconds: t0.elapsed().as_secs_f64(),
            epsilon,
        })
    }

    /// The sparse metric `K` (without the regularization), `K ⊗ I₃` being the metric.
    pub fn metric_matrix(&self) -> &SparseMatrix {
        &self.k
    }

    /// The regularization `ε` used in the factorization (0 for L²).
    pub fn epsilon(&self) -> f64 {
        self.epsilon
    }

    /// `y = A′ x` (exact).
    fn apply(&self, x: &[f64], y: &mut [f64]) -> Result<(), RsError> {
        let nd = 3 * self.n_vertices;
        let (xp, xl) = x.split_at(nd);
        let (yp, yl) = y.split_at_mut(nd);
        yp.iter_mut().for_each(|v| *v = 0.0);
        let (cp, ri, val) = (self.k.col_ptr(), self.k.row_idx(), self.k.values());
        for j in 0..self.n_vertices {
            for p in cp[j]..cp[j + 1] {
                let i = ri[p] as usize;
                for a in 0..3 {
                    yp[3 * i + a] += val[p] * xp[3 * j + a];
                }
            }
        }
        let err = |e: ddg_core::SparseError| RsError::InvalidOption {
            name: "constraint rows",
            value: e.to_string(),
            requirement: "matching sizes",
        };
        if !xl.is_empty() {
            let mut ctl = vec![0.0; nd];
            self.ct.mul_vec(xl, &mut ctl).map_err(err)?;
            for (a, b) in yp.iter_mut().zip(&ctl) {
                *a += b;
            }
            self.c.mul_vec(xp, yl).map_err(err)?;
        }
        Ok(())
    }
}

fn norm(x: &[f64]) -> f64 {
    x.iter().map(|v| v * v).sum::<f64>().sqrt()
}

impl SaddleSolver for SparseSaddle {
    fn n_dofs(&self) -> usize {
        3 * self.n_vertices
    }

    fn n_constraints(&self) -> usize {
        self.c.nrows()
    }

    fn solve(&mut self, a: &[f64], b: &[f64]) -> Result<SaddleSolution, RsError> {
        if a.len() != self.n_dofs() || b.len() != self.n_constraints() {
            return Err(RsError::InvalidOption {
                name: "saddle right-hand side",
                value: format!("({}, {})", a.len(), b.len()),
                requirement: "lengths (3|V|, k_Q) of the factored system",
            });
        }
        let rhs: Vec<f64> = a.iter().chain(b).copied().collect();
        let nrhs = norm(&rhs);
        let scale = if nrhs > 0.0 { nrhs } else { 1.0 };
        let residual = |x: &[f64]| -> Result<(Vec<f64>, f64), RsError> {
            let mut r = vec![0.0; rhs.len()];
            self.apply(x, &mut r)?;
            for (ri, bi) in r.iter_mut().zip(&rhs) {
                *ri = bi - *ri;
            }
            let rel = norm(&r) / scale;
            Ok((r, rel))
        };
        let mut x = self.factor.solve(&rhs)?;
        let (mut r, mut rel) = residual(&x)?;
        let mut steps = 0;
        // Refine while it helps, up to the tolerance or the cap.
        while rel > REFINEMENT_TOLERANCE && steps < REFINEMENT_STEPS {
            let dx = self.factor.solve(&r)?;
            let xn: Vec<f64> = x.iter().zip(&dx).map(|(a, b)| a + b).collect();
            let (rn, reln) = residual(&xn)?;
            if !(reln < rel) {
                break;
            }
            x = xn;
            r = rn;
            rel = reln;
            steps += 1;
        }
        if x.iter().any(|v| !v.is_finite()) {
            return Err(RsError::SingularSaddle { pivot: 0 });
        }
        let nd = self.n_dofs();
        Ok(SaddleSolution {
            lambda: x[nd..].to_vec(),
            x: {
                x.truncate(nd);
                x
            },
            record: SolveRecord {
                iterations: steps,
                relative_residual: rel,
                converged: rel <= 1e-10,
                gmres: None,
            },
        })
    }

    fn base_report(&self) -> LinearReport {
        LinearReport {
            solver: "sparse LDLT (comparison metric)",
            dimension: 3 * self.n_vertices + self.c.nrows(),
            schur_pivot: 1.0,
            factor_seconds: self.factor_seconds,
            ..Default::default()
        }
    }
}

/// The solver of one step.
pub(crate) enum StepSolver {
    Dense(DenseSaddle),
    Hierarchical(Box<RobustHierarchical>),
    Sparse(SparseSaddle),
}

/// The largest vertex count for which a failed GMRES solve falls back to a dense LU
/// (`|V|² · 8` bytes for the metric plus the same for the factor: 6000 vertices ≈ 600 MB).
pub const DENSE_FALLBACK_MAX_VERTICES: usize = 6000;

/// Everything needed to rebuild the step's saddle system with other options.
struct FallbackContext {
    pool: Arc<rayon::ThreadPool>,
    surface: Surface,
    positions: Vec<[f64; 3]>,
    exponents: Exponents,
    weights: SobolevWeights,
    extras: MetricExtras,
    c_q: SparseMatrix,
    options: HierarchicalOptions,
}

enum Active {
    Primary,
    Retry(Box<HierarchicalSaddle>),
    Dense(Box<DenseSaddle>),
}

/// The hierarchical saddle solver with the recovery chain of `SPEC.md` §6.3: when GMRES ends
/// above the divergence threshold ([`RsError::SolverDiverged`]) the solve is repeated by
///
/// 1. right GMRES with the augmented preconditioner on vertex pairs (the defaults, replacing a
///    reference-style left/printed/face-pair solver) with at least four times the iteration cap
///    and restart 60, on a rebuilt preconditioner; then
/// 2. a dense LU of the metric (at most [`DENSE_FALLBACK_MAX_VERTICES`] vertices).
///
/// The first stage that works stays active for the rest of the step. Every escalation is
/// recorded in [`LinearReport::fallbacks`]. Only when no stage works is the error returned.
pub(crate) struct RobustHierarchical {
    primary: HierarchicalSaddle,
    ctx: FallbackContext,
    active: Active,
    events: Vec<String>,
}

impl RobustHierarchical {
    fn escalate(&mut self, err: RsError, a: &[f64], b: &[f64]) -> Result<SaddleSolution, RsError> {
        let RsError::SolverDiverged { residual } = err else {
            return Err(err);
        };
        if matches!(self.active, Active::Primary) {
            let c = &self.ctx;
            let mut o = HierarchicalOptions::default();
            o.cluster = c.options.cluster;
            o.gmres.tolerance = c.options.gmres.tolerance;
            o.gmres.max_iterations = (4 * c.options.gmres.max_iterations).max(800);
            o.gmres.restart = 60;
            let built = hierarchical_saddle(
                c.pool.clone(),
                &c.surface,
                &c.positions,
                &c.exponents,
                &c.weights,
                &c.extras,
                &c.c_q,
                &o,
                None,
            );
            if let Ok(mut g) = built {
                match g.solve(a, b) {
                    Ok(sol) => {
                        self.events.push(format!(
                            "gmres diverged (relative residual {residual:.3e}): retried with right GMRES, \
                             augmented preconditioner on vertex pairs, restart 60, at most {} iterations \
                             (final residual {:.3e})",
                            o.gmres.max_iterations, sol.record.relative_residual
                        ));
                        self.active = Active::Retry(Box::new(g));
                        return Ok(sol);
                    }
                    Err(RsError::SolverDiverged { .. }) | Err(RsError::SingularSaddle { .. }) => {}
                    Err(e) => return Err(e),
                }
            }
        }
        let nv = self.ctx.surface.n_vertices();
        if nv <= DENSE_FALLBACK_MAX_VERTICES && !matches!(self.active, Active::Dense(_)) {
            let c = &self.ctx;
            let dense = dense_metric_with_extras(&c.pool, &c.surface, &c.positions, &c.exponents, &c.weights, &c.extras)
                .and_then(|m| DenseSaddle::new(&c.pool, &m, &c.c_q));
            if let Ok(mut d) = dense {
                let sol = d.solve(a, b)?;
                self.events.push(format!(
                    "gmres diverged (relative residual {residual:.3e}): dense LU of the {nv} x {nv} metric \
                     (residual {:.3e})",
                    sol.record.relative_residual
                ));
                self.active = Active::Dense(Box::new(d));
                return Ok(sol);
            }
        }
        Err(RsError::SolverDiverged { residual })
    }
}

impl SaddleSolver for RobustHierarchical {
    fn n_dofs(&self) -> usize {
        self.primary.n_dofs()
    }

    fn n_constraints(&self) -> usize {
        self.primary.n_constraints()
    }

    fn solve(&mut self, a: &[f64], b: &[f64]) -> Result<SaddleSolution, RsError> {
        let r = match &mut self.active {
            Active::Primary => self.primary.solve(a, b),
            Active::Retry(g) => g.solve(a, b),
            Active::Dense(d) => d.solve(a, b),
        };
        match r {
            Err(e @ RsError::SolverDiverged { .. }) => self.escalate(e, a, b),
            r => r,
        }
    }

    fn base_report(&self) -> LinearReport {
        let mut r = self.primary.base_report();
        r.fallbacks = self.events.clone();
        r
    }
}

impl StepSolver {
    pub fn as_dyn(&mut self) -> &mut dyn SaddleSolver {
        match self {
            StepSolver::Dense(s) => s,
            StepSolver::Hierarchical(s) => s.as_mut(),
            StepSolver::Sparse(s) => s,
        }
    }

    /// The sparse comparison-metric solver (comparison optimizers, WP9).
    pub fn as_sparse(&mut self) -> Option<&mut SparseSaddle> {
        match self {
            StepSolver::Sparse(s) => Some(s),
            _ => None,
        }
    }

    /// Operator applications and their seconds (hierarchical only).
    pub fn products(&self) -> Option<(usize, f64)> {
        match self {
            StepSolver::Hierarchical(s) => Some(s.primary.operator().usage()),
            _ => None,
        }
    }

    /// Seconds of the hierarchical back end: `[metric tree build, fractional tree build, L′
    /// factorization, preconditioner applications]`.
    pub fn build_timings(&self) -> Option<[f64; 4]> {
        match self {
            StepSolver::Hierarchical(s) => {
                let p = s.preconditioner();
                Some([
                    s.operator().tree().build_seconds(),
                    p.fractional_seconds(),
                    p.factor_seconds(),
                    p.usage().1,
                ])
            }
            _ => None,
        }
    }

    /// Builds the solver of one step.
    #[allow(clippy::too_many_arguments)]
    pub fn build(
        pool: &Arc<rayon::ThreadPool>,
        surface: &Surface,
        positions: &[[f64; 3]],
        exponents: &Exponents,
        metric: Metric,
        weights: &SobolevWeights,
        extras: &MetricExtras,
        c_q: &SparseMatrix,
        backend: &Backend,
        face_tree: Option<&FaceTree>,
        bh_leaf_size: Option<usize>,
    ) -> Result<StepSolver, RsError> {
        match metric {
            Metric::Sobolev { .. } => match backend.resolve(surface.n_faces()) {
                None => {
                    let a = dense_metric_with_extras(pool, surface, positions, exponents, weights, extras)?;
                    Ok(StepSolver::Dense(DenseSaddle::new(pool, &a, c_q)?))
                }
                Some(o) => {
                    let donor = match (face_tree, bh_leaf_size) {
                        (Some(t), Some(l)) if l == o.cluster.leaf_size => Some(t.clone()),
                        _ => None,
                    };
                    let primary = hierarchical_saddle(
                        pool.clone(),
                        surface,
                        positions,
                        exponents,
                        weights,
                        extras,
                        c_q,
                        &o,
                        donor,
                    )?;
                    Ok(StepSolver::Hierarchical(Box::new(RobustHierarchical {
                        primary,
                        ctx: FallbackContext {
                            pool: pool.clone(),
                            surface: surface.clone(),
                            positions: positions.to_vec(),
                            exponents: *exponents,
                            weights: *weights,
                            extras: extras.clone(),
                            c_q: c_q.clone(),
                            options: o,
                        },
                        active: Active::Primary,
                        events: Vec::new(),
                    })))
                }
            },
            _ => Ok(StepSolver::Sparse(SparseSaddle::new(surface, positions, metric, c_q)?)),
        }
    }
}
