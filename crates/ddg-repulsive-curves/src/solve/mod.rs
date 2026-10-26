//! Saddle-point solves of `SPEC.md` §6: the projected Sobolev gradient (Eq. 22), the chord
//! projection onto the constraints (Eq. 23, D10), the null-space projector (B.3.3, D4) and the
//! [`SaddleSolver`] trait that the back ends implement.
//!
//! Every solve is of the form
//!
//! ```text
//! [ Ā  Cᵀ ] [ x ]   [ a ]
//! [ C  0  ] [ λ ] = [ b ]
//! ```
//!
//! with `Ā = A ⊗ I_d` in the interleaved layout of [`DofLayout`] (D29). Back ends:
//! [`dense::DenseSaddle`] (one LU with partial pivoting of the whole matrix, `faer`, reused for
//! Eq. 22 and every chord iteration; Sec. 7.4), [`sparse::SparseSaddle`] (sparse LU,
//! `ddg_core::sparse`, for the L², H¹, H² metrics of §4.6), and later the multigrid of §8.4
//! (WP6), which uses [`projector::Projector`] for Eq. 30 and 33.

use std::sync::Arc;

use ddg_core::sparse::SparseError;

use crate::constraints::{ConstraintSet, rank::check_component_translations};
use crate::error::RcError;
use crate::kernel::Exponents;
use crate::metric::{DofLayout, GramQuadrature, Metric, classic_metric, sobolev_gram};
use crate::network::CurveNetwork;

pub mod dense;
pub mod projector;
pub mod sparse;

pub use dense::{DenseMethod, DenseSaddle};
pub use projector::Projector;
pub use sparse::SparseSaddle;

/// Maps a sparse linear-algebra error to the crate's error.
///
/// Pivot failures become [`RcError::SingularSaddle`]; resource errors (out of memory, index
/// overflow), which have no dedicated variant, become [`RcError::InvalidOption`] naming the
/// matrix.
pub(crate) fn sparse_error(e: SparseError, what: &'static str) -> RcError {
    match e {
        SparseError::NotPositiveDefinite { step }
        | SparseError::ZeroPivot { step }
        | SparseError::StructurallySingular { step } => RcError::SingularSaddle { pivot: step },
        SparseError::NonFiniteSolution { index } => RcError::SingularSaddle { pivot: index },
        other => RcError::InvalidOption {
            name: what,
            value: other.to_string(),
            requirement: "the sparse linear algebra must succeed (resource limits)",
        },
    }
}

/// How a linear solve went (`SPEC.md` §9.1, `StepReport::linear`).
#[derive(Debug, Clone, PartialEq, Default)]
pub struct LinearReport {
    /// The back end: `"dense Schur"`, `"dense LU"`, `"sparse LU"`, `"multigrid"`.
    pub solver: &'static str,
    /// Dimension of the saddle system `d|V| + k`.
    pub dimension: usize,
    /// Iterations (0 for direct solvers).
    pub iterations: usize,
    /// `‖r‖ / ‖rhs‖` of the full saddle system (absolute `‖r‖` if the right side is 0).
    pub relative_residual: f64,
    /// Whether a fallback was used (multigrid: flexible CG).
    pub fallback_used: bool,
    /// Multigrid levels (0 for direct solvers).
    pub levels: usize,
    /// Dense solvers: `min |pivot| / max |pivot|` of the LU, or the smallest `min L_ii² / max L_ii²`
    /// of the Schur method's Cholesky factors ([`DenseSaddle::pivot_ratio`]).
    pub pivot_ratio: Option<f64>,
    /// Seconds spent in the solve (factorization excluded).
    pub seconds: f64,
}

/// Solution of one saddle system.
#[derive(Debug, Clone, PartialEq)]
pub struct SaddleSolution {
    /// `x ∈ R^{d|V|}`, interleaved.
    pub x: Vec<f64>,
    /// `λ ∈ R^k`.
    pub lambda: Vec<f64>,
    /// How the solve went.
    pub report: LinearReport,
}

/// A solver for `[Ā Cᵀ; C 0] [x; λ] = [a; b]` with `Ā` and `C` fixed at construction (the
/// metric and the constraint Jacobian of one step). Implemented by the dense and sparse
/// direct solvers here and by the multigrid of WP6.
pub trait SaddleSolver {
    /// `d|V|`.
    fn n_dofs(&self) -> usize;
    /// `k`, the number of constraint rows.
    fn n_constraints(&self) -> usize;
    /// Solves the system for one right-hand side (`a.len() = n_dofs`, `b.len() =
    /// n_constraints`).
    ///
    /// # Errors
    /// [`RcError::InvalidOption`] for wrong lengths; back-end failures.
    fn solve(&mut self, a: &[f64], b: &[f64]) -> Result<SaddleSolution, RcError>;
}

pub(crate) fn check_rhs(solver: &dyn SaddleSolver, a: &[f64], b: &[f64]) -> Result<(), RcError> {
    if a.len() != solver.n_dofs() || b.len() != solver.n_constraints() {
        return Err(RcError::InvalidOption {
            name: "saddle right-hand side",
            value: format!("({}, {})", a.len(), b.len()),
            requirement: "lengths (d|V|, k) of the factored system",
        });
    }
    Ok(())
}

/// The length tolerance of the constraint projection (`SPEC.md` §9.1, D10).
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Tolerance {
    /// An absolute length (the paper's `10⁻⁴`, used by the evaluation).
    Absolute(f64),
    /// `t · ℓ̄₀`, a fraction of the mean edge length at session start (the default `1e-4`).
    MeanEdgeLength(f64),
}

/// Options of the constraint projection (§6.2).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ProjectionOptions {
    /// The length tolerance `τ_len`; the per-row tolerances follow from it (§5.1).
    pub tolerance: Tolerance,
    /// Maximal number of chord iterations.
    pub max_iterations: u32,
}

impl Default for ProjectionOptions {
    /// `MeanEdgeLength(1e-4)`, 3 iterations (`SPEC.md` §10).
    fn default() -> Self {
        ProjectionOptions {
            tolerance: Tolerance::MeanEdgeLength(1e-4),
            max_iterations: 3,
        }
    }
}

impl ProjectionOptions {
    /// Checks the options.
    ///
    /// # Errors
    /// [`RcError::InvalidOption`] for a non-positive or non-finite tolerance.
    pub fn validate(&self) -> Result<(), RcError> {
        let t = match self.tolerance {
            Tolerance::Absolute(t) | Tolerance::MeanEdgeLength(t) => t,
        };
        if !(t.is_finite() && t > 0.0) {
            return Err(RcError::InvalidOption {
                name: "projection.tolerance",
                value: t.to_string(),
                requirement: "a finite positive tolerance",
            });
        }
        Ok(())
    }
}

/// The projected Sobolev gradient of Eq. 22.
#[derive(Debug, Clone, PartialEq)]
pub struct ProjectedGradient {
    /// `g̃` per vertex (0 along the held axis in planar mode). The descent direction is `−g̃`.
    pub direction: Vec<[f64; 3]>,
    /// The multipliers `λ`.
    pub lambda: Vec<f64>,
    /// `‖g̃‖₂`, Euclidean over all coordinates (the stop rule of D14).
    pub gradient_norm: f64,
    /// `dF · g̃ = g̃ᵀ Ā g̃` (the reference's `soboDot`).
    pub sobolev_norm_sq: f64,
    /// The linear solve.
    pub linear: LinearReport,
}

/// Solves Eq. 22, `[Ā Cᵀ; C 0] [g̃; λ] = [dFᵀ; 0]`. If `Ā` is positive definite on `ker C`
/// the solution is unique and `g̃` is the `Ā`-closest point of `ker C` to the unconstrained
/// Sobolev gradient; `dF · g̃ = g̃ᵀ Ā g̃ ≥ 0` with equality only for `g̃ = 0`, so `−g̃` is a
/// descent direction tangent to the constraint set.
///
/// Under `γ → cγ` the direction scales as `c^{2σ+2+α−β}` (`dÊ` as `c^{1+α−β}`, `Ā` as
/// `c^{−(2σ+1)}`), when the constraint rows are homogeneous of degree 0 (all of §5.1 except
/// `Surface` with a non-homogeneous `f`): rescaling changes the direction's length, not its
/// direction (paper A.4.2, with the exponents corrected, D3).
///
/// [formula: rc-saddle-gradient]
/// [formula: rc-descent-direction]
/// [formula: rc-gradient-scaling]
///
/// # Errors
/// The solver's errors; [`RcError::InvalidOption`] if `differential` has the wrong length.
pub fn projected_gradient(
    solver: &mut dyn SaddleSolver,
    layout: DofLayout,
    differential: &[[f64; 3]],
) -> Result<ProjectedGradient, RcError> {
    let a = layout.gather(differential);
    let b = vec![0.0; solver.n_constraints()];
    let sol = solver.solve(&a, &b)?;
    let gradient_norm = sol.x.iter().map(|x| x * x).sum::<f64>().sqrt();
    let sobolev_norm_sq = a.iter().zip(&sol.x).map(|(p, q)| p * q).sum();
    Ok(ProjectedGradient {
        direction: layout.scatter(&sol.x),
        lambda: sol.lambda,
        gradient_norm,
        sobolev_norm_sq,
        linear: sol.report,
    })
}

/// How a constraint projection ended.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ProjectionStatus {
    /// Every row within its tolerance and every tangent site oriented along its direction.
    Converged,
    /// `max_iterations` chord steps without reaching the tolerances.
    MaxIterations,
    /// The largest relative violation grew by more than 10× or became non-finite.
    Diverged,
    /// The rows hold but a tangent site ended with `s·X ≤ 0` (D7).
    TangentFlipped {
        /// Index of the tangent constraint.
        constraint: usize,
    },
}

/// Result of [`project_constraints`]. A status other than `Converged` rejects the trial step
/// (§7); it is not an error.
#[derive(Debug, Clone, PartialEq)]
pub struct Projection {
    /// The projected positions (the last iterate when not converged).
    pub positions: Vec<[f64; 3]>,
    /// How it ended.
    pub status: ProjectionStatus,
    /// Number of chord iterations (solves).
    pub iterations: u32,
    /// `max_r |Φ_r|` at the end (in the rows' own units).
    pub max_violation: f64,
    /// `max_r |Φ_r| / tol_r` per iterate, starting with the trial point (≤ 1 at success).
    pub history: Vec<f64>,
    /// The linear solves, one per iteration.
    pub linear: Vec<LinearReport>,
}

/// The chord projection of §6.2 (Eq. 23): from `trial`, repeat
/// `[Ā Cᵀ; C 0] [x; µ] = [0; −Φ(γ̃)]`, `γ̃ ← γ̃ + x` with `Ā` and `C` frozen at the step's
/// base point (the `solver`), until every active row satisfies `|Φ_r| ≤ tol_r` (§5.1), at least
/// once (unless `Φ = 0` exactly; the paper's "update γ̃ ← γ̃ + x and repeat", §17l) and at most
/// `options.max_iterations` times. Each iteration gives the `Ā`-smallest `x` solving the
/// linearized constraints. Failure (a status, not an error): the iteration limit, growth of the
/// largest relative violation by more than 10× over the trial's, or a tangent site with
/// `s·X ≤ 0` at the end. Local convergence of the chord method near a regular point is cited
/// (`rc-chord-convergence`, D10).
///
/// [formula: rc-chord-projection]
///
/// # Errors
/// The solver's errors; [`RcError::PositionCountMismatch`].
pub fn project_constraints(
    solver: &mut dyn SaddleSolver,
    constraints: &ConstraintSet,
    network: &CurveNetwork,
    trial: &[[f64; 3]],
    options: &ProjectionOptions,
    mean_edge_length0: f64,
) -> Result<Projection, RcError> {
    options.validate()?;
    let layout = constraints.layout();
    let mut pos = trial.to_vec();
    let mut history = Vec::new();
    let mut linear = Vec::new();
    let zero = vec![0.0; solver.n_dofs()];
    let mut it = 0u32;
    loop {
        let phi = constraints.values(network, &pos)?;
        let tol =
            constraints.row_tolerances(network, &pos, options.tolerance, mean_edge_length0)?;
        let mut rel: f64 = 0.0;
        let mut max_violation: f64 = 0.0;
        for (p, t) in phi.iter().zip(&tol) {
            let r = p.abs() / t;
            rel = if r.is_nan() || rel.is_nan() {
                f64::NAN
            } else {
                rel.max(r)
            };
            max_violation = max_violation.max(p.abs());
        }
        history.push(rel);
        let finish = |status, pos: Vec<[f64; 3]>, history, linear| {
            Ok(Projection {
                positions: pos,
                status,
                iterations: it,
                max_violation,
                history,
                linear,
            })
        };
        // At least one step, also for a trial already within the tolerance (`SPEC.md` §17l):
        // otherwise the curve drifts unprojected up to the tolerance (second-order edge
        // growth lowers Ê) and every larger trial pays for the way back.
        if rel <= 1.0 && (it > 0 || rel == 0.0 || options.max_iterations == 0) {
            let status = match constraints.tangent_orientation_failure(network, &pos) {
                Some((constraint, _)) => ProjectionStatus::TangentFlipped { constraint },
                None => ProjectionStatus::Converged,
            };
            return finish(status, pos, history, linear);
        }
        if !rel.is_finite() || (it > 0 && rel > 10.0 * history[0]) {
            return finish(ProjectionStatus::Diverged, pos, history, linear);
        }
        if it >= options.max_iterations {
            return finish(ProjectionStatus::MaxIterations, pos, history, linear);
        }
        let b: Vec<f64> = phi.iter().map(|p| -p).collect();
        let sol = solver.solve(&zero, &b)?;
        layout.add_scaled(&mut pos, 1.0, &sol.x);
        linear.push(sol.report);
        it += 1;
    }
}

/// The second-order correction of the line search's trials (`SPEC.md` §17l).
#[derive(Debug, Clone, PartialEq)]
pub struct SecondOrderCorrection {
    /// `q` per vertex: trials are `γ − τ p + τ² q`.
    pub q: Vec<[f64; 3]>,
    /// `max_r |Φ₂_r|`, the curvature term `Φ₂ = ½ Φ''(γ)[p, p]` (rows' own units per length²
    /// of `p`).
    pub curvature: f64,
    /// The linear solve (`None` when `Φ₂ = 0`, e.g. linear constraints only).
    pub linear: Option<LinearReport>,
}

/// Relative step of the central difference of [`second_order_correction`]: `h · max_v |p_v|`
/// is this fraction of the mean edge length `ℓ̄₀`.
pub const SECOND_ORDER_DIFFERENCE_STEP: f64 = 1e-3;

/// The second-order correction of the trials `γ − τ p` (`SPEC.md` §17l; the "second-order
/// correction" of SQP methods, Nocedal & Wright, *Numerical Optimization*, 2nd ed., §18.3).
///
/// With `C p = 0` (the projected gradient of Eq. 22 satisfies it exactly), the constraint values
/// along the straight trial are `Φ(γ − τ p) = Φ(γ) + τ² Φ₂ + O(τ³)`, `Φ₂ = ½ Φ''(γ)[p, p]`. The
/// `Ā`-smallest `q` with `C q = −Φ₂` (the saddle system with right side `(0, −Φ₂)`, one solve with
/// the step's `solver`) gives `Φ(γ − τ p + τ² q) = Φ(γ) + O(τ³)` for every `τ`, so one solve per
/// step replaces the `O(τ²)` part of every trial's chord projection (§6.2). `τ² q` is the leading
/// term of the projection's own first correction from `γ − τ p`, so the projected curve is the
/// same to that order; the projection that follows keeps every tolerance of §6.2.
///
/// `Φ₂` is the central second difference
/// `(Φ(γ + h p) + Φ(γ − h p) − 2 Φ(γ)) / (2 h²)` with `h max_v |p_v| =`
/// [`SECOND_ORDER_DIFFERENCE_STEP`]` · ℓ̄₀` (truncation `O(h²)`, rounding `ε ℓ / h²`; both far
/// below the projection tolerance after the factor `τ²`), generic over every constraint type.
///
/// [formula: rc-chord-projection]
///
/// # Errors
/// The solver's errors; the errors of [`ConstraintSet::values`].
pub fn second_order_correction(
    solver: &mut dyn SaddleSolver,
    constraints: &ConstraintSet,
    network: &CurveNetwork,
    positions: &[[f64; 3]],
    direction: &[[f64; 3]],
    mean_edge_length0: f64,
) -> Result<SecondOrderCorrection, RcError> {
    let layout = constraints.layout();
    let zero_q = || vec![[0.0; 3]; positions.len()];
    let pmax = direction
        .iter()
        .map(|d| (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt())
        .fold(0.0f64, f64::max);
    if !(pmax > 0.0 && pmax.is_finite() && mean_edge_length0 > 0.0) || solver.n_constraints() == 0 {
        return Ok(SecondOrderCorrection {
            q: zero_q(),
            curvature: 0.0,
            linear: None,
        });
    }
    let h = SECOND_ORDER_DIFFERENCE_STEP * mean_edge_length0 / pmax;
    let shifted = |s: f64| -> Vec<[f64; 3]> {
        positions
            .iter()
            .zip(direction)
            .map(|(x, d)| [x[0] + s * d[0], x[1] + s * d[1], x[2] + s * d[2]])
            .collect()
    };
    let fp = constraints.values(network, &shifted(h))?;
    let fm = constraints.values(network, &shifted(-h))?;
    let f0 = constraints.values(network, positions)?;
    let b: Vec<f64> = (0..f0.len())
        .map(|r| -(fp[r] + fm[r] - 2.0 * f0[r]) / (2.0 * h * h))
        .collect();
    let curvature = b.iter().fold(0.0f64, |m, v| m.max(v.abs()));
    if curvature == 0.0 || !curvature.is_finite() {
        return Ok(SecondOrderCorrection {
            q: zero_q(),
            curvature,
            linear: None,
        });
    }
    let sol = solver.solve(&vec![0.0; solver.n_dofs()], &b)?;
    Ok(SecondOrderCorrection {
        q: layout.scatter(&sol.x),
        curvature,
        linear: Some(sol.report),
    })
}

/// Maximal number of Newton steps of the polish ([`newton_correction`]) in one projection.
pub const NEWTON_POLISH_MAX_STEPS: u32 = 8;

/// One Newton step of the projection's polish (`SPEC.md` §17l): the smallest correction `x` in
/// the lumped-mass metric `M = diag(m_v) ⊗ I_d` (`m_v` the dual lengths at `positions`) with
/// `C(γ̃) x = −Φ(γ̃)`, the Jacobian taken at the iterate `γ̃ = positions` itself:
/// `x = M⁻¹ Cᵀ (C M⁻¹ Cᵀ)⁻¹ (−Φ)`, by a sparse Cholesky factorization of `C M⁻¹ Cᵀ` (`O(k)` for
/// the constraints of §5.1). Newton's method for the constraints converges quadratically near
/// a regular point, where the chord iteration of Eq. 23 (Jacobian frozen at the base curve)
/// converges linearly at a rate proportional to the step; the metric only selects the
/// right-inverse, and the polish is used for small violations only.
///
/// [formula: rc-chord-projection]
///
/// # Errors
/// The errors of [`ConstraintSet::jacobian`] and [`Projector::new`] (a rank-deficient Jacobian
/// at the iterate).
pub fn newton_correction(
    constraints: &ConstraintSet,
    network: &CurveNetwork,
    positions: &[[f64; 3]],
    phi: &[f64],
) -> Result<Vec<f64>, RcError> {
    let d = constraints.layout().dim();
    let mut c = constraints.jacobian(network, positions)?;
    let m = network.dual_lengths(positions);
    let w: Vec<f64> = m
        .iter()
        .map(|x| 1.0 / x.max(f64::MIN_POSITIVE).sqrt())
        .collect();
    let cp = c.col_ptr().to_vec();
    let vals = c.values_mut();
    for j in 0..cp.len() - 1 {
        for v in &mut vals[cp[j]..cp[j + 1]] {
            *v *= w[j / d];
        }
    }
    let pr = Projector::new(&c, None)?;
    let b: Vec<f64> = phi.iter().map(|p| -p).collect();
    let mut x = pr.pseudo_inverse(&b)?;
    for (j, v) in x.iter_mut().enumerate() {
        *v *= w[j / d];
    }
    Ok(x)
}

/// The dense Sobolev saddle solver of one step: assembles `A` (Eq. 31) on `pool`, the
/// constraint Jacobian at `positions`, and factors the saddle matrix. A singular saddle (the
/// rank checks of §5.2 having passed) means `Ā` is singular on `ker C`, which is reported as
/// [`RcError::SingularMetric`] (e.g. a collinear network, §4.5).
///
/// # Errors
/// [`RcError::SingularMetric`]; the errors of [`ConstraintSet::jacobian`] and
/// [`DenseSaddle::new`].
pub fn dense_sobolev_solver(
    pool: Arc<rayon::ThreadPool>,
    network: &CurveNetwork,
    positions: &[[f64; 3]],
    exponents: Exponents,
    quadrature: GramQuadrature,
    constraints: &ConstraintSet,
) -> Result<DenseSaddle, RcError> {
    let gram = sobolev_gram(&pool, network, positions, exponents, quadrature);
    let c = constraints.jacobian(network, positions)?;
    DenseSaddle::new(pool, &gram.a, &c, constraints.layout()).map_err(|e| match e {
        RcError::SingularSaddle { .. } => RcError::SingularMetric {
            components: network.components().len(),
        },
        e => e,
    })
}

/// The sparse saddle solver for a classic metric (`L2`, `H1`, `H2`, §4.6). Checks first that
/// the constraints fix every component's translation (the kernel of `L` and `L M⁻¹ L`).
///
/// # Errors
/// [`RcError::SingularMetric`]; [`RcError::InvalidOption`] for [`Metric::Sobolev`]; the errors
/// of [`SparseSaddle::new`].
pub fn classic_solver(
    network: &CurveNetwork,
    positions: &[[f64; 3]],
    metric: &Metric,
    constraints: &ConstraintSet,
) -> Result<SparseSaddle, RcError> {
    let m = classic_metric(network, positions, metric)?;
    let c = constraints.jacobian(network, positions)?;
    if !matches!(metric, Metric::L2) {
        check_component_translations(&c, network, constraints.layout())?;
    }
    SparseSaddle::new(&m, &c, constraints.layout()).map_err(|e| match e {
        RcError::SingularSaddle { .. } => RcError::SingularMetric {
            components: network.components().len(),
        },
        e => e,
    })
}
