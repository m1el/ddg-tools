//! Linear solves for the constrained direction and the projection (`SPEC.md` §6).
//!
//! The positional rows `C_Q` (barycenters, pins) live inside the saddle matrix
//!
//! ```text
//! A′ = [ A₃   C_Qᵀ ]        A₃ = A ⊗ I₃ (interleaved: unknown 3 v + axis)
//!      [ C_Q  0    ]
//! ```
//!
//! and a [`SaddleSolver`] applies `A′⁻¹`. The nonlinear rows (area, volume) go through the
//! Schur complement of [`schur`], generic over the solver: [`dense::DenseSaddle`] (LU, small
//! meshes, the oracle) or [`gmres::GmresSaddle`] (GMRES with any operator and preconditioner;
//! WP6 plugs in the hierarchical product and the sparse preconditioner of `precond`).
//!
//! Filled by WP2 (`dense`, `schur`, `gmres`) and WP6 (`iterative`, `precond`, SPEC §15).

pub mod coarse;
pub mod dense;
pub mod gmres;
pub mod iterative;
pub mod ldlt;
pub mod precond;
pub mod schur;

pub use dense::{DensePreconditioner, DenseSaddle, DenseSaddleOperator};
pub use gmres::{
    GmresOptions, GmresReport, GmresSaddle, GmresSide, IdentityPreconditioner, LinearOperator,
    Preconditioner, gmres, gmres_batch, gmres_batch_on,
};
pub use schur::{
    Direction, NewtonSteps, Projection, ProjectionOptions, ProjectionStatus, SchurComplement,
    project,
};

use crate::error::RsError;

/// One solve with `A′` (`SPEC.md` §9.1 `StepReport::linear`).
#[derive(Debug, Clone, PartialEq)]
pub struct SolveRecord {
    /// Iterations (0 for a direct solve).
    pub iterations: usize,
    /// `‖rhs − A′ x‖ / ‖rhs‖` of the full saddle system (absolute if the right side is 0).
    pub relative_residual: f64,
    /// Whether the solve reached its tolerance (always for direct solves).
    pub converged: bool,
    /// The GMRES details (iterative solves only).
    pub gmres: Option<GmresReport>,
}

/// How the linear solves of one step went (`SPEC.md` §9.1).
#[derive(Debug, Clone, PartialEq, Default)]
pub struct LinearReport {
    /// The back end: `"dense LU"`, `"dense LU (per axis)"`, `"gmres"`, …
    pub solver: &'static str,
    /// Dimension of `A′`: `3|V| + k_Q`.
    pub dimension: usize,
    /// Every solve with `A′`, in order (`k_N` for the Schur columns, then the direction).
    pub solves: Vec<SolveRecord>,
    /// Dense solvers: `min |pivot| / max |pivot|` of the equilibrated LU.
    pub pivot_ratio: Option<f64>,
    /// Schur rows dropped for this step (indices into the constraint set's rows, §6.2).
    pub dropped_rows: Vec<usize>,
    /// The smallest pivot of the equilibrated Schur complement (1 without Schur rows).
    pub schur_pivot: f64,
    /// Seconds spent in factorizations.
    pub factor_seconds: f64,
    /// Seconds spent in solves.
    pub solve_seconds: f64,
    /// Recoverable solver failures of this step and what was done about them (empty normally):
    /// e.g. `"gmres diverged (relative residual 1.5e3): retried with right GMRES, …"` or
    /// `"…: dense LU of the 900 x 900 metric"` (`SPEC.md` §6.3).
    pub fallbacks: Vec<String>,
}

impl LinearReport {
    /// Total iterations over all solves.
    pub fn total_iterations(&self) -> usize {
        self.solves.iter().map(|s| s.iterations).sum()
    }

    /// The largest relative residual over all solves (0 without solves).
    pub fn max_relative_residual(&self) -> f64 {
        self.solves
            .iter()
            .fold(0.0f64, |m, s| m.max(s.relative_residual))
    }

    /// Whether every solve converged.
    pub fn converged(&self) -> bool {
        self.solves.iter().all(|s| s.converged)
    }
}

/// Solution of one saddle system.
#[derive(Debug, Clone, PartialEq)]
pub struct SaddleSolution {
    /// `x ∈ R^{3|V|}`, interleaved.
    pub x: Vec<f64>,
    /// `λ ∈ R^{k_Q}`.
    pub lambda: Vec<f64>,
    /// How the solve went.
    pub record: SolveRecord,
}

/// A solver for `A′ [x; λ] = [a; b]` with `A′ = [[A₃, C_Qᵀ], [C_Q, 0]]` fixed at construction
/// (the metric and the positional Jacobian of one step).
pub trait SaddleSolver {
    /// `3|V|`.
    fn n_dofs(&self) -> usize;
    /// `k_Q`, the number of positional rows.
    fn n_constraints(&self) -> usize;
    /// Solves for one right-hand side (`a.len() = n_dofs`, `b.len() = n_constraints`).
    ///
    /// # Errors
    /// [`RsError::InvalidOption`] for wrong lengths; back-end failures
    /// ([`RsError::SolverDiverged`], …).
    fn solve(&mut self, a: &[f64], b: &[f64]) -> Result<SaddleSolution, RsError>;
    /// Solves for several right-hand sides `(a, b)` (default: one [`SaddleSolver::solve`]
    /// each, in order). Iterative back ends advance the solves together and batch their
    /// products ([`gmres::gmres_batch`]); every solution is the one `solve` would return.
    ///
    /// # Errors
    /// As [`SaddleSolver::solve`] (the first failing solve).
    fn solve_many(&mut self, rhs: &[(&[f64], &[f64])]) -> Result<Vec<SaddleSolution>, RsError> {
        rhs.iter().map(|(a, b)| self.solve(a, b)).collect()
    }
    /// The report of the back end so far (factorization data; the solves are collected by the
    /// caller from [`SaddleSolution::record`]).
    fn base_report(&self) -> LinearReport;
}

pub(crate) fn check_rhs(solver: &dyn SaddleSolver, a: &[f64], b: &[f64]) -> Result<(), RsError> {
    if a.len() != solver.n_dofs() || b.len() != solver.n_constraints() {
        return Err(RsError::InvalidOption {
            name: "saddle right-hand side",
            value: format!("({}, {})", a.len(), b.len()),
            requirement: "lengths (3|V|, k_Q) of the factored system",
        });
    }
    Ok(())
}

/// Euclidean norm.
pub(crate) fn norm2(x: &[f64]) -> f64 {
    x.iter().map(|v| v * v).sum::<f64>().sqrt()
}
