//! The dense fallback of the iterative solvers (`SPEC.md` §17h).
//!
//! The multigrid solver (and any other iterative [`SolverFactory`]) can fail on a step:
//! [`RcError::MultigridDiverged`] (the residual grew) or [`RcError::MultigridNotConverged`] (it
//! stagnated above `failure_residual`). When the network has at most
//! [`FlowOptions::dense_fallback_max_vertices`](super::FlowOptions::dense_fallback_max_vertices)
//! vertices, `FallbackSaddle` builds the dense saddle solver of the same step (same metric,
//! same constraint rows, same pools), repeats the failed solve with it, and answers every later
//! solve of the step with it. The event is recorded in [`SolverFallback`] and ends up in
//! [`StepReport::solver_fallback`](super::StepReport::solver_fallback). Larger networks, and
//! every other error, are returned unchanged.

use std::sync::{Arc, Mutex};

use ddg_core::sparse::SparseMatrix;

use super::{DirectSolvers, SolverFactory, SolverRequest};
use crate::error::RcError;
use crate::kernel::Exponents;
use crate::metric::{DofLayout, Metric};
use crate::network::CurveNetwork;
use crate::solve::{SaddleSolution, SaddleSolver};

/// A step whose iterative solver failed and was replaced by the dense solver
/// ([`StepReport::solver_fallback`](super::StepReport::solver_fallback)).
#[derive(Debug, Clone, PartialEq)]
pub struct SolverFallback {
    /// The solver that failed ([`SolverFactory::name`]).
    pub from: &'static str,
    /// The solver that took over (`"direct"`).
    pub to: &'static str,
    /// The error of the first failed solve ([`RcError::MultigridDiverged`] or
    /// [`RcError::MultigridNotConverged`]).
    pub cause: RcError,
    /// Number of solves answered by the dense solver in this step (the failed one included).
    pub dense_solves: usize,
    /// Seconds spent building the dense solver.
    pub build_seconds: f64,
}

/// Shared slot through which the session reads the events of a step.
pub(crate) type FallbackSlot = Arc<Mutex<Option<SolverFallback>>>;

/// Owned copy of a [`SolverRequest`] (the dense solver is built lazily, after the borrowed
/// request is gone).
struct OwnedRequest {
    network: CurveNetwork,
    positions: Vec<[f64; 3]>,
    exponents: Exponents,
    metric: Metric,
    jacobian: SparseMatrix,
    layout: DofLayout,
    pool: Arc<rayon::ThreadPool>,
    dense_pool: Arc<rayon::ThreadPool>,
}

/// A saddle solver that falls back to the dense solver (module docs).
pub(crate) struct FallbackSaddle {
    primary: Box<dyn SaddleSolver>,
    from: &'static str,
    request: OwnedRequest,
    dense: Option<Box<dyn SaddleSolver>>,
    slot: FallbackSlot,
}

impl FallbackSaddle {
    /// Wraps `primary`, built from `request` by the factory named `from`.
    pub(crate) fn new(
        primary: Box<dyn SaddleSolver>,
        from: &'static str,
        request: &SolverRequest<'_>,
        slot: FallbackSlot,
    ) -> FallbackSaddle {
        FallbackSaddle {
            primary,
            from,
            request: OwnedRequest {
                network: request.network.clone(),
                positions: request.positions.to_vec(),
                exponents: request.exponents,
                metric: request.metric,
                jacobian: request.jacobian.clone(),
                layout: request.layout,
                pool: request.pool.clone(),
                dense_pool: request.dense_pool.clone(),
            },
            dense: None,
            slot,
        }
    }

    fn dense_solve(
        &mut self,
        cause: Option<RcError>,
        a: &[f64],
        b: &[f64],
    ) -> Result<SaddleSolution, RcError> {
        if self.dense.is_none() {
            let t = std::time::Instant::now();
            let r = &self.request;
            let built = DirectSolvers::new().build(&SolverRequest {
                network: &r.network,
                positions: &r.positions,
                exponents: r.exponents,
                metric: r.metric,
                jacobian: &r.jacobian,
                layout: r.layout,
                constraints: None,
                pool: &r.pool,
                dense_pool: &r.dense_pool,
            })?;
            self.dense = Some(built);
            let mut slot = self.slot.lock().unwrap_or_else(|e| e.into_inner());
            if let Some(cause) = cause {
                *slot = Some(SolverFallback {
                    from: self.from,
                    to: "direct",
                    cause,
                    dense_solves: 0,
                    build_seconds: t.elapsed().as_secs_f64(),
                });
            }
        }
        let sol = self
            .dense
            .as_mut()
            .expect("dense solver was just built")
            .solve(a, b)?;
        if let Some(f) = self.slot.lock().unwrap_or_else(|e| e.into_inner()).as_mut() {
            f.dense_solves += 1;
        }
        Ok(sol)
    }
}

impl SaddleSolver for FallbackSaddle {
    fn n_dofs(&self) -> usize {
        self.primary.n_dofs()
    }

    fn n_constraints(&self) -> usize {
        self.primary.n_constraints()
    }

    fn solve(&mut self, a: &[f64], b: &[f64]) -> Result<SaddleSolution, RcError> {
        if self.dense.is_some() {
            return self.dense_solve(None, a, b);
        }
        match self.primary.solve(a, b) {
            Err(
                e @ (RcError::MultigridDiverged { .. } | RcError::MultigridNotConverged { .. }),
            ) => self.dense_solve(Some(e), a, b),
            r => r,
        }
    }
}
