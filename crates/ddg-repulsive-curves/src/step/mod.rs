//! The flow: options, the step-at-a-time [`Session`], step anatomy and the [`flow`] driver
//! (`SPEC.md` §7, §9.1).
//!
//! # One step (§7.1)
//!
//! 0. **Re-target** (when a length schedule moves the targets, or a constraint is violated
//!    beyond its tolerance, e.g. new surface rows after a subdivision): new targets, projection
//!    of `γ` onto them (§6.2, matrices at `γ`), the projection path checked by the CCD when the
//!    line search is collision-guarded. On failure the increment is halved, up to 5 times, else
//!    the old targets are kept ([`RetargetStatus::Stalled`]). No Armijo test (D13).
//! 1. `F`, `dF` (energy back end + penalties); the metric and its saddle solver (back end).
//! 2. The direction: the projected gradient `g̃` of Eq. 22 ([`Method::ProjectedGradient`]) or,
//!    on soft problems, the L-BFGS or NCG direction with the mean-free preconditioner (§6.4).
//!    Gradient stop test (§7.5).
//! 3. Line search ([`line_search::LineSearch`], §7.2 or §7.3).
//! 4. Accept, or [`StepStatus::Stuck`].
//! 5. Subdivision check (§5.4): when the mean edge length exceeds twice `ℓ̄₀`, every edge is
//!    split (at most `max_subdivisions` times); constraints, penalties and the schedule follow
//!    the mapping of [`crate::subdivide`]. The vertex count of the session changes.
//!
//! Every part is timed in the [`StepReport`].
//!
//! # Back ends (hooks for later work packages)
//!
//! * Energy: [`EnergyModel`], chosen by [`energy_model`] from [`EnergyEval`].
//!   [`EnergyEval::Exact`] is [`ExactEnergy`]; [`EnergyEval::BarnesHut`] is [`BarnesHutEnergy`]
//!   `prepare` is called once per step at the base curve
//!   (Barnes–Hut: build the tree), the line search's trials only refit (D26).
//! * Metric and solver: [`SolverFactory`], one [`SaddleSolver`] per step and per re-target.
//!   [`Solver::Dense`] is [`DirectSolvers`] (dense LU for the Sobolev metric, sparse LU for
//!   L², H¹, H²); the hierarchical/multigrid solver of WP5/WP6 plugs in as
//!   [`Solver::Custom`].
//! * Penalties and obstacles: [`PenaltyTerm`] (WP7 implements its penalties against it; the
//!   Hencky penalty of the evaluation is [`crate::optimize::Hencky`]) and
//!   [`CollisionObstacle`] for the CCD (from [`FlowOptions::obstacles`] and from
//!   [`PenaltyTerm::obstacle`]).
//!
//! # Determinism (D25)
//!
//! The energy, the Gram matrix, the LU factorization and the CCD are bitwise independent of the
//! thread count, the triangular solves and every reduction here run in a fixed order: a run gives
//! bitwise identical curves for any `threads` and `dense_threads`.
//!
//! [formula: rc-gradient-scaling]
//! [formula: rc-isotopy]

use std::fmt;
use std::sync::Arc;
use std::sync::atomic::{AtomicBool, Ordering};
use std::time::{Duration, Instant};

use ddg_core::sparse::SparseMatrix;

use crate::ccd::{CcdOptions, CcdReport, CollisionObstacle, PreparedObstacle, first_contact};
use crate::constraints::{Constraint, ConstraintSet, default_constraints};
use crate::energy::barnes_hut::BhTree;
use crate::energy::exact::{energy_and_differential, energy_only, make_pool};
use crate::energy::{BhOptions, EnergyEval};
use crate::error::RcError;
use crate::kernel::Exponents;
use crate::metric::{DofLayout, Metric, classic_metric, sobolev_gram};
use crate::network::{Axis, CurveNetwork, check_planar};
use crate::optimize::{Lbfgs, Method, Ncg, dotv, mean_free_rows};
use crate::solve::{
    DenseSaddle, LinearReport, Projection, ProjectionOptions, ProjectionStatus, SaddleSolver,
    SparseSaddle, projected_gradient,
};
use crate::subdivide::{SubdivisionMap, subdivide, subdivide_constraints};

pub mod line_search;
pub mod schedule;
pub mod stop;

use line_search::{
    LineSearch, LineSearchReport, SearchEnv, SearchInput, Trial, line_search, project_with_path,
};
use schedule::{LengthSchedule, ScheduleState};
use stop::{Status, StepStatus, StopRule, StopState, mass_norm};

/// Default size of the dense solver's pool: `faer`'s LU stops scaling past 16–32 threads (WP2,
/// `SPEC.md` §17b; measured again by WP3, see the crate notes).
pub const DEFAULT_DENSE_THREADS: usize = 16;

// =============================================================================================
// Hooks
// =============================================================================================

/// A penalty term of the objective (`SPEC.md` §5.3). Work package WP7 implements its penalties
/// against this trait; the Hencky penalty of the evaluation is [`crate::optimize::Hencky`].
pub trait PenaltyTerm: Send + Sync + fmt::Debug {
    /// A short name for reports.
    fn name(&self) -> &str;
    /// Checks that the term fits `network` (e.g. one rest length per edge).
    ///
    /// # Errors
    /// [`RcError::InvalidOption`] naming the mismatch.
    fn check(&self, _network: &CurveNetwork) -> Result<(), RcError> {
        Ok(())
    }
    /// The value at `positions`; adds the gradient into `gradient` if given.
    ///
    /// # Errors
    /// Evaluation failures.
    fn evaluate(
        &self,
        network: &CurveNetwork,
        positions: &[[f64; 3]],
        gradient: Option<&mut [[f64; 3]]>,
    ) -> Result<f64, RcError>;
    /// The obstacle the CCD must avoid, if the term has one (surface potential).
    fn obstacle(&self) -> Option<Arc<dyn CollisionObstacle>> {
        None
    }
    /// The term on the subdivided network (§5.4).
    ///
    /// # Errors
    /// [`RcError::Unsupported`] if the term cannot follow a subdivision.
    fn subdivided(&self, map: &SubdivisionMap) -> Result<Arc<dyn PenaltyTerm>, RcError>;
}

/// The tangent-point energy at the base curve of a step.
#[derive(Debug, Clone, PartialEq)]
pub struct EnergyEvaluation {
    /// `Ê`.
    pub energy: f64,
    /// `∂Ê/∂γ_v`.
    pub differential: Vec<[f64; 3]>,
    /// Ordered edge pairs evaluated exactly.
    pub pairs_exact: u64,
    /// Lumped nodes used (Barnes–Hut).
    pub nodes_lumped: u64,
}

/// An evaluation back end of the tangent-point energy (exact here, Barnes–Hut in WP4).
pub trait EnergyModel: Send {
    /// A short name for reports.
    fn name(&self) -> &'static str;
    /// Called once per step at the base curve, before [`EnergyModel::energy_and_differential`]
    /// (Barnes–Hut: build the tree; the line search's trials only refit, D26).
    ///
    /// # Errors
    /// Back-end failures.
    fn prepare(&mut self, _network: &CurveNetwork, _positions: &[[f64; 3]]) -> Result<(), RcError> {
        Ok(())
    }
    /// `Ê` at a trial curve.
    ///
    /// # Errors
    /// Back-end failures.
    fn energy(&mut self, network: &CurveNetwork, positions: &[[f64; 3]]) -> Result<f64, RcError>;
    /// `Ê` and its differential at the base curve.
    ///
    /// # Errors
    /// Back-end failures.
    fn energy_and_differential(
        &mut self,
        network: &CurveNetwork,
        positions: &[[f64; 3]],
    ) -> Result<EnergyEvaluation, RcError>;
}

/// The exact `O(|E|²)` energy (`SPEC.md` §3) on an explicit pool.
#[derive(Debug, Clone)]
pub struct ExactEnergy {
    pool: Arc<rayon::ThreadPool>,
    exponents: Exponents,
}

impl ExactEnergy {
    /// The exact energy with `exponents` on `pool`.
    pub fn new(pool: Arc<rayon::ThreadPool>, exponents: Exponents) -> ExactEnergy {
        ExactEnergy { pool, exponents }
    }
}

fn exact_pairs(network: &CurveNetwork) -> u64 {
    let m = network.n_edges() as u64;
    let adjacent: u64 = (0..network.n_vertices() as u32)
        .map(|v| {
            let d = network.degree(v) as u64;
            d * (d - 1)
        })
        .sum();
    m * (m - 1) - adjacent
}

impl EnergyModel for ExactEnergy {
    fn name(&self) -> &'static str {
        "exact"
    }

    fn energy(&mut self, network: &CurveNetwork, positions: &[[f64; 3]]) -> Result<f64, RcError> {
        Ok(energy_only(&self.pool, network, positions, self.exponents))
    }

    fn energy_and_differential(
        &mut self,
        network: &CurveNetwork,
        positions: &[[f64; 3]],
    ) -> Result<EnergyEvaluation, RcError> {
        let (energy, differential) =
            energy_and_differential(&self.pool, network, positions, self.exponents);
        Ok(EnergyEvaluation {
            energy,
            differential,
            pairs_exact: exact_pairs(network),
            nodes_lumped: 0,
        })
    }
}

/// The Barnes–Hut energy (`SPEC.md` §8) on an explicit pool.
///
/// [`EnergyModel::prepare`] builds the tree at the base curve of the step (so it is rebuilt after
/// every accepted step and after a subdivision, the network of `prepare` being the current one);
/// the trial energies of the line search only [`refit`](BhTree::refit) it (D26). Both sides of
/// every Armijo comparison are Barnes–Hut energies (§17d).
#[derive(Debug)]
pub struct BarnesHutEnergy {
    pool: Arc<rayon::ThreadPool>,
    exponents: Exponents,
    options: BhOptions,
    tree: Option<BhTree>,
}

impl BarnesHutEnergy {
    /// The Barnes–Hut energy with `exponents` and `options` on `pool`.
    pub fn new(
        pool: Arc<rayon::ThreadPool>,
        exponents: Exponents,
        options: BhOptions,
    ) -> BarnesHutEnergy {
        BarnesHutEnergy {
            pool,
            exponents,
            options,
            tree: None,
        }
    }

    /// The current tree (after [`EnergyModel::prepare`]).
    pub fn tree(&self) -> Option<&BhTree> {
        self.tree.as_ref()
    }

    /// Refits the tree to `positions`, or builds it if `prepare` was not called. The topology
    /// only changes through `prepare` (a new step), never inside a line search.
    fn fit(&mut self, network: &CurveNetwork, positions: &[[f64; 3]]) -> Result<&BhTree, RcError> {
        match &mut self.tree {
            Some(t) => t.refit(positions)?,
            slot @ None => *slot = Some(BhTree::build(network, positions, self.options)?),
        }
        Ok(self.tree.as_ref().expect("tree was just set"))
    }
}

impl EnergyModel for BarnesHutEnergy {
    fn name(&self) -> &'static str {
        "barnes-hut"
    }

    fn prepare(&mut self, network: &CurveNetwork, positions: &[[f64; 3]]) -> Result<(), RcError> {
        self.tree = Some(BhTree::build(network, positions, self.options)?);
        Ok(())
    }

    fn energy(&mut self, network: &CurveNetwork, positions: &[[f64; 3]]) -> Result<f64, RcError> {
        let (pool, exponents) = (self.pool.clone(), self.exponents);
        Ok(self.fit(network, positions)?.energy(&pool, exponents).0)
    }

    fn energy_and_differential(
        &mut self,
        network: &CurveNetwork,
        positions: &[[f64; 3]],
    ) -> Result<EnergyEvaluation, RcError> {
        let (pool, exponents) = (self.pool.clone(), self.exponents);
        let (energy, differential, stats) = self
            .fit(network, positions)?
            .energy_and_differential(&pool, exponents);
        Ok(EnergyEvaluation {
            energy,
            differential,
            pairs_exact: stats.pairs_exact,
            nodes_lumped: stats.nodes_lumped,
        })
    }
}

/// The energy back end for `eval`.
///
/// # Errors
/// [`RcError::InvalidOption`] for bad Barnes–Hut options.
pub fn energy_model(
    eval: &EnergyEval,
    exponents: Exponents,
    pool: Arc<rayon::ThreadPool>,
) -> Result<Box<dyn EnergyModel>, RcError> {
    match eval {
        EnergyEval::Exact => Ok(Box::new(ExactEnergy::new(pool, exponents))),
        EnergyEval::BarnesHut(opts) => {
            opts.validate()?;
            Ok(Box::new(BarnesHutEnergy::new(pool, exponents, *opts)))
        }
    }
}

/// What a solver back end gets to build the saddle solver of one step.
#[derive(Clone, Copy)]
pub struct SolverRequest<'a> {
    /// The network.
    pub network: &'a CurveNetwork,
    /// The base curve of the step (where the metric and `C` are taken).
    pub positions: &'a [[f64; 3]],
    /// The exponents (Sobolev metric).
    pub exponents: Exponents,
    /// The metric.
    pub metric: Metric,
    /// The constraint rows `C` (`k × d|V|`): the Jacobian of the hard constraints, or the
    /// mean-free rows `M_c` of a soft problem.
    pub jacobian: &'a SparseMatrix,
    /// The layout of the unknowns.
    pub layout: DofLayout,
    /// The hard constraints (`None` on soft problems), for back ends that coarsen them (WP6).
    pub constraints: Option<&'a ConstraintSet>,
}

/// Builds the saddle solver of a step (`[Ā Cᵀ; C 0]`).
pub trait SolverFactory: Send + Sync {
    /// A short name for reports.
    fn name(&self) -> &'static str;
    /// Assembles the metric at `request.positions` and factors (or prepares) the saddle solve.
    ///
    /// # Errors
    /// [`RcError::SingularMetric`] and back-end failures.
    fn build(&self, request: &SolverRequest<'_>) -> Result<Box<dyn SaddleSolver>, RcError>;
}

/// The direct solvers: dense LU for the Sobolev metric (`faer`, on the dense pool; the Gram
/// matrix is assembled on the main pool), sparse LU for L², H¹, H² (§8.5).
#[derive(Debug, Clone)]
pub struct DirectSolvers {
    pool: Arc<rayon::ThreadPool>,
    dense_pool: Arc<rayon::ThreadPool>,
}

impl DirectSolvers {
    /// The direct solvers on `pool` (assembly) and `dense_pool` (dense LU).
    pub fn new(pool: Arc<rayon::ThreadPool>, dense_pool: Arc<rayon::ThreadPool>) -> Self {
        DirectSolvers { pool, dense_pool }
    }
}

impl SolverFactory for DirectSolvers {
    fn name(&self) -> &'static str {
        "direct"
    }

    fn build(&self, r: &SolverRequest<'_>) -> Result<Box<dyn SaddleSolver>, RcError> {
        let singular = |e| match e {
            RcError::SingularSaddle { .. } => RcError::SingularMetric {
                components: r.network.components().len(),
            },
            e => e,
        };
        match r.metric {
            Metric::Sobolev { quadrature } => {
                let gram =
                    sobolev_gram(&self.pool, r.network, r.positions, r.exponents, quadrature);
                let s = DenseSaddle::new(self.dense_pool.clone(), &gram.a, r.jacobian, r.layout)
                    .map_err(singular)?;
                Ok(Box::new(s))
            }
            Metric::L2 | Metric::H1 | Metric::H2 => {
                let m = classic_metric(r.network, r.positions, &r.metric)?;
                if !matches!(r.metric, Metric::L2) {
                    crate::constraints::check_component_translations(
                        r.jacobian, r.network, r.layout,
                    )?;
                }
                let s = SparseSaddle::new(&m, r.jacobian, r.layout).map_err(singular)?;
                Ok(Box::new(s))
            }
        }
    }
}

/// The solver back end (`SPEC.md` §8.5).
#[derive(Clone)]
pub enum Solver {
    /// [`DirectSolvers`].
    Dense,
    /// A caller-supplied back end (the multigrid of WP6 plugs in here).
    Custom(Arc<dyn SolverFactory>),
}

impl fmt::Debug for Solver {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Solver::Dense => write!(f, "Dense"),
            Solver::Custom(s) => write!(f, "Custom({})", s.name()),
        }
    }
}

/// Energy and solver back ends (`SPEC.md` §8.5).
#[derive(Debug, Clone)]
pub enum Backend {
    /// Both chosen explicitly.
    Explicit {
        /// The energy evaluation.
        energy: EnergyEval,
        /// The solver.
        solver: Solver,
    },
    /// Exact + dense up to `dense_max_vertices` vertices, Barnes–Hut + multigrid above. Until
    /// the accelerated back ends exist (WP4–WP6), exact + dense at every size (the report names
    /// the back end used).
    Auto {
        /// The crossover (placeholder 2000, measured in §14.1).
        dense_max_vertices: usize,
    },
}

impl Default for Backend {
    /// `Auto { dense_max_vertices: 2000 }` (§10).
    fn default() -> Self {
        Backend::Auto {
            dense_max_vertices: 2000,
        }
    }
}

impl Backend {
    /// The energy evaluation and solver used for a network of `n_vertices` vertices.
    pub fn resolve(&self, _n_vertices: usize) -> (EnergyEval, Solver) {
        match self {
            Backend::Explicit { energy, solver } => (*energy, solver.clone()),
            // The accelerated back ends are not available yet: exact + dense everywhere.
            Backend::Auto { .. } => (EnergyEval::Exact, Solver::Dense),
        }
    }
}

// =============================================================================================
// Options
// =============================================================================================

/// Options of a flow (`SPEC.md` §9.1, defaults §10).
#[derive(Clone)]
pub struct FlowOptions {
    /// `(α, β)`, default (3, 6).
    pub exponents: Exponents,
    /// The metric, default Sobolev with the midpoint quadrature.
    pub metric: Metric,
    /// Energy and solver back ends.
    pub backend: Backend,
    /// Projected gradient (hard constraints) or L-BFGS / NCG (soft problems).
    pub method: Method,
    /// Hard constraints, default `[Barycenter, EdgeLengths]`; must be empty for soft methods.
    pub constraints: Vec<Constraint>,
    /// Penalty terms added to the objective.
    pub penalties: Vec<Arc<dyn PenaltyTerm>>,
    /// Rigid obstacles for the CCD (besides those of the penalties).
    pub obstacles: Vec<Arc<dyn CollisionObstacle>>,
    /// Planar mode: the held axis (`SPEC.md` §2.2).
    pub planar: Option<Axis>,
    /// The line search.
    pub line_search: LineSearch,
    /// The constraint projection.
    pub projection: ProjectionOptions,
    /// The stop rule.
    pub stop: StopRule,
    /// Maximal number of steps of the session.
    pub max_iterations: usize,
    /// Wall-time limit of one [`Session::run`].
    pub time_limit: Option<Duration>,
    /// A length schedule (`None`: targets fixed).
    pub schedule: Option<LengthSchedule>,
    /// Maximal number of subdivisions (0: off, the default).
    pub max_subdivisions: u32,
    /// Size of the main pool (energy, Gram assembly, CCD); `None`: all cores.
    pub threads: Option<usize>,
    /// Size of the dense LU's pool; `None`: `min(DEFAULT_DENSE_THREADS, threads)`.
    pub dense_threads: Option<usize>,
    /// Keep the path of every accepted step in the report (for audits; memory heavy).
    pub record_path: bool,
}

impl fmt::Debug for FlowOptions {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("FlowOptions")
            .field("exponents", &self.exponents)
            .field("metric", &self.metric)
            .field("backend", &self.backend)
            .field("method", &self.method)
            .field("constraints", &self.constraints)
            .field("penalties", &self.penalties)
            .field("obstacles", &self.obstacles.len())
            .field("planar", &self.planar)
            .field("line_search", &self.line_search)
            .field("projection", &self.projection)
            .field("stop", &self.stop)
            .field("max_iterations", &self.max_iterations)
            .field("time_limit", &self.time_limit)
            .field("schedule", &self.schedule)
            .field("max_subdivisions", &self.max_subdivisions)
            .field("threads", &self.threads)
            .field("dense_threads", &self.dense_threads)
            .field("record_path", &self.record_path)
            .finish()
    }
}

impl Default for FlowOptions {
    /// The defaults of `SPEC.md` §10.
    fn default() -> Self {
        FlowOptions {
            exponents: Exponents::default(),
            metric: Metric::default(),
            backend: Backend::default(),
            method: Method::ProjectedGradient,
            constraints: default_constraints(),
            penalties: Vec::new(),
            obstacles: Vec::new(),
            planar: None,
            line_search: LineSearch::default(),
            projection: ProjectionOptions::default(),
            stop: StopRule::default(),
            max_iterations: 1000,
            time_limit: None,
            schedule: None,
            max_subdivisions: 0,
            threads: None,
            dense_threads: None,
            record_path: false,
        }
    }
}

impl FlowOptions {
    /// Checks the options that do not need the curve.
    ///
    /// # Errors
    /// [`RcError::InvalidOption`] naming the option; [`RcError::MethodNeedsSoftProblem`] for
    /// L-BFGS or NCG with hard constraints.
    pub fn validate(&self) -> Result<(), RcError> {
        self.line_search.validate()?;
        self.projection.validate()?;
        self.stop.validate()?;
        self.method.validate()?;
        if let Some(s) = &self.schedule {
            s.validate()?;
        }
        for (name, t) in [
            ("threads", self.threads),
            ("dense_threads", self.dense_threads),
        ] {
            if t == Some(0) {
                return Err(RcError::InvalidOption {
                    name,
                    value: "0".into(),
                    requirement: "at least 1 thread",
                });
            }
        }
        if self.method.is_soft() && !self.constraints.is_empty() {
            return Err(RcError::MethodNeedsSoftProblem);
        }
        Ok(())
    }
}

// =============================================================================================
// Reports
// =============================================================================================

/// Seconds per phase of a step.
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct Timings {
    /// The re-target sub-step (its metric, projection and CCD included).
    pub retarget: f64,
    /// `F` and `dF` at the base curve (energy back end and penalties).
    pub energy: f64,
    /// Metric assembly and factorization (the solver back end's build).
    pub metric: f64,
    /// The direction's solve(s).
    pub solve: f64,
    /// The whole line search.
    pub line_search: f64,
    /// Objective evaluations of the line search.
    pub line_search_energy: f64,
    /// Projections of the line search.
    pub projection: f64,
    /// CCD of the line search.
    pub ccd: f64,
    /// Subdivision (with the new constraint set and its rank checks).
    pub subdivision: f64,
    /// The whole step.
    pub total: f64,
}

/// How a re-target sub-step went.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum RetargetStatus {
    /// The curve was projected onto the new targets.
    Moved,
    /// Every attempt failed; the old targets are kept.
    Stalled,
}

/// The re-target sub-step of a step (§7.1 step 0).
#[derive(Debug, Clone, PartialEq)]
pub struct RetargetReport {
    /// The outcome.
    pub status: RetargetStatus,
    /// Schedule progress before.
    pub progress_from: f64,
    /// Schedule progress after (unchanged if stalled).
    pub progress_to: f64,
    /// Number of attempts (1 + halvings of the increment).
    pub attempts: u32,
    /// Chord iterations of the last attempt.
    pub projection_iterations: u32,
    /// Largest violation after the last attempt.
    pub max_violation: f64,
    /// Status of the last attempt's projection.
    pub projection_status: ProjectionStatus,
    /// Whether the last attempt's path had a contact.
    pub collision: bool,
}

/// L-BFGS / NCG details of a step.
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct SoftReport {
    /// Stored L-BFGS pairs.
    pub pairs: usize,
    /// L-BFGS pairs rejected by the curvature test so far.
    pub skipped_pairs: usize,
    /// NCG `β`.
    pub beta: f64,
    /// The direction fell back to the preconditioned gradient (restart).
    pub restarted: bool,
}

/// The report of one step (`SPEC.md` §9.1).
#[derive(Debug, Clone, PartialEq)]
pub struct StepReport {
    /// 1-based step number of the session.
    pub iteration: usize,
    /// The outcome.
    pub status: StepStatus,
    /// Vertices and edges at the start of the step.
    pub n_vertices: usize,
    /// Edges at the start of the step.
    pub n_edges: usize,
    /// `F` at the base curve (after a re-target).
    pub energy_before: f64,
    /// `F` after the step (`energy_before` if nothing was accepted).
    pub energy_after: f64,
    /// `Ê` after the step.
    pub tangent_point_after: f64,
    /// The penalty values after the step, in order.
    pub penalties_after: Vec<f64>,
    /// `‖g̃‖₂` (soft methods: of the preconditioned gradient `z`).
    pub gradient_norm: f64,
    /// `√(Σ m_v |g̃_v|²)`.
    pub gradient_mass_norm: f64,
    /// `dF · g̃` (projected gradient), or `dF · p` of the soft direction.
    pub sobolev_norm_sq: f64,
    /// Accepted step size (0 if none).
    pub tau: f64,
    /// First trial step size.
    pub tau_start: f64,
    /// `τ_hi` of the guarded search.
    pub tau_hi: Option<f64>,
    /// `τ_max` of the guarded search (`∞`: no contact up to `τ_hi`).
    pub tau_max: Option<f64>,
    /// Number of shrinks.
    pub halvings: u32,
    /// Every trial of the line search.
    pub trials: Vec<Trial>,
    /// Chord iterations of the accepted trial.
    pub projection_iterations: u32,
    /// Largest constraint violation after the step (its row's own unit).
    pub max_violation: f64,
    /// Constraint rows dropped as redundant (D9).
    pub dropped_rows: Vec<usize>,
    /// The direction's linear solve.
    pub linear: LinearReport,
    /// Number of linear solves of the line search's projections.
    pub projection_solves: usize,
    /// Largest relative residual of those solves.
    pub projection_residual: f64,
    /// CCD counters of the line search.
    pub ccd: CcdReport,
    /// Number of CCD queries of the line search.
    pub ccd_checks: u32,
    /// The re-target sub-step, if any.
    pub retarget: Option<RetargetReport>,
    /// The network was subdivided at the end of this step.
    pub subdivided: bool,
    /// Energy back end used.
    pub energy_backend: &'static str,
    /// Solver back end used.
    pub solver_backend: &'static str,
    /// Soft methods: L-BFGS / NCG details.
    pub soft: Option<SoftReport>,
    /// Seconds per phase.
    pub seconds: Timings,
    /// With `record_path`: the path of the accepted step (base, move, chord iterates).
    pub path: Option<Vec<Vec<[f64; 3]>>>,
}

/// What the progress callback of [`Session::run`] returns.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Control {
    /// Keep going.
    Continue,
    /// Stop the run ([`Status::Cancelled`]).
    Stop,
}

/// The report of a run.
#[derive(Debug, Clone, PartialEq)]
pub struct FlowReport {
    /// Why the run ended.
    pub status: Status,
    /// Every step of the run.
    pub steps: Vec<StepReport>,
    /// `F` at the end.
    pub final_energy: f64,
    /// Steps of the session so far.
    pub iterations: usize,
    /// Wall time of the run.
    pub seconds: f64,
}

/// The result of [`flow`].
#[derive(Debug, Clone)]
pub struct FlowResult {
    /// The final network (it changes with subdivision).
    pub network: CurveNetwork,
    /// The final positions.
    pub positions: Vec<[f64; 3]>,
    /// The run's report.
    pub report: FlowReport,
}

// =============================================================================================
// Session
// =============================================================================================

#[derive(Debug, Default)]
struct SoftState {
    lbfgs: Option<Lbfgs>,
    ncg: Option<Ncg>,
    previous: Option<(Vec<f64>, Vec<f64>)>,
}

impl SoftState {
    fn new(method: Method) -> SoftState {
        SoftState {
            lbfgs: match method {
                Method::Lbfgs { memory } => Some(Lbfgs::new(memory)),
                _ => None,
            },
            ncg: matches!(method, Method::Ncg).then(Ncg::new),
            previous: None,
        }
    }
}

/// A step-at-a-time flow with its state kept between steps (Blender's live mode,
/// `SPEC.md` §9.1). The session owns copies of the caller's data.
pub struct Session {
    network: CurveNetwork,
    positions: Vec<[f64; 3]>,
    options: FlowOptions,
    layout: DofLayout,
    constraints: ConstraintSet,
    penalties: Vec<Arc<dyn PenaltyTerm>>,
    obstacles: Vec<PreparedObstacle>,
    pool: Arc<rayon::ThreadPool>,
    energy: Box<dyn EnergyModel>,
    solver: Arc<dyn SolverFactory>,
    mean_edge_length0: f64,
    iteration: usize,
    status: Status,
    cancel: Arc<AtomicBool>,
    schedule: Option<ScheduleState>,
    subdivisions: u32,
    previous_tau: Option<f64>,
    stop: StopState,
    soft: SoftState,
    last_energy: f64,
    stalled_in_a_row: u32,
}

impl fmt::Debug for Session {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("Session")
            .field("n_vertices", &self.network.n_vertices())
            .field("iteration", &self.iteration)
            .field("status", &self.status)
            .field("energy", &self.last_energy)
            .finish_non_exhaustive()
    }
}

/// Evaluates `F` (and `dF`).
fn objective(
    energy: &mut dyn EnergyModel,
    penalties: &[Arc<dyn PenaltyTerm>],
    network: &CurveNetwork,
    positions: &[[f64; 3]],
) -> Result<(f64, f64, Vec<f64>), RcError> {
    let e = energy.energy(network, positions)?;
    let mut total = e;
    let mut values = Vec::with_capacity(penalties.len());
    for p in penalties {
        let v = p.evaluate(network, positions, None)?;
        total += v;
        values.push(v);
    }
    Ok((total, e, values))
}

struct StepEnv<'a> {
    network: &'a CurveNetwork,
    energy: &'a mut dyn EnergyModel,
    penalties: &'a [Arc<dyn PenaltyTerm>],
    solver: &'a mut dyn SaddleSolver,
    constraints: Option<&'a ConstraintSet>,
    projection: &'a ProjectionOptions,
    mean_edge_length0: f64,
    pool: &'a rayon::ThreadPool,
    obstacles: &'a [PreparedObstacle],
    ccd: CcdOptions,
    cancel: &'a AtomicBool,
}

impl SearchEnv for StepEnv<'_> {
    fn objective(&mut self, positions: &[[f64; 3]]) -> Result<f64, RcError> {
        Ok(objective(self.energy, self.penalties, self.network, positions)?.0)
    }

    fn project(
        &mut self,
        trial: Vec<[f64; 3]>,
    ) -> Result<line_search::ProjectionWithPath, RcError> {
        match self.constraints {
            Some(c) => project_with_path(
                self.solver,
                c,
                self.network,
                trial,
                self.projection,
                self.mean_edge_length0,
            ),
            None => Ok((
                Projection {
                    positions: trial.clone(),
                    status: ProjectionStatus::Converged,
                    iterations: 0,
                    max_violation: 0.0,
                    history: Vec::new(),
                    linear: Vec::new(),
                },
                vec![trial],
            )),
        }
    }

    fn ccd(
        &mut self,
        start: &[[f64; 3]],
        velocity: &[[f64; 3]],
        tau_end: f64,
    ) -> Result<CcdReport, RcError> {
        first_contact(
            self.pool,
            self.network,
            start,
            velocity,
            tau_end,
            self.obstacles,
            &self.ccd,
        )
    }

    fn cancelled(&self) -> bool {
        self.cancel.load(Ordering::Relaxed)
    }
}

fn prepare_obstacles(
    options: &FlowOptions,
    penalties: &[Arc<dyn PenaltyTerm>],
) -> Result<Vec<PreparedObstacle>, RcError> {
    options
        .obstacles
        .iter()
        .cloned()
        .chain(penalties.iter().filter_map(|p| p.obstacle()))
        .map(PreparedObstacle::new)
        .collect()
}

impl Session {
    /// Starts a session: validates the options, the curve (§2.3, planar mode), the constraints
    /// (feasibility at the start, §5.1; translations and rank, §5.2) and the penalties, and
    /// sets up the pools and the back ends.
    ///
    /// # Errors
    /// Every validation error of the crate; [`RcError::NonFiniteEnergy`] if `F` is not finite
    /// at the start.
    pub fn new(
        network: CurveNetwork,
        positions: Vec<[f64; 3]>,
        options: FlowOptions,
    ) -> Result<Session, RcError> {
        options.validate()?;
        network.validate_positions(&positions)?;
        if let Some(axis) = options.planar {
            check_planar(&positions, axis)?;
        }
        let pool = Arc::new(make_pool(options.threads)?);
        let n_threads = pool.current_num_threads();
        let dense_pool = Arc::new(make_pool(Some(
            options
                .dense_threads
                .unwrap_or(DEFAULT_DENSE_THREADS.min(n_threads)),
        ))?);
        let layout = DofLayout::new(options.planar);
        let mean_edge_length0 = network.mean_edge_length(&positions);
        let constraints =
            ConstraintSet::new(&network, &positions, &options.constraints, options.planar)?;
        if !options.method.is_soft() {
            constraints.check_feasible(
                &network,
                &positions,
                options.projection.tolerance,
                mean_edge_length0,
            )?;
            constraints.analyze(&network, &positions)?;
        }
        let schedule = match options.schedule {
            Some(s) => Some(ScheduleState::new(s, &constraints)?),
            None => None,
        };
        let penalties = options.penalties.clone();
        for p in &penalties {
            p.check(&network)?;
        }
        let obstacles = prepare_obstacles(&options, &penalties)?;
        let (eval, solver) = options.backend.resolve(network.n_vertices());
        let energy = energy_model(&eval, options.exponents, pool.clone())?;
        let solver: Arc<dyn SolverFactory> = match solver {
            Solver::Dense => Arc::new(DirectSolvers::new(pool.clone(), dense_pool)),
            Solver::Custom(f) => f,
        };
        let mut s = Session {
            soft: SoftState::new(options.method),
            network,
            positions,
            layout,
            constraints,
            penalties,
            obstacles,
            pool,
            energy,
            solver,
            mean_edge_length0,
            iteration: 0,
            status: Status::Running,
            cancel: Arc::new(AtomicBool::new(false)),
            schedule,
            subdivisions: 0,
            previous_tau: None,
            stop: StopState::default(),
            last_energy: f64::NAN,
            stalled_in_a_row: 0,
            options,
        };
        let (f, _, _) = objective(s.energy.as_mut(), &s.penalties, &s.network, &s.positions)?;
        if !f.is_finite() {
            return Err(RcError::NonFiniteEnergy { step: 0 });
        }
        s.last_energy = f;
        Ok(s)
    }

    /// The current positions.
    pub fn positions(&self) -> &[[f64; 3]] {
        &self.positions
    }

    /// The current network (changes with subdivision).
    pub fn network(&self) -> &CurveNetwork {
        &self.network
    }

    /// The current constraint set (targets move with a schedule).
    pub fn constraints(&self) -> &ConstraintSet {
        &self.constraints
    }

    /// The current penalty terms (they follow subdivisions).
    pub fn penalties(&self) -> &[Arc<dyn PenaltyTerm>] {
        &self.penalties
    }

    /// The options.
    pub fn options(&self) -> &FlowOptions {
        &self.options
    }

    /// The status after the last step.
    pub fn status(&self) -> Status {
        self.status.clone()
    }

    /// Steps taken so far.
    pub fn iteration(&self) -> usize {
        self.iteration
    }

    /// `F` at the current curve.
    pub fn energy(&self) -> f64 {
        self.last_energy
    }

    /// `ℓ̄₀`, the mean edge length at the start (unit of the relative tolerances).
    pub fn mean_edge_length0(&self) -> f64 {
        self.mean_edge_length0
    }

    /// Number of subdivisions so far.
    pub fn subdivisions(&self) -> u32 {
        self.subdivisions
    }

    /// The schedule's progress (`None` without a schedule).
    pub fn schedule_progress(&self) -> Option<f64> {
        self.schedule.as_ref().map(|s| s.progress())
    }

    /// The cancel flag, checked before every step and every line-search trial.
    pub fn cancel_flag(&self) -> Arc<AtomicBool> {
        self.cancel.clone()
    }

    /// Sets (or removes) the length schedule; its bases are the current targets.
    ///
    /// # Errors
    /// [`RcError::ScheduleWithoutLengthConstraint`]; [`RcError::InvalidOption`].
    pub fn set_schedule(&mut self, schedule: Option<LengthSchedule>) -> Result<(), RcError> {
        self.schedule = match schedule {
            Some(s) => Some(ScheduleState::new(s, &self.constraints)?),
            None => None,
        };
        self.options.schedule = schedule;
        Ok(())
    }

    fn schedule_moving(&self) -> bool {
        self.schedule.as_ref().is_some_and(|s| !s.finished())
    }

    fn build_solver(
        &self,
        jacobian: &SparseMatrix,
        hard: bool,
    ) -> Result<Box<dyn SaddleSolver>, RcError> {
        self.solver.build(&SolverRequest {
            network: &self.network,
            positions: &self.positions,
            exponents: self.options.exponents,
            metric: self.options.metric,
            jacobian,
            layout: self.layout,
            constraints: hard.then_some(&self.constraints),
        })
    }

    fn ccd_options(&self) -> CcdOptions {
        self.options.line_search.ccd_options()
    }

    /// Step 0 of §7.1.
    fn retarget(&mut self) -> Result<Option<RetargetReport>, RcError> {
        if self.options.method.is_soft() || self.constraints.n_active_rows() == 0 {
            return Ok(None);
        }
        let moving = self.schedule_moving();
        let violated = {
            let phi = self.constraints.values(&self.network, &self.positions)?;
            let tol = self.constraints.row_tolerances(
                &self.network,
                &self.positions,
                self.options.projection.tolerance,
                self.mean_edge_length0,
            )?;
            phi.iter()
                .zip(&tol)
                .any(|(p, t)| p.is_nan() || p.abs() > *t)
        };
        if !moving && !violated {
            return Ok(None);
        }
        let jac = self.constraints.jacobian(&self.network, &self.positions)?;
        let mut solver = self.build_solver(&jac, true)?;
        let progress_from = self.schedule.as_ref().map_or(0.0, |s| s.progress());
        let attempts_max = if moving { 6 } else { 1 };
        let guarded = self.options.line_search.is_guarded();
        let mut report = RetargetReport {
            status: RetargetStatus::Stalled,
            progress_from,
            progress_to: progress_from,
            attempts: 0,
            projection_iterations: 0,
            max_violation: 0.0,
            projection_status: ProjectionStatus::MaxIterations,
            collision: false,
        };
        for h in 0..attempts_max {
            if self.cancel.load(Ordering::Relaxed) {
                return Err(RcError::Cancelled);
            }
            report.attempts += 1;
            let mut cs = self.constraints.clone();
            let progress = match &self.schedule {
                Some(s) if moving => {
                    let p = s.next_progress(h);
                    s.apply(&mut cs, p)?;
                    p
                }
                _ => progress_from,
            };
            let (proj, path) = project_with_path(
                solver.as_mut(),
                &cs,
                &self.network,
                self.positions.clone(),
                &self.options.projection,
                self.mean_edge_length0,
            )?;
            report.projection_iterations = proj.iterations;
            report.max_violation = proj.max_violation;
            report.projection_status = proj.status;
            report.collision = false;
            if proj.status != ProjectionStatus::Converged {
                continue;
            }
            if guarded {
                let mut ccd = CcdReport::default();
                let hit = crate::ccd::path_contact(
                    &self.pool,
                    &self.network,
                    &path,
                    &self.obstacles,
                    &self.ccd_options(),
                    &mut ccd,
                )?;
                if hit.is_some() {
                    report.collision = true;
                    continue;
                }
            }
            self.constraints = cs;
            self.positions = proj.positions;
            if let Some(s) = self.schedule.as_mut() {
                s.commit(progress);
            }
            report.status = RetargetStatus::Moved;
            report.progress_to = progress;
            return Ok(Some(report));
        }
        Ok(Some(report))
    }

    /// One step of §7.1. The session's status is updated (converged, stuck, max iterations).
    ///
    /// # Errors
    /// [`RcError::Cancelled`] (the curve is that of the last completed sub-step);
    /// [`RcError::NonFiniteEnergy`]; back-end and constraint errors.
    pub fn step(&mut self) -> Result<StepReport, RcError> {
        let r = self.step_inner();
        match &r {
            Err(RcError::Cancelled) => self.status = Status::Cancelled,
            Err(e) => self.status = Status::Error(e.clone()),
            Ok(_) => {}
        }
        r
    }

    fn step_inner(&mut self) -> Result<StepReport, RcError> {
        let t_step = Instant::now();
        if self.cancel.load(Ordering::Relaxed) {
            return Err(RcError::Cancelled);
        }
        let iteration = self.iteration + 1;
        let mut seconds = Timings::default();

        // 0. Re-target.
        let t = Instant::now();
        let retarget = self.retarget()?;
        seconds.retarget = t.elapsed().as_secs_f64();
        let retarget_moved = retarget
            .as_ref()
            .is_some_and(|r| r.status == RetargetStatus::Moved);
        let retarget_stalled = retarget
            .as_ref()
            .is_some_and(|r| r.status == RetargetStatus::Stalled);

        // 1. F, dF at the base curve.
        let t = Instant::now();
        let (n_vertices, n_edges) = (self.network.n_vertices(), self.network.n_edges());
        self.energy.prepare(&self.network, &self.positions)?;
        let ev = self
            .energy
            .energy_and_differential(&self.network, &self.positions)?;
        let mut df = ev.differential;
        let mut f0 = ev.energy;
        let mut penalties0 = Vec::with_capacity(self.penalties.len());
        for p in &self.penalties {
            let v = p.evaluate(&self.network, &self.positions, Some(&mut df))?;
            f0 += v;
            penalties0.push(v);
        }
        seconds.energy = t.elapsed().as_secs_f64();
        if !f0.is_finite() || df.iter().flatten().any(|x| !x.is_finite()) {
            return Err(RcError::NonFiniteEnergy { step: iteration });
        }

        // 2. Metric, solver, direction.
        let soft = self.options.method.is_soft();
        let t = Instant::now();
        let jac = if soft {
            mean_free_rows(&self.network, &self.positions, self.layout)?
        } else {
            self.constraints.jacobian(&self.network, &self.positions)?
        };
        let mut solver = self.build_solver(&jac, !soft)?;
        seconds.metric = t.elapsed().as_secs_f64();
        let t = Instant::now();
        let masses = self.network.dual_lengths(&self.positions);
        let (direction, slope, gradient_norm, gradient_mass_norm, linear, soft_report) = if !soft {
            let pg = projected_gradient(solver.as_mut(), self.layout, &df)?;
            let mn = mass_norm(&pg.direction, &masses);
            (
                pg.direction,
                pg.sobolev_norm_sq,
                pg.gradient_norm,
                mn,
                pg.linear,
                None,
            )
        } else {
            self.soft_direction(solver.as_mut(), &df, &masses)?
        };
        seconds.solve = t.elapsed().as_secs_f64();

        let mut report = StepReport {
            iteration,
            status: StepStatus::Converged,
            n_vertices,
            n_edges,
            energy_before: f0,
            energy_after: f0,
            tangent_point_after: ev.energy,
            penalties_after: penalties0,
            gradient_norm,
            gradient_mass_norm,
            sobolev_norm_sq: slope,
            tau: 0.0,
            tau_start: 0.0,
            tau_hi: None,
            tau_max: None,
            halvings: 0,
            trials: Vec::new(),
            projection_iterations: 0,
            max_violation: 0.0,
            dropped_rows: self.constraints.dropped_rows().to_vec(),
            linear,
            projection_solves: 0,
            projection_residual: 0.0,
            ccd: CcdReport::default(),
            ccd_checks: 0,
            retarget,
            subdivided: false,
            energy_backend: self.energy.name(),
            solver_backend: self.solver.name(),
            soft: soft_report,
            seconds,
            path: None,
        };

        if !self.schedule_moving() {
            self.stop.seed_energy(f0);
        }
        let converged =
            self.stop
                .gradient_converged(&self.options.stop, gradient_norm, gradient_mass_norm);
        let mut accepted = false;
        if !converged {
            // 3. Line search.
            let t = Instant::now();
            let max_start = match (&self.options.method, &self.soft.lbfgs) {
                (Method::Lbfgs { .. }, Some(l)) if !l.is_empty() => Some(1.0),
                _ => None,
            };
            let mut env = StepEnv {
                network: &self.network,
                energy: self.energy.as_mut(),
                penalties: &self.penalties,
                solver: solver.as_mut(),
                constraints: (!soft).then_some(&self.constraints),
                projection: &self.options.projection,
                mean_edge_length0: self.mean_edge_length0,
                pool: &self.pool,
                obstacles: &self.obstacles,
                ccd: self.options.line_search.ccd_options(),
                cancel: &self.cancel,
            };
            let ls: LineSearchReport = line_search(
                &self.options.line_search,
                SearchInput {
                    base: &self.positions,
                    direction: &direction,
                    slope,
                    f0,
                    previous_tau: self.previous_tau,
                    max_start,
                },
                &mut env,
            )?;
            report.seconds.line_search = t.elapsed().as_secs_f64();
            report.seconds.line_search_energy = ls.seconds.energy;
            report.seconds.projection = ls.seconds.projection;
            report.seconds.ccd = ls.seconds.ccd;
            report.tau_start = ls.tau_start;
            report.tau_hi = ls.tau_hi;
            report.tau_max = ls.tau_max;
            report.halvings = ls.halvings;
            report.projection_solves = ls.projection_linear.len();
            report.projection_residual = ls
                .projection_linear
                .iter()
                .map(|l| l.relative_residual)
                .fold(0.0, f64::max);
            report.ccd = ls.ccd.clone();
            report.ccd_checks = ls.ccd_checks;
            report.max_violation = ls.max_violation;
            report.trials = ls.trials;
            // 4. Accept or stuck.
            if let (true, Some(new_pos)) = (ls.accepted, ls.positions) {
                accepted = true;
                let x_old = self.layout.gather(&self.positions);
                self.positions = new_pos;
                self.previous_tau = Some(ls.tau);
                report.tau = ls.tau;
                report.projection_iterations = ls.projection_iterations;
                if self.options.record_path {
                    report.path = Some(ls.path);
                }
                if soft {
                    self.soft.previous = Some((x_old, self.layout.gather(&df)));
                }
                let (f, e, pen) = objective(
                    self.energy.as_mut(),
                    &self.penalties,
                    &self.network,
                    &self.positions,
                )?;
                report.energy_after = f;
                report.tangent_point_after = e;
                report.penalties_after = pen;
                report.status = StepStatus::Accepted;
            } else {
                report.status = StepStatus::Stuck;
                self.soft.previous = None;
            }
        }
        if !accepted && retarget_stalled {
            report.status = StepStatus::RetargetStalled;
        }
        if !soft && !accepted {
            // The curve did not move in the line search: report its own violation.
            let phi = self.constraints.values(&self.network, &self.positions)?;
            report.max_violation = phi.iter().fold(0.0, |m: f64, p| m.max(p.abs()));
        }
        self.last_energy = report.energy_after;

        // Status.
        self.status = Status::Running;
        let moving = self.schedule_moving();
        match report.status {
            StepStatus::Converged if !moving => self.status = Status::Converged,
            StepStatus::Stuck if !moving => self.status = Status::Stuck,
            _ => {}
        }
        if report.status == StepStatus::RetargetStalled {
            self.stalled_in_a_row += 1;
            if self.stalled_in_a_row >= 5 {
                self.status = Status::Stuck;
            }
        } else {
            self.stalled_in_a_row = 0;
        }
        if accepted
            && !moving
            && !retarget_moved
            && self
                .stop
                .energy_stalled(&self.options.stop, report.energy_after)
        {
            self.status = Status::Converged;
        }

        // 5. Subdivision.
        let t = Instant::now();
        if self.subdivisions < self.options.max_subdivisions
            && self.network.mean_edge_length(&self.positions) > 2.0 * self.mean_edge_length0
        {
            self.subdivide_now()?;
            report.subdivided = true;
            self.last_energy = objective(
                self.energy.as_mut(),
                &self.penalties,
                &self.network,
                &self.positions,
            )?
            .0;
        }
        report.seconds.subdivision = t.elapsed().as_secs_f64();

        self.iteration = iteration;
        if self.status == Status::Running && self.iteration >= self.options.max_iterations {
            self.status = Status::MaxIterations;
        }
        report.seconds.total = t_step.elapsed().as_secs_f64();
        Ok(report)
    }

    /// The soft direction (§6.4): `z = M⁻¹ r` (mean-free solve), then L-BFGS or NCG.
    #[allow(clippy::type_complexity)]
    fn soft_direction(
        &mut self,
        solver: &mut dyn SaddleSolver,
        df: &[[f64; 3]],
        masses: &[f64],
    ) -> Result<
        (
            Vec<[f64; 3]>,
            f64,
            f64,
            f64,
            LinearReport,
            Option<SoftReport>,
        ),
        RcError,
    > {
        let layout = self.layout;
        let r = layout.gather(df);
        let x = layout.gather(&self.positions);
        let k = solver.n_constraints();
        let zero = vec![0.0; k];
        let sol = solver.solve(&r, &zero)?;
        let z = sol.x;
        let linear = sol.report;
        let zn = dotv(&z, &z).sqrt();
        let zmass = mass_norm(&layout.scatter(&z), masses);
        let mut rep = SoftReport::default();
        let mut p = match self.options.method {
            Method::Lbfgs { .. } => {
                let l = self.soft.lbfgs.get_or_insert_with(|| Lbfgs::new(30));
                if let Some((px, pr)) = self.soft.previous.take() {
                    let s: Vec<f64> = x.iter().zip(&px).map(|(a, b)| a - b).collect();
                    let y: Vec<f64> = r.iter().zip(&pr).map(|(a, b)| a - b).collect();
                    l.push(s, y);
                }
                let mut pre = |q: &[f64]| solver.solve(q, &zero).map(|s| s.x);
                let p = l.direction(&r, &mut pre)?;
                rep.pairs = l.len();
                rep.skipped_pairs = l.skipped();
                p
            }
            _ => {
                let n = self.soft.ncg.get_or_insert_with(Ncg::new);
                self.soft.previous = None;
                let d = n.direction(&r, &z);
                rep.beta = d.beta;
                rep.restarted = d.restarted;
                d.p
            }
        };
        let mut slope = dotv(&r, &p);
        if !(slope > 0.0 && slope.is_finite()) {
            if let Some(l) = self.soft.lbfgs.as_mut() {
                l.clear();
            }
            if let Some(n) = self.soft.ncg.as_mut() {
                n.clear();
            }
            p = z.clone();
            slope = dotv(&r, &p);
            rep.restarted = true;
        }
        Ok((layout.scatter(&p), slope, zn, zmass, linear, Some(rep)))
    }

    /// Splits every edge (§5.4) and maps the constraints, penalties and schedule.
    fn subdivide_now(&mut self) -> Result<(), RcError> {
        let (net, pos, map) = subdivide(&self.network, &self.positions)?;
        let (cons, index) =
            subdivide_constraints(self.constraints.constraints(), &self.network, &map);
        let cs = ConstraintSet::new(&net, &pos, &cons, self.options.planar)?;
        if !self.options.method.is_soft() {
            cs.analyze(&net, &pos)?;
        }
        let penalties = self
            .penalties
            .iter()
            .map(|p| p.subdivided(&map))
            .collect::<Result<Vec<_>, _>>()?;
        for p in &penalties {
            p.check(&net)?;
        }
        let schedule = self.schedule.as_ref().map(|s| s.subdivided(&map, &index));
        self.network = net;
        self.positions = pos;
        self.constraints = cs;
        self.penalties = penalties;
        self.schedule = schedule;
        self.subdivisions += 1;
        self.soft = SoftState::new(self.options.method);
        self.stop.reset();
        Ok(())
    }

    /// Steps until a stop rule, `max_iterations`, the time limit, a stuck line search (without
    /// a moving schedule), cancellation (flag or `progress` returning [`Control::Stop`]) or an
    /// error. Errors end the run with [`Status::Error`]; the curve is that of the last
    /// completed step.
    ///
    /// # Errors
    /// None at present (errors of steps are reported in the status); the `Result` is kept for
    /// the API of `SPEC.md` §9.1.
    pub fn run(
        &mut self,
        mut progress: Option<&mut dyn FnMut(&StepReport) -> Control>,
    ) -> Result<FlowReport, RcError> {
        let t0 = Instant::now();
        let mut steps = Vec::new();
        let status = loop {
            if self.cancel.load(Ordering::Relaxed) {
                self.status = Status::Cancelled;
                break Status::Cancelled;
            }
            if self.iteration >= self.options.max_iterations {
                self.status = Status::MaxIterations;
                break Status::MaxIterations;
            }
            if self.options.time_limit.is_some_and(|l| t0.elapsed() >= l) {
                self.status = Status::TimeLimit;
                break Status::TimeLimit;
            }
            match self.step() {
                Ok(rep) => {
                    let stop = progress.as_mut().is_some_and(|p| p(&rep) == Control::Stop);
                    steps.push(rep);
                    if self.status.is_terminal() {
                        break self.status.clone();
                    }
                    if stop {
                        self.status = Status::Cancelled;
                        break Status::Cancelled;
                    }
                }
                Err(_) => break self.status.clone(),
            }
        };
        Ok(FlowReport {
            status,
            steps,
            final_energy: self.last_energy,
            iterations: self.iteration,
            seconds: t0.elapsed().as_secs_f64(),
        })
    }
}

/// [`Session::new`] and [`Session::run`] in one call.
///
/// # Errors
/// The errors of [`Session::new`].
pub fn flow(
    network: CurveNetwork,
    positions: Vec<[f64; 3]>,
    options: FlowOptions,
    progress: Option<&mut dyn FnMut(&StepReport) -> Control>,
) -> Result<FlowResult, RcError> {
    let mut s = Session::new(network, positions, options)?;
    let report = s.run(progress)?;
    Ok(FlowResult {
        network: s.network,
        positions: s.positions,
        report,
    })
}

/// Builds [`FlowOptions`] for the soft problem of the evaluation (§6.4, D21): `method`,
/// no hard constraints, the Hencky penalty with weight `hencky_weight` at the current lengths.
///
/// # Errors
/// The errors of [`crate::optimize::Hencky::at_rest`].
pub fn soft_options(
    method: Method,
    hencky_weight: f64,
    network: &CurveNetwork,
    positions: &[[f64; 3]],
) -> Result<FlowOptions, RcError> {
    let h = crate::optimize::Hencky::at_rest(hencky_weight, network, positions)?;
    Ok(FlowOptions {
        method,
        constraints: Vec::new(),
        penalties: vec![Arc::new(h)],
        ..FlowOptions::default()
    })
}
