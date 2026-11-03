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
//!    the mapping of [`mod@crate::subdivide`]. The vertex count of the session changes.
//!
//! Every part is timed in the [`StepReport`].
//!
//! # Back ends
//!
//! * Energy: [`EnergyModel`], chosen by [`energy_model`] from [`EnergyEval`].
//!   [`EnergyEval::Exact`] is [`ExactEnergy`]; [`EnergyEval::BarnesHut`] is [`BarnesHutEnergy`]
//!   `prepare` is called once per step at the base curve
//!   (Barnes–Hut: build the tree), the line search's trials only refit (D26).
//! * Metric and solver: [`SolverFactory`], one [`SaddleSolver`] per step and per re-target.
//!   [`Solver::Dense`] is [`DirectSolvers`] (dense LU for the Sobolev metric, sparse LU for
//!   L², H¹, H²); [`Solver::Multigrid`] is the multigrid solver of [`crate::multigrid`] (WP6);
//!   [`Solver::Custom`] takes any factory. Factories run on the pools of the
//!   [`SolverRequest`]. An iterative solver that fails ([`RcError::MultigridDiverged`],
//!   [`RcError::MultigridNotConverged`]) falls back to the dense solver for the rest of the step
//!   when the network is small enough ([`fallback`], [`StepReport::solver_fallback`]).
//! * [`Backend::Auto`] picks exact + dense below a measured crossover and Barnes–Hut +
//!   multigrid above (`SPEC.md` §17h).
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
use crate::constraints::{Constraint, ConstraintSet, ResolvedConstraint, default_constraints};
use crate::energy::barnes_hut::{BhTree, InteractionLists};
use crate::energy::exact::{energy_and_differential, energy_only, make_pool};
use crate::energy::near_pair::{self, NearPairStats};
use crate::energy::{BhOptions, EnergyEval, NearPairOptions};
use crate::error::RcError;
use crate::kernel::Exponents;
use crate::metric::{DofLayout, Metric, classic_metric, sobolev_gram};
use crate::multigrid::{MgOptions, MultigridSolvers};
use crate::network::{Axis, CurveNetwork, check_planar};
use crate::optimize::{Lbfgs, Method, Ncg, dotv, mean_free_rows};
use crate::solve::{
    DenseSaddle, LinearReport, Projection, ProjectionOptions, ProjectionStatus, SaddleSolver,
    SparseSaddle, projected_gradient, second_order_correction,
};
use crate::subdivide::{SubdivisionMap, subdivide, subdivide_constraints};

pub mod fallback;
pub mod line_search;
pub mod schedule;
pub mod stop;

pub use fallback::SolverFallback;
use fallback::{FallbackSaddle, FallbackSlot};
use line_search::{
    LineSearch, LineSearchReport, SearchEnv, SearchInput, StartRule, Trial, TrialOutcome,
    line_search, project_with_path_polished,
};
pub use schedule::ScheduleProgress;
use schedule::{LengthSchedule, ScheduleState};
use stop::{Status, StepStatus, StopRule, StopState, mass_norm};

/// Default size of the dense solver's pool: `faer`'s LU stops scaling past 16–32 threads (WP2,
/// `SPEC.md` §17b; measured again by WP3, see the crate notes).
pub const DEFAULT_DENSE_THREADS: usize = 16;

/// Default of [`FlowOptions::min_edge_length_ratio`]: 1e-3.
pub const DEFAULT_MIN_EDGE_LENGTH_RATIO: f64 = 1e-3;

/// Default of [`FlowOptions::dense_fallback_max_vertices`].
pub const DEFAULT_DENSE_FALLBACK_MAX_VERTICES: usize = 5000;

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
    /// Edge pairs in the band of the near-pair quadrature (`SPEC.md` §3.4; 0 without it).
    pub near_pairs: u64,
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
    /// Seconds spent in the near-pair quadrature (`SPEC.md` §3.4) since the last call, which
    /// resets the count. The default: 0.
    fn take_near_pair_seconds(&mut self) -> f64 {
        0.0
    }
}

/// Adds the near-pair correction (`SPEC.md` §3.4) of the ordered pairs `include` accepts to
/// `energy` (and `differential`); unchanged, bit for bit, when no pair is near. Returns the
/// statistics and the seconds spent.
#[allow(clippy::too_many_arguments)]
fn add_near_pair(
    pool: &rayon::ThreadPool,
    network: &CurveNetwork,
    positions: &[[f64; 3]],
    exponents: Exponents,
    options: &NearPairOptions,
    include: &(dyn Fn(usize, usize) -> bool + Sync),
    energy: &mut f64,
    differential: Option<&mut [[f64; 3]]>,
) -> (NearPairStats, f64) {
    let t = Instant::now();
    let with_gradient = differential.is_some();
    let Some(near_pair::NearCorrection {
        energy: c,
        differential: g,
        stats,
    }) = near_pair::correction(
        pool,
        network,
        positions,
        exponents,
        options,
        include,
        with_gradient,
    )
    else {
        return (NearPairStats::default(), t.elapsed().as_secs_f64());
    };
    *energy += c;
    if let (Some(d), Some(g)) = (differential, g) {
        for (dv, gv) in d.iter_mut().zip(&g) {
            *dv = [dv[0] + gv[0], dv[1] + gv[1], dv[2] + gv[2]];
        }
    }
    (stats, t.elapsed().as_secs_f64())
}

/// The exact `O(|E|²)` energy (`SPEC.md` §3) on an explicit pool, with the near-pair
/// quadrature of §3.4 if set ([`ExactEnergy::with_near_pair`]).
#[derive(Debug, Clone)]
pub struct ExactEnergy {
    pool: Arc<rayon::ThreadPool>,
    exponents: Exponents,
    near_pair: Option<NearPairOptions>,
    near_seconds: f64,
}

impl ExactEnergy {
    /// The exact energy with `exponents` on `pool` (the paper's: no near-pair quadrature).
    pub fn new(pool: Arc<rayon::ThreadPool>, exponents: Exponents) -> ExactEnergy {
        ExactEnergy {
            pool,
            exponents,
            near_pair: None,
            near_seconds: 0.0,
        }
    }

    /// The same energy with the near-pair quadrature `near_pair` (`SPEC.md` §3.4; `None`: the
    /// paper's energy).
    pub fn with_near_pair(mut self, near_pair: Option<NearPairOptions>) -> ExactEnergy {
        self.near_pair = near_pair;
        self
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
        let mut e = energy_only(&self.pool, network, positions, self.exponents);
        if let Some(o) = &self.near_pair {
            let all = |_: usize, _: usize| true;
            let (_, sec) = add_near_pair(
                &self.pool,
                network,
                positions,
                self.exponents,
                o,
                &all,
                &mut e,
                None,
            );
            self.near_seconds += sec;
        }
        Ok(e)
    }

    fn energy_and_differential(
        &mut self,
        network: &CurveNetwork,
        positions: &[[f64; 3]],
    ) -> Result<EnergyEvaluation, RcError> {
        let (mut energy, mut differential) =
            energy_and_differential(&self.pool, network, positions, self.exponents);
        let mut near_pairs = 0;
        if let Some(o) = &self.near_pair {
            let all = |_: usize, _: usize| true;
            let (stats, sec) = add_near_pair(
                &self.pool,
                network,
                positions,
                self.exponents,
                o,
                &all,
                &mut energy,
                Some(&mut differential),
            );
            self.near_seconds += sec;
            near_pairs = stats.pairs;
        }
        Ok(EnergyEvaluation {
            energy,
            differential,
            pairs_exact: exact_pairs(network),
            nodes_lumped: 0,
            near_pairs,
        })
    }

    fn take_near_pair_seconds(&mut self) -> f64 {
        std::mem::take(&mut self.near_seconds)
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
    /// The interaction lists of the base curve, used for the trial energies of the line search.
    lists: Option<InteractionLists>,
    /// Vertex and edge count the tree was built for (a different network needs a rebuild).
    shape: (usize, usize),
    near_pair: Option<NearPairOptions>,
    near_seconds: f64,
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
            lists: None,
            shape: (0, 0),
            near_pair: None,
            near_seconds: 0.0,
        }
    }

    /// The same energy with the near-pair quadrature `near_pair` (`SPEC.md` §3.4; `None`: the
    /// paper's). The correction applies to the ordered pairs `(I, J)` with `J` among the exactly
    /// evaluated edges of row `I` in the frozen interaction lists (near pairs are not lumped by
    /// the default parameters; if one were, its lumped term stays), so the energy of a line
    /// search stays one smooth function and the differential its exact derivative.
    pub fn with_near_pair(mut self, near_pair: Option<NearPairOptions>) -> BarnesHutEnergy {
        self.near_pair = near_pair;
        self
    }

    /// The current tree (after [`EnergyModel::prepare`]).
    pub fn tree(&self) -> Option<&BhTree> {
        self.tree.as_ref()
    }

    fn build(&mut self, network: &CurveNetwork, positions: &[[f64; 3]]) -> Result<(), RcError> {
        let tree = BhTree::build(network, positions, self.options)?;
        self.lists = Some(tree.interaction_lists(&self.pool));
        self.tree = Some(tree);
        self.shape = (network.n_vertices(), network.n_edges());
        Ok(())
    }

    /// Refits the tree to `positions`; builds it if `prepare` was not called, or if the network
    /// is not the one the tree was built for (after a subdivision, when the energy is asked for
    /// before the next step's `prepare`).
    fn fit(&mut self, network: &CurveNetwork, positions: &[[f64; 3]]) -> Result<(), RcError> {
        let same = self.shape == (network.n_vertices(), network.n_edges());
        match &mut self.tree {
            Some(t) if same => t.refit(positions)?,
            _ => self.build(network, positions)?,
        }
        Ok(())
    }
}

impl EnergyModel for BarnesHutEnergy {
    fn name(&self) -> &'static str {
        "barnes-hut"
    }

    fn prepare(&mut self, network: &CurveNetwork, positions: &[[f64; 3]]) -> Result<(), RcError> {
        self.build(network, positions)
    }

    fn energy(&mut self, network: &CurveNetwork, positions: &[[f64; 3]]) -> Result<f64, RcError> {
        let (pool, exponents) = (self.pool.clone(), self.exponents);
        self.fit(network, positions)?;
        // The base curve's interaction lists (rebuilt together with the tree): a smooth energy
        // for the Armijo comparison, see `BhTree::interaction_lists`.
        let (tree, lists) = (self.tree.as_ref(), self.lists.as_ref());
        let (tree, lists) = (tree.expect("fitted"), lists.expect("fitted"));
        let mut e = tree.energy_with_lists(lists, &pool, exponents);
        if let Some(o) = &self.near_pair {
            let include = |i: usize, j: usize| lists.contains_exact(i, j);
            let (_, sec) = add_near_pair(
                &pool, network, positions, exponents, o, &include, &mut e, None,
            );
            self.near_seconds += sec;
        }
        Ok(e)
    }

    fn energy_and_differential(
        &mut self,
        network: &CurveNetwork,
        positions: &[[f64; 3]],
    ) -> Result<EnergyEvaluation, RcError> {
        let (pool, exponents) = (self.pool.clone(), self.exponents);
        self.fit(network, positions)?;
        let tree = self.tree.as_ref().expect("fitted");
        // The differential is taken at the base curve of a step: its interaction lists are the
        // ones the trial energies of the line search use.
        let lists = tree.interaction_lists(&pool);
        let (mut energy, mut differential, stats) =
            tree.energy_and_differential_with_lists(&lists, &pool, exponents)?;
        let mut near_pairs = 0;
        if let Some(o) = &self.near_pair {
            let include = |i: usize, j: usize| lists.contains_exact(i, j);
            let (st, sec) = add_near_pair(
                &pool,
                network,
                positions,
                exponents,
                o,
                &include,
                &mut energy,
                Some(&mut differential),
            );
            self.near_seconds += sec;
            near_pairs = st.pairs;
        }
        self.lists = Some(lists);
        Ok(EnergyEvaluation {
            energy,
            differential,
            pairs_exact: stats.pairs_exact,
            nodes_lumped: stats.nodes_lumped,
            near_pairs,
        })
    }

    fn take_near_pair_seconds(&mut self) -> f64 {
        std::mem::take(&mut self.near_seconds)
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
    energy_model_with_near_pair(eval, exponents, None, pool)
}

/// The energy back end for `eval` with the near-pair quadrature `near_pair` (`SPEC.md` §3.4;
/// `None`: [`energy_model`], the paper's energy).
///
/// # Errors
/// [`RcError::InvalidOption`] for bad Barnes–Hut or near-pair options.
pub fn energy_model_with_near_pair(
    eval: &EnergyEval,
    exponents: Exponents,
    near_pair: Option<NearPairOptions>,
    pool: Arc<rayon::ThreadPool>,
) -> Result<Box<dyn EnergyModel>, RcError> {
    if let Some(o) = &near_pair {
        o.validate()?;
    }
    match eval {
        EnergyEval::Exact => Ok(Box::new(
            ExactEnergy::new(pool, exponents).with_near_pair(near_pair),
        )),
        EnergyEval::BarnesHut(opts) => {
            opts.validate()?;
            Ok(Box::new(
                BarnesHutEnergy::new(pool, exponents, *opts).with_near_pair(near_pair),
            ))
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
    /// The session's main pool: every parallel loop of the solver (Gram assembly, block cluster
    /// tree products, hierarchy construction) runs on it, never on rayon's global pool, so that
    /// the thread count is the caller's and results are bitwise independent of it.
    pub pool: &'a Arc<rayon::ThreadPool>,
    /// The session's pool of the dense solver (`faer`), see [`FlowOptions::dense_threads`].
    pub dense_pool: &'a Arc<rayon::ThreadPool>,
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

/// The direct solvers: dense LU for the Sobolev metric (`faer`, on the request's dense pool;
/// the Gram matrix is assembled on its main pool), sparse LU for L², H¹, H² (§8.5).
#[derive(Debug, Clone, Copy, Default)]
pub struct DirectSolvers;

impl DirectSolvers {
    /// The direct solvers (the pools come with every [`SolverRequest`]).
    pub fn new() -> Self {
        DirectSolvers
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
                let gram = sobolev_gram(r.pool, r.network, r.positions, r.exponents, quadrature);
                let s = DenseSaddle::new(r.dense_pool.clone(), &gram.a, r.jacobian, r.layout)
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
    /// The multigrid saddle solver ([`MultigridSolvers`]) with these options, on the session's
    /// pools. The hierarchy is rebuilt at every step (and so after every subdivision); when it
    /// fails ([`RcError::MultigridDiverged`], [`RcError::MultigridNotConverged`]) the step falls
    /// back to the dense solver if the network has at most
    /// [`FlowOptions::dense_fallback_max_vertices`] vertices
    /// ([`StepReport::solver_fallback`]).
    Multigrid(MgOptions),
    /// A caller-supplied back end. It builds on the pools of the [`SolverRequest`], and gets
    /// the same dense fallback as [`Solver::Multigrid`].
    Custom(Arc<dyn SolverFactory>),
}

impl fmt::Debug for Solver {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Solver::Dense => write!(f, "Dense"),
            Solver::Multigrid(o) => write!(f, "Multigrid({o:?})"),
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
    /// Exact energy and dense solver up to `dense_max_vertices` vertices, Barnes–Hut energy
    /// ([`BhOptions::default`]) and multigrid ([`MgOptions::default`]) above. The session
    /// chooses at its start and again after every subdivision, with the vertex count of that
    /// moment. The crossover is measured in `SPEC.md` §17l (after the Schur-complement dense solve
    /// and the cheaper projections), on Freedman–He–Wang with the machine alone: the two choices
    /// cost the same per step at about 2000 vertices at 8 threads and 3000 at 64, so the default
    /// 2500 loses at most about 30 % on either side. The exact energy and the dense solver are
    /// the reference computation, so they win a tie.
    Auto {
        /// The crossover: networks with at most this many vertices use exact + dense.
        dense_max_vertices: usize,
    },
}

/// Default of [`FlowOptions::newton_polish`]: Newton steps once the largest violation is at most
/// 100 times its tolerance (`SPEC.md` §17l).
pub const DEFAULT_NEWTON_POLISH: f64 = 100.0;

/// Default of [`Backend::Auto::dense_max_vertices`](Backend::Auto): 2500 (`SPEC.md` §17l; 1500
/// in §17k, before the Schur-complement dense solve).
pub const DEFAULT_AUTO_DENSE_MAX_VERTICES: usize = 2500;

impl Default for Backend {
    /// `Auto { dense_max_vertices: 2500 }` (§10).
    fn default() -> Self {
        Backend::Auto {
            dense_max_vertices: DEFAULT_AUTO_DENSE_MAX_VERTICES,
        }
    }
}

impl Backend {
    /// The energy evaluation and solver used for a network of `n_vertices` vertices (the rule of
    /// [`Backend::Auto`] for it).
    pub fn resolve(&self, n_vertices: usize) -> (EnergyEval, Solver) {
        match self {
            Backend::Explicit { energy, solver } => (*energy, solver.clone()),
            Backend::Auto { dense_max_vertices } if n_vertices <= *dense_max_vertices => {
                (EnergyEval::Exact, Solver::Dense)
            }
            Backend::Auto { .. } => (
                EnergyEval::BarnesHut(BhOptions::default()),
                Solver::Multigrid(MgOptions::default()),
            ),
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
    /// Size of the dense solver's pool; `None`: `min(DEFAULT_DENSE_THREADS, threads)`.
    pub dense_threads: Option<usize>,
    /// Keep the path of every accepted step in the report (for audits; memory heavy).
    pub record_path: bool,
    /// When an iterative solver ([`Solver::Multigrid`], a custom factory) fails on a step, the
    /// step repeats the solve with the dense solver if the network has at most this many
    /// vertices (the Schur method holds about `4|V|² + k²` doubles, 1 GB at 5000); larger networks return
    /// the error. 0: never fall back. Default 5000 (`SPEC.md` §17h).
    pub dense_fallback_max_vertices: usize,
    /// Guard against collapsing edges (`SPEC.md` §17j): a step that would leave some edge
    /// shorter than `ratio · ℓ̄₀` (`ℓ̄₀` the mean edge length at the start of the session), and
    /// shorter than that edge was before the step, is not taken and fails with
    /// [`RcError::EdgeCollapse`], before the metric becomes singular. `None`: no guard. Not part
    /// of the paper or of the reference; it changes no energy and no step that is taken.
    /// Default `Some(1e-3)`.
    pub min_edge_length_ratio: Option<f64>,
    /// Project every trial of the line search back onto the constraints (§6.2, Eq. 23; the
    /// default). `false` is the paper's "Multigrid without projection" configuration of Fig. 20:
    /// the direction is still the projected gradient (tangent to the constraints, Eq. 22), but the
    /// trials `γ − τ p` are taken as they are, so nonlinear constraints (edge lengths, the
    /// length-weighted barycenter) drift by `O(τ²)` per step; the re-target sub-step then runs only
    /// for a moving length schedule, not for violated rows. Default `true` (`SPEC.md` §17l).
    pub back_projection: bool,
    /// Trials `γ − τ p + τ² q` with the second-order correction `q` of
    /// [`crate::solve::second_order_correction`] (one extra solve per
    /// step), so that the projection starts at an `O(τ³)` instead of an `O(τ²)` violation.
    /// Measured in `SPEC.md` §17l: the first trial of a step usually converges within the
    /// iteration limit instead of failing, which saves the failed trial's solves. Used only with
    /// hard constraints and [`FlowOptions::back_projection`]. Default `false`: on 24 perturbed
    /// Knot128 knots (`SPEC.md` §17l) it saves 8 % of the projection solves and 4 % of the time
    /// to 1.1 E* alone, but with the Newton polish (the default) it is 6 % slower than the polish
    /// alone.
    pub second_order_correction: bool,
    /// The Newton polish of the projection (`SPEC.md` §17l): once the largest violation is at
    /// most this factor times its tolerance, the projection continues with Newton steps (Jacobian
    /// at the iterate, lumped-mass metric, a sparse factorization each; see
    /// [`newton_correction`](crate::solve::newton_correction)) instead of chord steps with `Ā`.
    /// `None`: the chord method of §6.2 throughout (the paper's). Default
    /// `Some(`[`DEFAULT_NEWTON_POLISH`]`)`: measured on 24 perturbed Knot128 knots in `SPEC.md`
    /// §17l (0.69× the time to 1.1 E* of the chord method, no stuck step, no isotopy violation).
    /// It departs from the paper's projection only below 100 times the tolerance.
    pub newton_polish: Option<f64>,
    /// The near-pair quadrature of the energy (`SPEC.md` §3.4, §17o): edge pairs closer than
    /// `c₂ ℓ_IJ` get an accurate evaluation of the pair integral that the paper's endpoint rule
    /// (Eq. 18) approximates, blended in smoothly, so that two edges that come close in their
    /// interiors repel (the paper's rule does not see them). **Not part of the paper**; when no
    /// pair is that close, the energy and the flow are bitwise the paper's. `None`: the paper's
    /// energy. Default `Some(NearPairOptions::default())`.
    pub near_pair: Option<NearPairOptions>,
    /// Stop a chord projection early when its observed convergence rate cannot reach the target
    /// within the remaining iterations (the rate test of simplified-Newton iterations, Hairer &
    /// Wanner §IV.8; [`project_with_path_polished`]):
    /// the trial is rejected (`MaxIterations`) after fewer solves. Default `true`: on 24 perturbed
    /// Knot128 knots (`SPEC.md` §17l) the same iterations to 1.1 E* in 0.76× the time, no stuck
    /// step, no isotopy violation. `false`: every chord projection runs to its limit (the paper).
    pub projection_early_abort: bool,
    /// When a collision-guarded search with back-projection is stuck (not on the contact floor),
    /// project the base curve and, if that raises the objective, search again with the projected
    /// base's objective as the Armijo reference (WP16b, `SPEC.md` §17p;
    /// [`line_search::retry_with_projected_reference`]). **Not in the paper.** Default `false`;
    /// steps that are not stuck are unchanged either way.
    pub stuck_projected_reference: bool,
    /// When a search with back-projection and the Newton polish is stuck (not on the contact
    /// floor), repeat it with the chord projection of §6.2 (no polish) for this step (WP16b,
    /// `SPEC.md` §17p). The polish corrects in the lumped-mass metric, which favours the light
    /// vertices between short edges; its correction of the base's residual violation can raise
    /// the objective by more than any trial can gain (the chord step, in the metric `Ā`, is
    /// orthogonal to the direction to first order). **Not in the paper.** Default `true`
    /// (`SPEC.md` §18 Q10); steps that are not stuck are unchanged either way.
    pub stuck_retry_chord: bool,
    /// Back-off of the grow-previous start ([`StartRule::GrowPrevious`], `SPEC.md` §17t): when the
    /// first trial of a step was grown beyond the previous accepted step size and rejected, the
    /// next `w` steps start at the previous accepted step size (no growth; `min(⅔ τ_max,
    /// τ_prev)`), `w = 1` after the first such rejection and doubled after each further one
    /// (at most [`GROWTH_BACKOFF_MAX`]), back to 0 once a grown trial is accepted. It only lowers
    /// the first trial of some steps (every trial is still projected, CCD-checked and tested by
    /// Armijo), so the guarantees are those of the search. **Not in the paper** (whose start is
    /// `⅔ τ_max`). `false`: grow every step.
    pub growth_backoff: bool,
}

/// The longest wait of [`FlowOptions::growth_backoff`] (steps without growth).
pub const GROWTH_BACKOFF_MAX: u32 = 8;

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
            .field(
                "dense_fallback_max_vertices",
                &self.dense_fallback_max_vertices,
            )
            .field("min_edge_length_ratio", &self.min_edge_length_ratio)
            .field("back_projection", &self.back_projection)
            .field("second_order_correction", &self.second_order_correction)
            .field("newton_polish", &self.newton_polish)
            .field("near_pair", &self.near_pair)
            .field("projection_early_abort", &self.projection_early_abort)
            .field("stuck_projected_reference", &self.stuck_projected_reference)
            .field("stuck_retry_chord", &self.stuck_retry_chord)
            .field("growth_backoff", &self.growth_backoff)
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
            dense_fallback_max_vertices: DEFAULT_DENSE_FALLBACK_MAX_VERTICES,
            min_edge_length_ratio: Some(DEFAULT_MIN_EDGE_LENGTH_RATIO),
            back_projection: true,
            second_order_correction: false,
            newton_polish: Some(DEFAULT_NEWTON_POLISH),
            near_pair: Some(NearPairOptions::default()),
            projection_early_abort: true,
            stuck_projected_reference: false,
            stuck_retry_chord: true,
            growth_backoff: true,
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
        if let Some(o) = &self.near_pair {
            o.validate()?;
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
        if let Some(p) = self.newton_polish
            && !(p.is_finite() && p >= 1.0)
        {
            return Err(RcError::InvalidOption {
                name: "newton_polish",
                value: p.to_string(),
                requirement: "a finite factor ≥ 1 of the projection tolerance, or None",
            });
        }
        if let Some(r) = self.min_edge_length_ratio
            && !(r.is_finite() && r > 0.0 && r < 1.0)
        {
            return Err(RcError::InvalidOption {
                name: "min_edge_length_ratio",
                value: r.to_string(),
                requirement: "a ratio in (0, 1), or None",
            });
        }
        if let Backend::Explicit {
            solver: Solver::Multigrid(o),
            ..
        } = &self.backend
        {
            o.validate()?;
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
    /// The second-order correction of the trials (its constraint evaluations and solve; part of
    /// `line_search`).
    pub second_order: f64,
    /// CCD of the line search.
    pub ccd: f64,
    /// The near-pair quadrature (`SPEC.md` §3.4) at the base curve and in the line search's
    /// energies (part of `energy` and `line_search_energy`).
    pub near_pair: f64,
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
    /// Iterations (multigrid: V-cycles) of those solves, summed (0 for direct solvers).
    pub projection_cycles: usize,
    /// The solve of the second-order correction (`None`: not used, or `Φ₂ = 0`).
    pub second_order: Option<LinearReport>,
    /// The collision-guarded search took no step because two disjoint edges were closer than
    /// [`line_search::CONTACT_FLOOR`]` · ℓ̄₀` (`SPEC.md` §17l); the step is `Stuck`.
    pub contact_floor: bool,
    /// CCD counters of the line search.
    pub ccd: CcdReport,
    /// Edge pairs of the base curve in the band of the near-pair quadrature (`SPEC.md` §3.4).
    pub near_pairs: u64,
    /// Number of CCD queries of the line search.
    pub ccd_checks: u32,
    /// The re-target sub-step, if any.
    pub retarget: Option<RetargetReport>,
    /// The network was subdivided at the end of this step.
    pub subdivided: bool,
    /// The length schedule after this step: the final factor asked for, the factor imposed by
    /// the constraints, the factor the curve has (`None` without a schedule).
    pub schedule: Option<ScheduleProgress>,
    /// Energy back end used.
    pub energy_backend: &'static str,
    /// Solver back end used.
    pub solver_backend: &'static str,
    /// The iterative solver failed in this step and the dense solver answered instead
    /// ([`Solver::Multigrid`]); `None`: no fallback.
    pub solver_fallback: Option<SolverFallback>,
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

impl FlowReport {
    /// The error that ended the run ([`Status::Error`]), if any.
    pub fn error(&self) -> Option<&RcError> {
        match &self.status {
            Status::Error(e) => Some(e),
            _ => None,
        }
    }
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
    dense_pool: Arc<rayon::ThreadPool>,
    fallback: FallbackSlot,
    energy: Box<dyn EnergyModel>,
    solver: Arc<dyn SolverFactory>,
    mean_edge_length0: f64,
    iteration: usize,
    status: Status,
    cancel: Arc<AtomicBool>,
    schedule: Option<ScheduleState>,
    /// Total curve length when the schedule was set (the unit of `factor_measured`).
    schedule_length0: f64,
    subdivisions: u32,
    previous_tau: Option<f64>,
    /// [`FlowOptions::growth_backoff`]: steps left without growth, and the current wait.
    growth_wait: u32,
    growth_backoff: u32,
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
    second_order: Option<&'a [[f64; 3]]>,
    polish: Option<f64>,
    early_abort: bool,
}

impl SearchEnv for StepEnv<'_> {
    fn second_order_correction(&self) -> Option<&[[f64; 3]]> {
        self.second_order
    }

    fn below_contact_floor(&mut self, positions: &[[f64; 3]]) -> Result<bool, RcError> {
        let tol = line_search::CONTACT_FLOOR * self.mean_edge_length0;
        Ok(line_search::closest_pair_within(self.network, positions, tol).is_some())
    }

    fn objective(&mut self, positions: &[[f64; 3]]) -> Result<f64, RcError> {
        Ok(objective(self.energy, self.penalties, self.network, positions)?.0)
    }

    fn project(
        &mut self,
        trial: Vec<[f64; 3]>,
    ) -> Result<line_search::ProjectionWithPath, RcError> {
        match self.constraints {
            Some(c) => project_with_path_polished(
                self.solver,
                c,
                self.network,
                trial,
                self.projection,
                self.mean_edge_length0,
                line_search::ProjectionRules {
                    newton_polish: self.polish,
                    early_abort: self.early_abort,
                },
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

/// An energy model with a solver factory.
type Backends = (Box<dyn EnergyModel>, Arc<dyn SolverFactory>);

/// The energy model and solver factory for a network of `n_vertices` vertices
/// ([`Backend::resolve`]).
fn make_backends(
    options: &FlowOptions,
    n_vertices: usize,
    pool: &Arc<rayon::ThreadPool>,
) -> Result<Backends, RcError> {
    let (eval, solver) = options.backend.resolve(n_vertices);
    let energy =
        energy_model_with_near_pair(&eval, options.exponents, options.near_pair, pool.clone())?;
    let solver: Arc<dyn SolverFactory> = match solver {
        Solver::Dense => Arc::new(DirectSolvers::new()),
        Solver::Multigrid(o) => Arc::new(MultigridSolvers::new(o)?),
        Solver::Custom(f) => f,
    };
    Ok((energy, solver))
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
        let network_length: f64 = network.edge_geometry(&positions).length.iter().sum();
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
        let (energy, solver) = make_backends(&options, network.n_vertices(), &pool)?;
        let mut s = Session {
            soft: SoftState::new(options.method),
            network,
            positions,
            layout,
            constraints,
            penalties,
            obstacles,
            pool,
            dense_pool,
            fallback: FallbackSlot::default(),
            energy,
            solver,
            mean_edge_length0,
            iteration: 0,
            status: Status::Running,
            cancel: Arc::new(AtomicBool::new(false)),
            schedule,
            schedule_length0: network_length,
            subdivisions: 0,
            previous_tau: None,
            growth_wait: 0,
            growth_backoff: 0,
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
        self.schedule_length0 = self
            .network
            .edge_geometry(&self.positions)
            .length
            .iter()
            .sum();
        Ok(())
    }

    /// The schedule's progress as reported (`None` without a schedule).
    fn schedule_report(&self) -> Option<ScheduleProgress> {
        let s = self.schedule.as_ref()?;
        let length: f64 = self
            .network
            .edge_geometry(&self.positions)
            .length
            .iter()
            .sum();
        Some(s.report(length, self.schedule_length0))
    }

    fn schedule_moving(&self) -> bool {
        self.schedule.as_ref().is_some_and(|s| !s.finished())
    }

    fn build_solver(
        &self,
        jacobian: &SparseMatrix,
        hard: bool,
    ) -> Result<Box<dyn SaddleSolver>, RcError> {
        let request = SolverRequest {
            network: &self.network,
            positions: &self.positions,
            exponents: self.options.exponents,
            metric: self.options.metric,
            jacobian,
            layout: self.layout,
            constraints: hard.then_some(&self.constraints),
            pool: &self.pool,
            dense_pool: &self.dense_pool,
        };
        let solver = self.solver.build(&request)?;
        let iterative = !matches!(self.solver.name(), "direct");
        let n = self.network.n_vertices();
        if iterative && n <= self.options.dense_fallback_max_vertices {
            return Ok(Box::new(FallbackSaddle::new(
                solver,
                self.solver.name(),
                &request,
                self.fallback.clone(),
            )));
        }
        Ok(solver)
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
        // Without back-projection the rows drift by design: only a moving schedule re-targets.
        if !moving && (!violated || !self.options.back_projection) {
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
            let (proj, path) = project_with_path_polished(
                solver.as_mut(),
                &cs,
                &self.network,
                self.positions.clone(),
                &self.options.projection,
                self.mean_edge_length0,
                line_search::ProjectionRules {
                    newton_polish: self.options.newton_polish,
                    early_abort: self.options.projection_early_abort,
                },
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

    /// One step of §7.1. The session's status is updated (converged, stuck, max iterations; a
    /// failing step sets [`Status::Error`] and returns the same error, [`RcError::Cancelled`]
    /// sets [`Status::Cancelled`]). [`Session::run`] turns these errors into its status instead
    /// of returning them.
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
        let _ = self.energy.take_near_pair_seconds();
        if self.cancel.load(Ordering::Relaxed) {
            return Err(RcError::Cancelled);
        }
        let iteration = self.iteration + 1;
        let mut seconds = Timings::default();
        *self.fallback.lock().unwrap_or_else(|e| e.into_inner()) = None;

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
            let mut pg = projected_gradient(solver.as_mut(), self.layout, &df)?;
            // A `Point`-pinned vertex does not move: its exact direction is 0, the solve leaves
            // rounding noise (~1e-17) there. Without this, whether a trial is "exactly" on its
            // point target (SPEC §6.2: no projection step at Φ = 0) depended on that noise.
            let mut pinned = false;
            for c in self.constraints.constraints() {
                if let ResolvedConstraint::Point { vertex, .. } = c {
                    for &ax in self.layout.coords() {
                        pg.direction[*vertex as usize][ax] = 0.0;
                    }
                    pinned = true;
                }
            }
            if pinned {
                let x = self.layout.gather(&pg.direction);
                pg.gradient_norm = x.iter().map(|v| v * v).sum::<f64>().sqrt();
                pg.sobolev_norm_sq = self
                    .layout
                    .gather(&df)
                    .iter()
                    .zip(&x)
                    .map(|(p, q)| p * q)
                    .sum();
            }
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
            projection_cycles: 0,
            second_order: None,
            contact_floor: false,
            ccd: CcdReport::default(),
            ccd_checks: 0,
            near_pairs: ev.near_pairs,
            retarget,
            subdivided: false,
            schedule: None,
            energy_backend: self.energy.name(),
            solver_backend: self.solver.name(),
            solver_fallback: None,
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
            let mut max_start = match (&self.options.method, &self.soft.lbfgs) {
                (Method::Lbfgs { .. }, Some(l)) if !l.is_empty() => Some(1.0),
                _ => None,
            };
            // The growth back-off (`SPEC.md` §17t): hold the start at τ_prev while waiting.
            let backoff = self.options.growth_backoff
                && matches!(
                    self.options.line_search,
                    LineSearch::CollisionGuarded {
                        start: StartRule::GrowPrevious { .. },
                        ..
                    }
                );
            let holding = backoff && self.growth_wait > 0;
            if holding && let Some(prev) = self.previous_tau {
                max_start = Some(max_start.map_or(prev, |m: f64| m.min(prev)));
            }
            let project = !soft && self.options.back_projection;
            let correction = if project && self.options.second_order_correction {
                let t = Instant::now();
                let c = second_order_correction(
                    solver.as_mut(),
                    &self.constraints,
                    &self.network,
                    &self.positions,
                    &direction,
                    self.mean_edge_length0,
                )?;
                report.seconds.second_order = t.elapsed().as_secs_f64();
                report.second_order = c.linear.clone();
                c.linear.is_some().then_some(c.q)
            } else {
                None
            };
            let mut env = StepEnv {
                network: &self.network,
                energy: self.energy.as_mut(),
                penalties: &self.penalties,
                solver: solver.as_mut(),
                constraints: project.then_some(&self.constraints),
                projection: &self.options.projection,
                mean_edge_length0: self.mean_edge_length0,
                pool: &self.pool,
                obstacles: &self.obstacles,
                ccd: self.options.line_search.ccd_options(),
                cancel: &self.cancel,
                second_order: correction.as_deref(),
                polish: self.options.newton_polish,
                early_abort: self.options.projection_early_abort,
            };
            let input = SearchInput {
                base: &self.positions,
                direction: &direction,
                slope,
                f0,
                previous_tau: self.previous_tau,
                max_start,
            };
            let mut ls: LineSearchReport = line_search(&self.options.line_search, input, &mut env)?;
            if !ls.accepted
                && !ls.contact_floor
                && project
                && self.options.stuck_retry_chord
                && env.polish.is_some()
            {
                env.polish = None;
                let retry = line_search(&self.options.line_search, input, &mut env)?;
                ls = line_search::merge_after(&ls, retry);
            }
            if !ls.accepted
                && !ls.contact_floor
                && project
                && self.options.stuck_projected_reference
                && self.options.line_search.is_guarded()
                && let Some(retry) = line_search::retry_with_projected_reference(
                    &self.options.line_search,
                    input,
                    &mut env,
                    &ls,
                )?
            {
                ls = retry;
            }
            report.seconds.line_search = t.elapsed().as_secs_f64();
            report.seconds.line_search_energy = ls.seconds.energy;
            report.seconds.projection = ls.seconds.projection;
            report.seconds.ccd = ls.seconds.ccd;
            report.tau_start = ls.tau_start;
            report.tau_hi = ls.tau_hi;
            report.tau_max = ls.tau_max;
            report.halvings = ls.halvings;
            report.projection_solves = ls.projection_linear.len();
            report.projection_cycles = ls.projection_linear.iter().map(|l| l.iterations).sum();
            report.projection_residual = ls
                .projection_linear
                .iter()
                .map(|l| l.relative_residual)
                .fold(0.0, f64::max);
            report.ccd = ls.ccd.clone();
            report.ccd_checks = ls.ccd_checks;
            if backoff {
                if holding {
                    self.growth_wait -= 1;
                } else if let Some(prev) = self.previous_tau
                    && ls.tau_start > prev
                {
                    // The first trial was grown: accepted (reset) or rejected (wait longer).
                    let first_ok = ls
                        .trials
                        .first()
                        .is_some_and(|t| t.outcome == TrialOutcome::Accepted);
                    if first_ok {
                        self.growth_backoff = 0;
                    } else {
                        self.growth_backoff =
                            (2 * self.growth_backoff).clamp(1, GROWTH_BACKOFF_MAX);
                        self.growth_wait = self.growth_backoff;
                    }
                }
            }
            report.max_violation = ls.max_violation;
            report.trials = ls.trials;
            report.contact_floor = ls.contact_floor;
            // 4. Accept or stuck.
            if let (true, Some(new_pos)) = (ls.accepted, ls.positions) {
                self.check_edge_collapse(&new_pos, iteration)?;
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
        if !soft && (!accepted || !self.options.back_projection) {
            // The curve did not move in the line search, or was not projected: report its own
            // violation.
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
        report.schedule = self.schedule_report();
        report.solver_fallback = self
            .fallback
            .lock()
            .unwrap_or_else(|e| e.into_inner())
            .take();
        report.seconds.near_pair = self.energy.take_near_pair_seconds();
        report.seconds.total = t_step.elapsed().as_secs_f64();
        Ok(report)
    }

    /// The collapse guard of [`FlowOptions::min_edge_length_ratio`] on the positions of an
    /// accepted line search, before they are taken.
    fn check_edge_collapse(&self, new_pos: &[[f64; 3]], step: usize) -> Result<(), RcError> {
        let Some(ratio) = self.options.min_edge_length_ratio else {
            return Ok(());
        };
        let threshold = ratio * self.mean_edge_length0;
        let before = self.network.edge_geometry(&self.positions).length;
        let after = self.network.edge_geometry(new_pos).length;
        let worst = after
            .iter()
            .zip(&before)
            .enumerate()
            .filter(|(_, (a, b))| **a < threshold && **a < **b)
            .min_by(|x, y| x.1.0.total_cmp(y.1.0));
        match worst {
            Some((edge, (&length, &length_before))) => Err(RcError::EdgeCollapse {
                edge,
                length,
                length_before,
                threshold,
                step,
            }),
            None => Ok(()),
        }
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
        if matches!(self.options.backend, Backend::Auto { .. }) {
            // The vertex count changed: the rule of `Backend::Auto` may pick the other back ends.
            let (energy, solver) =
                make_backends(&self.options, self.network.n_vertices(), &self.pool)?;
            self.energy = energy;
            self.solver = solver;
        }
        Ok(())
    }

    /// Steps until a stop rule, `max_iterations`, the time limit, a stuck line search (without
    /// a moving schedule), cancellation (flag or `progress` returning [`Control::Stop`]) or a
    /// failing step.
    ///
    /// **Errors of steps do not make `run` fail.** [`Session::step`] returns the error of a
    /// failing step; `run` catches it and ends with [`Status::Error`] carrying that error
    /// ([`FlowReport::status`], [`FlowReport::error`]), or with [`Status::Cancelled`] for
    /// [`RcError::Cancelled`]. The report keeps the steps completed before it, and the curve is
    /// that of the last completed step (or sub-step: a re-target that finished before the error
    /// is kept). Python raises the error from the status.
    ///
    /// # Errors
    /// None at present: the `Result` is kept for the API of `SPEC.md` §9.1 (a future start-up
    /// check can fail here without changing the signature).
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
