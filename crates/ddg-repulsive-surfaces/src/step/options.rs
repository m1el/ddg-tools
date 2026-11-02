//! The options of a flow (`SPEC.md` §7, §9.1, §11).

// `!(x > tol)` is deliberate throughout: a NaN must fail the test.
#![allow(clippy::neg_cmp_op_on_partial_ord)]

use std::sync::Arc;
use std::time::Duration;

use ddg_repulsive_curves::obstacle::Obstacle;

use crate::ccd::CcdOptions;
use crate::constraints::{BarycenterMode, Constraint};
use crate::energy::{BhOptions, Energy, EnergyEval};
use crate::error::RsError;
use crate::kernel::Exponents;
use crate::metric::{FractionalForm, L2Kind, Metric};
use crate::penalties::{Penalty, Quadrature};
use crate::remesh::{RemeshOptions, RemeshOrder};
use crate::solve::iterative::{Backend, HierarchicalOptions};
use crate::solve::precond::PreconditionerForm;
use crate::solve::{GmresSide, NewtonSteps};

/// The sufficient-decrease test of the line search (`SPEC.md` §7.2, D15).
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum ArmijoRule {
    /// Accept `τ` if `F(f − τx) ≤ F(f) − c τ ⟨g, x⟩` (default, `c = 1e-2`).
    Standard {
        /// The constant `c ∈ (0, 1)`.
        c: f64,
    },
    /// The reference's test `F(f) − F(f − τx) ≥ c τ ‖x‖ cos∠(g, x)` (`ref: line_search.cpp:86,103`;
    /// the `‖g‖` factor is missing, so the test is not scale-invariant: it is the standard test
    /// with `c' = c/‖g‖`).
    ///
    /// [formula: rs-armijo-reference]
    Reference {
        /// The constant `c`.
        c: f64,
    },
}

impl Default for ArmijoRule {
    fn default() -> Self {
        ArmijoRule::Standard { c: 1e-2 }
    }
}

/// The first step size `τ₀` of a line search (`SPEC.md` §7.2).
#[derive(Debug, Clone, PartialEq)]
pub enum StartRule {
    /// `τ₀ = 2/‖x‖_F` (absolute units; `ref: surface_flow.cpp:60–66`).
    Reference,
    /// `τ₀` moves the farthest vertex by `lengths · ℓ̄₀` (`ℓ̄₀` the mean edge length of the
    /// initial mesh; scale-invariant).
    MaxDisplacement {
        /// The displacement in mean edge lengths.
        lengths: f64,
    },
    /// `factor ×` the last accepted `τ`; `first` for the first step (and until a step is
    /// accepted).
    GrowPrevious {
        /// The growth factor.
        factor: f64,
        /// The rule without a previous step.
        first: Box<StartRule>,
    },
    /// A fixed `τ₀` (absolute; the oracle's `start abs`, tests).
    Absolute {
        /// The step size.
        tau: f64,
    },
}

impl Default for StartRule {
    /// The provisional default `GrowPrevious { factor: 2, first: MaxDisplacement { lengths: 1 } }`
    /// (`SPEC.md` §7.2, to be confirmed by the measurement of §14.3).
    fn default() -> Self {
        StartRule::GrowPrevious {
            factor: 2.0,
            first: Box::new(StartRule::MaxDisplacement { lengths: 1.0 }),
        }
    }
}

impl StartRule {
    fn validate(&self) -> Result<(), RsError> {
        let pos = |name: &'static str, x: f64| {
            if x.is_finite() && x > 0.0 {
                Ok(())
            } else {
                Err(RsError::InvalidOption {
                    name,
                    value: x.to_string(),
                    requirement: "finite and positive",
                })
            }
        };
        match self {
            StartRule::Reference => Ok(()),
            StartRule::MaxDisplacement { lengths } => pos("line_search.start.lengths", *lengths),
            StartRule::GrowPrevious { factor, first } => {
                pos("line_search.start.factor", *factor)?;
                first.validate()
            }
            StartRule::Absolute { tau } => pos("line_search.start.tau", *tau),
        }
    }
}

/// The backtracking line search (`SPEC.md` §7.2, §7.3).
#[derive(Debug, Clone, PartialEq)]
pub struct LineSearch {
    /// The sufficient-decrease test.
    pub armijo: ArmijoRule,
    /// The first step size.
    pub start: StartRule,
    /// The search stops (the step is `Stuck`) below `τ_min = tau_min_rel · ℓ̄₀ / ‖x‖_∞`
    /// (default `1e-12`).
    pub tau_min_rel: f64,
    /// An optional cap on the number of halvings (`None`: only `τ_min`; the oracle's `maxhalv`).
    pub max_halvings: Option<u32>,
    /// Continuous collision detection (`SPEC.md` §7.3, D20; WP8): the start becomes
    /// `min(τ₀, ⅔ τ_max)` when a contact closer than the gap is found along `f − t x`. The
    /// `reference_length` of the options defaults to `ℓ̄₀` of the run. Obstacles are
    /// [`FlowOptions::ccd_obstacles`]; remeshing collapses are guarded too.
    pub ccd: Option<CcdOptions>,
}

impl Default for LineSearch {
    fn default() -> Self {
        LineSearch {
            armijo: ArmijoRule::default(),
            start: StartRule::default(),
            tau_min_rel: 1e-12,
            max_halvings: None,
            ccd: None,
        }
    }
}

/// When a run has converged (`SPEC.md` §7.4, D29).
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum StopRule {
    /// `F_{k−window} − F_k ≤ relative · |F_k|` over the last `window` accepted steps (default
    /// `1e-7`, 10).
    EnergyStall {
        /// The relative decrease.
        relative: f64,
        /// Number of accepted steps.
        window: usize,
    },
    /// `√(⟨g, x⟩_k / ⟨g, x⟩_1) < eps`: the `Hs` norm of the gradient relative to its value at the
    /// session's first step.
    SobolevNorm {
        /// The ratio.
        eps: f64,
    },
    /// No convergence test (only `max_iterations`, the time limit, stuck steps).
    Off,
}

impl Default for StopRule {
    fn default() -> Self {
        StopRule::EnergyStall {
            relative: 1e-7,
            window: 10,
        }
    }
}

/// When a growth schedule raises the targets (`SPEC.md` §5.1, D11).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum ScheduleTiming {
    /// Step `k` projects onto the target of step `k` (ours).
    #[default]
    BeforeProjection,
    /// Step `k` projects onto the target of step `k − 1` and raises it afterwards (the
    /// reference).
    AfterProjection,
}

/// A linear schedule of an area or volume target (`SPEC.md` §5.1, D11): the target moves from
/// its initial value `T₀` to `factor · T₀ + offset` over `steps` steps (the reference's
/// `constrain area 3 200 1.5` is `factor 3, offset 1.5, steps 200`).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Schedule {
    /// Index of an `Area` or `Volume` constraint in [`FlowOptions::constraints`].
    pub constraint: usize,
    /// The factor on the initial target.
    pub factor: f64,
    /// Added to the final target.
    pub offset: f64,
    /// Number of steps, at least 1.
    pub steps: usize,
}

/// Moving pins (`SPEC.md` §5.1; the reference's `constrain vertices … move`): the targets of a
/// `Pins` constraint move linearly by `displacement` (one per pin) over `steps` steps.
#[derive(Debug, Clone, PartialEq)]
pub struct PinSchedule {
    /// Index of a `Pins` constraint in [`FlowOptions::constraints`].
    pub constraint: usize,
    /// The total displacement of every pin of the constraint, in its order.
    pub displacement: Vec<[f64; 3]>,
    /// Number of steps, at least 1.
    pub steps: usize,
}

/// The L-BFGS variants of the evaluation (`SPEC.md` D26, D27; WP9).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum LbfgsVariant {
    /// The standard two-loop recursion.
    #[default]
    Standard,
    /// The reference's variant.
    Reference,
}

/// The optimizer (`SPEC.md` §14.2). Only [`Method::ProjectedGradient`] is implemented by the
/// session; the comparison optimizers are WP9's.
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub enum Method {
    /// Projected gradient descent in the chosen metric (the paper).
    #[default]
    ProjectedGradient,
    /// Accelerated quadratic proxy (Kovalsky et al. 2016).
    Aqp {
        /// The momentum parameter.
        kappa: f64,
    },
    /// L-BFGS with the H¹ initial Hessian.
    H1Lbfgs {
        /// Number of stored pairs.
        memory: usize,
        /// The variant.
        variant: LbfgsVariant,
    },
    /// Blended quasi-Newton (Zhu et al. 2018).
    Bcqn {
        /// Number of stored pairs.
        memory: usize,
    },
}

impl Method {
    /// A short name.
    pub fn name(&self) -> &'static str {
        match self {
            Method::ProjectedGradient => "projected gradient",
            Method::Aqp { .. } => "AQP",
            Method::H1Lbfgs { .. } => "H1 L-BFGS",
            Method::Bcqn { .. } => "BCQN",
        }
    }
}

/// The options of a flow (`SPEC.md` §9.1). Plain public fields; [`Default`] gives the paper's
/// flow with our decisions (§11); [`FlowOptions::validate`] checks the values (the
/// [`Session`](super::Session) also checks them against the surface).
#[derive(Debug, Clone)]
pub struct FlowOptions {
    /// `(α, β)` of the tangent-point energy and the metric (default `(6, 12)`).
    pub exponents: Exponents,
    /// The energy (`SPEC.md` §3.5): tangent-point (default), Coulomb, Willmore.
    pub energy: Energy,
    /// How the tangent-point energy (and the tangent-point obstacles) are evaluated: Barnes–Hut
    /// with `θ = 0.5`, leaf size 8 (default) or exact. The Barnes–Hut lists are built once per
    /// step and frozen during the line search (§3.4).
    pub energy_eval: EnergyEval,
    /// `w_TP`: the objective is `F = w_TP ℰ + Σ penalties` and the Sobolev metric is scaled by
    /// it (D6). Default 1.
    pub tp_weight: f64,
    /// The metric (default `Sobolev { low_order: true }`); `L2`, `H1`, `H2` are the comparison
    /// metrics of §4.5, solved directly (sparse `LDLᵀ`) with all constraints.
    pub metric: Metric,
    /// The solver of the Sobolev saddle (default `Auto`, dense up to 1000 faces).
    pub backend: Backend,
    /// The optimizer (default projected gradient; the others are WP9's).
    pub method: Method,
    /// Hard constraints (default: the barycenter of every component, `SPEC.md` §5.1).
    pub constraints: Vec<Constraint>,
    /// `None` = `Auto` (D8 as revised): one global barycenter when the low-order term `B⁰` is
    /// on, per component otherwise (and for the comparison metrics).
    pub barycenter_mode: Option<BarycenterMode>,
    /// The paper's barycenter shift (Sec. 7.1.1, Fig. 8; documented as unit-inconsistent): the
    /// `a_v`-weighted mean of the L² differential of each component is added to the direction
    /// and the barycenter is not projected (its target follows the surface).
    pub barycenter_shift: bool,
    /// Growth schedules of area/volume targets.
    pub schedules: Vec<Schedule>,
    /// Moving pins.
    pub pin_schedules: Vec<PinSchedule>,
    /// When schedules raise their targets (D11).
    pub schedule_timing: ScheduleTiming,
    /// Penalties `w_i P_i` (§5.2).
    pub penalties: Vec<Penalty>,
    /// The line search.
    pub line_search: LineSearch,
    /// The Newton rounds of the projection (default until tolerance, at most 5; D12).
    pub newton_steps: NewtonSteps,
    /// Remeshing after every step (default on, `SPEC.md` §10); `None` keeps the topology.
    pub remesh: Option<RemeshOptions>,
    /// The convergence test.
    pub stop: StopRule,
    /// Maximal number of steps of [`Session::run`](super::Session::run) (default 1000).
    pub max_iterations: usize,
    /// Wall-clock limit of a run.
    pub time_limit: Option<Duration>,
    /// Consecutive stuck steps that end a run with `Status::Stuck` (default 5, at least 1).
    pub stuck_limit: usize,
    /// Size of the explicit thread pool (`None`: all cores); results do not depend on it.
    pub threads: Option<usize>,
    /// Static triangle obstacles for the CCD of the line search (the penalties' obstacles are
    /// point sets without connectivity).
    pub ccd_obstacles: Vec<Arc<Obstacle>>,
    /// Sets every reference mode of `SPEC.md` §11 at once ([`FlowOptions::effective`]); for the
    /// black-box comparison only.
    pub reference_compat: bool,
}

impl Default for FlowOptions {
    fn default() -> Self {
        FlowOptions {
            exponents: Exponents::default(),
            energy: Energy::TangentPoint,
            energy_eval: EnergyEval::BarnesHut(BhOptions::default()),
            tp_weight: 1.0,
            metric: Metric::default(),
            backend: Backend::default(),
            method: Method::ProjectedGradient,
            constraints: vec![Constraint::Barycenter {
                component: None,
                target: None,
            }],
            barycenter_mode: None,
            barycenter_shift: false,
            schedules: Vec::new(),
            pin_schedules: Vec::new(),
            schedule_timing: ScheduleTiming::default(),
            penalties: Vec::new(),
            line_search: LineSearch::default(),
            newton_steps: NewtonSteps::default(),
            remesh: Some(RemeshOptions::default()),
            stop: StopRule::default(),
            max_iterations: 1000,
            time_limit: None,
            stuck_limit: 5,
            threads: None,
            ccd_obstacles: Vec::new(),
            reference_compat: false,
        }
    }
}

fn invalid(name: &'static str, value: impl ToString, requirement: &'static str) -> RsError {
    RsError::InvalidOption {
        name,
        value: value.to_string(),
        requirement,
    }
}

impl FlowOptions {
    /// Checks every value that does not depend on the surface.
    ///
    /// # Errors
    /// [`RsError::InvalidOption`] naming the option; [`RsError::Unsupported`] for a comparison
    /// optimizer (WP9).
    pub fn validate(&self) -> Result<(), RsError> {
        if !(self.tp_weight.is_finite() && self.tp_weight > 0.0) {
            return Err(invalid("tp_weight", self.tp_weight, "finite and positive"));
        }
        if let Energy::Coulomb { alpha } = self.energy
            && !(alpha.is_finite() && alpha > 0.0)
        {
            return Err(invalid("energy.alpha", alpha, "finite and positive"));
        }
        if let EnergyEval::BarnesHut(o) = &self.energy_eval {
            o.validate()?;
        }
        match self.backend {
            Backend::Dense => {}
            Backend::Hierarchical(o) | Backend::Auto { hierarchical: o, .. } => {
                o.gmres.validate()?;
                o.cluster.validate()?;
                o.preconditioner.validate()?;
            }
        }
        let c = match self.line_search.armijo {
            ArmijoRule::Standard { c } | ArmijoRule::Reference { c } => c,
        };
        if !(c > 0.0 && c < 1.0) {
            return Err(invalid("line_search.armijo.c", c, "in (0, 1)"));
        }
        self.line_search.start.validate()?;
        if !(self.line_search.tau_min_rel.is_finite() && self.line_search.tau_min_rel > 0.0) {
            return Err(invalid(
                "line_search.tau_min_rel",
                self.line_search.tau_min_rel,
                "finite and positive",
            ));
        }
        if let Some(o) = &self.line_search.ccd {
            o.validate()?;
        }
        if let Some(r) = &self.remesh {
            r.validate()?;
        }
        match self.stop {
            StopRule::EnergyStall { relative, window } => {
                if window == 0 {
                    return Err(invalid("stop.window", 0, "at least one step"));
                }
                if !(relative.is_finite() && relative > 0.0) {
                    return Err(invalid("stop.relative", relative, "finite and positive"));
                }
            }
            StopRule::SobolevNorm { eps } => {
                if !(eps.is_finite() && eps > 0.0) {
                    return Err(invalid("stop.eps", eps, "finite and positive"));
                }
            }
            StopRule::Off => {}
        }
        if self.stuck_limit == 0 {
            return Err(invalid("stuck_limit", 0, "at least 1"));
        }
        if self.threads == Some(0) {
            return Err(invalid("threads", 0, "at least 1 thread"));
        }
        if let NewtonSteps::UntilTolerance { max: 0 } = self.newton_steps {
            return Err(invalid("newton_steps.max", 0, "at least 1 (use Fixed(0) for none)"));
        }
        for s in &self.schedules {
            if !(s.factor.is_finite() && s.offset.is_finite()) {
                return Err(invalid(
                    "schedule",
                    format!("factor {}, offset {}", s.factor, s.offset),
                    "finite factor and offset",
                ));
            }
            if s.steps == 0 {
                return Err(invalid("schedule.steps", 0, "at least 1"));
            }
        }
        for s in &self.pin_schedules {
            if s.steps == 0 {
                return Err(invalid("pin_schedule.steps", 0, "at least 1"));
            }
            if s.displacement.iter().flatten().any(|x| !x.is_finite()) {
                return Err(invalid(
                    "pin_schedule.displacement",
                    format!("{:?}", s.displacement),
                    "finite displacements",
                ));
            }
        }
        if self.method != Method::ProjectedGradient {
            if !self.constraints.is_empty() {
                return Err(RsError::MethodNeedsSoftProblem {
                    method: self.method.name(),
                });
            }
            return Err(RsError::Unsupported {
                what: "the comparison optimizers (AQP, H1 L-BFGS, BCQN) are WP9's",
            });
        }
        Ok(())
    }

    /// The options actually used: with [`FlowOptions::reference_compat`] every reference mode of
    /// `SPEC.md` §11 that the crate implements is set (D8 global barycenter, D11 schedules after
    /// the projection, D12 one Newton step, D13 left GMRES with the reference's preconditioner —
    /// face pairs D32 and the printed form D33 —, D15 the reference's Armijo test and start, D16
    /// smooth before flip, D21 per-vertex quadrature of the implicit potentials, D25 the identity
    /// L² metric); otherwise a copy.
    ///
    /// Not reproduced (no option exists): the reference's 3× volume Jacobian (D9), its absolute
    /// preconditioner `ε` (D14; ours is relative), its absolute `τ_min = 1e-20`, its penalty value
    /// conventions (D24), and the comparison optimizers' variants (WP9).
    pub fn effective(&self) -> FlowOptions {
        let mut o = self.clone();
        if !self.reference_compat {
            return o;
        }
        o.barycenter_mode = Some(BarycenterMode::Global);
        o.schedule_timing = ScheduleTiming::AfterProjection;
        o.newton_steps = NewtonSteps::Fixed(1);
        let reference_hier = |h: HierarchicalOptions| {
            let mut r = HierarchicalOptions::reference();
            r.cluster = h.cluster;
            r.gmres.restart = h.gmres.restart;
            r.gmres.tolerance = h.gmres.tolerance;
            r.gmres.max_iterations = h.gmres.max_iterations;
            r.gmres.side = GmresSide::Left;
            r.preconditioner.fractional = FractionalForm::FacePairs;
            r.preconditioner.form = PreconditionerForm::Congruence;
            r.preconditioner.eps_multi_rel = 1e-2;
            r
        };
        o.backend = match self.backend {
            Backend::Dense => Backend::Dense,
            Backend::Hierarchical(h) => Backend::Hierarchical(reference_hier(h)),
            Backend::Auto {
                dense_max_faces,
                hierarchical,
            } => Backend::Auto {
                dense_max_faces,
                hierarchical: reference_hier(hierarchical),
            },
        };
        let c = match self.line_search.armijo {
            ArmijoRule::Standard { c } | ArmijoRule::Reference { c } => c,
        };
        o.line_search.armijo = ArmijoRule::Reference { c };
        o.line_search.start = StartRule::Reference;
        if let Some(r) = o.remesh.as_mut() {
            r.order = RemeshOrder::SmoothThenFlip;
        }
        for p in &mut o.penalties {
            match p {
                Penalty::ImplicitObstacle { quadrature, .. }
                | Penalty::ImplicitAttractor { quadrature, .. } => *quadrature = Quadrature::PerVertex,
                _ => {}
            }
        }
        if let Metric::L2(_) = o.metric {
            o.metric = Metric::L2(L2Kind::Identity);
        }
        o
    }

    /// The barycenter mode this run uses (`Auto` resolved, D8 as revised): global when the
    /// metric is the Sobolev metric with `B⁰`, per component otherwise.
    pub fn resolved_barycenter_mode(&self) -> BarycenterMode {
        self.barycenter_mode.unwrap_or(match self.metric {
            Metric::Sobolev { low_order: true } => BarycenterMode::Global,
            _ => BarycenterMode::PerComponent,
        })
    }
}
