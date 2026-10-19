//! Option and report types of the optimizer (SPEC §4, §7.1).

use std::fmt;

/// An option value violates its requirement (converted to `DevError::InvalidOption` by the crate).
#[derive(Debug, Clone, PartialEq)]
pub struct InvalidOption {
    /// Option name.
    pub name: &'static str,
    /// Offending value.
    pub value: f64,
    /// What is required.
    pub requirement: &'static str,
}

impl fmt::Display for InvalidOption {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "invalid option {} = {}: must be {}",
            self.name, self.value, self.requirement
        )
    }
}

impl std::error::Error for InvalidOption {}

/// Optimization method (SPEC §4).
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Method {
    /// L-BFGS with memory `memory` (§4.2).
    Lbfgs {
        /// Number of stored correction pairs (>= 1).
        memory: usize,
    },
    /// Steepest descent, direction scaled as in §4.1, with a line search (§4.4).
    GradientDescent,
    /// `x <- x - step * ℓ̄² * g`, no line search (§4.5).
    FixedStep {
        /// Step in units of `ℓ̄²`.
        step: f64,
    },
}

/// Line search (§4.3, §4.4).
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum LineSearch {
    /// Lewis–Overton weak Wolfe bracketing/bisection (Algorithm 2.6).
    WeakWolfe {
        /// Armijo constant, `0 <= c1 < c2`.
        c1: f64,
        /// Curvature constant, `c1 < c2 < 1`.
        c2: f64,
        /// Maximum number of doublings.
        max_doublings: u32,
        /// Maximum number of bisections.
        max_bisections: u32,
    },
    /// Armijo backtracking (N&W Alg. 3.1).
    Backtracking {
        /// Armijo constant in `(0, 1)`.
        c1: f64,
        /// Shrink factor in `(0, 1)`.
        shrink: f64,
        /// Maximum number of shrinks.
        max_shrinks: u32,
    },
}

/// Options of [`minimize`](super::minimize) and `optimize` (SPEC §7.1).
#[derive(Debug, Clone, PartialEq)]
pub struct OptimizeOptions {
    /// Method.
    pub method: Method,
    /// Line search (ignored by `FixedStep`).
    pub line_search: LineSearch,
    /// Maximum number of accepted iterations.
    pub max_iterations: usize,
    /// Stop when `ℓ̄ ‖g‖_∞ <= gradient_tolerance` (dimensionless).
    pub gradient_tolerance: f64,
    /// Relative energy decrease below which an iteration counts as stalled.
    pub energy_tolerance: f64,
    /// Number of consecutive stalled iterations that stop the run.
    pub stall_iterations: usize,
    /// Steepest-descent trial step: the fastest vertex moves by `initial_step · ℓ̄` at `t = 1`.
    pub initial_step: f64,
    /// Trial steps are capped so the fastest vertex moves at most `max_displacement · ℓ̄` (D19).
    pub max_displacement: f64,
    /// Pinned vertices (groups): zero gradient, never move (D18).
    pub pinned: Vec<u32>,
    /// Admissibility floor of the mesh optimizer (D31): a trial point with a face whose smallest
    /// corner angle is below this (radians) is rejected like an undefined energy (`+∞` in the line
    /// search); a starting point below it is an error. `0` switches the check off. Ignored by the
    /// generic [`minimize`](super::minimize).
    pub min_face_angle: f64,
}

/// Default of [`OptimizeOptions::min_face_angle`] (radians).
pub const DEFAULT_MIN_FACE_ANGLE: f64 = 1e-4;

impl Default for OptimizeOptions {
    fn default() -> Self {
        OptimizeOptions {
            method: Method::Lbfgs { memory: 10 },
            line_search: LineSearch::WeakWolfe {
                c1: 0.0,
                c2: 0.5,
                max_doublings: 30,
                max_bisections: 30,
            },
            max_iterations: 1000,
            gradient_tolerance: 1e-8,
            energy_tolerance: 1e-10,
            stall_iterations: 10,
            initial_step: 0.1,
            max_displacement: 0.5,
            pinned: Vec::new(),
            min_face_angle: DEFAULT_MIN_FACE_ANGLE,
        }
    }
}

fn check(ok: bool, name: &'static str, value: f64, requirement: &'static str) -> Result<(), InvalidOption> {
    if ok { Ok(()) } else { Err(InvalidOption { name, value, requirement }) }
}

impl OptimizeOptions {
    /// Checks every option; the first violation is returned.
    pub fn validate(&self) -> Result<(), InvalidOption> {
        match self.method {
            Method::Lbfgs { memory } => check(memory >= 1, "memory", memory as f64, ">= 1")?,
            Method::GradientDescent => {}
            Method::FixedStep { step } => {
                check(step.is_finite() && step > 0.0, "step", step, "finite and > 0")?
            }
        }
        match self.line_search {
            LineSearch::WeakWolfe { c1, c2, .. } => {
                check(c1.is_finite() && c1 >= 0.0, "c1", c1, "finite and >= 0")?;
                check(c2.is_finite() && c2 < 1.0, "c2", c2, "< 1")?;
                check(c1 < c2, "c2", c2, "> c1")?;
            }
            LineSearch::Backtracking { c1, shrink, .. } => {
                check(c1 > 0.0 && c1 < 1.0, "c1", c1, "in (0, 1)")?;
                check(shrink > 0.0 && shrink < 1.0, "shrink", shrink, "in (0, 1)")?;
            }
        }
        check(
            self.gradient_tolerance.is_finite() && self.gradient_tolerance >= 0.0,
            "gradient_tolerance",
            self.gradient_tolerance,
            "finite and >= 0",
        )?;
        check(
            self.energy_tolerance.is_finite() && self.energy_tolerance >= 0.0,
            "energy_tolerance",
            self.energy_tolerance,
            "finite and >= 0",
        )?;
        check(self.stall_iterations >= 1, "stall_iterations", self.stall_iterations as f64, ">= 1")?;
        check(
            self.initial_step.is_finite() && self.initial_step > 0.0,
            "initial_step",
            self.initial_step,
            "finite and > 0",
        )?;
        check(
            self.max_displacement.is_finite() && self.max_displacement > 0.0,
            "max_displacement",
            self.max_displacement,
            "finite and > 0",
        )?;
        check(
            self.min_face_angle.is_finite() && (0.0..std::f64::consts::FRAC_PI_3).contains(&self.min_face_angle),
            "min_face_angle",
            self.min_face_angle,
            "in [0, pi/3)",
        )?;
        Ok(())
    }
}

/// Callback return value.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Control {
    /// Keep going.
    Continue,
    /// Stop; the last accepted state is returned (`Cancelled`).
    Break,
}

/// Progress of the generic optimizer (after each accepted iteration).
#[derive(Debug, Clone, Copy)]
pub struct FlatProgress<'a> {
    /// Accepted iterations so far (>= 1).
    pub iteration: usize,
    /// Current energy.
    pub energy: f64,
    /// `ℓ̄ ‖g‖_∞`.
    pub gradient_norm: f64,
    /// Current variables.
    pub x: &'a [f64],
}

/// Why a run ended (SPEC §4.6).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum OptimizeStatus {
    /// `ℓ̄ ‖g‖_∞ <= gradient_tolerance`.
    GradientTolerance,
    /// Energy decrease below tolerance for `stall_iterations` consecutive iterations.
    EnergyStalled,
    /// `max_iterations` accepted iterations done.
    MaxIterations,
    /// No acceptable step was found (also from steepest descent), and at least one trial point
    /// of the last line search was undefined or inadmissible (a wall, D12/D31): the run stopped
    /// at an obstacle, not at a stationary point.
    LineSearchFailed,
    /// The line search from steepest descent found no lower energy at any trial step (down to
    /// `2^-max_bisections` of the first one), and every trial point was defined and admissible:
    /// no descent is possible from here, a stationary point to working precision (or a nonsmooth
    /// point where `−g` is not a descent direction; Lewis–Overton's "breakdown", their normal
    /// end of a nonsmooth BFGS run). SPEC §4.6, D34.
    Stationary,
    /// The progress callback returned `Break`.
    Cancelled,
    /// A steepest descent direction was not a descent direction (rounding).
    NonDescent,
}

/// Report of one optimization run (SPEC §4.7).
#[derive(Debug, Clone, PartialEq)]
pub struct OptimizeReport {
    /// Why the run ended.
    pub status: OptimizeStatus,
    /// Accepted iterations.
    pub iterations: usize,
    /// Energy+gradient evaluations (including undefined ones).
    pub evaluations: usize,
    /// Line search trial points.
    pub line_search_trials: usize,
    /// Trials where the energy was undefined (treated as `+∞`).
    pub undefined_trials: usize,
    /// Trials rejected by the admissibility floors (`min_face_angle`, `min_fold_angle`; treated
    /// as `+∞`, D31/D32).
    pub inadmissible_trials: usize,
    /// Steps accepted at the cap `max_displacement`.
    pub capped_steps: usize,
    /// Steps accepted after a limit was hit (satisfy Armijo only).
    pub weak_steps: usize,
    /// Pairs not stored because `sᵀy <= 0` or not finite.
    pub skipped_pairs: usize,
    /// Memory resets (non-descent, line search failure).
    pub memory_resets: usize,
    /// Energy at the start.
    pub energy_initial: f64,
    /// Energy at the end.
    pub energy_final: f64,
    /// Energy at each accepted iterate, including the start.
    pub energy_history: Vec<f64>,
    /// `ℓ̄ ‖g‖_∞` at each accepted iterate, including the start.
    pub gradient_norm_history: Vec<f64>,
    /// The length unit `ℓ̄`.
    pub mean_edge_length: f64,
    /// Largest total vertex displacement from the start, in units of `ℓ̄`.
    pub max_vertex_displacement: f64,
    /// Fixed step only: iterations where the energy increased.
    pub energy_increases: usize,
    /// Near ties of the last evaluation (objective-defined).
    pub near_ties: usize,
    /// Wall time.
    pub seconds: f64,
}
