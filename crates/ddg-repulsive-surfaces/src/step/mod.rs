//! The session, one step, line searches, projection, schedules and stop rules (`SPEC.md` §7;
//! WP5).
//!
//! One step ([`Session::step`], `SPEC.md` §7.1):
//!
//! 1. the constraint set of the step (schedules move area/volume targets and pins, D11; the
//!    rank check drops dependent rows, D10);
//! 2. the Barnes–Hut tree and lists of the energy and of every tangent-point obstacle, built at
//!    the base positions and frozen for the step; `F` and `dF` (energy + penalties);
//! 3. the solver of `A′` (dense or hierarchical for the Sobolev metric, with the Willmore term
//!    and the obstacle diagonal of §4.4; sparse `LDLᵀ` for the comparison metrics of §4.5), the
//!    Schur complement of the nonlinear rows and the direction `x` (§6.1); the barycenter
//!    shift (optional);
//! 4. the backtracking line search on `F(f − τx)` (§7.2: standard or reference Armijo, start
//!    rules, `τ_min`; every trial refits the frozen trees; a non-finite trial is rejected;
//!    cancellation between trials), its start capped by continuous collision detection when on
//!    (§7.3, WP8);
//! 5. the projection (§6.2): positional rows exactly, Newton rounds with the safeguard;
//! 6. remeshing (§10) with `L₀` of the initial mesh; pins, squared-error targets and attractor
//!    sets follow the vertex map, components and pins are renumbered;
//! 7. the stop rules (§7.4).
//!
//! [`Session::run`] loops until a stop rule, `max_iterations`, the time limit, `stuck_limit`
//! consecutive stuck steps or cancellation; [`flow`] is both in one call. [`run_scene`] runs a
//! reference scene file.
//!
//! **Determinism.** Every part is bitwise independent of the thread count (explicit pool), so
//! trajectories are identical at any `threads`, and stepping a session by hand equals [`flow`].

mod direction;
mod objective;
mod options;
mod scene_run;
mod session;

pub use direction::{COMPARISON_EPS_REL, REFINEMENT_STEPS, REFINEMENT_TOLERANCE, SparseSaddle};
pub use options::{
    ArmijoRule, FlowOptions, LbfgsVariant, LineSearch, Method, PinSchedule, Schedule,
    ScheduleTiming, StartRule, StopRule,
};
pub use scene_run::{RunReport, RunSceneOptions, run_scene, scene_options};
pub use session::{
    Control, FlowReport, FlowResult, Session, Status, StepReport, StepStatus, Timings, Trial,
    TrialOutcome, flow,
};
