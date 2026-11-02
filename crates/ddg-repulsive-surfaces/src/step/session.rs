//! The session: one step of `SPEC.md` §7.1 and the run loop (§7.4).

// `!(x > tol)` is deliberate throughout: a NaN must fail the test.
#![allow(clippy::neg_cmp_op_on_partial_ord)]

use std::collections::{HashMap, VecDeque};
use std::sync::Arc;
use std::sync::atomic::{AtomicBool, Ordering};
use std::time::Instant;

use crate::ccd::{
    CcdCollapseGuard, CcdObstacle, CcdOptions, CcdReport, IntersectionBand, first_contact_with_band,
};
use crate::constraints::{BarycenterMode, Constraint, ConstraintSet, RowKind};
use crate::energy::exact::make_pool;
use crate::energy::{Energy, EnergyEval};
use crate::error::RsError;
use crate::geometry::mean_edge_length;
use crate::metric::{Metric, SobolevWeights};
use crate::penalties::{Penalty, willmore_matrix};
use crate::remesh::{CollapseGuard, RemeshReport, Remeshed, remesh};
use crate::solve::iterative::MetricExtras;
use crate::solve::{
    LinearReport, Projection, ProjectionOptions, ProjectionStatus, SchurComplement, project,
};
use crate::surface::Surface;

use crate::optimize::{MethodDiagnostics, MethodState};
use super::direction::StepSolver;
use super::objective::{Evaluation, Objective, component_means, is_geometric};
use super::options::{ArmijoRule, FlowOptions, Method, ScheduleTiming, StartRule, StopRule};

// =============================================================================================
// Reports
// =============================================================================================

/// The state of a session or a run (`SPEC.md` §7.4).
#[derive(Debug, Clone, PartialEq)]
pub enum Status {
    /// Not finished.
    Running,
    /// The stop rule holds (and no schedule is moving its targets).
    Converged,
    /// `max_iterations` steps were taken.
    MaxIterations,
    /// The run's time limit was reached.
    TimeLimit,
    /// `stuck_limit` consecutive steps found no acceptable step size.
    Stuck,
    /// Cancelled through the cancel flag or the progress callback.
    Cancelled,
    /// A step failed; the surface is that of the last completed step.
    Error(RsError),
}

impl Status {
    /// `true` for every status but [`Status::Running`].
    pub fn is_terminal(&self) -> bool {
        !matches!(self, Status::Running)
    }
}

/// The outcome of one step.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum StepStatus {
    /// A step size was accepted.
    Accepted,
    /// The line search reached `τ_min` (or the halving cap) or the direction was not a descent
    /// direction; the positions were not moved by the descent (projection, schedules and
    /// remeshing still ran, as the reference).
    Stuck,
}

/// What the progress callback of [`Session::run`] returns.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Control {
    /// Keep going.
    Continue,
    /// Stop the run ([`Status::Cancelled`]).
    Stop,
}

/// How one line-search trial ended.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum TrialOutcome {
    /// The Armijo test held.
    Accepted,
    /// The Armijo test failed; the step was halved.
    ArmijoFailed,
    /// The trial's objective was not finite (a degenerate face, coincident barycenters, an
    /// implicit obstacle crossed); the step was halved.
    NonFinite,
    /// The Armijo test held but the projection of the trial was not admissible (`F` not
    /// finite after it, or its motion crosses the CCD gap: `LineSearch::guard_projection`);
    /// the step was halved.
    ProjectionBlocked,
}

/// One line-search trial.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Trial {
    /// The step size.
    pub tau: f64,
    /// `F(f − τx)` (NaN if not finite).
    pub energy: f64,
    /// How it ended.
    pub outcome: TrialOutcome,
}

/// Wall times of one step, seconds (`SPEC.md` §9.1).
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct Timings {
    /// Barnes–Hut tree, lists and the obstacles' cross lists.
    pub bvh: f64,
    /// `F` and `dF` at the base positions.
    pub energy: f64,
    /// The solver of `A′`: dense assembly + LU, or the block cluster tree, `L′` factorization
    /// and the fractional tree; the constraint set.
    pub factor: f64,
    /// Part of `factor`: the block cluster tree of `B`, `B⁰` (hierarchical back end).
    pub metric_tree: f64,
    /// Part of `factor`: the tree of the preconditioner's `L_{2−s}` (hierarchical back end).
    pub fractional_tree: f64,
    /// Part of `factor`: assembly and sparse `LDLᵀ` of the preconditioner's `L′`.
    pub lprime: f64,
    /// The Schur columns and the direction (GMRES or back substitutions).
    pub solve: f64,
    /// Part of `solve`: the applications of the preconditioner (hierarchical back end).
    pub preconditioner: f64,
    /// The line search (all trials).
    pub line_search: f64,
    /// Continuous collision detection.
    pub ccd: f64,
    /// The projection and the evaluation after it.
    pub projection: f64,
    /// Remeshing and the transport of the per-vertex data.
    pub remesh: f64,
    /// The whole step.
    pub total: f64,
}

/// The report of one step (`SPEC.md` §9.1).
#[derive(Debug, Clone, PartialEq)]
pub struct StepReport {
    /// The step number of the session (1-based).
    pub iteration: usize,
    /// Accepted or stuck.
    pub status: StepStatus,
    /// `F` at the base positions (this step's Barnes–Hut lists).
    pub energy_before: f64,
    /// `F` at the accepted trial (`energy_before` when stuck).
    pub energy_trial: f64,
    /// `F` after the projection (same lists; before remeshing).
    pub energy_after: f64,
    /// `ℰ` (unweighted) after the projection.
    pub energy_term_after: f64,
    /// `w_i P_i` after the projection.
    pub penalties_after: Vec<f64>,
    /// `‖dF‖₂` over the `3|V|` coordinates.
    pub gradient_norm: f64,
    /// `⟨g, x⟩`, the squared norm of the gradient in the metric.
    pub sobolev_norm_sq: f64,
    /// `‖x‖_∞` (largest vertex coordinate of the direction).
    pub direction_max: f64,
    /// `‖x‖_F` over the `3|V|` coordinates.
    pub direction_norm: f64,
    /// The accepted step size (0 when stuck).
    pub tau: f64,
    /// The first trial's step size (after the CCD cap).
    pub tau0: f64,
    /// The CCD's `τ_max` (when CCD is on).
    pub tau_ccd: Option<f64>,
    /// The CCD report.
    pub ccd: Option<CcdReport>,
    /// Number of halvings.
    pub halvings: u32,
    /// Every trial.
    pub trials: Vec<Trial>,
    /// Newton rounds kept by the projection.
    pub newton_steps: u32,
    /// How the projection ended.
    pub projection_status: ProjectionStatus,
    /// Projections that were not admissible (`LineSearch::guard_projection`): rejected trials,
    /// plus one if the stuck step's projection had to be limited.
    pub projection_blocked: u32,
    /// The fraction of the projection's motion applied (1 unless `ProjectionStatus::Blocked`).
    pub projection_fraction: f64,
    /// `max |Φ_r|` over the active rows after the projection.
    pub max_violation: f64,
    /// Constraint rows dropped for this step: by the rank check and from the Schur complement
    /// (indices into the step's constraint rows).
    pub dropped_rows: Vec<usize>,
    /// The linear solves.
    pub linear: LinearReport,
    /// Hierarchical back end: number of `A′` products and their seconds.
    pub products: Option<(usize, f64)>,
    /// The remeshing (when on).
    pub remesh: Option<RemeshReport>,
    /// The remeshing pass was computed but not applied: its moves met a CCD obstacle
    /// (`LineSearch::guard_projection`).
    pub remesh_rejected: bool,
    /// `|V|` after the step (and the remeshing).
    pub vertices: usize,
    /// `|F|` after the step.
    pub faces: usize,
    /// Wall times.
    pub seconds: Timings,
}

/// The report of a run.
#[derive(Debug, Clone, PartialEq)]
pub struct FlowReport {
    /// Why the run ended.
    pub status: Status,
    /// Every step of the run.
    pub steps: Vec<StepReport>,
    /// `F` at the end (after the last step's projection; the initial `F` without steps).
    pub final_energy: f64,
    /// Steps of the session so far.
    pub iterations: usize,
    /// Wall time of the run.
    pub seconds: f64,
}

impl FlowReport {
    /// The error that ended the run, if any.
    pub fn error(&self) -> Option<&RsError> {
        match &self.status {
            Status::Error(e) => Some(e),
            _ => None,
        }
    }
}

/// The result of [`flow`].
#[derive(Debug, Clone)]
pub struct FlowResult {
    /// The final surface (it changes with remeshing).
    pub surface: Surface,
    /// The final positions.
    pub positions: Vec<[f64; 3]>,
    /// The run's report.
    pub report: FlowReport,
}

// =============================================================================================
// Session
// =============================================================================================

/// A constraint with explicit targets (and explicit components for per-component
/// barycenters), remembering which of the caller's constraints it came from.
#[derive(Debug, Clone)]
struct Resolved {
    origin: usize,
    constraint: Constraint,
}

/// The initial targets of a scheduled constraint.
#[derive(Debug, Clone)]
enum Scheduled {
    Scalar {
        resolved: usize,
        initial: f64,
        last: f64,
        steps: usize,
    },
    Pins {
        resolved: usize,
        initial: Vec<[f64; 3]>,
        displacement: Vec<[f64; 3]>,
        steps: usize,
    },
}

/// A flow in progress (`SPEC.md` §7.1, §9.1): the surface, the positions and the per-vertex
/// data (pins, squared-error targets, attractor sets) between steps. Cancellable between steps
/// and between line-search trials ([`Session::cancel_flag`]).
pub struct Session {
    options: FlowOptions,
    eff: FlowOptions,
    mode: BarycenterMode,
    pool: Arc<rayon::ThreadPool>,
    surface: Surface,
    positions: Vec<[f64; 3]>,
    penalties: Vec<Penalty>,
    resolved: Vec<Resolved>,
    scheduled: Vec<Scheduled>,
    ccd_obstacles: Vec<CcdObstacle>,
    face_origin: Vec<u32>,
    pinned: Vec<u32>,
    l0: f64,
    iteration: usize,
    last_tau: Option<f64>,
    last_energy: f64,
    first_sobolev: Option<f64>,
    accepted_energies: VecDeque<f64>,
    stuck_count: usize,
    status: Status,
    cancel: Arc<AtomicBool>,
    method_state: std::sync::Mutex<Option<MethodState>>,
}

impl std::fmt::Debug for Session {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("Session")
            .field("vertices", &self.surface.n_vertices())
            .field("faces", &self.surface.n_faces())
            .field("iteration", &self.iteration)
            .field("status", &self.status)
            .finish_non_exhaustive()
    }
}

/// Maps the constraint index of a constraint-set error back to the caller's list.
fn map_constraint_error(e: RsError, resolved: &[Resolved]) -> RsError {
    let origin = |i: usize| resolved.get(i).map_or(i, |r| r.origin);
    match e {
        RsError::InfeasibleConstraint {
            constraint,
            violation,
        } => RsError::InfeasibleConstraint {
            constraint: origin(constraint),
            violation,
        },
        RsError::ConstraintOutOfRange { constraint, index } => RsError::ConstraintOutOfRange {
            constraint: origin(constraint),
            index,
        },
        e => e,
    }
}

fn pinned_of(resolved: &[Resolved]) -> Vec<u32> {
    let mut seen = std::collections::HashSet::new();
    let mut out = Vec::new();
    for r in resolved {
        if let Constraint::Pins { vertices, .. } = &r.constraint {
            for &v in vertices {
                if seen.insert(v) {
                    out.push(v);
                }
            }
        }
    }
    out
}

fn flatten(v: &[[f64; 3]]) -> Vec<f64> {
    v.iter().flatten().copied().collect()
}

fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b).map(|(x, y)| x * y).sum()
}

impl Session {
    /// Validates the options against the surface, resolves the constraints' targets (`None` =
    /// the initial values), and evaluates `F` once (a non-finite energy at the start is an
    /// error, `SPEC.md` §2.3).
    ///
    /// # Errors
    /// [`FlowOptions::validate`]; the errors of [`Surface::validate_positions`], of
    /// [`ConstraintSet::new`] (indices, orientation, translations), of the penalties'
    /// validation; [`RsError::ConstraintOutOfRange`] / [`RsError::InvalidOption`] for a schedule
    /// that does not fit its constraint; [`RsError::NonFiniteEnergy`] or
    /// [`RsError::InvalidOption`] (a penalty) for a non-finite `F`.
    pub fn new(
        surface: Surface,
        positions: Vec<[f64; 3]>,
        options: FlowOptions,
    ) -> Result<Session, RsError> {
        options.validate()?;
        let eff = options.effective();
        surface.validate_positions(&positions)?;
        for p in &eff.penalties {
            p.validate(&surface)?;
        }
        let pool = Arc::new(make_pool(eff.threads)?);
        let mode = eff.resolved_barycenter_mode();
        let set = ConstraintSet::new(&surface, &positions, &eff.constraints, mode)?;
        // Resolve: explicit targets per row group.
        let mut resolved: Vec<Resolved> = Vec::new();
        let rows = set.rows();
        let mut r = 0;
        while r < rows.len() {
            let row = rows[r];
            match row.kind {
                RowKind::Barycenter { component, .. } => {
                    resolved.push(Resolved {
                        origin: row.constraint,
                        constraint: Constraint::Barycenter {
                            component,
                            target: Some([rows[r].target, rows[r + 1].target, rows[r + 2].target]),
                        },
                    });
                    r += 3;
                }
                RowKind::Pin { vertex, .. } => {
                    let t = [rows[r].target, rows[r + 1].target, rows[r + 2].target];
                    let append = matches!(resolved.last(), Some(Resolved { origin, constraint: Constraint::Pins { .. } }) if *origin == row.constraint);
                    if !append {
                        resolved.push(Resolved {
                            origin: row.constraint,
                            constraint: Constraint::Pins {
                                vertices: Vec::new(),
                                targets: Some(Vec::new()),
                            },
                        });
                    }
                    if let Some(Resolved {
                        constraint:
                            Constraint::Pins {
                                vertices,
                                targets: Some(targets),
                            },
                        ..
                    }) = resolved.last_mut()
                    {
                        vertices.push(vertex);
                        targets.push(t);
                    }
                    r += 3;
                }
                RowKind::Area { component } => {
                    resolved.push(Resolved {
                        origin: row.constraint,
                        constraint: Constraint::Area {
                            component,
                            target: Some(row.target),
                        },
                    });
                    r += 1;
                }
                RowKind::Volume { component } => {
                    resolved.push(Resolved {
                        origin: row.constraint,
                        constraint: Constraint::Volume {
                            component,
                            target: Some(row.target),
                        },
                    });
                    r += 1;
                }
            }
        }
        // Schedules.
        let mut scheduled = Vec::new();
        for s in &eff.schedules {
            if s.constraint >= eff.constraints.len() {
                return Err(RsError::ConstraintOutOfRange {
                    constraint: s.constraint,
                    index: s.constraint as u32,
                });
            }
            let j = resolved.iter().position(|r| r.origin == s.constraint);
            let (j, t0) = match j.map(|j| (j, &resolved[j].constraint)) {
                Some((j, Constraint::Area { target: Some(t), .. }))
                | Some((j, Constraint::Volume { target: Some(t), .. })) => (j, *t),
                _ => {
                    return Err(RsError::InvalidOption {
                        name: "schedule.constraint",
                        value: s.constraint.to_string(),
                        requirement: "the index of an Area or Volume constraint",
                    });
                }
            };
            scheduled.push(Scheduled::Scalar {
                resolved: j,
                initial: t0,
                last: s.factor * t0 + s.offset,
                steps: s.steps,
            });
        }
        for s in &eff.pin_schedules {
            if s.constraint >= eff.constraints.len() {
                return Err(RsError::ConstraintOutOfRange {
                    constraint: s.constraint,
                    index: s.constraint as u32,
                });
            }
            let Constraint::Pins { vertices, .. } = &eff.constraints[s.constraint] else {
                return Err(RsError::InvalidOption {
                    name: "pin_schedule.constraint",
                    value: s.constraint.to_string(),
                    requirement: "the index of a Pins constraint",
                });
            };
            if s.displacement.len() != vertices.len() {
                return Err(RsError::InvalidOption {
                    name: "pin_schedule.displacement",
                    value: format!("{} displacements for {} pins", s.displacement.len(), vertices.len()),
                    requirement: "one displacement per pin",
                });
            }
            if let Some(j) = resolved.iter().position(|r| r.origin == s.constraint)
                && let Constraint::Pins {
                    targets: Some(t), ..
                } = &resolved[j].constraint
            {
                scheduled.push(Scheduled::Pins {
                    resolved: j,
                    initial: t.clone(),
                    displacement: s.displacement.clone(),
                    steps: s.steps,
                });
            }
        }
        let mut ccd_obstacles = Vec::new();
        if let Some(o) = &eff.line_search.ccd {
            for ob in &eff.ccd_obstacles {
                ccd_obstacles.push(CcdObstacle::new(ob, o.leaf_size)?);
            }
        }
        let l0 = mean_edge_length(&surface, &positions);
        let pinned = pinned_of(&resolved);
        let method_state = (eff.method != Method::ProjectedGradient).then(|| MethodState::new(eff.method));
        let mut s = Session {
            penalties: eff.penalties.clone(),
            options,
            eff,
            mode,
            pool,
            face_origin: (0..surface.n_faces() as u32).collect(),
            surface,
            positions,
            resolved,
            scheduled,
            ccd_obstacles,
            pinned,
            l0,
            iteration: 0,
            last_tau: None,
            last_energy: f64::NAN,
            first_sobolev: None,
            accepted_energies: VecDeque::new(),
            stuck_count: 0,
            method_state: std::sync::Mutex::new(method_state),
            status: Status::Running,
            cancel: Arc::new(AtomicBool::new(false)),
        };
        let f0 = {
            let mut obj = s.objective()?;
            let ev = obj.evaluate(&s.positions, false)?;
            if let Some(i) = ev.penalties.iter().position(|v| !v.is_finite()) {
                return Err(RsError::InvalidOption {
                    name: "penalties",
                    value: format!("penalty {i} ({}) = {}", s.penalties[i].name(), ev.penalties[i]),
                    requirement: "a finite value at the initial positions",
                });
            }
            ev.total
        };
        if !f0.is_finite() {
            return Err(RsError::NonFiniteEnergy {
                faces: [u32::MAX, u32::MAX],
            });
        }
        s.last_energy = f0;
        Ok(s)
    }

    fn objective(&self) -> Result<Objective<'_>, RsError> {
        Objective::build(
            &self.pool,
            &self.surface,
            &self.positions,
            self.eff.exponents,
            self.eff.energy,
            self.eff.energy_eval,
            self.eff.tp_weight,
            &self.penalties,
        )
    }

    /// The current surface (it changes with remeshing).
    pub fn surface(&self) -> &Surface {
        &self.surface
    }

    /// The current positions.
    pub fn positions(&self) -> &[[f64; 3]] {
        &self.positions
    }

    /// The pinned vertices (current ids, in the order of the constraints, each once).
    pub fn pinned(&self) -> &[u32] {
        &self.pinned
    }

    /// The flag checked before every step and every line-search trial; set it to cancel.
    pub fn cancel_flag(&self) -> Arc<AtomicBool> {
        self.cancel.clone()
    }

    /// The status.
    pub fn status(&self) -> Status {
        self.status.clone()
    }

    /// The options as given.
    pub fn options(&self) -> &FlowOptions {
        &self.options
    }

    /// The options in use ([`FlowOptions::effective`]).
    pub fn effective_options(&self) -> &FlowOptions {
        &self.eff
    }

    /// The resolved barycenter mode.
    pub fn barycenter_mode(&self) -> BarycenterMode {
        self.mode
    }

    /// What the comparison optimizer did so far (`None` for projected gradient).
    pub fn method_diagnostics(&self) -> Option<MethodDiagnostics> {
        self.method_state.lock().unwrap_or_else(std::sync::PoisonError::into_inner).as_ref().map(MethodState::diagnostics)
    }

    /// Steps taken.
    pub fn iteration(&self) -> usize {
        self.iteration
    }

    /// `F` after the last step's projection (the initial `F` before the first step).
    pub fn energy(&self) -> f64 {
        self.last_energy
    }

    /// `ℓ̄₀`, the mean edge length of the initial mesh (remeshing's `L₀`, the line search's unit).
    pub fn mean_edge_length0(&self) -> f64 {
        self.l0
    }

    /// For every current face, the face of the initial surface it descends from.
    pub fn face_origin(&self) -> &[u32] {
        &self.face_origin
    }

    /// The current penalties (squared-error targets and attractor sets follow remeshing).
    pub fn penalties(&self) -> &[Penalty] {
        &self.penalties
    }

    /// The constraints with their current targets (barycenters per component when resolved
    /// so; pins with their current ids).
    pub fn constraints(&self) -> Vec<Constraint> {
        self.resolved.iter().map(|r| r.constraint.clone()).collect()
    }

    fn schedule_moving(&self, k: usize) -> bool {
        self.scheduled.iter().any(|s| match s {
            Scheduled::Scalar { steps, .. } | Scheduled::Pins { steps, .. } => k < *steps,
        })
    }

    /// The constraints with the schedules' targets of step `k`.
    fn constraints_at(&self, k: usize) -> Vec<Resolved> {
        let mut out = self.resolved.clone();
        for s in &self.scheduled {
            match s {
                Scheduled::Scalar {
                    resolved,
                    initial,
                    last,
                    steps,
                } => {
                    let p = k.min(*steps) as f64 / *steps as f64;
                    let t = initial + p * (last - initial);
                    match &mut out[*resolved].constraint {
                        Constraint::Area { target, .. } | Constraint::Volume { target, .. } => {
                            *target = Some(t);
                        }
                        _ => {}
                    }
                }
                Scheduled::Pins {
                    resolved,
                    initial,
                    displacement,
                    steps,
                } => {
                    let p = k.min(*steps) as f64 / *steps as f64;
                    if let Constraint::Pins {
                        targets: Some(t), ..
                    } = &mut out[*resolved].constraint
                    {
                        for (i, ti) in t.iter_mut().enumerate() {
                            for a in 0..3 {
                                ti[a] = initial[i][a] + p * displacement[i][a];
                            }
                        }
                    }
                }
            }
        }
        out
    }

    /// One step of `SPEC.md` §7.1. The session's status is updated (converged, stuck, max
    /// iterations); a failing step sets [`Status::Error`] and returns the error,
    /// [`RsError::Cancelled`] sets [`Status::Cancelled`]. The surface stays that of the last
    /// completed step on any error.
    ///
    /// # Errors
    /// [`RsError::Cancelled`]; [`RsError::NonFiniteEnergy`]; constraint, solver and remeshing
    /// errors ([`RsError::Internal`] if a remeshed surface fails its revalidation).
    pub fn step(&mut self) -> Result<StepReport, RsError> {
        let e_before = self.last_energy;
        if let Some(m) = self.method_state.get_mut().unwrap_or_else(std::sync::PoisonError::into_inner).as_mut() {
            m.begin_step(&mut self.positions);
        }
        let r = self.step_inner();
        if let Some(m) = self.method_state.get_mut().unwrap_or_else(std::sync::PoisonError::into_inner).as_mut() {
            match &r {
                Ok(rep) => m.end_step(
                    rep.remesh.as_ref().is_some_and(|x| x.splits + x.collapses > 0),
                    e_before,
                    self.last_energy,
                ),
                Err(_) => m.abort_step(&mut self.positions),
            }
        }
        match &r {
            Err(RsError::Cancelled) => self.status = Status::Cancelled,
            Err(e) => self.status = Status::Error(e.clone()),
            Ok(_) => {}
        }
        r
    }

    #[allow(clippy::too_many_lines)]
    fn step_inner(&mut self) -> Result<StepReport, RsError> {
        let t_step = Instant::now();
        if self.cancel.load(Ordering::Relaxed) {
            return Err(RsError::Cancelled);
        }
        let iteration = self.iteration + 1;
        let mut seconds = Timings::default();
        let eff = &self.eff;
        let nv = self.surface.n_vertices();

        // 0. The constraints of this step (schedules, D11).
        let target_step = match eff.schedule_timing {
            ScheduleTiming::BeforeProjection => iteration,
            ScheduleTiming::AfterProjection => iteration - 1,
        };
        let mut resolved = self.constraints_at(target_step);
        let t = Instant::now();
        let list: Vec<Constraint> = resolved.iter().map(|r| r.constraint.clone()).collect();
        let mut set = ConstraintSet::new(&self.surface, &self.positions, &list, self.mode)
            .map_err(|e| map_constraint_error(e, &resolved))?;
        seconds.factor += t.elapsed().as_secs_f64();

        // 1. Trees and lists; F, dF.
        let mut obj = self.objective()?;
        seconds.bvh = obj.build_seconds;
        let t = Instant::now();
        let ev0 = obj.evaluate(&self.positions, true)?;
        seconds.energy = t.elapsed().as_secs_f64();
        let f0 = ev0.total;
        let gv = ev0.differential.expect("asked for the differential");
        if !f0.is_finite() || gv.iter().flatten().any(|x| !x.is_finite()) {
            if let Some(i) = ev0.penalties.iter().position(|v| !v.is_finite()) {
                return Err(RsError::InvalidOption {
                    name: "penalties",
                    value: format!("penalty {i} ({}) = {}", self.penalties[i].name(), ev0.penalties[i]),
                    requirement: "a finite value at the step's positions",
                });
            }
            return Err(RsError::NonFiniteEnergy {
                faces: [u32::MAX, u32::MAX],
            });
        }
        let g = flatten(&gv);
        let gradient_norm = dot(&g, &g).sqrt();

        // 2. The solver of A′ and the direction (§6.1).
        let t = Instant::now();
        let (weights, extras) = match eff.metric {
            Metric::Sobolev { low_order } => {
                let ww = obj.willmore_weight();
                let willmore = if ww > 0.0 {
                    Some(willmore_matrix(&self.surface, &self.positions, true)?)
                } else {
                    None
                };
                (
                    SobolevWeights {
                        tp_weight: eff.tp_weight,
                        low_order,
                        willmore_weight: ww,
                    },
                    MetricExtras {
                        willmore,
                        obstacle: obj.obstacle_diagonal(&self.positions, low_order, self.options.reference_compat)?,
                    },
                )
            }
            _ => (SobolevWeights::default(), MetricExtras::default()),
        };
        let qrows = set.positional_rows();
        let nrows = set.nonlinear_rows();
        let c_q = set.jacobian(&self.positions, &qrows)?;
        let c_n = set.jacobian(&self.positions, &nrows)?;
        let bh_leaf = match eff.energy_eval {
            EnergyEval::BarnesHut(o) if eff.energy == Energy::TangentPoint => Some(o.leaf_size),
            _ => None,
        };
        let mut solver = StepSolver::build(
            &self.pool,
            &self.surface,
            &self.positions,
            &eff.exponents,
            eff.metric,
            &weights,
            &extras,
            &c_q,
            &eff.backend,
            obj.bh_tree().map(|t| t.face_tree()),
            bh_leaf,
        )?;
        seconds.factor += t.elapsed().as_secs_f64();
        let t = Instant::now();
        let mut guard = self.method_state.lock().unwrap_or_else(std::sync::PoisonError::into_inner);
        let quasi_newton = guard.is_some() && solver.as_sparse().is_some();
        let (schur, dir_res) = if quasi_newton {
            let schur = SchurComplement::new(solver.as_dyn(), nrows.clone(), c_n)?;
            let (Some(m), Some(sp)) = (guard.as_mut(), solver.as_sparse()) else {
                unreachable!("checked above")
            };
            let geom = crate::geometry::FaceGeometry::new(&self.surface, &self.positions)?;
            let area: f64 = geom.areas.iter().sum();
            let d = m.direction(sp, &self.positions, &g, &c_q, area);
            (schur, d)
        } else {
            // The Schur columns and the direction's solve as one batch (same results).
            let (schur, d) = SchurComplement::with_direction(solver.as_dyn(), nrows.clone(), c_n, &g)?;
            (schur, Ok((d.x, d.record)))
        };
        drop(guard);
        let (dir_x, dir_record) = dir_res?;
        seconds.solve = t.elapsed().as_secs_f64();
        let mut linear = solver.as_dyn().base_report();
        if linear.factor_seconds == 0.0 {
            linear.factor_seconds = seconds.factor;
        }
        schur.fill_report(&mut linear);
        linear.solves.push(dir_record);
        linear.solve_seconds = seconds.solve;
        let mut x = dir_x;
        if eff.barycenter_shift {
            match self.mode {
                BarycenterMode::Global => {
                    let geom = crate::geometry::FaceGeometry::new(&self.surface, &self.positions)?;
                    let va = crate::geometry::vertex_areas(&self.surface, &geom.areas);
                    let tot: f64 = va.iter().sum();
                    let mut m = [0.0; 3];
                    for v in 0..nv {
                        for a in 0..3 {
                            m[a] += va[v] * gv[v][a];
                        }
                    }
                    for v in 0..nv {
                        for a in 0..3 {
                            x[3 * v + a] += m[a] / tot;
                        }
                    }
                }
                BarycenterMode::PerComponent => {
                    let means = component_means(&self.surface, &self.positions, &gv)?;
                    let comp = self.surface.vertex_components();
                    for v in 0..nv {
                        for a in 0..3 {
                            x[3 * v + a] += means[comp[v] as usize][a];
                        }
                    }
                }
            }
        }
        let slope = dot(&g, &x);
        let xinf = x.iter().fold(0.0f64, |m, v| m.max(v.abs()));
        let xnorm = dot(&x, &x).sqrt();
        let max_disp = (0..nv)
            .map(|v| (x[3 * v] * x[3 * v] + x[3 * v + 1] * x[3 * v + 1] + x[3 * v + 2] * x[3 * v + 2]).sqrt())
            .fold(0.0f64, f64::max);

        // 3. Line search (§7.2) with the CCD cap (§7.3).
        let mut tau0 = start_tau(&eff.line_search.start, self.last_tau, xnorm, max_disp, self.l0);
        // Quasi-Newton directions with stored pairs have a natural unit step (WP9).
        if self
            .method_state
            .lock()
            .unwrap_or_else(std::sync::PoisonError::into_inner)
            .as_ref()
            .is_some_and(MethodState::unit_step)
        {
            tau0 = 1.0;
        }
        let mut ccd_report = None;
        // The faces that intersect at the start of the step exempt their neighbourhood from the
        // gap rule (SPEC 7.3); `None` for an embedded mesh.
        let mut band: Option<IntersectionBand> = None;
        if eff.line_search.ccd.is_some() {
            let t = Instant::now();
            let b = IntersectionBand::new(&self.pool, &self.surface, &self.positions)?;
            if !b.is_empty() {
                band = Some(b);
            }
            seconds.ccd = t.elapsed().as_secs_f64();
        }
        if let Some(copts) = &eff.line_search.ccd
            && tau0.is_finite()
            && tau0 > 0.0
        {
            let t = Instant::now();
            let mut o: CcdOptions = *copts;
            if o.reference_length.is_none() {
                o.reference_length = Some(self.l0);
            }
            let velocity: Vec<[f64; 3]> = (0..nv)
                .map(|v| [-x[3 * v], -x[3 * v + 1], -x[3 * v + 2]])
                .collect();
            let rep = first_contact_with_band(
                &self.pool,
                &self.surface,
                &self.positions,
                &velocity,
                tau0,
                &self.ccd_obstacles,
                &o,
                band.as_ref(),
            )?;
            if rep.contact.is_some() {
                tau0 = tau0.min(2.0 / 3.0 * rep.tau_max);
            }
            ccd_report = Some(rep);
            seconds.ccd += t.elapsed().as_secs_f64();
        }
        let t = Instant::now();
        let tau_min = eff.line_search.tau_min_rel * self.l0 / xinf;
        let gnorm = gradient_norm;
        let mut trials = Vec::new();
        let mut halvings = 0u32;
        // The projection guard (§6.2, §7.3): an Armijo-accepted trial is kept only if its
        // projection is admissible too (finite F, no contact along the projection's motion).
        let ccd_opts = eff.line_search.ccd.map(|mut o| {
            if o.reference_length.is_none() {
                o.reference_length = Some(self.l0);
            }
            o
        });
        let mut pguard = ProjectionGuard {
            enabled: eff.line_search.guard_projection,
            pool: &self.pool,
            surface: &self.surface,
            ccd: ccd_opts,
            obstacles: &self.ccd_obstacles,
            band: band.as_ref(),
            schur: (!nrows.is_empty()).then_some(&schur),
            newton: eff.newton_steps,
            barycenter_shift: eff.barycenter_shift,
            n_penalties: self.penalties.len(),
            base_measure: 0.0,
        };
        pguard.base_measure = pguard.measure(&set, &set.values(&self.positions)?);
        let mut projection_blocked = 0u32;
        let mut projection_fraction = 1.0;
        // Whether the projection of the unmoved positions already meets a contact (computed
        // when a trial's projection does): then halving cannot help.
        let mut base_contact: Option<bool> = None;
        let mut accepted: Option<AcceptedTrial> = None;
        let mut t_projection = 0.0;
        if slope > 0.0 && tau0.is_finite() && tau0 > 0.0 {
            let mut tau = tau0;
            loop {
                if !(tau >= tau_min) || eff.line_search.max_halvings.is_some_and(|m| halvings > m) {
                    break;
                }
                if self.cancel.load(Ordering::Relaxed) {
                    return Err(RsError::Cancelled);
                }
                let trial: Vec<[f64; 3]> = (0..nv)
                    .map(|v| {
                        let p = self.positions[v];
                        [
                            p[0] - tau * x[3 * v],
                            p[1] - tau * x[3 * v + 1],
                            p[2] - tau * x[3 * v + 2],
                        ]
                    })
                    .collect();
                let ft = match obj.evaluate(&trial, false) {
                    Ok(e) => e.total,
                    Err(e) if is_geometric(&e) => f64::NAN,
                    Err(e) => return Err(e),
                };
                let ok = ft.is_finite()
                    && match eff.line_search.armijo {
                        ArmijoRule::Standard { c } => ft <= f0 - c * tau * slope,
                        ArmijoRule::Reference { c } => f0 - ft >= c * tau * slope / gnorm,
                    };
                let mut outcome = if !ft.is_finite() {
                    TrialOutcome::NonFinite
                } else if ok {
                    TrialOutcome::Accepted
                } else {
                    TrialOutcome::ArmijoFailed
                };
                if ok {
                    let tp = Instant::now();
                    let checked = pguard.project(&mut set, &mut resolved, &mut obj, &trial)?;
                    if checked.admissible || !pguard.enabled {
                        t_projection += tp.elapsed().as_secs_f64();
                        trials.push(Trial { tau, energy: ft, outcome });
                        accepted = Some((tau, ft, trial, checked.projection, checked.after));
                        break;
                    }
                    projection_blocked += 1;
                    if checked.contact {
                        let base = match base_contact {
                            Some(b) => b,
                            None => {
                                let c = pguard.project(&mut set, &mut resolved, &mut obj, &self.positions)?;
                                let b = !c.admissible && c.contact;
                                base_contact = Some(b);
                                b
                            }
                        };
                        if base {
                            // The contact blocks the projection whatever the step: keep the
                            // descent (CCD-capped above) and apply the admissible part of the
                            // projection.
                            let (p, a, frac) = pguard.limit(&set, &mut obj, &trial, checked)?;
                            projection_fraction = frac;
                            t_projection += tp.elapsed().as_secs_f64();
                            trials.push(Trial { tau, energy: ft, outcome });
                            accepted = Some((tau, ft, trial, p, a));
                            break;
                        }
                    }
                    t_projection += tp.elapsed().as_secs_f64();
                    outcome = TrialOutcome::ProjectionBlocked;
                }
                trials.push(Trial {
                    tau,
                    energy: if ft.is_finite() { ft } else { f64::NAN },
                    outcome,
                });
                tau *= 0.5;
                halvings += 1;
            }
        }
        seconds.line_search = t.elapsed().as_secs_f64() - t_projection;

        // 4. Projection (§6.2), done above for the accepted trial. A stuck step projects the
        // unmoved positions; if that projection is not admissible either, only the largest
        // admissible fraction of its motion is applied (`ProjectionStatus::Blocked`).
        let t = Instant::now();
        let (status, tau, energy_trial, proj, after) = match accepted {
            Some((tau, ft, _, p, a)) => (StepStatus::Accepted, tau, ft, p, a),
            None => {
                let moved = self.positions.clone();
                let checked = pguard.project(&mut set, &mut resolved, &mut obj, &moved)?;
                if checked.admissible || !pguard.enabled {
                    (StepStatus::Stuck, 0.0, f0, checked.projection, checked.after)
                } else {
                    projection_blocked += 1;
                    let (p, a, frac) = pguard.limit(&set, &mut obj, &moved, checked)?;
                    projection_fraction = frac;
                    (StepStatus::Stuck, 0.0, f0, p, a)
                }
            }
        };
        seconds.projection = t.elapsed().as_secs_f64() + t_projection;
        let products = solver.products();
        if let Some(t) = solver.build_timings() {
            seconds.metric_tree = t[0];
            seconds.fractional_tree = t[1];
            seconds.lprime = t[2];
            seconds.preconditioner = t[3];
        }
        let mut dropped_rows: Vec<usize> = set.dropped_rows().to_vec();
        dropped_rows.extend(proj.dropped_rows.iter().copied());
        dropped_rows.sort_unstable();
        dropped_rows.dedup();
        drop(obj);

        // 5. Remeshing (§10) of the projected positions, with the per-vertex data.
        let t = Instant::now();
        let mut remesh_report = None;
        let mut remesh_rejected = false;
        let mut new_state = None;
        if let Some(ropts) = &eff.remesh {
            let pinned = pinned_of(&resolved);
            let mut guard = match &eff.line_search.ccd {
                Some(o) => {
                    let mut g = CcdCollapseGuard::new(o, self.l0)?;
                    if let Some(b) = &band {
                        g.set_band(b);
                    }
                    Some(g)
                }
                // Without CCD the collapses are still checked (default gap, with the band of
                // the faces that intersect now): unguarded collapses on a self-intersecting
                // subcritical flow raised F by up to 1e6 in one pass (SOLVERFIX §4).
                None if ropts.guard_collapses => {
                    let mut g = CcdCollapseGuard::new(&CcdOptions::default(), self.l0)?;
                    let b = IntersectionBand::new(&self.pool, &self.surface, &proj.positions)?;
                    g.set_band(&b);
                    Some(g)
                }
                None => None,
            };
            let (out, rep) = remesh(
                &self.surface,
                &proj.positions,
                &pinned,
                self.l0,
                ropts,
                guard.as_mut().map(|g| g as &mut dyn CollapseGuard),
            )
            .map_err(|e| match e {
                e @ (RsError::InvalidOption { .. } | RsError::Cancelled) => e,
                e => RsError::Internal {
                    what: "remeshing",
                    detail: format!(
                        "step {iteration}: the remeshed surface failed its revalidation: {e} \
                         ({} vertices, {} faces before)",
                        self.surface.n_vertices(),
                        self.surface.n_faces()
                    ),
                },
            })?;
            // The remeshing's moves (split midpoints, collapses, smoothing) must not take the
            // surface into a CCD obstacle either: the motion from the old surface (the new
            // vertices at their vertex-map positions on it) to the new positions is checked
            // against the obstacles; a pass that meets one is not applied this step.
            let keep = match &pguard.ccd {
                Some(o) if pguard.enabled && !self.ccd_obstacles.is_empty() => {
                    let start = out.vertex_map.transport_vec3(&proj.positions)?;
                    let velocity: Vec<[f64; 3]> = start
                        .iter()
                        .zip(&out.positions)
                        .map(|(a, b)| [b[0] - a[0], b[1] - a[1], b[2] - a[2]])
                        .collect();
                    let all = IntersectionBand::all(out.surface.n_vertices());
                    let r = first_contact_with_band(
                        &self.pool,
                        &out.surface,
                        &start,
                        &velocity,
                        1.0,
                        &self.ccd_obstacles,
                        o,
                        Some(&all),
                    )?;
                    r.contact.is_none()
                }
                _ => true,
            };
            if keep {
                new_state = Some(self.transport(out, &pinned, &mut resolved)?);
            } else {
                remesh_rejected = true;
            }
            remesh_report = Some(rep);
        }
        seconds.remesh = t.elapsed().as_secs_f64();

        // 6. Commit.
        match new_state {
            Some((surface, positions, penalties, face_origin)) => {
                self.surface = surface;
                self.positions = positions;
                self.penalties = penalties;
                self.face_origin = face_origin;
            }
            None => self.positions = proj.positions,
        }
        // Keep the non-scheduled targets (barycenter shift) and the transported ids.
        self.resolved = resolved;
        self.pinned = pinned_of(&self.resolved);
        self.iteration = iteration;
        self.last_energy = after.total;
        if status == StepStatus::Accepted {
            self.last_tau = Some(tau);
            self.stuck_count = 0;
        } else {
            self.stuck_count += 1;
        }
        if self.first_sobolev.is_none() {
            self.first_sobolev = Some(slope);
        }

        // 7. Stop rules (§7.4).
        let moving = self.schedule_moving(iteration);
        let converged = !moving
            && match self.eff.stop {
                StopRule::EnergyStall { relative, window } => {
                    if status == StepStatus::Accepted {
                        self.accepted_energies.push_back(after.total);
                        while self.accepted_energies.len() > window + 1 {
                            self.accepted_energies.pop_front();
                        }
                    }
                    self.accepted_energies.len() == window + 1 && {
                        let first = self.accepted_energies[0];
                        let last = self.accepted_energies[window];
                        first - last <= relative * last.abs()
                    }
                }
                StopRule::SobolevNorm { eps } => {
                    let s0 = self.first_sobolev.unwrap_or(slope);
                    s0 > 0.0 && (slope / s0).max(0.0).sqrt() < eps
                }
                StopRule::Off => false,
            };
        self.status = if self.stuck_count >= self.eff.stuck_limit {
            Status::Stuck
        } else if converged {
            Status::Converged
        } else if self.iteration >= self.eff.max_iterations {
            Status::MaxIterations
        } else {
            Status::Running
        };
        seconds.total = t_step.elapsed().as_secs_f64();
        Ok(StepReport {
            iteration,
            status,
            energy_before: f0,
            energy_trial,
            energy_after: after.total,
            energy_term_after: after.energy,
            penalties_after: after.penalties,
            gradient_norm,
            sobolev_norm_sq: slope,
            direction_max: xinf,
            direction_norm: xnorm,
            tau,
            tau0,
            tau_ccd: ccd_report.as_ref().map(|r| r.tau_max),
            ccd: ccd_report,
            halvings,
            trials,
            newton_steps: proj.rounds,
            projection_status: proj.status,
            projection_blocked,
            projection_fraction,
            max_violation: proj.max_violation,
            dropped_rows,
            linear,
            products,
            remesh: remesh_report,
            remesh_rejected,
            vertices: self.surface.n_vertices(),
            faces: self.surface.n_faces(),
            seconds,
        })
    }

    /// The new state after remeshing: the surface, positions, penalties with transported
    /// per-vertex data and the face origins; `resolved` gets the new component and pin ids.
    #[allow(clippy::type_complexity)]
    fn transport(
        &self,
        out: Remeshed,
        pinned: &[u32],
        resolved: &mut [Resolved],
    ) -> Result<(Surface, Vec<[f64; 3]>, Vec<Penalty>, Vec<u32>), RsError> {
        let old = &self.surface;
        let new = &out.surface;
        let internal = |detail: String| RsError::Internal {
            what: "remeshing transport",
            detail,
        };
        // Components: old → new through the face parents.
        let nc_old = old.n_components();
        let mut comp_map: Vec<Option<u32>> = vec![None; nc_old];
        for (f, &parent) in out.face_parents.iter().enumerate() {
            let cn = new.face_components()[f];
            let co = old.face_components()[parent as usize] as usize;
            match comp_map[co] {
                None => comp_map[co] = Some(cn),
                Some(c) if c == cn => {}
                Some(c) => {
                    return Err(internal(format!(
                        "old component {co} maps to new components {c} and {cn}"
                    )));
                }
            }
        }
        if new.n_components() != nc_old || comp_map.iter().any(Option::is_none) {
            return Err(internal(format!(
                "{nc_old} components before remeshing, {} after",
                new.n_components()
            )));
        }
        let map_c = |c: &mut Option<u32>| {
            if let Some(x) = c.as_mut() {
                *x = comp_map[*x as usize].expect("checked");
            }
        };
        let pin_map: HashMap<u32, u32> = pinned.iter().copied().zip(out.pinned.iter().copied()).collect();
        for r in resolved.iter_mut() {
            match &mut r.constraint {
                Constraint::Barycenter { component, .. }
                | Constraint::Area { component, .. }
                | Constraint::Volume { component, .. } => map_c(component),
                Constraint::Pins { vertices, .. } => {
                    for v in vertices.iter_mut() {
                        *v = *pin_map
                            .get(v)
                            .ok_or_else(|| internal(format!("pin {v} was not passed to the remesher")))?;
                    }
                }
            }
        }
        let mut penalties = self.penalties.clone();
        for p in penalties.iter_mut() {
            match p {
                Penalty::SquaredError { targets, .. } => {
                    *targets = out.vertex_map.transport_vec3(targets)?;
                }
                Penalty::ImplicitAttractor {
                    vertices: Some(vs), ..
                } => {
                    *vs = out.vertex_map.transport_set(vs, 0.5)?;
                }
                _ => {}
            }
        }
        let face_origin = out
            .face_parents
            .iter()
            .map(|&p| self.face_origin[p as usize])
            .collect();
        Ok((out.surface, out.positions, penalties, face_origin))
    }

    /// Steps until a stop rule, `max_iterations`, the time limit, `stuck_limit` stuck steps,
    /// cancellation (flag or `progress` returning [`Control::Stop`]) or a failing step.
    ///
    /// **Errors of steps do not make `run` fail** (as the curves crate): the run ends with
    /// [`Status::Error`] carrying the error (or [`Status::Cancelled`]), the report keeps the
    /// completed steps and the surface is that of the last completed step.
    ///
    /// # Errors
    /// None at present (the `Result` is the API of `SPEC.md` §9.1).
    pub fn run(
        &mut self,
        mut progress: Option<&mut dyn FnMut(&StepReport) -> Control>,
    ) -> Result<FlowReport, RsError> {
        let t0 = Instant::now();
        let mut steps = Vec::new();
        if self.status.is_terminal() && !matches!(self.status, Status::Error(_)) {
            self.status = Status::Running;
        }
        let status = loop {
            if self.cancel.load(Ordering::Relaxed) {
                self.status = Status::Cancelled;
                break Status::Cancelled;
            }
            if self.iteration >= self.eff.max_iterations {
                self.status = Status::MaxIterations;
                break Status::MaxIterations;
            }
            if self.eff.time_limit.is_some_and(|l| t0.elapsed() >= l) {
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

/// `τ₀` of a start rule (`SPEC.md` §7.2).
/// An accepted line-search trial: `τ`, `F` at the trial, the trial, its projection and `F`
/// after the projection.
type AcceptedTrial = (f64, f64, Vec<[f64; 3]>, Projection, Evaluation);

/// The projection of one line-search trial and whether it is admissible.
struct Checked {
    projection: Projection,
    after: Evaluation,
    admissible: bool,
    /// Not admissible because of a contact (CCD) or a non-finite `F`, rather than the
    /// constraint test.
    contact: bool,
    /// The CCD's `τ_max` along the projection's motion when it found a contact.
    tau_ccd: Option<f64>,
}

/// The projection with its admissibility test (`SPEC.md` §6.2, §7.3): after the projection,
/// `F` must be finite (implicit obstacles, degenerate faces) and, with CCD on, the linear motion
/// from the trial to the projected positions must keep the gap rule of the line search against
/// the surface itself and the CCD obstacles. The Newton projection onto area/volume targets
/// and the positional resets are not descent steps: without the test a growing area target
/// pushes a surface jammed against an obstacle through it (cobblestones, WP15).
struct ProjectionGuard<'a> {
    enabled: bool,
    pool: &'a rayon::ThreadPool,
    surface: &'a Surface,
    ccd: Option<CcdOptions>,
    obstacles: &'a [CcdObstacle],
    band: Option<&'a IntersectionBand>,
    schur: Option<&'a SchurComplement>,
    newton: crate::solve::NewtonSteps,
    barycenter_shift: bool,
    n_penalties: usize,
    /// `max_r |Φ_r| / tol_r` over the nonlinear rows at the step's base positions.
    base_measure: f64,
}

impl ProjectionGuard<'_> {
    /// `max_r |Φ_r| / tol_r` over all nonlinear rows of the step, including those the Schur
    /// complement dropped for this step (a dropped row is not corrected by the projection, so
    /// the line search must not let it drift; 0 without nonlinear rows).
    fn measure(&self, set: &ConstraintSet, phi: &[f64]) -> f64 {
        let Some(s) = self.schur else { return 0.0 };
        s.rows().iter().fold(0.0f64, |m, &r| {
            let v = phi[r].abs() / set.tolerance(r);
            if v.is_nan() || m.is_nan() { f64::NAN } else { m.max(v) }
        })
    }

    fn evaluate(&self, obj: &mut Objective<'_>, positions: &[[f64; 3]]) -> Result<Evaluation, RsError> {
        match obj.evaluate(positions, false) {
            Ok(e) => Ok(e),
            Err(e) if is_geometric(&e) => Ok(Evaluation {
                energy: f64::NAN,
                penalties: vec![f64::NAN; self.n_penalties],
                total: f64::NAN,
                differential: None,
            }),
            Err(e) => Err(e),
        }
    }

    /// Projects `moved` (the barycenter shift first moves the barycenter targets to `moved`)
    /// and tests the result.
    fn project(
        &self,
        set: &mut ConstraintSet,
        resolved: &mut [Resolved],
        obj: &mut Objective<'_>,
        moved: &[[f64; 3]],
    ) -> Result<Checked, RsError> {
        if self.barycenter_shift {
            let phi = set.values(moved)?;
            let rows = set.rows().to_vec();
            for (r, row) in rows.iter().enumerate() {
                if let RowKind::Barycenter { .. } = row.kind {
                    set.set_target(r, row.target + phi[r])?;
                }
            }
            for (r, row) in set.rows().iter().enumerate() {
                if let (RowKind::Barycenter { axis, .. }, Constraint::Barycenter { target: Some(t), .. }) =
                    (row.kind, &mut resolved[row.constraint].constraint)
                {
                    t[axis as usize] = set.rows()[r].target;
                }
            }
        }
        let projection = project(set, self.schur, moved, &ProjectionOptions { newton: self.newton })?;
        let after = self.evaluate(obj, &projection.positions)?;
        // The projection must not leave the kept nonlinear rows further from their targets
        // than the base positions are (the curves crate's `ProjectionFailed`): a lost area or
        // volume constraint lets the surface grow without bound (WP15, subcritical + remeshing).
        let m_after = self.measure(set, &projection.phi);
        let finite = after.total.is_finite();
        let mut admissible = finite && (!self.enabled || m_after <= self.base_measure.max(1.0));
        let mut contact = !finite;
        let mut tau_ccd = None;
        if self.enabled
            && finite
            && let Some(o) = &self.ccd
        {
            let velocity: Vec<[f64; 3]> = moved
                .iter()
                .zip(&projection.positions)
                .map(|(a, b)| [b[0] - a[0], b[1] - a[1], b[2] - a[2]])
                .collect();
            if velocity.iter().any(|v| v.iter().any(|&c| c != 0.0)) {
                let rep = first_contact_with_band(
                    self.pool,
                    self.surface,
                    moved,
                    &velocity,
                    1.0,
                    self.obstacles,
                    o,
                    self.band,
                )?;
                if rep.contact.is_some() {
                    admissible = false;
                    contact = true;
                    tau_ccd = Some(rep.tau_max);
                }
            }
        }
        Ok(Checked {
            projection,
            after,
            admissible,
            contact,
            tau_ccd,
        })
    }

    /// The largest admissible fraction of an inadmissible projection of `moved`: `⅔ τ_max` of
    /// its motion (as the line search's CCD cap), halved while `F` is not finite (at most 12
    /// times, then 0). The result has `ProjectionStatus::Blocked` and the values at the
    /// positions actually used.
    fn limit(
        &self,
        set: &ConstraintSet,
        obj: &mut Objective<'_>,
        moved: &[[f64; 3]],
        checked: Checked,
    ) -> Result<(Projection, Evaluation, f64), RsError> {
        let target = &checked.projection.positions;
        let mut frac = checked.tau_ccd.map_or(0.5, |t| (2.0 / 3.0 * t).min(1.0));
        let lerp = |f: f64| -> Vec<[f64; 3]> {
            moved
                .iter()
                .zip(target)
                .map(|(a, b)| [0, 1, 2].map(|k| a[k] + f * (b[k] - a[k])))
                .collect()
        };
        let mut chosen = None;
        for _ in 0..12 {
            if !(frac > 0.0) {
                break;
            }
            let pos = lerp(frac);
            let ev = self.evaluate(obj, &pos)?;
            if ev.total.is_finite() {
                chosen = Some((pos, ev));
                break;
            }
            frac *= 0.5;
        }
        let (positions, after) = match chosen {
            Some(x) => x,
            None => {
                frac = 0.0;
                let ev = self.evaluate(obj, moved)?;
                (moved.to_vec(), ev)
            }
        };
        let phi = set.values(&positions)?;
        let max_violation = (0..set.n_rows())
            .filter(|&r| set.is_active(r))
            .fold(0.0f64, |acc, r| acc.max(phi[r].abs()));
        let mut p = checked.projection;
        p.positions = positions;
        p.phi = phi;
        p.max_violation = max_violation;
        p.status = ProjectionStatus::Blocked;
        Ok((p, after, frac))
    }
}

fn start_tau(rule: &StartRule, previous: Option<f64>, xnorm: f64, max_disp: f64, l0: f64) -> f64 {
    match rule {
        StartRule::Reference => 2.0 / xnorm,
        StartRule::MaxDisplacement { lengths } => lengths * l0 / max_disp,
        StartRule::GrowPrevious { factor, first } => match previous {
            Some(t) => factor * t,
            None => start_tau(first, None, xnorm, max_disp, l0),
        },
        StartRule::Absolute { tau } => *tau,
    }
}

/// [`Session::new`] and [`Session::run`] in one call.
///
/// # Errors
/// The errors of [`Session::new`].
pub fn flow(
    surface: Surface,
    positions: Vec<[f64; 3]>,
    options: FlowOptions,
    progress: Option<&mut dyn FnMut(&StepReport) -> Control>,
) -> Result<FlowResult, RsError> {
    let mut s = Session::new(surface, positions, options)?;
    let report = s.run(progress)?;
    Ok(FlowResult {
        surface: s.surface,
        positions: s.positions,
        report,
    })
}
