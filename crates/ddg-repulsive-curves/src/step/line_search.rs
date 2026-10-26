//! Line searches (`SPEC.md` §7.2–7.3, D11, D12).
//!
//! Both are backtracking searches (Boyd & Vandenberghe, Alg. 9.2) along a descent direction:
//! the trial for step size `τ` is `γ_τ = γ − τ p`, projected onto the constraints (§6.2) to
//! `γ̂_τ`, and accepted when the projection succeeds and the Armijo condition holds on the
//! projected curve,
//!
//! ```text
//! F(γ̂_τ) ≤ F(γ) − c τ dF·p,
//! ```
//!
//! otherwise `τ ← shrink · τ`. After `max_halvings` shrinks, or below `1e-12 τ₀`, the search is
//! stuck. For a descent direction (`dF·p > 0`) and a projection that succeeds for small `τ`, the
//! search terminates (cited, `rc-armijo`).
//!
//! * [`LineSearch::CollisionGuarded`] (the default, D11): `τ_hi = diam(bbox γ) / max_v |p_v|`,
//!   `τ_max` = the first contact time of the linear motion `γ − τ p` on `[0, τ_hi]` by the CCD
//!   of §7.4 (∞ if none), start `τ₀` by the [`StartRule`] (the paper's `⅔ τ_max`, or `τ_hi` without contact; or
//!   grow-previous `min(⅔ τ_max, f τ_prev)`).
//!   A trial is accepted only if, in addition, every chord step of its projection (each a linear
//!   motion) is free of contacts. The move itself is free because `τ ≤ τ_max`. So the whole path
//!   of every accepted step is a piecewise linear path of embeddings, and no two edges of `𝒫`
//!   and no edge and obstacle ever cross: the isotopy class is kept (`rc-isotopy`, cited).
//! * [`LineSearch::NormalizedArmijo`] (Sec. 5.4, D12): the reference's start rule
//!   (`ref: src/tpe_flow_sc.cpp`, `LineSearchStep`): `n = ‖p‖₂`, `τ₀ = 1/n` if `n > 1` else
//!   `1/√n`, and if the previous accepted step `τ_prev > 1e-5`, `τ₀ = min(1.5 τ_prev, 4 τ₀)`.
//!   Standard Armijo test (the reference divides by `‖dE‖`). No collision guarantee.
//!
//! [formula: rc-armijo]
//! [formula: rc-isotopy]

use std::time::Instant;

use crate::ccd::{CcdOptions, CcdReport, Contact};
use crate::constraints::ConstraintSet;
use crate::error::RcError;
use crate::network::CurveNetwork;
use crate::solve::{LinearReport, Projection, ProjectionOptions, ProjectionStatus, SaddleSolver};

/// A projection and its path (the trial, then every chord iterate).
pub type ProjectionWithPath = (Projection, Vec<Vec<[f64; 3]>>);

/// Where the collision-guarded search starts (D11, §17f).
///
/// Every trial is bounded by `τ_max` and checked by the CCD whatever the rule: the rule only
/// moves the first trial, so the isotopy guarantee is the same.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum StartRule {
    /// `τ₀ = ⅔ τ_max` (`τ_hi` if there is no contact): the paper's rule.
    TwoThirdsTauMax,
    /// `τ₀ = min(⅔ τ_max, factor · τ_prev)`, with `τ_prev` the previous accepted step size of the
    /// session (the paper's rule for the first step). Needs far fewer trials per step, because the
    /// accepted `τ` is usually much smaller than `⅔ τ_max`.
    GrowPrevious {
        /// Growth of the previous accepted step size, `≥ 1` (2 in the study of §17f).
        factor: f64,
    },
}

impl StartRule {
    /// The grow-previous rule with factor 2, the measured one (§17f).
    pub const fn grow_previous() -> StartRule {
        StartRule::GrowPrevious { factor: 2.0 }
    }

    fn start(&self, tau_max: Option<f64>, tau_hi: f64, previous_tau: Option<f64>) -> f64 {
        let paper = match tau_max {
            Some(t) => 2.0 / 3.0 * t,
            None => tau_hi,
        };
        match (*self, previous_tau) {
            (StartRule::GrowPrevious { factor }, Some(prev)) if prev > 0.0 => {
                paper.min(factor * prev)
            }
            _ => paper,
        }
    }
}

/// The line search (`SPEC.md` §9.1).
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum LineSearch {
    /// CCD-guarded backtracking (D11, the default).
    CollisionGuarded {
        /// Armijo constant `c ∈ (0, 1)`.
        armijo: f64,
        /// Shrink factor `∈ (0, 1)`.
        shrink: f64,
        /// The start rule (the paper's ⅔ τ_max, or grow-previous).
        start: StartRule,
        /// The CCD gap `s ∈ (0, 1)`.
        gap: f64,
        /// Maximal number of shrinks.
        max_halvings: u32,
    },
    /// Backtracking from the reference's normalized start (D12); no collision check.
    NormalizedArmijo {
        /// Armijo constant `c ∈ (0, 1)`.
        armijo: f64,
        /// Shrink factor `∈ (0, 1)`.
        shrink: f64,
        /// Maximal number of shrinks.
        max_halvings: u32,
    },
}

impl Default for LineSearch {
    /// `CollisionGuarded { armijo 1e-4, shrink ½, start ⅔, gap 0.1, 40 halvings }` (§10).
    fn default() -> Self {
        LineSearch::CollisionGuarded {
            armijo: 1e-4,
            shrink: 0.5,
            start: StartRule::TwoThirdsTauMax,
            gap: 0.1,
            max_halvings: 40,
        }
    }
}

impl LineSearch {
    /// The default guarded search with the start rule `start`.
    pub fn guarded_with_start(start: StartRule) -> LineSearch {
        LineSearch::CollisionGuarded {
            armijo: 1e-4,
            shrink: 0.5,
            start,
            gap: 0.1,
            max_halvings: 40,
        }
    }

    /// `NormalizedArmijo { armijo 1e-4, shrink ½, 40 halvings }` (§10).
    pub fn normalized_armijo() -> LineSearch {
        LineSearch::NormalizedArmijo {
            armijo: 1e-4,
            shrink: 0.5,
            max_halvings: 40,
        }
    }

    /// `true` for [`LineSearch::CollisionGuarded`].
    pub fn is_guarded(&self) -> bool {
        matches!(self, LineSearch::CollisionGuarded { .. })
    }

    /// The CCD options of a guarded search (the default options otherwise).
    pub fn ccd_options(&self) -> CcdOptions {
        match *self {
            LineSearch::CollisionGuarded { gap, .. } => CcdOptions {
                gap,
                ..CcdOptions::default()
            },
            LineSearch::NormalizedArmijo { .. } => CcdOptions::default(),
        }
    }

    fn params(&self) -> (f64, f64, u32) {
        match *self {
            LineSearch::CollisionGuarded {
                armijo,
                shrink,
                max_halvings,
                ..
            }
            | LineSearch::NormalizedArmijo {
                armijo,
                shrink,
                max_halvings,
            } => (armijo, shrink, max_halvings),
        }
    }

    /// Checks the parameters.
    ///
    /// # Errors
    /// [`RcError::InvalidOption`] naming the parameter.
    pub fn validate(&self) -> Result<(), RcError> {
        let (armijo, shrink, _) = self.params();
        let open = |name: &'static str, v: f64| {
            if v > 0.0 && v < 1.0 {
                Ok(())
            } else {
                Err(RcError::InvalidOption {
                    name,
                    value: v.to_string(),
                    requirement: "a number in (0, 1)",
                })
            }
        };
        open("line_search.armijo", armijo)?;
        open("line_search.shrink", shrink)?;
        if let LineSearch::CollisionGuarded { start, gap, .. } = *self {
            if let StartRule::GrowPrevious { factor } = start
                && !(factor >= 1.0 && factor.is_finite())
            {
                return Err(RcError::InvalidOption {
                    name: "line_search.start.factor",
                    value: factor.to_string(),
                    requirement: "a finite number ≥ 1",
                });
            }
            open("line_search.gap", gap)?;
        }
        Ok(())
    }
}

/// Why a trial step was rejected, or that it was accepted.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum TrialOutcome {
    /// Accepted.
    Accepted,
    /// The projection failed (§6.2).
    ProjectionFailed(ProjectionStatus),
    /// A chord step of the projection has a contact (guarded search).
    PathCollision {
        /// Index of the piece of the path (0 = first chord step).
        piece: usize,
        /// The contact.
        contact: Contact,
    },
    /// The objective at the projected curve is not finite.
    NonFinite,
    /// The Armijo condition fails.
    ArmijoFailed,
}

/// One trial of a line search.
#[derive(Debug, Clone, PartialEq)]
pub struct Trial {
    /// The step size.
    pub tau: f64,
    /// The objective at the projected curve (`None` if not evaluated).
    pub energy: Option<f64>,
    /// Chord iterations of the projection.
    pub projection_iterations: u32,
    /// Largest constraint violation after the projection.
    pub max_violation: f64,
    /// The outcome.
    pub outcome: TrialOutcome,
}

/// Seconds per phase of a line search.
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct LineSearchTimings {
    /// Objective evaluations.
    pub energy: f64,
    /// Projections (linear solves and constraint evaluations).
    pub projection: f64,
    /// CCD (the move and the projection paths).
    pub ccd: f64,
    /// The whole search.
    pub total: f64,
}

/// Result of [`line_search`].
#[derive(Debug, Clone, PartialEq)]
pub struct LineSearchReport {
    /// A step was accepted.
    pub accepted: bool,
    /// The accepted step size (0 if none).
    pub tau: f64,
    /// The first trial step size `τ₀`.
    pub tau_start: f64,
    /// `τ_hi` (guarded search).
    pub tau_hi: Option<f64>,
    /// `τ_max`, the CCD's first contact time on `[0, τ_hi]` (guarded search; `∞` if none).
    pub tau_max: Option<f64>,
    /// Number of shrinks.
    pub halvings: u32,
    /// Every trial in order.
    pub trials: Vec<Trial>,
    /// The objective of the accepted curve (`F(γ)` if none).
    pub energy: f64,
    /// The accepted curve.
    pub positions: Option<Vec<[f64; 3]>>,
    /// The path of the accepted step: `γ`, `γ − τ p`, then every chord iterate (empty if none).
    pub path: Vec<Vec<[f64; 3]>>,
    /// Chord iterations of the accepted trial.
    pub projection_iterations: u32,
    /// Largest constraint violation of the accepted trial (of the last trial if none).
    pub max_violation: f64,
    /// Linear solves of all projections.
    pub projection_linear: Vec<LinearReport>,
    /// CCD counters (all queries).
    pub ccd: CcdReport,
    /// Number of CCD queries.
    pub ccd_checks: u32,
    /// Seconds per phase.
    pub seconds: LineSearchTimings,
}

/// What a line search needs from its caller: the objective, the projection, the CCD and the
/// cancel flag.
pub trait SearchEnv {
    /// The objective `F` at `positions` (may be non-finite; such trials are rejected).
    ///
    /// # Errors
    /// Evaluation failures (they abort the search).
    fn objective(&mut self, positions: &[[f64; 3]]) -> Result<f64, RcError>;
    /// Projects `trial` onto the constraints (§6.2) and returns the projection with its path
    /// (`trial`, then every chord iterate).
    ///
    /// # Errors
    /// Solver failures.
    fn project(&mut self, trial: Vec<[f64; 3]>) -> Result<ProjectionWithPath, RcError>;
    /// The first contact of `start + τ velocity`, `τ ∈ [0, tau_end]`.
    ///
    /// # Errors
    /// CCD input errors.
    fn ccd(
        &mut self,
        start: &[[f64; 3]],
        velocity: &[[f64; 3]],
        tau_end: f64,
    ) -> Result<CcdReport, RcError>;
    /// `true` once the run was cancelled (checked before every trial).
    fn cancelled(&self) -> bool;
}

/// The input of one line search.
#[derive(Debug, Clone, Copy)]
pub struct SearchInput<'a> {
    /// The base curve `γ`.
    pub base: &'a [[f64; 3]],
    /// The direction `p`; trials are `γ − τ p` (0 along a held axis).
    pub direction: &'a [[f64; 3]],
    /// `dF · p` (must be positive).
    pub slope: f64,
    /// `F(γ)`.
    pub f0: f64,
    /// The previous accepted step size (normalized search, grow-previous start).
    pub previous_tau: Option<f64>,
    /// An upper bound on `τ₀` (L-BFGS: 1 once curvature pairs exist), `None` for none.
    pub max_start: Option<f64>,
}

fn bbox_diameter(p: &[[f64; 3]]) -> f64 {
    let mut lo = [f64::INFINITY; 3];
    let mut hi = [f64::NEG_INFINITY; 3];
    for x in p {
        for k in 0..3 {
            lo[k] = lo[k].min(x[k]);
            hi[k] = hi[k].max(x[k]);
        }
    }
    let d = [hi[0] - lo[0], hi[1] - lo[1], hi[2] - lo[2]];
    (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt()
}

/// The start of the normalized search (D12, the reference's rule).
pub fn normalized_start(direction_norm: f64, previous_tau: Option<f64>) -> f64 {
    let n = direction_norm;
    let mut t = if n > 1.0 { 1.0 / n } else { 1.0 / n.sqrt() };
    if let Some(prev) = previous_tau
        && prev > 1e-5
    {
        t = (1.5 * prev).min(4.0 * t);
    }
    t
}

/// Runs the line search `ls` (§7.2 or §7.3) from `input` in `env`. Returns a report; a stuck
/// search is `accepted = false`, not an error.
///
/// [formula: rc-armijo]
///
/// # Errors
/// [`RcError::Cancelled`] if `env` is cancelled before a trial; errors of `env`.
pub fn line_search(
    ls: &LineSearch,
    input: SearchInput<'_>,
    env: &mut dyn SearchEnv,
) -> Result<LineSearchReport, RcError> {
    let t_all = Instant::now();
    let (armijo, shrink, max_halvings) = ls.params();
    let mut rep = LineSearchReport {
        accepted: false,
        tau: 0.0,
        tau_start: 0.0,
        tau_hi: None,
        tau_max: None,
        halvings: 0,
        trials: Vec::new(),
        energy: input.f0,
        positions: None,
        path: Vec::new(),
        projection_iterations: 0,
        max_violation: 0.0,
        projection_linear: Vec::new(),
        ccd: CcdReport::default(),
        ccd_checks: 0,
        seconds: LineSearchTimings::default(),
    };
    let max_speed = input
        .direction
        .iter()
        .map(|p| (p[0] * p[0] + p[1] * p[1] + p[2] * p[2]).sqrt())
        .fold(0.0f64, f64::max);
    if !(input.slope > 0.0 && input.slope.is_finite() && max_speed > 0.0 && max_speed.is_finite()) {
        rep.seconds.total = t_all.elapsed().as_secs_f64();
        return Ok(rep);
    }
    let velocity: Vec<[f64; 3]> = input
        .direction
        .iter()
        .map(|p| [-p[0], -p[1], -p[2]])
        .collect();
    let mut tau0 = match *ls {
        LineSearch::CollisionGuarded { start, .. } => {
            let tau_hi = bbox_diameter(input.base) / max_speed;
            rep.tau_hi = Some(tau_hi);
            let t = Instant::now();
            let r = env.ccd(input.base, &velocity, tau_hi)?;
            rep.seconds.ccd += t.elapsed().as_secs_f64();
            rep.ccd.accumulate(&r);
            rep.ccd_checks += 1;
            let tau_max = r.contact.map(|c| c.tau);
            rep.tau_max = Some(tau_max.unwrap_or(f64::INFINITY));
            start.start(tau_max, tau_hi, input.previous_tau)
        }
        LineSearch::NormalizedArmijo { .. } => {
            let n = input
                .direction
                .iter()
                .map(|p| p[0] * p[0] + p[1] * p[1] + p[2] * p[2])
                .sum::<f64>()
                .sqrt();
            normalized_start(n, input.previous_tau)
        }
    };
    if let Some(m) = input.max_start {
        tau0 = tau0.min(m);
    }
    rep.tau_start = tau0;
    let mut tau = tau0;
    if !(tau0 > 0.0 && tau0.is_finite()) {
        rep.seconds.total = t_all.elapsed().as_secs_f64();
        return Ok(rep);
    }
    loop {
        if env.cancelled() {
            return Err(RcError::Cancelled);
        }
        let trial_pos: Vec<[f64; 3]> = input
            .base
            .iter()
            .zip(input.direction)
            .map(|(b, p)| [b[0] - tau * p[0], b[1] - tau * p[1], b[2] - tau * p[2]])
            .collect();
        let t = Instant::now();
        let (proj, chord_path) = env.project(trial_pos)?;
        rep.seconds.projection += t.elapsed().as_secs_f64();
        rep.projection_linear.extend(proj.linear.iter().cloned());
        rep.max_violation = proj.max_violation;
        let mut trial = Trial {
            tau,
            energy: None,
            projection_iterations: proj.iterations,
            max_violation: proj.max_violation,
            outcome: TrialOutcome::ArmijoFailed,
        };
        let mut ok = proj.status == ProjectionStatus::Converged;
        if !ok {
            trial.outcome = TrialOutcome::ProjectionFailed(proj.status);
        }
        if ok && ls.is_guarded() && chord_path.len() > 1 {
            let t = Instant::now();
            for (k, w) in chord_path.windows(2).enumerate() {
                let vel: Vec<[f64; 3]> = w[0]
                    .iter()
                    .zip(&w[1])
                    .map(|(a, b)| [b[0] - a[0], b[1] - a[1], b[2] - a[2]])
                    .collect();
                let r = env.ccd(&w[0], &vel, 1.0)?;
                rep.ccd.accumulate(&r);
                rep.ccd_checks += 1;
                if let Some(contact) = r.contact {
                    trial.outcome = TrialOutcome::PathCollision { piece: k, contact };
                    ok = false;
                    break;
                }
            }
            rep.seconds.ccd += t.elapsed().as_secs_f64();
        }
        if ok {
            let t = Instant::now();
            let f = env.objective(&proj.positions)?;
            rep.seconds.energy += t.elapsed().as_secs_f64();
            trial.energy = Some(f);
            if !f.is_finite() {
                trial.outcome = TrialOutcome::NonFinite;
            } else if f <= input.f0 - armijo * tau * input.slope {
                trial.outcome = TrialOutcome::Accepted;
                rep.trials.push(trial);
                rep.accepted = true;
                rep.tau = tau;
                rep.energy = f;
                rep.projection_iterations = proj.iterations;
                let mut path = Vec::with_capacity(chord_path.len() + 1);
                path.push(input.base.to_vec());
                path.extend(chord_path);
                rep.path = path;
                rep.positions = Some(proj.positions);
                rep.seconds.total = t_all.elapsed().as_secs_f64();
                return Ok(rep);
            }
        }
        rep.trials.push(trial);
        tau *= shrink;
        rep.halvings += 1;
        if rep.halvings > max_halvings || tau < 1e-12 * tau0 {
            rep.seconds.total = t_all.elapsed().as_secs_f64();
            return Ok(rep);
        }
    }
}

/// The chord projection of §6.2 with its path: the same iteration and the same success and
/// failure rules as [`crate::solve::project_constraints`], and in addition the list of iterates
/// (`trial` first), so that the guarded search can check every chord step with the CCD.
/// Local convergence of the chord method near a regular point is cited (D10).
///
/// [formula: rc-chord-projection]
/// [formula: rc-chord-convergence]
///
/// # Errors
/// The solver's errors; [`RcError::PositionCountMismatch`].
pub fn project_with_path(
    solver: &mut dyn SaddleSolver,
    constraints: &ConstraintSet,
    network: &CurveNetwork,
    trial: Vec<[f64; 3]>,
    options: &ProjectionOptions,
    mean_edge_length0: f64,
) -> Result<ProjectionWithPath, RcError> {
    options.validate()?;
    let layout = constraints.layout();
    let mut path = vec![trial.clone()];
    let mut pos = trial;
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
        let status = if rel <= 1.0 {
            Some(
                match constraints.tangent_orientation_failure(network, &pos) {
                    Some((constraint, _)) => ProjectionStatus::TangentFlipped { constraint },
                    None => ProjectionStatus::Converged,
                },
            )
        } else if !rel.is_finite() || (it > 0 && rel > 10.0 * history[0]) {
            Some(ProjectionStatus::Diverged)
        } else if it >= options.max_iterations {
            Some(ProjectionStatus::MaxIterations)
        } else {
            None
        };
        if let Some(status) = status {
            return Ok((
                Projection {
                    positions: pos,
                    status,
                    iterations: it,
                    max_violation,
                    history,
                    linear,
                },
                path,
            ));
        }
        let b: Vec<f64> = phi.iter().map(|p| -p).collect();
        let sol = solver.solve(&zero, &b)?;
        layout.add_scaled(&mut pos, 1.0, &sol.x);
        path.push(pos.clone());
        linear.push(sol.report);
        it += 1;
    }
}
