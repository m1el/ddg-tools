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
//! otherwise `τ ← shrink · τ`. With a second-order correction `q` ([`SearchEnv::second_order_correction`],
//! `SPEC.md` §17l) the trial is `γ − τ p + τ² q` before the projection: `C q = −½ Φ''(γ)[p, p]`
//! cancels the `O(τ²)` violation of the straight trial (the leading term of the projection's
//! own correction), so the projection starts at an `O(τ³)` violation; the piece from `γ − τ p`
//! to the predicted trial is checked by the CCD like a chord step.
//! After `max_halvings` shrinks, or below `1e-12 τ₀`, the search is
//! stuck. For a descent direction (`dF·p > 0`) and a projection that succeeds for small `τ`, the
//! search terminates (cited, `rc-armijo`).
//!
//! * [`LineSearch::CollisionGuarded`] (the default, D11): `τ_hi = diam(bbox γ) / max_v |p_v|`,
//!   `τ_max` = the first contact time of the linear motion `γ − τ p` on `[0, τ_hi]` by the CCD
//!   of §7.4 (∞ if none; 0 if two disjoint edges are already closer than [`CONTACT_FLOOR`]
//!   times `ℓ̄₀`, §17l), start `τ₀` by the [`StartRule`] (the paper's `⅔ τ_max`, or `τ_hi` without contact; or
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
use crate::solve::{
    LinearReport, NEWTON_POLISH_MAX_STEPS, Projection, ProjectionOptions, ProjectionStatus,
    SaddleSolver, newton_correction,
};

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
    /// CCD-guarded backtracking (D11, the default with the grow-previous start).
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
    /// `CollisionGuarded { armijo 1e-4, shrink ½, start GrowPrevious 2, gap 0.1, 40 halvings }`
    /// (§10, §17f).
    fn default() -> Self {
        LineSearch::guarded_with_start(StartRule::grow_previous())
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
    /// The path of the accepted step: `γ`, `γ − τ p`, the predicted trial `γ − τ p + τ² q` (with a
    /// second-order correction), then every chord iterate (empty if none).
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
    /// The guarded search found two disjoint edges closer than the contact floor at the base
    /// curve ([`CONTACT_FLOOR`]) and took no step (`τ_max = 0`).
    pub contact_floor: bool,
}

/// The contact floor of the collision-guarded search (`SPEC.md` §17l), relative to the mean
/// edge length `ℓ̄₀` at the start of the session: two disjoint edges of the base curve closer
/// than `CONTACT_FLOOR · ℓ̄₀` count as a contact, so `τ_max = 0` and the step is
/// [`StepStatus::Stuck`](super::StepStatus::Stuck). Without it the CCD bounds each step by a
/// fraction (the gap) of the current distance, which a metric blind to the repulsion (L², H¹,
/// H²) lets shrink geometrically down to rounding, until the curve is no longer valid. The floor
/// is 10⁶ times the validation tolerance of §2.3 (`INTERSECTION_TOLERANCE = 1e-12`), so the
/// last accepted curve keeps a distance of at least `gap · CONTACT_FLOOR · ℓ̄₀`; the
/// tangent-point energy of such a pair is enormous (`ℓ² d^{−(β−α)}`, about 10¹⁸ times its value
/// at `d = ℓ` for (3, 6)), so a Sobolev flow never gets there and its steps are unchanged.
pub const CONTACT_FLOOR: f64 = 1e-6;

/// The pair `(I, J)`, `I < J`, of disjoint curve edges at the smallest distance `≤ tol`, if any
/// (a uniform grid of cell size `max ℓ + 2 tol` as broad phase, exact segment distances;
/// deterministic: candidates are sorted).
pub fn closest_pair_within(
    network: &CurveNetwork,
    positions: &[[f64; 3]],
    tol: f64,
) -> Option<(u32, u32, f64)> {
    let edges = network.edges();
    if edges.len() < 2 || tol.is_nan() || tol < 0.0 {
        return None;
    }
    let len = |&[a, b]: &[u32; 2]| {
        let (p, q) = (positions[a as usize], positions[b as usize]);
        ((p[0] - q[0]).powi(2) + (p[1] - q[1]).powi(2) + (p[2] - q[2]).powi(2)).sqrt()
    };
    let h = edges.iter().map(len).fold(0.0f64, f64::max) + 2.0 * tol;
    if !(h > 0.0 && h.is_finite()) {
        return None;
    }
    let cell = |x: f64| (x / h).floor() as i64;
    let mut grid: std::collections::HashMap<[i64; 3], Vec<u32>> = std::collections::HashMap::new();
    for (e, &[a, b]) in edges.iter().enumerate() {
        let (p, q) = (positions[a as usize], positions[b as usize]);
        let lo = [0, 1, 2].map(|c| cell(p[c].min(q[c]) - tol));
        let hi = [0, 1, 2].map(|c| cell(p[c].max(q[c]) + tol));
        for x in lo[0]..=hi[0] {
            for y in lo[1]..=hi[1] {
                for z in lo[2]..=hi[2] {
                    grid.entry([x, y, z]).or_default().push(e as u32);
                }
            }
        }
    }
    let mut pairs: Vec<(u32, u32)> = Vec::new();
    for list in grid.values() {
        for (k, &i) in list.iter().enumerate() {
            for &j in &list[k + 1..] {
                let (i, j) = (i.min(j), i.max(j));
                if !network.edges_adjacent(i as usize, j as usize) {
                    pairs.push((i, j));
                }
            }
        }
    }
    pairs.sort_unstable();
    pairs.dedup();
    let mut best: Option<(u32, u32, f64)> = None;
    for (i, j) in pairs {
        let [a, b] = edges[i as usize];
        let [c, d] = edges[j as usize];
        let dist = crate::ccd::segment_segment_distance(
            positions[a as usize],
            positions[b as usize],
            positions[c as usize],
            positions[d as usize],
        );
        if dist <= tol && best.is_none_or(|(_, _, bd)| dist < bd) {
            best = Some((i, j, dist));
        }
    }
    best
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
    /// `true` if two disjoint edges of `positions` are closer than the contact floor
    /// ([`CONTACT_FLOOR`]); the guarded search then takes no step. The default: never.
    ///
    /// # Errors
    /// Evaluation failures.
    fn below_contact_floor(&mut self, positions: &[[f64; 3]]) -> Result<bool, RcError> {
        let _ = positions;
        Ok(false)
    }
    /// The second-order correction `q` of the trials (`SPEC.md` §17l): trials are
    /// `γ − τ p + τ² q` (then projected); `None` (the default) for the straight trials `γ − τ p`.
    fn second_order_correction(&self) -> Option<&[[f64; 3]]> {
        None
    }
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
        contact_floor: false,
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
            let floor = env.below_contact_floor(input.base)?;
            rep.seconds.ccd += t.elapsed().as_secs_f64();
            if floor {
                rep.contact_floor = true;
                rep.tau_max = Some(0.0);
                rep.seconds.total = t_all.elapsed().as_secs_f64();
                return Ok(rep);
            }
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
        let correction = env.second_order_correction().map(<[[f64; 3]]>::to_vec);
        let (proj, chord_path) = match correction.as_deref() {
            // The predicted trial γ − τp + τ²q; the piece from γ − τp to it is part of the path
            // (checked by the CCD like a chord step), so the path stays piecewise linear.
            Some(q) => {
                let t2 = tau * tau;
                let predicted: Vec<[f64; 3]> = trial_pos
                    .iter()
                    .zip(q)
                    .map(|(x, q)| [x[0] + t2 * q[0], x[1] + t2 * q[1], x[2] + t2 * q[2]])
                    .collect();
                let (proj, chord) = env.project(predicted)?;
                let mut path = Vec::with_capacity(chord.len() + 1);
                path.push(trial_pos);
                path.extend(chord);
                (proj, path)
            }
            None => env.project(trial_pos)?,
        };
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

/// The retry of a stuck search with the projected base as the Armijo reference (WP16b,
/// `SPEC.md` §17p; [`FlowOptions::stuck_projected_reference`](super::FlowOptions::stuck_projected_reference)).
///
/// With back-projection every trial is projected with at least one step, so as `τ → 0` the
/// trials tend to `P(γ)`, the projection of the base, not to `γ`: the base satisfies the
/// constraints only within the tolerance, and `F(P(γ)) − F(γ) ≈ −λ·Φ(γ)` (the multipliers of
/// Eq. 22 times the residual violation). When that exceeds the decrease the Armijo test can
/// find, no trial passes and the search is stuck although `γ` is far from a minimum. This
/// projects the base (`τ = 0`), and if that succeeds and `F(P(γ)) > F(γ)`, runs the same search
/// again with `F(P(γ))` as the reference: `F(γ̂_τ) ≤ F(P(γ)) − c τ dF·p`. An accepted step can
/// raise `F` by at most `F(P(γ)) − F(γ)`. Returns `None` when the projection fails or does not
/// raise the objective (then the first search's verdict stands), else the second search with
/// the first one's trials, counters and times prepended.
///
/// # Errors
/// As [`line_search`].
pub fn retry_with_projected_reference(
    ls: &LineSearch,
    input: SearchInput<'_>,
    env: &mut dyn SearchEnv,
    first: &LineSearchReport,
) -> Result<Option<LineSearchReport>, RcError> {
    let t = Instant::now();
    let (proj, _) = env.project(input.base.to_vec())?;
    let t_proj = t.elapsed().as_secs_f64();
    if proj.status != ProjectionStatus::Converged {
        return Ok(None);
    }
    let t = Instant::now();
    let f_ref = env.objective(&proj.positions)?;
    let t_energy = t.elapsed().as_secs_f64();
    if !(f_ref.is_finite() && f_ref > input.f0) {
        return Ok(None);
    }
    let mut second = line_search(ls, SearchInput { f0: f_ref, ..input }, env)?;
    second.projection_linear.splice(0..0, proj.linear);
    second.seconds.energy += t_energy;
    second.seconds.projection += t_proj;
    second.seconds.total += t_proj + t_energy;
    if !second.accepted {
        second.energy = input.f0;
    }
    Ok(Some(merge_after(first, second)))
}

/// `second`, a repeated search of the same step, with the trials, counters and times of
/// `first` prepended (the retries of WP16b, `SPEC.md` §17p).
pub fn merge_after(first: &LineSearchReport, mut second: LineSearchReport) -> LineSearchReport {
    let mut trials = first.trials.clone();
    trials.append(&mut second.trials);
    second.trials = trials;
    second.halvings += first.halvings;
    let mut linear = first.projection_linear.clone();
    linear.append(&mut second.projection_linear);
    second.projection_linear = linear;
    second.ccd.accumulate(&first.ccd);
    second.ccd_checks += first.ccd_checks;
    second.seconds.energy += first.seconds.energy;
    second.seconds.projection += first.seconds.projection;
    second.seconds.ccd += first.seconds.ccd;
    second.seconds.total += first.seconds.total;
    second
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
    project_with_path_polished(
        solver,
        constraints,
        network,
        trial,
        options,
        mean_edge_length0,
        ProjectionRules::default(),
    )
}

/// The departures from the paper's chord projection (`SPEC.md` §17l) used by
/// [`project_with_path_polished`]; the default (`None`, `false`) is the paper's projection.
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct ProjectionRules {
    /// The Newton polish threshold ([`FlowOptions::newton_polish`](super::FlowOptions::newton_polish)).
    pub newton_polish: Option<f64>,
    /// The rate test ([`FlowOptions::projection_early_abort`](super::FlowOptions::projection_early_abort)).
    pub early_abort: bool,
}

/// [`project_with_path`] with the Newton polish of `SPEC.md` §17l: once the largest relative
/// violation `max_r |Φ_r| / tol_r` is at most `rules.newton_polish` (`Some`), the remaining iterations are
/// Newton steps with the Jacobian at the iterate in the lumped-mass metric
/// ([`crate::solve::newton_correction`], no solve with `Ā`), at most
/// [`NEWTON_POLISH_MAX_STEPS`] of them; chord iterations (with `solver`) still count against
/// `options.max_iterations`. The success and failure rules are those of §6.2. `None`: the chord
/// method throughout.
///
/// With `rules.early_abort`, a chord iteration that cannot reach its target within the remaining
/// chord iterations is stopped early (`MaxIterations`, fewer iterations): after a chord step from
/// `h_{k−1}` to `h_k` (largest relative violations) the observed rate `θ = h_k / h_{k−1}` predicts
/// `h_k θ^m` after `m` more chord steps (the chord method converges linearly, at a rate
/// proportional to the step); the target is the polish threshold (or 1 without polish). The trial
/// is given up when `θ ≥ 1` or the predicted number of steps exceeds the remaining budget. This is
/// the convergence-rate test of simplified-Newton iterations in implicit integrators (Hairer &
/// Wanner, *Solving ODEs II*, 2nd ed., §IV.8, stopping rule (8.12)), `SPEC.md` §17l.
///
/// [formula: rc-chord-projection]
///
/// # Errors
/// As [`project_with_path`]; the errors of [`crate::solve::newton_correction`].
pub fn project_with_path_polished(
    solver: &mut dyn SaddleSolver,
    constraints: &ConstraintSet,
    network: &CurveNetwork,
    trial: Vec<[f64; 3]>,
    options: &ProjectionOptions,
    mean_edge_length0: f64,
    rules: ProjectionRules,
) -> Result<ProjectionWithPath, RcError> {
    options.validate()?;
    let mut polish = rules.newton_polish;
    let early_abort = rules.early_abort;
    let layout = constraints.layout();
    let mut path = vec![trial.clone()];
    let mut pos = trial;
    let mut history = Vec::new();
    let mut linear = Vec::new();
    let zero = vec![0.0; solver.n_dofs()];
    let mut it = 0u32;
    let (mut chord, mut newton) = (0u32, 0u32);
    let mut last_chord = false;
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
        let use_newton = polish.is_some_and(|p| rel <= p);
        // The rate test (Hairer & Wanner IV.8) after a chord step.
        let hopeless = early_abort && last_chord && !use_newton && rel > 1.0 && it > 0 && {
            let target = polish.unwrap_or(1.0).max(1.0);
            let theta = rel / history[history.len() - 2];
            let remaining = options.max_iterations.saturating_sub(chord);
            theta.is_nan()
                || theta >= 1.0
                || (rel / target).ln() / (1.0 / theta).ln() > f64::from(remaining)
        };
        // At least one step, also for a trial already within the tolerance (§17l, as
        // `project_constraints`).
        let status = if rel <= 1.0 && (it > 0 || rel == 0.0 || options.max_iterations == 0) {
            Some(
                match constraints.tangent_orientation_failure(network, &pos) {
                    Some((constraint, _)) => ProjectionStatus::TangentFlipped { constraint },
                    None => ProjectionStatus::Converged,
                },
            )
        } else if !rel.is_finite() || (it > 0 && rel > 10.0 * history[0]) {
            Some(ProjectionStatus::Diverged)
        } else if (use_newton && newton >= NEWTON_POLISH_MAX_STEPS)
            || (!use_newton && chord >= options.max_iterations)
            || hopeless
        {
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
        // A Jacobian that is rank deficient at the iterate (e.g. a collapsing edge) has no
        // Newton step: the iteration continues with chord steps, which use the base curve's.
        let newton_step = if use_newton {
            match newton_correction(constraints, network, &pos, &phi) {
                Ok(x) => Some(x),
                Err(RcError::RedundantConstraints { .. } | RcError::SingularSaddle { .. }) => {
                    polish = None;
                    None
                }
                Err(e) => return Err(e),
            }
        } else {
            None
        };
        last_chord = newton_step.is_none();
        if let Some(x) = newton_step {
            layout.add_scaled(&mut pos, 1.0, &x);
            newton += 1;
        } else if use_newton && chord >= options.max_iterations {
            return Ok((
                Projection {
                    positions: pos,
                    status: ProjectionStatus::MaxIterations,
                    iterations: it,
                    max_violation,
                    history,
                    linear,
                },
                path,
            ));
        } else {
            let b: Vec<f64> = phi.iter().map(|p| -p).collect();
            let sol = solver.solve(&zero, &b)?;
            layout.add_scaled(&mut pos, 1.0, &sol.x);
            linear.push(sol.report);
            chord += 1;
        }
        path.push(pos.clone());
        it += 1;
    }
}
