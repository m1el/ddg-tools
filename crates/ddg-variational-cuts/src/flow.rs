//! The flow (SPEC §3.8–§3.12, §4): one step moves the level-set function by the harmonic
//! extension of the speed on the cut, applies the implicit length step, and makes the function
//! a signed distance again.

use std::collections::HashMap;
use std::time::Instant;

use ddg_core::distance::upwind_distance;
use ddg_core::sparse::{Method, Parallelism, Symbolic};
use ddg_core::{Geometry, SparseMatrix, TriMesh};

use crate::cut::{Cut, region_of};
use crate::cut_mesh::CutTriangles;
use crate::error::CutError;
use crate::evaluate::{EnergyTerms, Evaluation, Evaluator};
use crate::problem::{CutProblem, Energy, EvalOptions, UnitAreaScale, check_length};
use crate::system::LaplaceSystem;

/// How the time step `τ` of a flow step is chosen (SPEC §4).
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum StepSize {
    /// The same `τ > 0` in every step, as in the paper (Sec. 6).
    Fixed(f64),
    /// `τ = min(tau_max, max_motion · h / ρ̂)` with `h` the shortest crossed edge and `ρ̂` the
    /// speed of the cut that is known before the step: the cut moves by at most the fraction
    /// `max_motion` of that edge per step (the paper's remark in Sec. 7.4). With
    /// [`LengthStep::OnCut`] `ρ = s + λ⁰`, distortion and length together; with
    /// [`LengthStep::OnMesh`] the length term never passes through the cut points and
    /// `ρ = s`. Which `h` and which statistic of `|ρ_c|` over the cut are taken is
    /// [`FlowOptions::step_statistics`] (SPEC D22). If the speed vanishes, `τ = tau_max`.
    Capped {
        /// The largest step, `> 0`.
        tau_max: f64,
        /// The largest motion of the cut in a step, as a fraction of the shortest crossed edge,
        /// `> 0`.
        max_motion: f64,
    },
    /// `τ = min(tau_max, time_scale / ω, max_motion · h / ρ̂)`: [`StepSize::Capped`] with a
    /// largest step that follows the state. `ω` is the roughening rate of the state (the same
    /// quantile of the local rates), which has the unit of a speed: a cut of the size of
    /// the surface (length 1 on the unit-area surface) reacts to a displacement at rates of
    /// the order `ω` (on the sphere of unit area, without a length weight, the equator returns
    /// at the rate `ω/R = 3.5 ω` and a ripple with `m` periods grows at the rate
    /// `3.5 (m − 1) ω`), and the distortion term is explicit. Unlike a fixed `tau_max` this does not depend on the
    /// energy or on the scale of its weights. The default (SPEC §4, D9f).
    Adaptive {
        /// The step in units of `1/ω`, `> 0`.
        time_scale: f64,
        /// The largest motion of the cut in a step, as a fraction of the shortest crossed edge,
        /// `> 0`.
        max_motion: f64,
        /// The largest step, `> 0`; it acts where `ω` and the speed vanish.
        tau_max: f64,
    },
}

impl Default for StepSize {
    /// `Adaptive { time_scale: 0.05, max_motion: 0.25, tau_max: 1 }` (SPEC §4, D9f).
    fn default() -> Self {
        StepSize::Adaptive {
            time_scale: 0.05,
            max_motion: 0.25,
            tau_max: 1.0,
        }
    }
}

impl StepSize {
    /// The step for the speed `ρ̂ ≥ 0` of the cut, the edge length `h` and the roughening rate
    /// `ω`, each of them the statistic over the cut that [`StepStatistics`] names (for the
    /// robust ones the speed is in crossed edges per unit time and `h = 1`). A speed or a rate
    /// that is not positive does not limit the step.
    ///
    /// [formula: cut-step-cap]
    pub fn step(&self, speed: f64, shortest_crossed_edge: f64, roughening_rate: f64) -> f64 {
        let by_motion = |max_motion: f64| {
            if speed > 0.0 {
                max_motion * shortest_crossed_edge / speed
            } else {
                f64::INFINITY
            }
        };
        match *self {
            StepSize::Fixed(tau) => tau,
            StepSize::Capped {
                tau_max,
                max_motion,
            } => tau_max.min(by_motion(max_motion)),
            StepSize::Adaptive {
                time_scale,
                max_motion,
                tau_max,
            } => {
                let by_rate = if roughening_rate > 0.0 {
                    time_scale / roughening_rate
                } else {
                    f64::INFINITY
                };
                tau_max.min(by_rate).min(by_motion(max_motion))
            }
        }
    }

    /// The local limit of SPEC §4 (D22, [`StepStatistics::Robust`]): the factor `≤ 1` by which
    /// the rate `rate` of a cut point on a crossed edge of length `shortest_crossed_edge` is
    /// reduced so that the point moves by at most `max_motion` of that edge in the step `tau`,
    /// where its speed is taken to be at least half of `excess`, the amount by which its own
    /// roughening rate exceeds that of the state ([`StepSize::Adaptive`] only; a cut reaches
    /// the speed `ω_c/2` within the face of a spike of `ω_c`). 1 for [`StepSize::Fixed`].
    ///
    /// `speed_statistic` is the `m̂` that the step rule used. A cut point whose speed in edges
    /// per unit time does not exceed it is not limited. In exact arithmetic that follows from
    /// `τ ≤ max_motion / m̂`; it is tested first because the cut point of the quantile lies on
    /// the limit exactly, and the rounding of `τ` would otherwise decide whether it counts as
    /// limited.
    ///
    /// [formula: cut-step-local]
    pub fn local_factor(
        &self,
        tau: f64,
        speed_statistic: f64,
        shortest_crossed_edge: f64,
        rate: f64,
        excess: f64,
    ) -> f64 {
        let (max_motion, speed) = match *self {
            StepSize::Fixed(_) => return 1.0,
            StepSize::Capped { max_motion, .. } => (max_motion, rate.abs()),
            StepSize::Adaptive { max_motion, .. } => {
                (max_motion, rate.abs().max(0.5 * excess.max(0.0)))
            }
        };
        if speed / shortest_crossed_edge <= speed_statistic {
            return 1.0;
        }
        let moved = tau * speed;
        let allowed = max_motion * shortest_crossed_edge;
        if moved > allowed {
            allowed / moved
        } else {
            1.0
        }
    }
}

/// The quantile of the local roughening rates that is the roughening rate `ω` of a state
/// (SPEC §1, D22).
pub const ROUGHENING_QUANTILE: f64 = 0.9;

/// The `q`-quantile of `values` weighted by `weights ≥ 0` (SPEC §1): with the values in
/// ascending order (equal values in the order of their index), the first value at which the
/// weights summed so far reach `q` times the sum of all weights (both sums taken in that
/// order). `q = 1` gives the largest value, an empty list 0.
///
/// [formula: cut-weighted-quantile]
pub fn weighted_quantile(values: &[f64], weights: &[f64], q: f64) -> f64 {
    let mut order: Vec<usize> = (0..values.len()).collect();
    order.sort_by(|a, b| values[*a].total_cmp(&values[*b]));
    let total: f64 = order.iter().map(|&i| weights[i]).sum();
    let mut seen = 0.0;
    for &i in &order {
        seen += weights[i];
        if seen >= q * total {
            return values[i];
        }
    }
    order.last().map_or(0.0, |&i| values[i])
}

/// The statistics of the cut that [`StepSize::Capped`] and [`StepSize::Adaptive`] use (SPEC §4,
/// D22). They do not matter for [`StepSize::Fixed`].
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum StepStatistics {
    /// Quantiles, and a local limit. The speed of a cut point is measured in units of the
    /// shortest crossed edge `h_c` at the two ends of its own edge, `m_c = |ρ_c| / h_c`; the step is
    /// `τ = min(tau_max, time_scale / ω, max_motion / m̂)` with `m̂` and `ω` the
    /// [`ROUGHENING_QUANTILE`] quantiles, weighted by the length of the cut, of `m_c` and of
    /// the local roughening rates `ω_c`. A cut point beyond them takes a smaller step of its
    /// own ([`StepSize::local_factor`]): no cut point moves by more than `max_motion · h_c`.
    /// One spike of curvature or one short edge then does not set the step of the whole cut.
    /// With the length step on the mesh the step is split, and `m̂` is the larger of the
    /// quantiles for `ρ = s` (the explicit part) and for `ρ = s + λ⁰` (the whole motion; the
    /// length force is computed on the cut for this purpose only). The default.
    #[default]
    Robust,
    /// The largest speed and roughening rate over the cut and the shortest crossed edge,
    /// `τ = min(tau_max, time_scale / max ω_c, max_motion · min h / max |ρ_c|)`, and no local
    /// limit: the rule of the paper's remark (Sec. 7.4) taken literally, and the rule of this
    /// crate before D22.
    Extremes,
}

/// The edge length `h_c` of every cut point for the step rule (SPEC §4): the shortest of the
/// crossed edges that share a vertex with its own edge (its own edge included). A change of
/// the level-set function at a vertex moves the crossings of all its edges, so this is the
/// length that the motion of the cut point has to be compared with.
///
/// [formula: cut-step-cap]
pub(crate) fn local_edges(mesh: &TriMesh, geometry: &Geometry, cut: &Cut) -> Vec<f64> {
    let lengths = geometry.edge_lengths();
    let mut shortest_at = vec![f64::INFINITY; mesh.n_vertices()];
    for point in &cut.points {
        let length = lengths[point.edge as usize];
        for v in mesh.edge_vertices(point.edge as usize) {
            shortest_at[v] = shortest_at[v].min(length);
        }
    }
    cut.points
        .iter()
        .map(|point| {
            let [a, b] = mesh.edge_vertices(point.edge as usize);
            shortest_at[a].min(shortest_at[b])
        })
        .collect()
}

/// For the local limit at the vertices (SPEC §4, D26): the edge length `h_i` of every vertex,
/// the shortest crossed edge at a cut-adjacent vertex and the shortest edge at any other one,
/// and the largest local roughening rate `ω_c` of the cut points on the edges at the vertex
/// (`−∞` where there is none).
///
/// [formula: cut-step-cap]
pub(crate) fn vertex_limits(
    mesh: &TriMesh,
    geometry: &Geometry,
    cut: &Cut,
    local_rates: &[f64],
) -> (Vec<f64>, Vec<f64>) {
    let lengths = geometry.edge_lengths();
    let nv = mesh.n_vertices();
    let mut crossed = vec![f64::INFINITY; nv];
    let mut rate = vec![f64::NEG_INFINITY; nv];
    for (point, omega) in cut.points.iter().zip(local_rates) {
        let length = lengths[point.edge as usize];
        for v in mesh.edge_vertices(point.edge as usize) {
            crossed[v] = crossed[v].min(length);
            rate[v] = rate[v].max(*omega);
        }
    }
    let mut any = vec![f64::INFINITY; nv];
    for (edge, length) in mesh.edges().iter().zip(lengths) {
        for v in edge {
            any[*v as usize] = any[*v as usize].min(*length);
        }
    }
    let edges = (0..nv)
        .map(|v| {
            if crossed[v].is_finite() {
                crossed[v]
            } else {
                any[v]
            }
        })
        .collect();
    (edges, rate)
}

/// The step of SPEC §4 for the speeds `speed` at the cut points (or boundary nodes) with the
/// lengths `weights` of the cut and the crossed edges `edges`, and the roughening rate
/// `roughening_rate` of the same statistics. Returns `τ` and the speed statistic that was
/// used: `m̂` (in crossed edges per unit time) or the largest `|ρ|`.
///
/// [formula: cut-step-cap]
pub(crate) fn rule_step(
    step_size: &StepSize,
    statistics: StepStatistics,
    speed: &[f64],
    edges: &[f64],
    weights: &[f64],
    roughening_rate: f64,
) -> (f64, f64) {
    match statistics {
        StepStatistics::Robust => {
            let in_edges: Vec<f64> = speed.iter().zip(edges).map(|(r, h)| r.abs() / h).collect();
            let fastest = weighted_quantile(&in_edges, weights, ROUGHENING_QUANTILE);
            (step_size.step(fastest, 1.0, roughening_rate), fastest)
        }
        StepStatistics::Extremes => {
            let fastest = speed.iter().fold(0.0f64, |m, r| m.max(r.abs()));
            let shortest = edges.iter().fold(f64::INFINITY, |m, h| m.min(*h));
            (step_size.step(fastest, shortest, roughening_rate), fastest)
        }
    }
}

/// How the penalty `η` enters the implicit step (SPEC §3.10).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum PenaltyRule {
    /// `(W_V + τ L_w) φ′ = W_V φ̃` with `w = α_L + η`: backward Euler for `φ_t = −Δ_w φ`, the
    /// full shape derivative `∂η/∂n + (α_L + η) κ` of the length terms (SPEC D8). The default.
    #[default]
    Gradient,
    /// `(W_V D⁻¹ + τ L) φ′ = W_V D⁻¹ φ̃` with `D = diag(α_L + η_i)`: the rule of the paper
    /// (Sec. 7.1), which keeps `(α_L + η) κ` and drops `∂η/∂n`. Needs `α_L + η_i > 0` at every
    /// vertex (or zero at every vertex, which skips the step).
    PaperDiagonal,
}

/// Where the length and penalty terms of the energy are integrated (SPEC §3.10, D14).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum LengthStep {
    /// On the cut, where the energy `Σ ½ (w_p + w_q) |pq|` is defined. The first variation of
    /// that sum with respect to the position of every cut point, per unit of area swept, is the
    /// discrete `∂w/∂n + w κ`; it is made implicit along the cut by a 1D system and added to the
    /// speed `s` before the harmonic extension. The default.
    #[default]
    OnCut,
    /// On the mesh, as in the paper (Fig. 9 step 6): `(W_V + τ L_w) φ′ = W_V φ̃`. The discrete
    /// Laplacian differentiates the redistanced function twice across the boundary between the
    /// cut-adjacent vertices, which keep their values, and their neighbours, which are
    /// recomputed to first order; the curvature speed it yields is wrong by a factor that does
    /// not go to 1 with the mesh width unless `τ` is much larger than its square.
    OnMesh,
}

/// How the speed on the cut becomes a rate of change of the level-set function (SPEC §3.9).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum SpeedScaling {
    /// `φ_t = s` on the cut, as in the paper (Eq. 21): exact while `|∇φ| = 1` there.
    Unit,
    /// `φ_t = s · ḡ` on the cut, with `ḡ` the mean of `|∇φ|` over the connected piece of the
    /// cut: the level-set equation for a function whose slope across the cut is not 1. The
    /// default (SPEC D16): vertices that became cut-adjacent in the last step carry first-order
    /// distances, so `|∇φ|` in crossed faces is a few percent below 1 and the cut would run
    /// that much too fast. The mean, not the local value, is used: multiplying by the local
    /// `|∇φ|` feeds a deviation of one vertex back into its own rate and is unstable for large
    /// steps.
    #[default]
    MeanGradient,
}

/// How the rate on the cut is extended to the vertices of the mesh (SPEC §3.8, D15).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum Extension {
    /// From the cut outward. A cut-adjacent vertex takes the mean of the rates at the cut points
    /// on its own crossed edges, each weighted by its dual length on the cut times the
    /// barycentric weight of the vertex at the cut point; every other vertex takes the mean,
    /// weighted by inverse edge length, of its neighbours in the same region that are one ring
    /// closer to the cut. No linear system. The values lie between the smallest and the largest
    /// rate, and the weights at the cut are the transpose of the interpolation from the
    /// vertices back to the cut points, so the mean displacement of the cut is the mean of the
    /// rates. The default.
    #[default]
    Local,
    /// The harmonic extension with the cotan Laplacian of the cut mesh, as in the paper
    /// (Sec. 5.3.1): `L_OO σ_O = −L_OC s_C`. On a closed mesh the matrix is the one the
    /// evaluation factorized. The triangles that cutting produces are arbitrarily thin, and
    /// through them the value at a cut-adjacent vertex is governed by cotangent weights that
    /// have nothing to do with which cut points are near it (measured: a vertex whose own
    /// crossings have the rates `+4`, `+16`, `+6` receives `−4`); the two ends of a crossed
    /// edge are then pushed apart and the level-set function stops being a distance.
    Harmonic,
}

/// What the redistancing does with the cut-adjacent vertices (SPEC §3.12).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum Redistance {
    /// The cut-adjacent vertices keep their values, every other vertex gets its distance from
    /// them by fast marching (paper Sec. 6.4). The cut does not move. The default.
    #[default]
    KeepAdjacent,
    /// No redistancing: the level-set function is left as the length step produced it.
    Off,
}

/// Options of a flow step: the options of the evaluation, and the controls of SPEC §4, §5 that
/// are not part of the energy.
#[derive(Debug, Clone)]
pub struct FlowOptions<'a> {
    /// Options of the evaluation of the cut (crossing offset, flux recovery, parallelism).
    pub eval: EvalOptions<'a>,
    /// The time step.
    pub step_size: StepSize,
    /// Where the length and penalty terms are integrated. Default: on the cut.
    pub length_step: LengthStep,
    /// How the penalty `η` enters the length step. Default: the gradient rule.
    pub penalty_rule: PenaltyRule,
    /// Symmetry classes (SPEC §3.11): one label per vertex; vertices with equal labels are kept
    /// at equal values of the level-set function. `None`: no symmetry.
    pub symmetry: Option<Vec<u32>>,
    /// How the speed becomes a rate of change of the level-set function.
    pub speed_scaling: SpeedScaling,
    /// How the rate on the cut is extended to the vertices. Default: by distance.
    pub extension: Extension,
    /// The redistancing.
    pub redistance: Redistance,
    /// Which statistics of the cut the step rule uses (SPEC §4, D22). Default: robust ones.
    pub step_statistics: StepStatistics,
}

impl FlowOptions<'_> {
    /// Default options with the given step size.
    pub fn new(step_size: StepSize) -> Self {
        FlowOptions {
            eval: EvalOptions::default(),
            step_size,
            length_step: LengthStep::default(),
            penalty_rule: PenaltyRule::default(),
            symmetry: None,
            speed_scaling: SpeedScaling::default(),
            extension: Extension::default(),
            redistance: Redistance::default(),
            step_statistics: StepStatistics::default(),
        }
    }

    /// The stages of the paper as printed, with the given step size: nodal normal derivative
    /// with one pass of sliver smoothing ([`EvalOptions::paper`]), harmonic extension, `φ_t = s`,
    /// the length step on the mesh. The penalty rule stays the gradient rule; set
    /// [`PenaltyRule::PaperDiagonal`] for the paper's.
    pub fn paper(step_size: StepSize) -> Self {
        FlowOptions {
            eval: EvalOptions::paper(),
            step_size,
            length_step: LengthStep::OnMesh,
            penalty_rule: PenaltyRule::default(),
            symmetry: None,
            speed_scaling: SpeedScaling::Unit,
            extension: Extension::Harmonic,
            redistance: Redistance::KeepAdjacent,
            step_statistics: StepStatistics::default(),
        }
    }

    /// The stages for more than two regions (SPEC §8.4, D20): those of the paper (length step
    /// on the mesh, harmonic extension per region, `φ_t = s`) with the flux recovery of SPEC
    /// §3.5, which is what [`FlowOptions::paper`] does not have.
    ///
    /// [`FlowOptions::new`] works on networks too: every region then moves its own boundary by
    /// the length force of SPEC §3.10 (§8.5). It keeps symmetric partitions of the sphere in
    /// place six to eight times better than the stages of the paper and follows single curves
    /// exactly, but at triple points it is no better than the paper's stages (junction angles
    /// settle 4° to 9° off 120° at the mesh widths tested, against 3° to 5°).
    ///
    /// The paper's stages were chosen here when, from a speckled start with a strong length
    /// weight, the default stages came to rest with two to three times the patches and the
    /// energy. Most of that was the step rule of the time: with [`StepStatistics::Robust`]
    /// (SPEC D22) the default stages end at the lower energy on the three models tested with
    /// the automatic length weight, and on one of three with the strong one (README, "Several
    /// regions"; SPEC §8.5). Which set to prefer for `N > 2` is open; this function is as it
    /// was.
    pub fn for_regions(step_size: StepSize) -> Self {
        FlowOptions {
            eval: EvalOptions::default(),
            ..FlowOptions::paper(step_size)
        }
    }

    pub(crate) fn validate(&self, n_vertices: usize) -> Result<(), CutError> {
        let positive = |x: f64| x.is_finite() && x > 0.0;
        match self.step_size {
            StepSize::Fixed(tau) if !positive(tau) => Err(CutError::InvalidOptions {
                reason: format!("the time step must be positive and finite, got {tau}"),
            }),
            StepSize::Capped {
                tau_max,
                max_motion,
            } if !(positive(tau_max) && positive(max_motion)) => Err(CutError::InvalidOptions {
                reason: format!(
                    "tau_max and max_motion must be positive and finite, got {tau_max} and \
                     {max_motion}"
                ),
            }),
            StepSize::Adaptive {
                time_scale,
                max_motion,
                tau_max,
            } if !(positive(time_scale) && positive(max_motion) && positive(tau_max)) => {
                Err(CutError::InvalidOptions {
                    reason: format!(
                        "time_scale, max_motion and tau_max must be positive and finite, got \
                         {time_scale}, {max_motion} and {tau_max}"
                    ),
                })
            }
            _ => Ok(()),
        }?;
        if let Some(classes) = &self.symmetry {
            check_length("symmetry classes", "vertices", classes.len(), n_vertices)?;
        }
        Ok(())
    }
}

/// When [`Flow::run`] stops (SPEC §4).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct StopRule {
    /// Stop after this many steps.
    pub max_steps: usize,
    /// Stop when the cut motion of a step ([`StepReport::cut_motion`], a length on the
    /// unit-area surface) is at most this. `0.0` never stops for this reason.
    pub motion_tolerance: f64,
    /// Stop when the energy has not fallen by more than this fraction of its size over the
    /// last `energy_window` steps: `E_{k − window} − E_k ≤ energy_tolerance · |E_k|`. On an
    /// irregular mesh the rate of the cut has a noise floor of the size of the speeds that
    /// balance each other, the cut keeps moving by a fraction of an edge per step, and the
    /// motion never falls below a small tolerance; the energy then stands still. `0.0` never
    /// stops for this reason.
    pub energy_tolerance: f64,
    /// The number of steps the energy is compared over, `≥ 1`.
    pub energy_window: usize,
}

impl Default for StopRule {
    /// At most 1000 steps; stop when the cut moves by at most `10⁻⁵` in a step, or when the
    /// energy has fallen by less than `10⁻⁴` of its size over 25 steps (SPEC §4, D9f).
    fn default() -> Self {
        StopRule {
            max_steps: 1000,
            motion_tolerance: 1e-5,
            energy_tolerance: 1e-4,
            energy_window: 25,
        }
    }
}

/// Why [`Flow::run`] stopped.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum StopReason {
    /// The limit on the number of steps was reached.
    MaxSteps,
    /// The cut motion of the last step was within the tolerance.
    Converged,
    /// The energy stood still over the window of the stop rule.
    EnergyStalled,
    /// One region is empty: there is no cut left.
    CutVanished,
}

/// `|∇φ|` of the linear interpolant in the crossed faces: 1 for a signed distance function.
/// With it one sees whether the level-set function drifts away from a distance near the cut.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct GradientNorm {
    /// The smallest value over the crossed faces.
    pub min: f64,
    /// The largest value.
    pub max: f64,
    /// The mean, weighted by the lengths of the cut segments.
    pub mean: f64,
}

/// Wall-clock time of the stages of one step, in seconds.
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct StepTimings {
    /// Stages 1–3: cut mesh, scale factor, adjoint, speed (one factorization).
    pub evaluate: f64,
    /// Stage 4: harmonic extension (a second factorization only on meshes with boundary).
    pub extension: f64,
    /// Stage 6: the implicit step (a factorization only when `τ` or the weights changed).
    pub length: f64,
    /// Stage 7: redistancing; and the symmetry projections.
    pub redistance: f64,
}

/// What one flow step did. Lengths are on the unit-area surface.
#[derive(Debug, Clone, PartialEq)]
pub struct StepReport {
    /// The flow time before the step.
    pub time: f64,
    /// The time step that was taken (0 if there was no cut).
    pub tau: f64,
    /// The energy of the state the step started from, with all its terms.
    pub energy: EnergyTerms,
    /// The largest `|s_c|` on the cut of the state the step started from.
    pub max_speed: f64,
    /// The largest `|λ⁰_c|`, the length force per unit length of the cut ([`LengthStep::OnCut`];
    /// 0 for [`LengthStep::OnMesh`], where the length term does not pass through the cut
    /// points).
    pub max_length_speed: f64,
    /// The roughening rate `ω` of the state the step started from: the
    /// [`ROUGHENING_QUANTILE`] quantile, weighted by the length of the cut, of the local rates
    /// `ω_c = 2 Σ_r Q_r,c` at the cut points, with `Q_r = G_r − β_r − w_r b(0)` the part of the
    /// gradient density that is quadratic in the normal derivatives (`w_r (∂u/∂n)²` for the
    /// Dirichlet energy). Where the local rate is `ω_c`, the distortion term amplifies a ripple
    /// of the cut of wavenumber `k` at the rate `ω_c k`, and the length term damps it at the
    /// rate `(α_L + η) k²` (SPEC §1, D18, D22): on nine tenths of the cut, ripples shorter than
    /// `2π (α_L + η)/ω` decay. If that length is below a few mesh edges, the flow amplifies
    /// the noise of the mesh; with `α_L + η = 0` it does so at every scale. For energies other
    /// than the Dirichlet energy `ω_c` can be negative.
    pub roughening_rate: f64,
    /// The largest local roughening rate `ω_c` on the cut. Far above
    /// [`StepReport::roughening_rate`] where the cut passes a spike of curvature.
    pub roughening_rate_max: f64,
    /// The speed of the cut that entered the step rule ([`StepStatistics`]): `m̂`, in crossed
    /// edges per unit time, or the largest `|ρ_c|`.
    pub step_speed: f64,
    /// `h`: the shortest crossed edge, on the unit-area surface.
    pub shortest_crossed_edge: f64,
    /// Number of cut points that were slowed down in this step by the local limit of the step
    /// rule ([`StepSize::local_factor`]). On a network: the cut points and the triple points
    /// with a copy that was slowed down, each counted once.
    pub n_limited: usize,
    /// The largest `|φ_new − φ_old|` over the vertices that are cut-adjacent before or after
    /// the step: how far the cut moved.
    pub cut_motion: f64,
    /// `|∇φ|` in the crossed faces of the state the step started from; `None` without a cut.
    pub gradient_norm: Option<GradientNorm>,
    /// Number of cut points before the step.
    pub n_cut_points: usize,
    /// Number of crossings that were moved away from a vertex in the evaluation.
    pub n_clamped_crossings: usize,
    /// Number of closed components the cut did not touch.
    pub n_uncut_closed_components: usize,
    /// Number of vertices the redistancing did not reach (no cut-adjacent vertex of their sign
    /// in their component); they keep their value.
    pub n_unreached: usize,
    /// Residual `max |(L σ)_i| / max |(L_OC s_C)_i|` of the harmonic extension.
    pub extension_residual: f64,
    /// Residual of the implicit step, relative to the largest entry of its right-hand side;
    /// `None` if the step was skipped (`α_L + η ≡ 0`).
    pub length_residual: Option<f64>,
    /// Whether the matrix of the implicit step was factorized in this step (it is reused while
    /// `τ`, `α_L`, `η` and the rule do not change).
    pub length_factorized: bool,
    /// Whether no cut is left after the step (or there was none to begin with).
    pub cut_vanished: bool,
    /// Wall-clock time of the stages.
    pub timings: StepTimings,
}

/// The intermediate functions of one step, on the unit-area surface, one value per vertex:
/// what the stages of SPEC §4 produced. For tests and for comparing implementations.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct StepStages {
    /// `ω_c` at every cut point: the local roughening rate (SPEC §1).
    pub roughening: Vec<f64>,
    /// `ℓ_c` at every cut point: the length of the cut that belongs to it, the weight of the
    /// quantiles of the step rule (SPEC §1).
    pub cut_length: Vec<f64>,
    /// `h_c` at every cut point: the shortest crossed edge that shares a vertex with its edge
    /// (SPEC §4).
    pub local_edge: Vec<f64>,
    /// `Θ` at every cut point: the integrated length force. Computed with
    /// [`LengthStep::OnCut`], and with [`LengthStep::OnMesh`] where the step rule needs it
    /// ([`StepStatistics::Robust`] with a step that is not fixed); empty otherwise.
    pub length_force: Vec<f64>,
    /// `λ⁰` at every cut point: the length force per unit length of the cut. Computed and
    /// empty as `length_force`.
    pub length_speed: Vec<f64>,
    /// `ḡ` at every cut point: the mean of `|∇φ|` over its piece of the cut (as computed; it
    /// is applied with [`SpeedScaling::MeanGradient`] only).
    pub gradient_scale: Vec<f64>,
    /// The rate at every cut point before the local limit of the step rule and before the
    /// scaling: `s`, or with [`LengthStep::OnCut`] the rate `r`.
    pub unscaled_rate: Vec<f64>,
    /// The factor of the local limit of the step rule at every cut point
    /// ([`StepSize::local_factor`]): 1 where the point was not slowed down. The rate that is
    /// scaled and extended is `unscaled_rate` times this factor.
    pub limit_factor: Vec<f64>,
    /// The rate at every cut point that was extended: `s` (with [`LengthStep::OnCut`] the rate
    /// `r` that is implicit in the length force), times the mean `|∇φ|` with
    /// [`SpeedScaling::MeanGradient`].
    pub rate: Vec<f64>,
    /// `σ`, the extension of the rate on the cut (stage 4).
    pub extension: Vec<f64>,
    /// `φ̃ = φ + τ σ` (stage 5).
    pub explicit: Vec<f64>,
    /// `φ′` after the length step (stage 6; equal to `explicit` if there is none).
    pub length: Vec<f64>,
    /// The factor of the local limit at every vertex (SPEC §4, D26): with
    /// [`LengthStep::OnMesh`] and [`StepStatistics::Robust`] the change `φ′ − φ` of the whole
    /// step is limited at the vertices. Empty where that limit is not in use.
    pub vertex_factor: Vec<f64>,
    /// `h_i` at every vertex for that limit: the shortest crossed edge at a cut-adjacent
    /// vertex, the shortest edge at any other one. Empty as `vertex_factor`.
    pub vertex_edge: Vec<f64>,
    /// `φ + vertex_factor · (φ′ − φ)`: after the local limit at the vertices (equal to `length`
    /// where there is none).
    pub limited: Vec<f64>,
    /// After the first symmetry projection (equal to `limited` without symmetry classes).
    pub symmetric: Vec<f64>,
    /// Whether the vertex is cut-adjacent in the function that was redistanced.
    pub adjacent: Vec<bool>,
}

/// The history of a run.
#[derive(Debug, Clone, PartialEq)]
pub struct RunReport {
    /// One report per step taken.
    pub steps: Vec<StepReport>,
    /// Why the run stopped.
    pub stop: StopReason,
    /// The energy of the final state.
    pub final_energy: EnergyTerms,
}

/// The state of the implicit step that its factorization belongs to.
#[derive(Debug, Clone, PartialEq)]
struct LengthKey {
    tau: f64,
    rule: PenaltyRule,
    weights: Vec<f64>,
}

/// The cut flow on one mesh: the evaluator, the matrices on the mesh and the level-set function.
///
/// The level-set function crosses the API in the caller's length unit; inside it is kept on the
/// unit-area surface. A flow can be created from, and reset to, any level-set function
/// ([`Flow::new`], [`Flow::set_level_set`]); the problem and the options are arguments of every
/// step, so every control of SPEC §5 can change between steps.
#[derive(Debug)]
pub struct Flow {
    pub(crate) evaluator: Evaluator,
    /// The geometry of the unit-area surface, for the redistancing.
    pub(crate) geometry: Geometry,
    /// The level-set function on the unit-area surface.
    phi: Vec<f64>,
    time: f64,
    n_steps: usize,
    /// The Laplacian of the cut mesh on all vertices, for the harmonic extension on a mesh with
    /// boundary. On a closed mesh the system of the evaluator is that matrix, and this one only
    /// carries the pattern (shared symbolic analysis) for the length step.
    extension: LaplaceSystem,
    /// For every vertex, where its neighbours start in `neighbours`.
    neighbour_start: Vec<u32>,
    /// The neighbours of every vertex with the edge that leads to them.
    neighbours: Vec<(u32, u32)>,
    /// Whether the mesh has no boundary.
    pub(crate) closed: bool,
    /// The matrix of the implicit step; it shares the symbolic analysis of the pattern on all
    /// vertices.
    length: LaplaceSystem,
    length_key: Option<LengthKey>,
}

impl Flow {
    /// Sets up the flow on `mesh` with the intrinsic `geometry`, starting from the level-set
    /// function `phi` (one value per vertex, in the length unit of `geometry`; region 0 is where
    /// it is `≤ 0`). `phi` is taken as it is: see [`crate::init`] for ways to make one.
    ///
    /// # Errors
    /// As [`Evaluator::new`]; [`CutError::WrongLength`], [`CutError::NonFiniteLevelSet`].
    pub fn new(mesh: &TriMesh, geometry: &Geometry, phi: &[f64]) -> Result<Flow, CutError> {
        let evaluator = Evaluator::new(mesh, geometry)?;
        let scale = evaluator.scale();
        let scaled = Geometry::from_edge_lengths(
            mesh,
            geometry
                .edge_lengths()
                .iter()
                .map(|l| scale.length_to_unit_area(*l))
                .collect(),
        )?;
        let closed = mesh.n_boundary_edges() == 0;
        let extension = if closed {
            evaluator.system().clone_pattern()
        } else {
            LaplaceSystem::new(mesh, |_| true)?
        };
        let length = extension.clone_pattern();
        // Vertex adjacency in compressed form.
        let nv = mesh.n_vertices();
        let mut neighbour_start = vec![0u32; nv + 1];
        for &[a, b] in mesh.edges() {
            neighbour_start[a as usize + 1] += 1;
            neighbour_start[b as usize + 1] += 1;
        }
        for v in 0..nv {
            neighbour_start[v + 1] += neighbour_start[v];
        }
        let mut next = neighbour_start.clone();
        let mut neighbours = vec![(0u32, 0u32); 2 * mesh.n_edges()];
        for (e, &[a, b]) in mesh.edges().iter().enumerate() {
            for (from, to) in [(a, b), (b, a)] {
                neighbours[next[from as usize] as usize] = (to, e as u32);
                next[from as usize] += 1;
            }
        }
        let mut flow = Flow {
            evaluator,
            geometry: scaled,
            phi: Vec::new(),
            time: 0.0,
            n_steps: 0,
            extension,
            neighbour_start,
            neighbours,
            closed,
            length,
            length_key: None,
        };
        flow.set_level_set(phi)?;
        Ok(flow)
    }

    /// The mesh.
    #[inline]
    pub fn mesh(&self) -> &TriMesh {
        self.evaluator.mesh()
    }

    /// The scaling of the input surface to unit area.
    #[inline]
    pub fn scale(&self) -> UnitAreaScale {
        self.evaluator.scale()
    }

    /// The evaluator of the flow (for its curvature, vertex areas, feature-size weight).
    #[inline]
    pub fn evaluator(&mut self) -> &mut Evaluator {
        &mut self.evaluator
    }

    /// The level-set function, in the caller's length unit.
    pub fn level_set(&self) -> Vec<f64> {
        let scale = self.scale();
        self.phi
            .iter()
            .map(|x| scale.length_from_unit_area(*x))
            .collect()
    }

    /// The level-set function on the unit-area surface, as the flow keeps it.
    #[inline]
    pub fn level_set_unit_area(&self) -> &[f64] {
        &self.phi
    }

    /// Replaces the level-set function (caller's length unit) and resets the flow time and the
    /// step count: the flow resumes from this state.
    ///
    /// # Errors
    /// [`CutError::WrongLength`], [`CutError::NonFiniteLevelSet`].
    pub fn set_level_set(&mut self, phi: &[f64]) -> Result<(), CutError> {
        check_length(
            "level-set function",
            "vertices",
            phi.len(),
            self.mesh().n_vertices(),
        )?;
        if let Some(v) = phi.iter().position(|x| !x.is_finite()) {
            return Err(CutError::NonFiniteLevelSet {
                vertex: v as u32,
                value: phi[v],
            });
        }
        let scale = self.scale();
        self.phi = phi.iter().map(|x| scale.length_to_unit_area(*x)).collect();
        self.time = 0.0;
        self.n_steps = 0;
        Ok(())
    }

    /// The flow time: the sum of the steps taken since the level-set function was set.
    #[inline]
    pub fn time(&self) -> f64 {
        self.time
    }

    /// The number of steps taken since the level-set function was set.
    #[inline]
    pub fn n_steps(&self) -> usize {
        self.n_steps
    }

    /// Evaluates the current state (SPEC §3.1–§3.7, §6) without changing it.
    ///
    /// # Errors
    /// As [`Evaluator::evaluate`].
    pub fn evaluate(
        &mut self,
        problem: &CutProblem,
        options: &EvalOptions<'_>,
    ) -> Result<Evaluation, CutError> {
        self.evaluator.evaluate(&self.phi, problem, options)
    }

    /// The roughening rate `ω` of the current state ([`StepReport::roughening_rate`]), without
    /// taking a step: the scale for the length weight ([`length_weight_for_ripples`]) and for
    /// the time step ([`StepSize::Adaptive`]). `0` if there is no cut.
    ///
    /// # Errors
    /// As [`Evaluator::evaluate`].
    pub fn roughening_rate(
        &mut self,
        problem: &CutProblem,
        options: &EvalOptions<'_>,
    ) -> Result<f64, CutError> {
        let e = self.evaluator.evaluate(&self.phi, problem, options)?;
        Ok(roughening_rate(&e, problem))
    }

    /// Makes the current level-set function a signed distance to its own cut, by the
    /// redistancing of SPEC §3.12 (cut-adjacent vertices keep their values). Returns the number
    /// of vertices no front reached.
    ///
    /// # Errors
    /// [`CutError::Distance`] (cannot happen for a finite function).
    pub fn redistance(&mut self) -> Result<usize, CutError> {
        let (phi, unreached, _) = redistance(self.evaluator.mesh(), &self.geometry, &self.phi)?;
        self.phi = phi;
        Ok(unreached)
    }

    /// One step of the flow (SPEC §4): evaluation, harmonic extension of the speed, explicit
    /// step, implicit length step, symmetry, redistancing.
    ///
    /// If the current state has no cut, nothing changes and the report says so
    /// ([`StepReport::cut_vanished`]).
    ///
    /// # Errors
    /// As [`Evaluator::evaluate`]; [`CutError::InvalidOptions`] for a step size that is not
    /// positive, symmetry classes of the wrong length, or the paper's penalty rule with
    /// `α_L + η_i = 0` at some but not all vertices.
    pub fn step(
        &mut self,
        problem: &CutProblem,
        options: &FlowOptions<'_>,
    ) -> Result<StepReport, CutError> {
        self.advance(problem, options, None)
    }

    /// [`Flow::step`], also returning the intermediate functions of the stages. If the state has
    /// no cut the stages are empty.
    ///
    /// # Errors
    /// As [`Flow::step`].
    pub fn step_with_stages(
        &mut self,
        problem: &CutProblem,
        options: &FlowOptions<'_>,
    ) -> Result<(StepReport, StepStages), CutError> {
        let mut stages = StepStages::default();
        let report = self.advance(problem, options, Some(&mut stages))?;
        Ok((report, stages))
    }

    fn advance(
        &mut self,
        problem: &CutProblem,
        options: &FlowOptions<'_>,
        mut stages: Option<&mut StepStages>,
    ) -> Result<StepReport, CutError> {
        let nv = self.mesh().n_vertices();
        options.validate(nv)?;
        let par = options.eval.parallelism;
        let mut timings = StepTimings::default();

        // Stages 1–3.
        let start = Instant::now();
        let e = self.evaluator.evaluate(&self.phi, problem, &options.eval)?;
        timings.evaluate = start.elapsed().as_secs_f64();
        let max_speed = e.speed.iter().fold(0.0f64, |m, s| m.max(s.abs()));
        let mut report = StepReport {
            time: self.time,
            tau: 0.0,
            energy: e.energy.clone(),
            max_speed,
            max_length_speed: 0.0,
            roughening_rate: 0.0,
            roughening_rate_max: 0.0,
            step_speed: 0.0,
            shortest_crossed_edge: 0.0,
            n_limited: 0,
            cut_motion: 0.0,
            gradient_norm: None,
            n_cut_points: e.cut.points.len(),
            n_clamped_crossings: e.report.n_clamped_crossings,
            n_uncut_closed_components: e.report.uncut_closed_components.len(),
            n_unreached: 0,
            extension_residual: 0.0,
            length_residual: None,
            length_factorized: false,
            cut_vanished: e.cut.points.is_empty(),
            timings,
        };
        if report.cut_vanished {
            return Ok(report);
        }
        let (gradient_norm, mean_gradient) =
            gradient_norms(&e.cut, &e.triangles, |_, a| self.phi[a].abs());
        report.gradient_norm = Some(gradient_norm);
        let local_rates = roughening_rates(&e, problem);
        let dual = cut_dual_lengths(&e.cut, &e.triangles);
        report.roughening_rate = weighted_quantile(&local_rates, &dual, ROUGHENING_QUANTILE);
        report.roughening_rate_max = weighted_quantile(&local_rates, &dual, 1.0);
        // The statistics of the step rule (SPEC §4).
        let statistics = options.step_statistics;
        let omega = match statistics {
            StepStatistics::Robust => report.roughening_rate,
            StepStatistics::Extremes => report.roughening_rate_max,
        };
        let edges = local_edges(self.evaluator.mesh(), &self.geometry, &e.cut);
        report.shortest_crossed_edge = edges.iter().fold(f64::INFINITY, |m, h| m.min(*h));
        // Stage 6 on the cut (SPEC §3.10): the length force, the step, and the rate that is
        // implicit in the length force.
        let start = Instant::now();
        let weights = length_weights(problem, nv);
        let active = weights.iter().any(|w| *w != 0.0);
        let mut rate = e.speed.clone();
        let (mut tau, speed) = rule_step(
            &options.step_size,
            statistics,
            &e.speed,
            &edges,
            &dual,
            omega,
        );
        report.step_speed = speed;
        // With the length step on the mesh the length force is still computed on the cut for
        // the robust step rule: it is part of the speed of the cut (SPEC §4). The step is
        // then the smaller of the two that the speeds `s` and `s + λ⁰` allow.
        let on_cut_step = options.length_step == LengthStep::OnCut;
        let for_step_only = !on_cut_step
            && statistics == StepStatistics::Robust
            && !matches!(options.step_size, StepSize::Fixed(_));
        if active && (on_cut_step || for_step_only) {
            let smoothing_length = match options.eval.flux_recovery {
                crate::problem::FluxRecovery::Smoothed { length } => length,
                crate::problem::FluxRecovery::Nodal => 0.0,
            };
            let on_cut = rate_on_cut(
                self.evaluator.mesh(),
                self.evaluator.base(),
                &e.cut,
                &e.triangles,
                &e.speed,
                &weights,
                options.penalty_rule,
                smoothing_length,
                &options.step_size,
                statistics,
                &edges,
                omega,
                par,
            )?;
            if on_cut_step {
                tau = on_cut.tau;
                report.step_speed = on_cut.step_speed;
            } else {
                // The step is split: neither the explicit part `τ s` nor the whole motion
                // `τ (s + λ⁰)` may exceed the limit.
                tau = tau.min(on_cut.tau);
                report.step_speed = report.step_speed.max(on_cut.step_speed);
            }
            if on_cut_step {
                report.length_residual = Some(on_cut.residual);
                report.max_length_speed = on_cut
                    .length_speed
                    .iter()
                    .fold(0.0f64, |m, x| m.max(x.abs()));
                rate = on_cut.rate;
            }
            if let Some(stages) = stages.as_deref_mut() {
                stages.length_force = on_cut.force;
                stages.length_speed = on_cut.length_speed;
            }
        }
        report.tau = tau;
        if let Some(stages) = stages.as_deref_mut() {
            stages.roughening = local_rates.clone();
            stages.cut_length = dual.clone();
            stages.local_edge = edges.clone();
            stages.gradient_scale = mean_gradient.clone();
            stages.unscaled_rate = rate.clone();
        }
        // The local limit of the step rule: cut points beyond the quantile take a smaller
        // step of their own. With the length step on the mesh the rate on the cut is only the
        // explicit part of the motion: limiting it alone would leave the length step, which
        // opposes it, in full, and the cut would move uphill. There the limit acts on the
        // whole change at the vertices, after the length step (SPEC §4, D26).
        let mut factors = vec![1.0; rate.len()];
        if statistics == StepStatistics::Robust && !for_step_only {
            for (p, r) in rate.iter_mut().enumerate() {
                let factor = options.step_size.local_factor(
                    tau,
                    report.step_speed,
                    edges[p],
                    *r,
                    local_rates[p] - omega,
                );
                if factor < 1.0 {
                    *r *= factor;
                    factors[p] = factor;
                    report.n_limited += 1;
                }
            }
        }
        if let Some(stages) = stages.as_deref_mut() {
            stages.limit_factor = factors;
        }
        report.timings.length = start.elapsed().as_secs_f64();

        // Stage 4: σ, the extension of the rate on the cut.
        let start = Instant::now();
        if options.speed_scaling == SpeedScaling::MeanGradient {
            for (r, g) in rate.iter_mut().zip(&mean_gradient) {
                *r *= g;
            }
        }
        let (sigma, residual) = self.extend(&e, &rate, options.extension, par)?;
        report.extension_residual = residual;
        report.timings.extension = start.elapsed().as_secs_f64();

        // Stage 5: the explicit step.
        let mut phi: Vec<f64> = self
            .phi
            .iter()
            .zip(&sigma)
            .map(|(x, s)| explicit_step(*x, tau, *s))
            .collect();

        if let Some(stages) = stages.as_deref_mut() {
            stages.rate = rate;
            stages.extension = sigma;
            stages.explicit = phi.clone();
        }

        // Stage 6 on the mesh (the paper's), then the symmetry projection.
        let start = Instant::now();
        if active && options.length_step == LengthStep::OnMesh {
            let (residual, factorized) =
                self.length_step(&weights, options.penalty_rule, tau, &mut phi, par)?;
            report.length_residual = Some(residual);
            report.length_factorized = factorized;
        }
        report.timings.length += start.elapsed().as_secs_f64();
        let start = Instant::now();
        if let Some(stages) = stages.as_deref_mut() {
            stages.length = phi.clone();
        }
        let adjacent_before = e.cut.cut_adjacent_vertices(self.evaluator.mesh());
        // The local limit at the vertices: the change of the whole step, both parts together.
        if for_step_only {
            let (vertex_edges, vertex_rates) =
                vertex_limits(self.evaluator.mesh(), &self.geometry, &e.cut, &local_rates);
            let mut vertex_factor = vec![1.0; nv];
            for (v, (new, old)) in phi.iter_mut().zip(&self.phi).enumerate() {
                let change = *new - old;
                let factor = options.step_size.local_factor(
                    tau,
                    0.0,
                    vertex_edges[v],
                    change / tau,
                    vertex_rates[v] - omega,
                );
                if factor < 1.0 {
                    *new = old + factor * change;
                    vertex_factor[v] = factor;
                    report.n_limited += usize::from(adjacent_before[v]);
                }
            }
            if let Some(stages) = stages.as_deref_mut() {
                stages.vertex_factor = vertex_factor;
                stages.vertex_edge = vertex_edges;
            }
        }
        if let Some(stages) = stages.as_deref_mut() {
            stages.limited = phi.clone();
        }
        if let Some(classes) = &options.symmetry {
            symmetry_projection(classes, &mut phi);
        }
        if let Some(stages) = stages.as_deref_mut() {
            stages.symmetric = phi.clone();
        }

        // Stage 7: redistancing, then the symmetry projection again.
        let (adjacent_after, vanished) = match options.redistance {
            Redistance::KeepAdjacent => {
                let (new, unreached, adjacent) =
                    redistance(self.evaluator.mesh(), &self.geometry, &phi)?;
                report.n_unreached = unreached;
                phi = new;
                let vanished = !adjacent.iter().any(|&a| a);
                (adjacent, vanished)
            }
            Redistance::Off => {
                let adjacent = cut_adjacent(self.evaluator.mesh(), &phi);
                let vanished = !adjacent.iter().any(|&a| a);
                (adjacent, vanished)
            }
        };
        if let Some(classes) = &options.symmetry {
            symmetry_projection(classes, &mut phi);
        }
        if let Some(stages) = stages {
            stages.adjacent = adjacent_after.clone();
        }
        report.timings.redistance = start.elapsed().as_secs_f64();
        report.cut_vanished = vanished;
        report.cut_motion = (0..nv)
            .filter(|&v| adjacent_before[v] || adjacent_after[v])
            .map(|v| (phi[v] - self.phi[v]).abs())
            .fold(0.0, f64::max);

        self.phi = phi;
        self.time += tau;
        self.n_steps += 1;
        Ok(report)
    }

    /// Steps until `stop` says so, with the same problem and options in every step, and returns
    /// the history. For schedules (SPEC §4) call [`Flow::step`] or this function repeatedly
    /// with changed arguments; the state carries over.
    ///
    /// # Errors
    /// As [`Flow::step`].
    pub fn run(
        &mut self,
        problem: &CutProblem,
        options: &FlowOptions<'_>,
        stop: &StopRule,
    ) -> Result<RunReport, CutError> {
        let mut steps = Vec::new();
        let reason = loop {
            if steps.len() >= stop.max_steps {
                break StopReason::MaxSteps;
            }
            let report = self.step(problem, options)?;
            let (vanished, motion) = (report.cut_vanished, report.cut_motion);
            steps.push(report);
            if vanished {
                break StopReason::CutVanished;
            }
            if motion <= stop.motion_tolerance {
                break StopReason::Converged;
            }
            let window = stop.energy_window.max(1);
            if stop.energy_tolerance > 0.0 && steps.len() > window {
                let now = steps[steps.len() - 1].energy.total;
                let before = steps[steps.len() - 1 - window].energy.total;
                if before - now <= stop.energy_tolerance * now.abs() {
                    break StopReason::EnergyStalled;
                }
            }
        };
        let final_energy = self.evaluate(problem, &options.eval)?.energy;
        Ok(RunReport {
            steps,
            stop: reason,
            final_energy,
        })
    }

    /// The extension of the rate `rate` on the cut points to every vertex of the mesh
    /// (SPEC §3.8), and the residual of the solve (0 if there is none).
    fn extend(
        &mut self,
        e: &Evaluation,
        rate: &[f64],
        rule: Extension,
        par: Parallelism<'_>,
    ) -> Result<(Vec<f64>, f64), CutError> {
        match rule {
            Extension::Local => Ok((self.extend_from_cut(&e.cut, &e.triangles, rate), 0.0)),
            Extension::Harmonic => self.extend_harmonically(e, rate, par),
        }
    }

    /// [`Extension::Local`]: with `r_p` the rate at cut point `p`, `ℓ_p` its dual length on the
    /// cut (half the lengths of its segments) and `β_ip` the barycentric weight of the
    /// cut-adjacent vertex `i` at `p` (`1 − t` or `t` on the edge of `p`),
    ///
    /// ```text
    /// σ_i = Σ_p ℓ_p β_ip r_p / Σ_p ℓ_p β_ip               over the crossed edges at i
    /// σ_i = Σ_j (σ_j / ℓ_ij) / Σ_j (1 / ℓ_ij)              over the neighbours j of i in its region
    ///                                                      that are one ring closer to the cut
    /// ```
    ///
    /// where the ring of a vertex is the number of uncrossed edges on the shortest path from it
    /// to a cut-adjacent vertex of its region. A vertex that no such path reaches gets `σ = 0`.
    ///
    /// The first line is the transpose of the interpolation `Σ_i β_ip σ_i` that the explicit
    /// step applies at the cut points, normalized: `Σ_p ℓ_p Σ_i β_ip σ_i = Σ_p ℓ_p r_p`, the
    /// mean rate is kept.
    ///
    /// [formula: cut-extension-local]
    pub(crate) fn extend_from_cut(&self, cut: &Cut, tri: &CutTriangles, rate: &[f64]) -> Vec<f64> {
        let mesh = self.evaluator.mesh();
        let region = &cut.vertex_region;
        let nv = region.len();
        let dual = cut_dual_lengths(cut, tri);
        let (mut sum, mut weight) = (vec![0.0; nv], vec![0.0; nv]);
        for (k, point) in cut.points.iter().enumerate() {
            let [a, b] = mesh.edge_vertices(point.edge as usize);
            for (v, barycentric) in [(a, 1.0 - point.t), (b, point.t)] {
                sum[v] += dual[k] * barycentric * rate[k];
                weight[v] += dual[k] * barycentric;
            }
        }
        self.extend_from_seeds(region, sum, weight)
    }

    /// The second line of [`Extension::Local`]: `sum[i] / weight[i]` at the vertices with a
    /// positive weight (the cut-adjacent ones), and from there ring by ring within each region.
    ///
    /// [formula: cut-extension-local]
    pub(crate) fn extend_from_seeds(
        &self,
        region: &[u32],
        sum: Vec<f64>,
        weight: Vec<f64>,
    ) -> Vec<f64> {
        let lengths = self.geometry.edge_lengths();
        let nv = region.len();
        const UNSEEN: u32 = u32::MAX;
        let mut ring = vec![UNSEEN; nv];
        let mut sigma = vec![0.0; nv];
        let mut front: Vec<u32> = Vec::new();
        for v in 0..nv {
            if weight[v] > 0.0 {
                sigma[v] = sum[v] / weight[v];
                ring[v] = 0;
                front.push(v as u32);
            }
        }
        let mut k = 0;
        while !front.is_empty() {
            k += 1;
            let mut next: Vec<u32> = Vec::new();
            for &v in &front {
                let range = self.neighbour_start[v as usize] as usize
                    ..self.neighbour_start[v as usize + 1] as usize;
                for &(u, _) in &self.neighbours[range] {
                    if ring[u as usize] == UNSEEN && region[u as usize] == region[v as usize] {
                        ring[u as usize] = k;
                        next.push(u);
                    }
                }
            }
            for &u in &next {
                let range = self.neighbour_start[u as usize] as usize
                    ..self.neighbour_start[u as usize + 1] as usize;
                let (mut sum, mut weight) = (0.0, 0.0);
                for &(j, edge) in &self.neighbours[range] {
                    if ring[j as usize] == k - 1 && region[j as usize] == region[u as usize] {
                        sum += sigma[j as usize] / lengths[edge as usize];
                        weight += 1.0 / lengths[edge as usize];
                    }
                }
                sigma[u as usize] = sum / weight;
            }
            front = next;
        }
        sigma
    }

    /// [`Extension::Harmonic`]: in each region `L_OO σ_O = −L_OC s_C` with the cotan Laplacian
    /// of the cut mesh, all vertices of the mesh unknown (zero flux through `∂M`) and `σ = 0`
    /// in a component without cut vertices. Returns `σ` and the residual of the solve.
    ///
    /// [formula: cut-harmonic-extension]
    fn extend_harmonically(
        &mut self,
        e: &Evaluation,
        rate: &[f64],
        par: Parallelism<'_>,
    ) -> Result<(Vec<f64>, f64), CutError> {
        let tri = &e.triangles;
        let mut sigma = vec![0.0; tri.n_vertices()];
        for (k, &s) in rate.iter().enumerate() {
            for v in tri.point_vertices(k) {
                sigma[v] = s;
            }
        }
        self.extend_harmonically_from(e, sigma, par)
    }

    /// The harmonic extension of the values that `sigma` has at the boundary vertices of the
    /// cut mesh (the copies of the cut points and of the triple points; its entries at the
    /// vertices of the input mesh are ignored and must be zero): the values at the vertices of
    /// the input mesh, and the residual of the solve.
    ///
    /// [formula: cut-harmonic-extension]
    pub(crate) fn extend_harmonically_from(
        &mut self,
        e: &Evaluation,
        mut sigma: Vec<f64>,
        par: Parallelism<'_>,
    ) -> Result<(Vec<f64>, f64), CutError> {
        let tri = &e.triangles;
        let nv = tri.n_original;
        let n = tri.n_vertices();
        // (L x)_O for x = s on the cut vertices and 0 elsewhere is L_OC s_C.
        let mut product = vec![0.0; n];
        tri.apply_laplacian(None, &sigma, &mut product)?;
        let rhs = |v: usize| -product[v];
        let largest = product[..nv].iter().fold(0.0f64, |m, x| m.max(x.abs()));
        let mesh = self.evaluator.mesh();
        // Vertices of components without a cut vertex: σ = 0 there.
        let pinned: Vec<bool> = if self.closed {
            let components = mesh.vertex_components();
            components
                .iter()
                .map(|c| e.report.uncut_closed_components.contains(c))
                .collect()
        } else {
            components_without_cut(mesh, &e.cut.vertex_region, &e.cut)
        };
        if self.closed {
            // The matrix is the one the evaluation factorized, with identity rows for the
            // closed components without a cut.
            self.evaluator.system().solve(rhs, &mut sigma, par)?;
        } else {
            self.extension.assemble(mesh, tri, &pinned);
            self.extension.factor(par)?;
            self.extension.solve(rhs, &mut sigma, par)?;
        }
        tri.apply_laplacian(None, &sigma, &mut product)?;
        let worst = (0..nv)
            .filter(|&v| !pinned[v])
            .map(|v| product[v].abs())
            .fold(0.0f64, f64::max);
        sigma.truncate(nv);
        let residual = if largest > 0.0 {
            worst / largest
        } else {
            worst
        };
        Ok((sigma, residual))
    }

    /// The implicit step on the mesh (SPEC §3.10, [`LengthStep::OnMesh`]) applied to `phi` in
    /// place, for vertex weights `w = α_L + η` that are not all zero. Returns the residual and
    /// whether the matrix was factorized.
    ///
    /// [formula: cut-length-step]
    /// [formula: cut-penalty-step-paper]
    pub(crate) fn length_step(
        &mut self,
        weights: &[f64],
        rule: PenaltyRule,
        tau: f64,
        phi: &mut [f64],
        par: Parallelism<'_>,
    ) -> Result<(f64, bool), CutError> {
        let nv = phi.len();
        if rule == PenaltyRule::PaperDiagonal
            && let Some(v) = weights.iter().position(|w| *w == 0.0)
        {
            return Err(CutError::InvalidOptions {
                reason: format!(
                    "the paper's penalty rule divides by alpha_L + eta, which is 0 at vertex {v} \
                     (and not at every vertex)"
                ),
            });
        }
        let mesh = self.evaluator.mesh();
        let base = self.evaluator.base();
        let area = self.evaluator.vertex_areas();
        // The diagonal mass and the scale of every face's stiffness matrix.
        let mass = |v: usize| match rule {
            PenaltyRule::Gradient => area[v],
            PenaltyRule::PaperDiagonal => area[v] / weights[v],
        };
        let face_scale = |f: usize| match rule {
            PenaltyRule::Gradient => {
                let [a, b, c] = mesh.faces()[f];
                tau * (weights[a as usize] + weights[b as usize] + weights[c as usize]) / 3.0
            }
            PenaltyRule::PaperDiagonal => tau,
        };
        let key = LengthKey {
            tau,
            rule,
            weights: weights.to_vec(),
        };
        let factorized = self.length_key.as_ref() != Some(&key);
        if factorized {
            self.length_key = None;
            self.length
                .assemble_on_mesh(mesh, &base.cot, face_scale, mass);
            self.length.factor(par)?;
            self.length_key = Some(key);
        }
        let rhs: Vec<f64> = (0..nv).map(|v| mass(v) * phi[v]).collect();
        let mut solution = vec![0.0; nv];
        self.length.solve(|v| rhs[v], &mut solution, par)?;
        // Residual: mass · φ′ + Σ_f scale_f L_f φ′ − rhs.
        let mut product: Vec<f64> = (0..nv).map(|v| mass(v) * solution[v] - rhs[v]).collect();
        for (f, (face, cot)) in mesh.faces().iter().zip(&base.cot).enumerate() {
            let s = 0.5 * face_scale(f);
            for k in 0..3 {
                let (i, j) = (face[(k + 1) % 3] as usize, face[(k + 2) % 3] as usize);
                let flow = s * cot[k] * (solution[i] - solution[j]);
                product[i] += flow;
                product[j] -= flow;
            }
        }
        let worst = product.iter().fold(0.0f64, |m, x| m.max(x.abs()));
        let largest = rhs.iter().fold(0.0f64, |m, x| m.max(x.abs()));
        phi.copy_from_slice(&solution);
        let residual = if largest > 0.0 {
            worst / largest
        } else {
            worst
        };
        Ok((residual, factorized))
    }
}

/// The vertex weights `w_i = α_L + η_i` of the length terms.
pub(crate) fn length_weights(problem: &CutProblem, n_vertices: usize) -> Vec<f64> {
    match &problem.penalty {
        Some(eta) => eta.iter().map(|e| problem.length_weight + e).collect(),
        None => vec![problem.length_weight; n_vertices],
    }
}

/// What [`rate_on_cut`] returns.
pub(crate) struct RateOnCut {
    /// `Θ` at every cut point: the integrated length force.
    pub(crate) force: Vec<f64>,
    /// `λ⁰` at every cut point: the length force per unit length of the cut, recovered like the
    /// normal derivative.
    pub(crate) length_speed: Vec<f64>,
    /// `r` at every cut point: the rate of change of the level-set function.
    pub(crate) rate: Vec<f64>,
    /// The step `τ`: the step size rule applied to `s + λ⁰`.
    pub(crate) tau: f64,
    /// `ρ̂`: the quantile of `|s + λ⁰|` that entered the step rule.
    pub(crate) step_speed: f64,
    /// The larger of the residuals of the two systems, each relative to the largest entry of
    /// its right-hand side.
    pub(crate) residual: f64,
}

/// The length `ℓ_p` of the cut that belongs to every cut point: half the length of each of its
/// segments to another cut point, and the whole length of each of its segments to a triple
/// point (SPEC §3.10, §8.5). On a network these sum to the length of the cut.
pub(crate) fn cut_dual_lengths(cut: &Cut, tri: &CutTriangles) -> Vec<f64> {
    let n = cut.points.len();
    let mut dual = vec![0.0; n];
    for (segment, length) in cut.segments.iter().zip(&tri.segment_length) {
        let [p, q] = segment.points;
        if (q as usize) < n {
            dual[p as usize] += 0.5 * length;
            dual[q as usize] += 0.5 * length;
        } else {
            dual[p as usize] += length;
        }
    }
    dual
}

/// The rate `r` of the level-set function at every cut point with the length terms taken on the
/// cut ([`LengthStep::OnCut`], SPEC §3.10, D14).
///
/// **Force.** The integrated length force at cut point `p`, for the vertex weight
/// `w = α_L + η` interpolated to the cut, is
///
/// ```text
/// Θ_p = w_p · θ_p  +  Σ_{σ = pq at p} ½ |pq| (∂w/∂n)_σ
/// ```
///
/// * `θ_p = π − β_p` is the turning angle of the cut at `p`, with `β_p` the interior angle of
///   the lower region at `p` (the sum of the corner angles of its triangles there; the two
///   faces at a crossed edge unfold into a plane). At a cut point on `∂M` the cut is reflected
///   in the boundary, which is the zero-flux condition of the paper's step:
///   `θ_p = π − 2 β_p`, and half of it acts. The turning angle itself is used, not
///   `2 sin(θ_p/2)`, the derivative of the length under a displacement along the bisector:
///   the angles add up along the cut to the change of its direction whatever happens in
///   between, so the force on a run of cut points crowded around a mesh vertex is the turning
///   of the cut across the run, and the area enclosed by a closed cut on a flat surface
///   shrinks at the exact rate `2π w`.
/// * `(∂w/∂n)_σ` is the derivative of the linear `w` of the face of `σ` across the segment, out
///   of the lower region: with `a` the lone vertex of the face, `h_a` its height over `pq` and
///   `f` the foot of that height, `(w_f − w_a)/h_a`, negated if `a` is in the higher region.
///   The paper's penalty rule drops this term.
///
/// **Recovery.** `Θ_p/ℓ_p` with the dual length `ℓ_p = ½ Σ |pq|` is as unusable as `v_p/ℓ_p`
/// is for the normal derivative (SPEC D12), for the same reason. The force per unit length
/// `λ⁰` is recovered by the same system, with the same smoothing length `c`:
///
/// ```text
/// ℓ_p λ⁰_p + Σ_{σ = pq at p} ((c h_σ)² / |pq|) (λ⁰_p − λ⁰_q) = Θ_p
/// ```
///
/// **Rate.** Moving `p` out by `δ` relative to its neighbours turns the cut by
/// `Σ (δ_p − δ_q)/|pq|` more, so after a step in which the function changes at the rate `r` on
/// the cut, the force per unit length is `λ⁰ − τ ℓ⁻¹ S_w r` to first order, with the stiffness
/// `S_w` of the polyline weighted by `w̄_σ = ½ (w_p + w_q)`. Taking the length force at the end
/// of the step and the distortion speed `s` at its start, `r = s + λ⁰ − τ ℓ⁻¹ S_w r`:
///
/// ```text
/// ℓ_p r_p + τ Σ_{σ = pq at p} (w̄_σ / |pq|) (r_p − r_q) = ℓ_p (s_p + λ⁰_p)
/// ```
///
/// Both matrices are symmetric positive definite M-matrices with row sums `ℓ_p`, so
/// `Σ ℓ_p λ⁰_p = Σ Θ_p`, `Σ ℓ_p r_p = Σ ℓ_p (s_p + λ⁰_p)`, and `r` lies between the smallest
/// and the largest `s + λ⁰`. A cut on which `s + λ⁰ = 0` does not move, whatever `τ`.
///
/// **Step.** `τ` is the step size rule applied to `s + λ⁰`, which bounds `|r|`: with
/// [`StepSize::Capped`] the cut moves by at most the given fraction of the shortest crossed
/// edge, by distortion and length together.
///
/// This is the rule for two regions. On a network every region has rates of its own along its
/// whole boundary (`regions::region_rates`, SPEC §8.5).
///
/// [formula: cut-length-step-on-cut]
#[allow(clippy::too_many_arguments)]
pub(crate) fn rate_on_cut(
    mesh: &TriMesh,
    base: &crate::cut_mesh::BaseGeometry,
    cut: &Cut,
    tri: &CutTriangles,
    speed: &[f64],
    weights: &[f64],
    rule: PenaltyRule,
    smoothing_length: f64,
    step_size: &StepSize,
    statistics: StepStatistics,
    edges: &[f64],
    roughening_rate: f64,
    par: Parallelism<'_>,
) -> Result<RateOnCut, CutError> {
    let nv = tri.n_original;
    let n = cut.points.len();
    let w_point = cut.interpolate(mesh, weights);
    // Interior angles at the copies of the cut points.
    let mut angle = vec![0.0; tri.n_vertices() - nv];
    for (face, cot) in tri.faces.iter().zip(&tri.cot) {
        for k in 0..3 {
            let v = face[k] as usize;
            if v >= nv {
                angle[v - nv] += 1f64.atan2(cot[k]);
            }
        }
    }
    let mut force = vec![0.0; n];
    for (p, value) in force.iter_mut().enumerate() {
        let lower = tri.point_vertices_by_region(p)[0];
        let interior = angle[lower - nv];
        *value = if mesh.is_boundary_edge(cut.points[p].edge as usize) {
            w_point[p] * (0.5 * std::f64::consts::PI - interior)
        } else {
            w_point[p] * (std::f64::consts::PI - interior)
        };
    }
    let mut dual = vec![0.0; n];
    // Per segment: the stiffness of the recovery and that of the implicit step.
    let mut recovery = Vec::with_capacity(cut.segments.len());
    let mut implicit = Vec::with_capacity(cut.segments.len());
    for (i, segment) in cut.segments.iter().enumerate() {
        let t = tri.segment_triangle[i] as usize;
        let length = tri.segment_length[i];
        // The corner triangle (a, p, q): the segment is its edge from corner 1 to corner 2.
        let a = tri.faces[t][0] as usize;
        let [p, q] = [1, 2].map(|k| (tri.faces[t][k] as usize - nv) / 2);
        if rule == PenaltyRule::Gradient {
            let (cot_p, cot_q) = (tri.cot[t][1], tri.cot[t][2]);
            let foot = w_point[p] + cot_p / (cot_p + cot_q) * (w_point[q] - w_point[p]);
            let height = 2.0 * tri.area[t] / length;
            let across = (foot - weights[a]) / height;
            // Out of the lower region: the segment separates a from the two other vertices.
            let outward = if cut.vertex_region[a] == cut.segment_regions(mesh, i)[0] {
                across
            } else {
                -across
            };
            force[p] += 0.5 * length * outward;
            force[q] += 0.5 * length * outward;
        }
        dual[p] += 0.5 * length;
        dual[q] += 0.5 * length;
        let len2 = base.len2[segment.face as usize];
        let mean_squared_edge = (len2[0] + len2[1] + len2[2]) / 3.0;
        recovery.push(smoothing_length * smoothing_length * mean_squared_edge / length);
        implicit.push(0.5 * (w_point[p] + w_point[q]) / length);
    }
    let pairs: Vec<[u32; 2]> = cut.segments.iter().map(|segment| segment.points).collect();
    let (length_speed, first) = solve_on_cut(&pairs, &dual, &recovery, &force, par)?;
    let explicit: Vec<f64> = (0..n).map(|p| speed[p] + length_speed[p]).collect();
    let (tau, step_speed) = rule_step(
        step_size,
        statistics,
        &explicit,
        edges,
        &dual,
        roughening_rate,
    );
    for k in &mut implicit {
        *k *= tau;
    }
    let rhs: Vec<f64> = (0..n).map(|p| dual[p] * explicit[p]).collect();
    let (rate, second) = solve_on_cut(&pairs, &dual, &implicit, &rhs, par)?;
    Ok(RateOnCut {
        force,
        length_speed,
        rate,
        tau,
        step_speed,
        residual: first.max(second),
    })
}

/// Solves `(diag(ℓ) + S) x = b`, with `S` the stiffness matrix of the polylines that join the
/// nodes `pairs[σ]` with the stiffness `stiffness[σ]`. Returns `x` and the largest entry of the
/// residual relative to the largest entry of `b`.
pub(crate) fn solve_on_cut(
    pairs: &[[u32; 2]],
    dual: &[f64],
    stiffness: &[f64],
    rhs: &[f64],
    par: Parallelism<'_>,
) -> Result<(Vec<f64>, f64), CutError> {
    let n = dual.len();
    if stiffness.iter().all(|k| *k == 0.0) {
        return Ok((rhs.iter().zip(dual).map(|(b, l)| b / l).collect(), 0.0));
    }
    let mut triplets: Vec<(u32, u32, f64)> = Vec::with_capacity(n + 3 * stiffness.len());
    for (p, &l) in dual.iter().enumerate() {
        triplets.push((p as u32, p as u32, l));
    }
    for (&[p, q], &k) in pairs.iter().zip(stiffness) {
        triplets.push((p, p, k));
        triplets.push((q, q, k));
        triplets.push((p.max(q), p.min(q), -k));
    }
    let matrix = SparseMatrix::from_triplets(n, n, &triplets)?;
    let factorization = Symbolic::analyze(&matrix, Method::Llt)?.factor(&matrix, par)?;
    let mut x = rhs.to_vec();
    factorization.solve_in_place(&mut x, par)?;
    let mut residual: Vec<f64> = (0..n).map(|p| dual[p] * x[p] - rhs[p]).collect();
    for (&[p, q], &k) in pairs.iter().zip(stiffness) {
        let flow = k * (x[p as usize] - x[q as usize]);
        residual[p as usize] += flow;
        residual[q as usize] -= flow;
    }
    let worst = residual.iter().fold(0.0f64, |m, r| m.max(r.abs()));
    let largest = rhs.iter().fold(0.0f64, |m, b| m.max(b.abs()));
    Ok((
        x,
        if largest > 0.0 {
            worst / largest
        } else {
            worst
        },
    ))
}

/// The length weight `α_L = ω ℓ / (2π)` at which ripples of the cut longer than `wavelength`
/// (`ℓ`, a length on the unit-area surface) grow and shorter ones decay, for the roughening
/// rate `ω` (SPEC §1, D18): a ripple of wavenumber `k` grows at the rate `ω k − α_L k²`.
/// On the sphere the equator is stable from a third of its length on (`ℓ = 2πR/3`).
///
/// [formula: cut-roughening-rate]
pub fn length_weight_for_ripples(roughening_rate: f64, wavelength: f64) -> f64 {
    roughening_rate * wavelength / (2.0 * std::f64::consts::PI)
}

/// The local roughening rate `ω_c = 2 Σ_r (G_r − β_r − w_r b(0))` at every cut point (SPEC §1):
/// see [`StepReport::roughening_rate`].
///
/// [formula: cut-roughening-rate]
pub(crate) fn roughening_rates(e: &Evaluation, problem: &CutProblem) -> Vec<f64> {
    let offset = |region: usize| {
        let weights = problem.regions[region];
        let at_zero = match &problem.energy {
            Energy::Pointwise(b) => b.value(0.0),
            _ => 0.0,
        };
        weights.area + weights.distortion * at_zero
    };
    let tri = &e.triangles;
    (0..e.cut.points.len())
        .map(|p| {
            2.0 * tri
                .point_vertices(p)
                .iter()
                .map(|&v| e.gradient[v] - offset(tri.vertex_region[v] as usize))
                .sum::<f64>()
        })
        .collect()
}

/// The roughening rate `ω` of a state: the [`ROUGHENING_QUANTILE`] quantile of the local rates,
/// weighted by the length of the cut at every cut point.
///
/// [formula: cut-roughening-rate]
pub(crate) fn roughening_rate(e: &Evaluation, problem: &CutProblem) -> f64 {
    weighted_quantile(
        &roughening_rates(e, problem),
        &cut_dual_lengths(&e.cut, &e.triangles),
        ROUGHENING_QUANTILE,
    )
}

/// The explicit step `φ̃ = φ + τ σ` (SPEC §3.9).
///
/// [formula: cut-explicit-step]
#[inline]
pub fn explicit_step(phi: f64, tau: f64, sigma: f64) -> f64 {
    phi + tau * sigma
}

/// The symmetry projection (SPEC §3.11): every vertex takes the arithmetic mean of `phi` over
/// its class (the sum in increasing vertex index, divided by the size of the class).
///
/// # Panics
/// If `classes` and `phi` differ in length.
///
/// [formula: cut-symmetry-projection]
pub fn symmetry_projection(classes: &[u32], phi: &mut [f64]) {
    assert_eq!(classes.len(), phi.len());
    let mut sums: HashMap<u32, (f64, usize)> = HashMap::new();
    for (class, value) in classes.iter().zip(phi.iter()) {
        let entry = sums.entry(*class).or_insert((0.0, 0));
        entry.0 += value;
        entry.1 += 1;
    }
    for (class, value) in classes.iter().zip(phi.iter_mut()) {
        let (sum, count) = sums[class];
        *value = sum / count as f64;
    }
}

/// For every vertex, whether it is an end of an edge whose ends lie in different regions of
/// `phi`.
fn cut_adjacent(mesh: &TriMesh, phi: &[f64]) -> Vec<bool> {
    let mut adjacent = vec![false; phi.len()];
    for &[a, b] in mesh.edges() {
        if region_of(phi[a as usize]) != region_of(phi[b as usize]) {
            adjacent[a as usize] = true;
            adjacent[b as usize] = true;
        }
    }
    adjacent
}

/// The redistancing of SPEC §3.12 on the geometry `geometry`: the cut-adjacent vertices of `phi`
/// keep their values; every other vertex gets `sign · T` with `T` the fast-marching distance
/// from them within its own sign class; a vertex no front reaches keeps its value. Returns the
/// new function, the number of unreached vertices, and the cut-adjacent vertices.
///
/// [formula: cut-redistance]
pub(crate) fn redistance(
    mesh: &TriMesh,
    geometry: &Geometry,
    phi: &[f64],
) -> Result<(Vec<f64>, usize, Vec<bool>), CutError> {
    let region: Vec<u32> = phi.iter().map(|x| region_of(*x)).collect();
    let adjacent = cut_adjacent(mesh, phi);
    let seeds: Vec<(u32, f64)> = (0..phi.len())
        .filter(|&v| adjacent[v])
        .map(|v| (v as u32, phi[v].abs()))
        .collect();
    if seeds.is_empty() {
        return Ok((phi.to_vec(), 0, adjacent));
    }
    let distance = upwind_distance(mesh, geometry, &seeds, Some(&region))?;
    let mut unreached = 0;
    let out = (0..phi.len())
        .map(|v| {
            if adjacent[v] {
                phi[v]
            } else if distance.values[v].is_finite() {
                if region[v] == 0 {
                    -distance.values[v]
                } else {
                    distance.values[v]
                }
            } else {
                unreached += 1;
                phi[v]
            }
        })
        .collect();
    Ok((out, unreached, adjacent))
}

/// The vertices of the components of the regions that contain no cut-adjacent vertex: the
/// components of the graph of uncrossed edges without an end of a crossed edge.
fn components_without_cut(mesh: &TriMesh, region: &[u32], cut: &crate::cut::Cut) -> Vec<bool> {
    let nv = region.len();
    let mut parent: Vec<u32> = (0..nv as u32).collect();
    fn find(parent: &mut [u32], mut x: u32) -> u32 {
        while parent[x as usize] != x {
            parent[x as usize] = parent[parent[x as usize] as usize];
            x = parent[x as usize];
        }
        x
    }
    for &[a, b] in mesh.edges() {
        if region[a as usize] == region[b as usize] {
            let (ra, rb) = (find(&mut parent, a), find(&mut parent, b));
            parent[ra as usize] = rb;
        }
    }
    let mut has_cut = vec![false; nv];
    for point in &cut.points {
        for v in mesh.edge_vertices(point.edge as usize) {
            let root = find(&mut parent, v as u32);
            has_cut[root as usize] = true;
        }
    }
    (0..nv as u32)
        .map(|v| !has_cut[find(&mut parent, v) as usize])
        .collect()
}

/// `|∇φ|` of the linear interpolant in every face with two regions (its statistics, weighted by
/// the lengths of the segments), and for every cut point the mean of `|∇φ|` over the connected
/// piece of the cut it lies on. In a face whose lone vertex `a` is cut off by the segment `pq`,
/// the interpolant vanishes on the line `pq`, so
/// `|∇φ| = |φ_a| / dist(a, pq) = |φ_a| · |pq| / (2 · area(a, p, q))`. The mean over a piece is
/// the harmonic one, `Σ |pq| / Σ (|pq| / |∇φ|)`: a rate `r · ḡ` of the function moves the
/// segment by `r ḡ / |∇φ|` per unit time, and the mean of that over the piece is `r`.
///
/// A face with `|∇φ| = 0` is left out of the mean (it is in the statistics), and a piece
/// without any face that counts takes the mean 1.
///
/// `apex(σ, a)` is `|φ_a|` for segment `σ`. On a network (SPEC §8.5) it is
/// `½ |φ_i(a) − φ_j(a)|` for the two regions `i, j` of the segment, faces with three regions
/// are left out, and a piece is a connected piece of the network, through triple points.
///
/// [formula: cut-gradient-norm]
pub(crate) fn gradient_norms(
    cut: &Cut,
    tri: &CutTriangles,
    apex: impl Fn(usize, usize) -> f64,
) -> (GradientNorm, Vec<f64>) {
    let n_points = cut.points.len();
    // Connected pieces of the cut: union-find over the segments.
    let mut parent: Vec<u32> = (0..n_points as u32).collect();
    fn find(parent: &mut [u32], mut x: u32) -> u32 {
        while parent[x as usize] != x {
            parent[x as usize] = parent[parent[x as usize] as usize];
            x = parent[x as usize];
        }
        x
    }
    let join = |parent: &mut [u32], p: u32, q: u32| {
        let (a, b) = (find(parent, p), find(parent, q));
        if a != b {
            parent[a.max(b) as usize] = a.min(b);
        }
    };
    for segment in &cut.segments {
        if (segment.points[1] as usize) < n_points {
            join(&mut parent, segment.points[0], segment.points[1]);
        }
    }
    for triple in &cut.triple_points {
        join(&mut parent, triple.points[0], triple.points[1]);
        join(&mut parent, triple.points[0], triple.points[2]);
    }
    let (mut min, mut max, mut sum, mut total) = (f64::INFINITY, 0.0f64, 0.0, 0.0);
    let mut length_of = vec![0.0; n_points];
    let mut time_of = vec![0.0; n_points];
    for (i, segment) in cut.segments.iter().enumerate() {
        if segment.points[1] as usize >= n_points {
            continue;
        }
        let t = tri.segment_triangle[i] as usize;
        let a = tri.faces[t][0] as usize;
        let length = tri.segment_length[i];
        let g = apex(i, a) * length / (2.0 * tri.area[t]);
        min = min.min(g);
        max = max.max(g);
        sum += g * length;
        total += length;
        if g > 0.0 {
            // A face in which the function has no slope (it vanishes at the lone vertex) is
            // left out of the mean: its term would be infinite.
            let root = find(&mut parent, segment.points[0]) as usize;
            length_of[root] += length;
            time_of[root] += length / g;
        }
    }
    let mean_of_piece: Vec<f64> = (0..n_points as u32)
        .map(|p| {
            let root = find(&mut parent, p) as usize;
            let mean = length_of[root] / time_of[root];
            if mean.is_finite() && mean > 0.0 {
                mean
            } else {
                1.0
            }
        })
        .collect();
    (
        GradientNorm {
            min,
            max,
            mean: sum / total,
        },
        mean_of_piece,
    )
}
