//! The flow (SPEC §3.8–§3.12, §4): one step moves the level-set function by the harmonic
//! extension of the speed on the cut, applies the implicit length step, and makes the function
//! a signed distance again.

use std::collections::HashMap;
use std::time::Instant;

use ddg_core::distance::upwind_distance;
use ddg_core::sparse::{Method, Parallelism, Symbolic};
use ddg_core::{Geometry, SparseMatrix, TriMesh};

use crate::cut::region_of;
use crate::error::CutError;
use crate::evaluate::{EnergyTerms, Evaluation, Evaluator};
use crate::problem::{CutProblem, EvalOptions, UnitAreaScale, check_length};
use crate::system::LaplaceSystem;

/// How the time step `τ` of a flow step is chosen (SPEC §4).
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum StepSize {
    /// The same `τ > 0` in every step, as in the paper (Sec. 6).
    Fixed(f64),
    /// `τ = min(tau_max, max_motion · h / max_c |s_c|)` with `h` the shortest crossed edge: the
    /// distortion term moves the cut by at most the fraction `max_motion` of that edge per step
    /// (the paper's remark in Sec. 7.4). If the speed vanishes, `τ = tau_max`.
    Capped {
        /// The largest step, `> 0`.
        tau_max: f64,
        /// The largest motion by the distortion term, as a fraction of the shortest crossed
        /// edge, `> 0`.
        max_motion: f64,
    },
}

impl StepSize {
    /// The step for the speeds `speed` on the cut and the shortest crossed edge `h`.
    ///
    /// [formula: cut-step-cap]
    pub fn step(&self, speed: &[f64], shortest_crossed_edge: f64) -> f64 {
        match *self {
            StepSize::Fixed(tau) => tau,
            StepSize::Capped {
                tau_max,
                max_motion,
            } => {
                let fastest = speed.iter().fold(0.0f64, |m, s| m.max(s.abs()));
                if fastest > 0.0 {
                    tau_max.min(max_motion * shortest_crossed_edge / fastest)
                } else {
                    tau_max
                }
            }
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
    #[default]
    Unit,
    /// `φ_t = s · |∇φ|` on the cut, the level-set equation for a function that is not a signed
    /// distance: the cut moves with the normal speed `−s` whatever `|∇φ|` is.
    GradientNorm,
}

/// The matrix of the harmonic extension of the rate on the cut (SPEC §3.8, D15).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum Extension {
    /// The cotan Laplacian of the cut mesh with negative corner cotangents replaced by zero: an
    /// M-matrix, so the extension stays within the range of the rates on the cut. The default.
    /// It needs a numeric factorization of its own in every step.
    #[default]
    Monotone,
    /// The cotan Laplacian of the cut mesh, as in the paper (Sec. 5.3.1). On a closed mesh this
    /// is the matrix the evaluation factorized, so the extension costs two triangular solves.
    /// On the obtuse triangles that cutting produces its weights are negative, and the
    /// extension at the vertices next to the cut can lie far outside the range of the rates
    /// (measured: a tenth and three times the neighbouring values on the two sides of one
    /// edge), which drives the level-set function away from a distance function.
    Cotan,
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
    /// The matrix of the harmonic extension. Default: the monotone one.
    pub extension: Extension,
    /// The redistancing.
    pub redistance: Redistance,
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
        }
    }

    fn validate(&self, n_vertices: usize) -> Result<(), CutError> {
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
}

/// Why [`Flow::run`] stopped.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum StopReason {
    /// The limit on the number of steps was reached.
    MaxSteps,
    /// The cut motion of the last step was within the tolerance.
    Converged,
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
    /// The largest `|λ_c|`, the implicit length speed on the cut ([`LengthStep::OnCut`]; 0 for
    /// [`LengthStep::OnMesh`], where the length term does not pass through the cut points).
    pub max_length_speed: f64,
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
    /// `σ`, the extension of the rate on the cut (stage 4).
    pub extension: Vec<f64>,
    /// `φ̃ = φ + τ σ` (stage 5).
    pub explicit: Vec<f64>,
    /// `φ′` after the length step on the mesh (stage 6; equal to `explicit` if there is none).
    pub length: Vec<f64>,
    /// After the first symmetry projection (equal to `length` without symmetry classes).
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
    evaluator: Evaluator,
    /// The geometry of the unit-area surface, for the redistancing.
    geometry: Geometry,
    /// The level-set function on the unit-area surface.
    phi: Vec<f64>,
    time: f64,
    n_steps: usize,
    /// The Laplacian of the cut mesh on all vertices, for the extension. On a closed mesh it
    /// shares the pattern (and the symbolic analysis) of the evaluator's system, and with
    /// [`Extension::Cotan`] that system's factorization is used instead.
    extension: LaplaceSystem,
    /// Whether the mesh has no boundary.
    closed: bool,
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
        let mut flow = Flow {
            evaluator,
            geometry: scaled,
            phi: Vec::new(),
            time: 0.0,
            n_steps: 0,
            extension,
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
        let (gradient_norm, gradient_at_point) = gradient_norms(&e, &self.phi);
        report.gradient_norm = Some(gradient_norm);
        let shortest = e
            .cut
            .points
            .iter()
            .map(|p| self.geometry.edge_lengths()[p.edge as usize])
            .fold(f64::INFINITY, f64::min);
        let tau = options.step_size.step(&e.speed, shortest);
        report.tau = tau;

        // The length term on the cut (SPEC §3.10): an implicit speed λ next to s.
        let start = Instant::now();
        let weights = length_weights(problem, nv);
        let active = weights.iter().any(|w| *w != 0.0);
        let mut lambda = vec![0.0; e.speed.len()];
        if active && options.length_step == LengthStep::OnCut {
            let smoothing_length = match options.eval.flux_recovery {
                crate::problem::FluxRecovery::Smoothed { length } => length,
                crate::problem::FluxRecovery::Nodal => 0.0,
            };
            let (solution, residual) = length_speed_on_cut(
                self.evaluator.mesh(),
                self.evaluator.base(),
                &e,
                &weights,
                options.penalty_rule,
                smoothing_length,
                tau,
                par,
            )?;
            lambda = solution;
            report.length_residual = Some(residual);
            report.max_length_speed = lambda.iter().fold(0.0f64, |m, x| m.max(x.abs()));
        }
        report.timings.length = start.elapsed().as_secs_f64();

        // Stage 4: σ, the harmonic extension of the rate on the cut.
        let start = Instant::now();
        let rate: Vec<f64> = (0..e.speed.len())
            .map(|k| {
                let total = e.speed[k] + lambda[k];
                match options.speed_scaling {
                    SpeedScaling::Unit => total,
                    SpeedScaling::GradientNorm => total * gradient_at_point[k],
                }
            })
            .collect();
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
            stages.extension = sigma;
            stages.explicit = phi.clone();
        }

        // Stage 6: the implicit step on the mesh (the paper's), then the symmetry projection.
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
        if let Some(classes) = &options.symmetry {
            symmetry_projection(classes, &mut phi);
        }
        if let Some(stages) = stages.as_deref_mut() {
            stages.symmetric = phi.clone();
        }

        // Stage 7: redistancing, then the symmetry projection again.
        let adjacent_before = e.cut.cut_adjacent_vertices(self.evaluator.mesh());
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
        };
        let final_energy = self.evaluate(problem, &options.eval)?.energy;
        Ok(RunReport {
            steps,
            stop: reason,
            final_energy,
        })
    }

    /// The harmonic extension of the rate `rate` on the cut points to every vertex of the mesh
    /// (SPEC §3.8): in each region `L_OO σ_O = −L_OC s_C`, with all vertices of the mesh unknown
    /// (zero flux through `∂M`) and `σ = 0` in a component without cut vertices. `L` is the
    /// cotan Laplacian of the cut mesh, or its monotone variant. Returns `σ` and the residual of
    /// the solve.
    ///
    /// [formula: cut-harmonic-extension]
    fn extend(
        &mut self,
        e: &Evaluation,
        rate: &[f64],
        matrix: Extension,
        par: Parallelism<'_>,
    ) -> Result<(Vec<f64>, f64), CutError> {
        let tri = &e.triangles;
        let nv = tri.n_original;
        let n = tri.n_vertices();
        let monotone = matrix == Extension::Monotone;
        let apply = |x: &[f64], y: &mut [f64]| {
            if monotone {
                tri.apply_monotone_laplacian(x, y)
            } else {
                tri.apply_laplacian(None, x, y)
            }
        };
        let mut sigma = vec![0.0; n];
        for (k, &s) in rate.iter().enumerate() {
            for v in tri.point_vertices(k) {
                sigma[v] = s;
            }
        }
        // (L x)_O for x = s on the cut vertices and 0 elsewhere is L_OC s_C.
        let mut product = vec![0.0; n];
        apply(&sigma, &mut product)?;
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
        if self.closed && !monotone {
            // The matrix is the one the evaluation factorized, with identity rows for the
            // closed components without a cut.
            self.evaluator.system().solve(rhs, &mut sigma, par)?;
        } else {
            self.extension.assemble(mesh, tri, &pinned, monotone);
            self.extension.factor(par)?;
            self.extension.solve(rhs, &mut sigma, par)?;
        }
        apply(&sigma, &mut product)?;
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
    fn length_step(
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
fn length_weights(problem: &CutProblem, n_vertices: usize) -> Vec<f64> {
    match &problem.penalty {
        Some(eta) => eta.iter().map(|e| problem.length_weight + e).collect(),
        None => vec![problem.length_weight; n_vertices],
    }
}

/// The implicit length speed `λ` at every cut point ([`LengthStep::OnCut`], SPEC §3.10), and
/// the residual of its 1D system.
///
/// The cut is a polyline; the first variation of its weighted length `∫_γ w ds` under a normal
/// displacement `δ` (positive out of the region with the lower id) is
/// `∫_γ (∂w/∂n + w κ) δ ds` (SPEC §1), and on a polyline the curvature sits in the vertices:
///
/// ```text
/// Θ_p = w_p · 2 sin(θ_p / 2)  +  Σ_{σ = pq at p} ½ |pq| (∂w/∂n)_σ
/// ```
///
/// * `θ_p = π − β_p` is the turning angle of the cut at `p`, with `β_p` the interior angle of
///   the lower region at `p` (the sum of the corner angles of its triangles there; the two
///   faces at a crossed edge unfold into a plane). At a cut point on `∂M` the cut is reflected
///   in the boundary, which is the zero-flux condition of the paper's step:
///   `θ_p = π − 2 β_p`, and half of the force acts.
/// * `(∂w/∂n)_σ` is the derivative of the linear `w` of the face of `σ` across the segment, out
///   of the lower region: with `a` the lone vertex of the face, `h_a` its height over `pq` and
///   `f` the foot of that height, `(w_f − w_a)/h_a`, negated if `a` is in the higher region.
///   The paper's penalty rule drops this term.
///
/// `λ = Θ/ℓ` taken explicitly would restrict the step to the square of the shortest segment.
/// It is made implicit along the cut, with the lumped mass `ℓ_p = ½ Σ |pq|` and the stiffness
/// of the polyline (moving `p` out by `δ` relative to its neighbours turns the cut by
/// `Σ (δ_p − δ_q)/|pq|` more):
///
/// ```text
/// ℓ_p λ_p + τ Σ_{σ = pq at p} (w̄_σ / |pq|) (λ_p − λ_q) = Θ_p,      w̄_σ = ½ (w_p + w_q)
/// ```
///
/// The matrix is symmetric positive definite, an M-matrix with row sums `ℓ_p`; `Σ ℓ_p λ_p =
/// Σ Θ_p`, and a vertex is never moved beyond the chord of its neighbours, however large `τ`.
///
/// [formula: cut-length-step-on-cut]
fn length_speed_on_cut(
    mesh: &TriMesh,
    base: &crate::cut_mesh::BaseGeometry,
    e: &Evaluation,
    weights: &[f64],
    rule: PenaltyRule,
    smoothing_length: f64,
    tau: f64,
    par: Parallelism<'_>,
) -> Result<(Vec<f64>, f64), CutError> {
    let (cut, tri) = (&e.cut, &e.triangles);
    let nv = tri.n_original;
    let n = cut.points.len();
    let w_point = cut.interpolate(mesh, weights);
    // Interior angles at the copies of the cut points.
    let mut angle = vec![0.0; 2 * n];
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
            w_point[p] * (0.5 * std::f64::consts::PI - interior).sin()
        } else {
            w_point[p] * 2.0 * (0.5 * (std::f64::consts::PI - interior)).sin()
        };
    }
    let mut dual = vec![0.0; n];
    let mut triplets: Vec<(u32, u32, f64)> = Vec::with_capacity(4 * n);
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
            let other = tri.vertex_region[tri.point_vertices(p)[0]]
                .max(tri.vertex_region[tri.point_vertices(p)[1]]);
            let outward = if cut.vertex_region[a] < other {
                across
            } else {
                -across
            };
            force[p] += 0.5 * length * outward;
            force[q] += 0.5 * length * outward;
        }
        dual[p] += 0.5 * length;
        dual[q] += 0.5 * length;
        let stiffness = segment_stiffness(
            base,
            segment.face,
            smoothing_length,
            tau,
            w_point[p],
            w_point[q],
        ) / length;
        let (p, q) = (p as u32, q as u32);
        triplets.push((p, p, stiffness));
        triplets.push((q, q, stiffness));
        triplets.push((p.max(q), p.min(q), -stiffness));
    }
    for (p, &l) in dual.iter().enumerate() {
        triplets.push((p as u32, p as u32, l));
    }
    let matrix = SparseMatrix::from_triplets(n, n, &triplets)?;
    let factorization = Symbolic::analyze(&matrix, Method::Llt)?.factor(&matrix, par)?;
    let mut lambda = force.clone();
    factorization.solve_in_place(&mut lambda, par)?;
    // Residual of the system, assembled again from the segments.
    let mut residual: Vec<f64> = (0..n).map(|p| dual[p] * lambda[p] - force[p]).collect();
    for (i, segment) in cut.segments.iter().enumerate() {
        let [p, q] = segment.points.map(|x| x as usize);
        let flow = segment_stiffness(
            base,
            segment.face,
            smoothing_length,
            tau,
            w_point[p],
            w_point[q],
        ) / tri.segment_length[i]
            * (lambda[p] - lambda[q]);
        residual[p] += flow;
        residual[q] -= flow;
    }
    let worst = residual.iter().fold(0.0f64, |m, x| m.max(x.abs()));
    let largest = force.iter().fold(0.0f64, |m, x| m.max(x.abs()));
    Ok((
        lambda,
        if largest > 0.0 {
            worst / largest
        } else {
            worst
        },
    ))
}

/// `(c h_f)² + τ w̄`: the stiffness of a cut segment in the system of the length speed, times
/// its length. `c` is the smoothing length of the flux recovery (SPEC §3.5), `h_f²` the mean
/// squared edge length of the face of the segment, `w̄` the mean weight of its two ends.
#[inline]
fn segment_stiffness(
    base: &crate::cut_mesh::BaseGeometry,
    face: u32,
    smoothing_length: f64,
    tau: f64,
    w_p: f64,
    w_q: f64,
) -> f64 {
    let len2 = base.len2[face as usize];
    let mean_squared_edge = (len2[0] + len2[1] + len2[2]) / 3.0;
    smoothing_length * smoothing_length * mean_squared_edge + tau * 0.5 * (w_p + w_q)
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

/// `|∇φ|` of the linear interpolant of `phi` in every crossed face (its statistics, weighted by
/// the lengths of the segments) and at every cut point (the length-weighted mean over its
/// segments). In a face whose lone vertex `a` is cut off by the segment `pq`, the interpolant
/// vanishes on the line `pq`, so `|∇φ| = |φ_a| / dist(a, pq) = |φ_a| · |pq| / (2 · area(a, p, q))`.
fn gradient_norms(e: &Evaluation, phi: &[f64]) -> (GradientNorm, Vec<f64>) {
    let tri = &e.triangles;
    let n_points = e.cut.points.len();
    let (mut min, mut max, mut sum, mut total) = (f64::INFINITY, 0.0f64, 0.0, 0.0);
    let mut at_point = vec![0.0; n_points];
    let mut weight = vec![0.0; n_points];
    for (i, segment) in e.cut.segments.iter().enumerate() {
        let t = tri.segment_triangle[i] as usize;
        let a = tri.faces[t][0] as usize;
        let length = tri.segment_length[i];
        let g = phi[a].abs() * length / (2.0 * tri.area[t]);
        min = min.min(g);
        max = max.max(g);
        sum += g * length;
        total += length;
        for p in segment.points {
            at_point[p as usize] += g * length;
            weight[p as usize] += length;
        }
    }
    for (g, w) in at_point.iter_mut().zip(&weight) {
        *g = if *w > 0.0 { *g / w } else { 1.0 };
    }
    (
        GradientNorm {
            min,
            max,
            mean: sum / total,
        },
        at_point,
    )
}
