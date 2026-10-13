//! The flow for `N` regions (SPEC §8): one level-set function per region, a network of cuts
//! with triple points.
//!
//! Everything that is not specific to networks is the code of the two-region flow: the
//! evaluation ([`Evaluator::evaluate_regions`](crate::Evaluator::evaluate_regions)), the rate on
//! the cut, the extensions, the implicit step and the fast marching. With two functions
//! `φ, −φ` every stage reproduces [`Flow`] (tested).

use std::time::Instant;

use ddg_core::distance::upwind_distance;
use ddg_core::{Geometry, TriMesh};

use crate::cut::regions_of;
use crate::error::CutError;
use crate::evaluate::{Evaluation, Evaluator};
use crate::flow::{
    Extension, Flow, FlowOptions, LengthStep, Redistance, RunReport, SpeedScaling, StepReport,
    StepTimings, StopReason, StopRule, explicit_step, gradient_norms, length_weights, rate_on_cut,
    roughening_rate, symmetry_projection,
};
use crate::problem::{CutProblem, EvalOptions, UnitAreaScale};

/// The cut flow for `N ≥ 2` regions (SPEC §8).
///
/// The state is one level-set function per region, `φ_k`, negative inside region `k`; a vertex
/// belongs to the region whose function is smallest there (ties: the lowest index). The
/// problem needs one [`RegionWeights`](crate::RegionWeights) per region
/// ([`CutProblem::for_regions`]); options, reports and the stop rule are those of [`Flow`].
///
/// A two-region function `φ` is the state `[φ, −φ]`, and on it this flow takes the same steps
/// as [`Flow`], to rounding.
#[derive(Debug)]
pub struct RegionFlow {
    /// The two-region flow of the same mesh: it owns the evaluator and the matrices, and its
    /// own level-set function is not used.
    core: Flow,
    /// The level-set functions on the unit-area surface.
    phis: Vec<Vec<f64>>,
    time: f64,
    n_steps: usize,
}

impl RegionFlow {
    /// Sets up the flow on `mesh` with the intrinsic `geometry`, starting from the level-set
    /// functions `phis` (`phis[k]` has one value per vertex, in the length unit of `geometry`,
    /// and is negative inside region `k`). The functions are taken as they are: see
    /// [`crate::init`] for ways to make them.
    ///
    /// # Errors
    /// [`CutError::UnsupportedRegionCount`] for fewer than two functions,
    /// [`CutError::WrongLength`], [`CutError::NonFiniteLevelSet`], and what
    /// [`Evaluator::new`] returns.
    pub fn new(
        mesh: &TriMesh,
        geometry: &Geometry,
        phis: &[Vec<f64>],
    ) -> Result<RegionFlow, CutError> {
        regions_of(mesh, phis)?;
        let core = Flow::new(mesh, geometry, &phis[0])?;
        let mut flow = RegionFlow {
            core,
            phis: Vec::new(),
            time: 0.0,
            n_steps: 0,
        };
        flow.set_level_sets(phis)?;
        Ok(flow)
    }

    /// The mesh.
    #[inline]
    pub fn mesh(&self) -> &TriMesh {
        self.core.mesh()
    }

    /// The scaling between the caller's unit and the unit-area surface.
    #[inline]
    pub fn scale(&self) -> UnitAreaScale {
        self.core.scale()
    }

    /// The evaluator (for its cached data and options that are not part of a step).
    #[inline]
    pub fn evaluator(&mut self) -> &mut Evaluator {
        self.core.evaluator()
    }

    /// The number of regions `N`.
    #[inline]
    pub fn n_regions(&self) -> usize {
        self.phis.len()
    }

    /// The level-set functions in the caller's length unit: the state, from which the flow can
    /// be resumed ([`RegionFlow::new`]).
    pub fn level_sets(&self) -> Vec<Vec<f64>> {
        let scale = self.scale();
        self.phis
            .iter()
            .map(|phi| {
                phi.iter()
                    .map(|x| scale.length_from_unit_area(*x))
                    .collect()
            })
            .collect()
    }

    /// The level-set functions on the unit-area surface.
    #[inline]
    pub fn level_sets_unit_area(&self) -> &[Vec<f64>] {
        &self.phis
    }

    /// The region of every vertex.
    pub fn regions(&self) -> Vec<u32> {
        regions_of(self.mesh(), &self.phis).expect("the state is valid")
    }

    /// Replaces the state (in the caller's length unit) and resets time and step count. The
    /// number of regions may change.
    ///
    /// # Errors
    /// As [`RegionFlow::new`] for the functions.
    pub fn set_level_sets(&mut self, phis: &[Vec<f64>]) -> Result<(), CutError> {
        regions_of(self.mesh(), phis)?;
        let scale = self.scale();
        self.phis = phis
            .iter()
            .map(|phi| phi.iter().map(|x| scale.length_to_unit_area(*x)).collect())
            .collect();
        self.time = 0.0;
        self.n_steps = 0;
        Ok(())
    }

    /// The flow time: the sum of the steps taken since the state was set.
    #[inline]
    pub fn time(&self) -> f64 {
        self.time
    }

    /// The number of steps taken since the state was set.
    #[inline]
    pub fn n_steps(&self) -> usize {
        self.n_steps
    }

    /// Evaluates the current state (SPEC §8.1–§8.3) without changing it.
    ///
    /// # Errors
    /// As [`Evaluator::evaluate_regions`].
    pub fn evaluate(
        &mut self,
        problem: &CutProblem,
        options: &EvalOptions<'_>,
    ) -> Result<Evaluation, CutError> {
        self.core
            .evaluator
            .evaluate_regions(&self.phis, problem, options)
    }

    /// The roughening rate `ω` of the current state ([`StepReport::roughening_rate`]).
    ///
    /// # Errors
    /// As [`Evaluator::evaluate_regions`].
    pub fn roughening_rate(
        &mut self,
        problem: &CutProblem,
        options: &EvalOptions<'_>,
    ) -> Result<f64, CutError> {
        let e = self.evaluate(problem, options)?;
        Ok(roughening_rate(&e, problem))
    }

    /// Makes every function a signed distance to the boundary of its region by the
    /// redistancing of SPEC §8.4 (the projection at the cut-adjacent vertices, then fast
    /// marching per function). Returns the number of vertices that a front did not reach.
    ///
    /// # Errors
    /// [`CutError::Distance`] (cannot happen for finite functions).
    pub fn redistance(&mut self) -> Result<usize, CutError> {
        let (phis, unreached, _) =
            redistance_regions(self.core.evaluator.mesh(), &self.core.geometry, &self.phis)?;
        self.phis = phis;
        Ok(unreached)
    }

    /// One step of the flow (SPEC §8.4, §8.5) with the given problem and options, which may
    /// differ from step to step. The stages are those of [`Flow::step`]; the report has the
    /// same meaning, with the cut motion taken over all functions.
    ///
    /// # Errors
    /// As [`Flow::step`]; the state is unchanged if a step fails.
    pub fn step(
        &mut self,
        problem: &CutProblem,
        options: &FlowOptions<'_>,
    ) -> Result<StepReport, CutError> {
        let nv = self.mesh().n_vertices();
        let n_regions = self.phis.len();
        options.validate(nv)?;
        let par = options.eval.parallelism;

        let start = Instant::now();
        let e = self
            .core
            .evaluator
            .evaluate_regions(&self.phis, problem, &options.eval)?;
        let timings = StepTimings {
            evaluate: start.elapsed().as_secs_f64(),
            ..StepTimings::default()
        };
        let max_speed = e.speed.iter().fold(0.0f64, |m, s| m.max(s.abs()));
        let mut report = StepReport {
            time: self.time,
            tau: 0.0,
            energy: e.energy.clone(),
            max_speed,
            max_length_speed: 0.0,
            roughening_rate: 0.0,
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
        let mesh = self.core.evaluator.mesh();
        let (cut, tri) = (&e.cut, &e.triangles);
        let phis = &self.phis;
        // Half the difference of the two functions of a segment plays the role of φ.
        let (gradient_norm, mean_gradient) = gradient_norms(cut, tri, |segment, a| {
            let [i, j] = cut.segment_regions(mesh, segment);
            0.5 * (phis[i as usize][a] - phis[j as usize][a]).abs()
        });
        report.gradient_norm = Some(gradient_norm);
        report.roughening_rate = roughening_rate(&e, problem);
        let shortest = cut
            .points
            .iter()
            .map(|p| self.core.geometry.edge_lengths()[p.edge as usize])
            .fold(f64::INFINITY, f64::min);

        // The length force on the cut and the rate that is implicit in it (SPEC §3.10, §8.5).
        let start = Instant::now();
        let weights = length_weights(problem, nv);
        let active = weights.iter().any(|w| *w != 0.0);
        let mut rate = e.speed.clone();
        let mut tau = options
            .step_size
            .step(&e.speed, shortest, report.roughening_rate);
        if active && options.length_step == LengthStep::OnCut {
            let smoothing_length = match options.eval.flux_recovery {
                crate::problem::FluxRecovery::Smoothed { length } => length,
                crate::problem::FluxRecovery::Nodal => 0.0,
            };
            let on_cut = rate_on_cut(
                mesh,
                self.core.evaluator.base(),
                cut,
                tri,
                &e.speed,
                &weights,
                options.penalty_rule,
                smoothing_length,
                &options.step_size,
                shortest,
                report.roughening_rate,
                par,
            )?;
            tau = on_cut.tau;
            report.length_residual = Some(on_cut.residual);
            report.max_length_speed = on_cut
                .length_speed
                .iter()
                .fold(0.0f64, |m, x| m.max(x.abs()));
            rate = on_cut.rate;
        }
        report.tau = tau;
        report.timings.length = start.elapsed().as_secs_f64();

        // The extension: at a vertex of region k, the extension of the rate of φ_k.
        let start = Instant::now();
        if options.speed_scaling == SpeedScaling::MeanGradient {
            for (r, g) in rate.iter_mut().zip(&mean_gradient) {
                *r *= g;
            }
        }
        let (extension, residual) = match options.extension {
            Extension::Local => (self.core.extend_from_cut(cut, tri, &rate, true), 0.0),
            Extension::Harmonic => {
                let data = boundary_rates(mesh, &e, &rate);
                self.core.extend_harmonically_from(&e, data, par)?
            }
        };
        report.extension_residual = residual;
        report.timings.extension = start.elapsed().as_secs_f64();

        // The explicit step with the substitution of App. A.3: the function of the vertex's own
        // region takes the rate, every other function its negative.
        let region = &cut.vertex_region;
        let mut new: Vec<Vec<f64>> = (0..n_regions)
            .map(|k| {
                (0..nv)
                    .map(|v| {
                        let sigma = if region[v] as usize == k {
                            extension[v]
                        } else {
                            -extension[v]
                        };
                        explicit_step(self.phis[k][v], tau, sigma)
                    })
                    .collect()
            })
            .collect();

        // The implicit step on the mesh (the paper's): one matrix, N right-hand sides.
        let start = Instant::now();
        if active && options.length_step == LengthStep::OnMesh {
            let mut worst = 0.0f64;
            for phi in &mut new {
                let (residual, factorized) =
                    self.core
                        .length_step(&weights, options.penalty_rule, tau, phi, par)?;
                worst = worst.max(residual);
                report.length_factorized |= factorized;
            }
            report.length_residual = Some(worst);
        }
        report.timings.length += start.elapsed().as_secs_f64();
        let start = Instant::now();
        if let Some(classes) = &options.symmetry {
            for phi in &mut new {
                symmetry_projection(classes, phi);
            }
        }

        // Redistancing (SPEC §8.4), then the symmetry projection again.
        let mesh = self.core.evaluator.mesh();
        let adjacent_before = cut.cut_adjacent_vertices(mesh);
        let adjacent_after = match options.redistance {
            Redistance::KeepAdjacent => {
                let (phis, unreached, adjacent) =
                    redistance_regions(mesh, &self.core.geometry, &new)?;
                report.n_unreached = unreached;
                new = phis;
                adjacent
            }
            Redistance::Off => {
                let region = regions_of(mesh, &new)?;
                let mut adjacent = vec![false; nv];
                for &[a, b] in mesh.edges() {
                    if region[a as usize] != region[b as usize] {
                        adjacent[a as usize] = true;
                        adjacent[b as usize] = true;
                    }
                }
                adjacent
            }
        };
        if let Some(classes) = &options.symmetry {
            for phi in &mut new {
                symmetry_projection(classes, phi);
            }
        }
        report.timings.redistance = start.elapsed().as_secs_f64();
        report.cut_vanished = !adjacent_after.iter().any(|&a| a);
        report.cut_motion = (0..nv)
            .filter(|&v| adjacent_before[v] || adjacent_after[v])
            .map(|v| {
                (0..n_regions)
                    .map(|k| (new[k][v] - self.phis[k][v]).abs())
                    .fold(0.0, f64::max)
            })
            .fold(0.0, f64::max);

        self.phis = new;
        self.time += tau;
        self.n_steps += 1;
        Ok(report)
    }

    /// Steps until `stop` says so, with the same problem and options in every step, and returns
    /// the history ([`Flow::run`]).
    ///
    /// # Errors
    /// As [`RegionFlow::step`].
    pub fn run(
        &mut self,
        problem: &CutProblem,
        options: &FlowOptions<'_>,
        stop: &StopRule,
    ) -> Result<RunReport, CutError> {
        let mut steps: Vec<StepReport> = Vec::new();
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
}

/// The boundary data of the harmonic extension on a network (SPEC §8.4): at the copy of a cut
/// point in region `k` the rate of `φ_k` there (`rate` is that of the lower region, the higher
/// one has its negative); at the copy of a triple point in region `k` the linear interpolation,
/// in arclength, between the two cut points next to it on the boundary of region `k`. Zero at
/// the vertices of the input mesh.
///
/// [formula: cut-region-extension]
fn boundary_rates(mesh: &TriMesh, e: &Evaluation, rate: &[f64]) -> Vec<f64> {
    let (cut, tri) = (&e.cut, &e.triangles);
    let n_points = cut.points.len();
    let mut data = vec![0.0; tri.n_vertices()];
    for (p, &r) in rate.iter().enumerate() {
        let [lower, higher] = tri.point_vertices_by_region(p);
        data[lower] = r;
        data[higher] = -r;
    }
    // The three segments of every triple point, in the order of the edges of its face.
    let mut segments = vec![[0usize; 3]; cut.triple_points.len()];
    let mut filled = vec![0usize; cut.triple_points.len()];
    for (i, segment) in cut.segments.iter().enumerate() {
        if segment.points[1] as usize >= n_points {
            let m = segment.points[1] as usize - n_points;
            segments[m][filled[m]] = i;
            filled[m] += 1;
        }
    }
    for (m, triple) in cut.triple_points.iter().enumerate() {
        let copies = tri.triple_vertices(m);
        for c in 0..3 {
            let region = tri.vertex_region[copies[c]];
            let previous = (c + 2) % 3;
            let value = |edge: usize| {
                let p = triple.points[edge] as usize;
                if cut.point_regions(mesh, p)[0] == region {
                    rate[p]
                } else {
                    -rate[p]
                }
            };
            let (to_p, to_q) = (
                tri.segment_length[segments[m][c]],
                tri.segment_length[segments[m][previous]],
            );
            data[copies[c]] = (to_q * value(c) + to_p * value(previous)) / (to_p + to_q);
        }
    }
    data
}

/// The redistancing of SPEC §8.4 on the geometry `geometry`. Regions are those of `phis`.
/// `Adj_k` are the ends of the edges with exactly one end in region `k`.
///
/// 1. At every vertex of `Adj = ∪ Adj_k` the mean of the two smallest values is subtracted from
///    all of them (Losasso, Shinar, Selle, Fedkiw: *Multiple Interacting Liquids*, ACM TOG
///    25(3), 2006, Sec. 3.1): the smallest becomes `−d`, the second `+d`.
/// 2. For every `k`, fast marching from the seeds `Adj_k` with the values `|φ_k|`, within the
///    two classes "in region `k`" and "not in region `k`"; `φ_k = −T` inside, `+T` outside. A
///    vertex that no front reaches keeps its value; a function without seeds is left alone.
///
/// Returns the new functions, the number of vertices some front did not reach (over the
/// functions that have seeds), and the vertices of `Adj`.
///
/// [formula: cut-region-projection]
pub(crate) fn redistance_regions(
    mesh: &TriMesh,
    geometry: &Geometry,
    phis: &[Vec<f64>],
) -> Result<(Vec<Vec<f64>>, usize, Vec<bool>), CutError> {
    let nv = mesh.n_vertices();
    let n_regions = phis.len();
    let region = regions_of(mesh, phis)?;
    let mut adjacent = vec![false; nv];
    for &[a, b] in mesh.edges() {
        if region[a as usize] != region[b as usize] {
            adjacent[a as usize] = true;
            adjacent[b as usize] = true;
        }
    }
    // Step 1: the projection.
    let mut out: Vec<Vec<f64>> = phis.to_vec();
    for v in (0..nv).filter(|&v| adjacent[v]) {
        let (mut first, mut second) = (f64::INFINITY, f64::INFINITY);
        for phi in phis {
            let x = phi[v];
            if x < first {
                second = first;
                first = x;
            } else if x < second {
                second = x;
            }
        }
        let mean = 0.5 * (first + second);
        for phi in &mut out {
            phi[v] -= mean;
        }
    }
    // Step 2: fast marching per function.
    let mut unreached = vec![false; nv];
    for k in 0..n_regions {
        let mut seed = vec![false; nv];
        for &[a, b] in mesh.edges() {
            let (ra, rb) = (region[a as usize] as usize, region[b as usize] as usize);
            if ra != rb && (ra == k || rb == k) {
                seed[a as usize] = true;
                seed[b as usize] = true;
            }
        }
        let seeds: Vec<(u32, f64)> = (0..nv)
            .filter(|&v| seed[v])
            .map(|v| (v as u32, out[k][v].abs()))
            .collect();
        if seeds.is_empty() {
            continue;
        }
        let class: Vec<u32> = region.iter().map(|&r| u32::from(r as usize != k)).collect();
        let distance = upwind_distance(mesh, geometry, &seeds, Some(&class))?;
        for v in (0..nv).filter(|&v| !seed[v]) {
            if distance.values[v].is_finite() {
                out[k][v] = if region[v] as usize == k {
                    -distance.values[v]
                } else {
                    distance.values[v]
                };
            } else {
                unreached[v] = true;
            }
        }
    }
    Ok((out, unreached.iter().filter(|&&u| u).count(), adjacent))
}
