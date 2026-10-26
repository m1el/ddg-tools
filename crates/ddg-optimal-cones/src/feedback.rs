//! The feedback loop: lower the distortion weight where the solver's linearized log scale
//! factor is far from the one a flattening realizes.
//!
//! # Why
//!
//! The canonical problem measures distortion with the solution `u` of the **linearized**
//! Yamabe equation `L u = µ − Ω`. Where that linearization is poor, the solver optimizes a
//! distortion that no flattening has. The typical case is a needle-like spike (README,
//! "Needle-like spikes"): at the tip of a needle of height `h` over a base of radius `ρ` the
//! linearized `u` is of order `−h/ρ` (`−1.5·10⁴` on the dataset's `camelshotnoise`), while a
//! flattening places the tip between its neighbours (about `−4` there). The 22 tips then carry
//! almost all of the distortion, and the solver spends its cones on them.
//!
//! # The loop
//!
//! 0. The cone-free probe ([`FeedbackOptions::cone_free_probe`], where every component has a
//!    boundary): steps 2 to 4 on the surface without cones, before any solve. A cone that the
//!    solver puts at a needle tip absorbs its curvature and hides the failure of the
//!    linearization there; without cones the tip stands out (`|δ| ≈ 1.7·10⁴` on
//!    `camelshotnoise`, at most 3.5 elsewhere on the dataset).
//! 1. Solve with the current weights `wE`.
//! 2. Flatten with the cones (the caller does this, e.g. with Boundary First Flattening; this
//!    crate does not depend on a flattener) and measure the log area factor of every face of the
//!    flattening, `r_t = ½ ln(|B_t|/A_t)`.
//! 3. Compare at every vertex the solver's `u_i` with the realized factor around it
//!    ([`scale_mismatch`]): `δ_i = u_i − R_i − c`, with `R_i` the area-weighted mean of `r_t`
//!    over the faces at `i` and `c` the mass-weighted median of `u − R` over the connected
//!    component (which removes the free constant of `u` on a closed mesh and the scale of the
//!    plane; `δ` is invariant under both, `lean/proofs`, `cone-feedback-mismatch`).
//! 4. Every interior vertex with `wE_i > 0`, not reweighted yet, with `|δ_i| > threshold` has its
//!    weight multiplied by [`FeedbackOptions::reduced_weight`] (default `0`, which excludes it
//!    from the distortion measure, as the README recommends by hand for needle tips).
//! 5. Repeat until no vertex is added ([`FeedbackStop::Stable`]), or a limit ends it.
//!
//! The rule is a design choice, not a formula (`lean/proofs/CITED.md`, `cone-feedback-rule`):
//! the threshold, the per-vertex reweighting and the limits were chosen by measurement (crate
//! README, "The feedback loop"). Each solve is certified for its own weights; the final
//! solution is the certified optimum of the canonical problem with the returned
//! [`FeedbackReport::weights`].
//!
//! Weights only ever decrease and every vertex is reweighted at most once, so the loop ends
//! after at most [`FeedbackOptions::max_rounds`] solves. If more than
//! [`FeedbackOptions::max_reweighted_fraction`] of the vertices would be reweighted, the loop
//! stops and keeps the last solution: a mismatch everywhere is not a spike, and lowering
//! weights everywhere would change the problem, not repair it.
//!
//! # Interfaces
//!
//! [`Feedback`] is the rule as a state machine that the caller drives with its own solve and
//! flattening (`ddg-uv` does that for the pieces of its pipeline). [`solve_with_feedback`] drives
//! it with [`solve_with_boundary`] and a closure that flattens.

use std::fmt;
use std::time::Instant;

use ddg_core::operators::lumped_mass;
use ddg_core::{Geometry, TriMesh};

use crate::error::ConeError;
use crate::evaluate::evaluate_cones;
use crate::problem::{BoundaryTreatment, ConeProblem};
use crate::solve::{Solution, SolveOptions, solve_with_boundary};
use crate::{Cone, Parallelism};

/// The options of the feedback loop. Plain values; [`FeedbackOptions::validate`] checks them.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct FeedbackOptions {
    /// A vertex is reweighted where `|δ_i|` exceeds this (natural log units of the length
    /// scale; `δ = 1` is a factor `e` in length, `e²` in area). Positive and finite.
    /// Default `20`: on the 23 meshes of the published dataset the vertices of the clean
    /// meshes have `|δ|` at most 3.5 without cones and 7.6 after a solve, the 22 needle tips of
    /// `camelshotnoise` about `10⁴` without cones (crate README, "The feedback loop").
    pub threshold: f64,
    /// The factor by which the distortion weight of a reweighted vertex is multiplied,
    /// `0 ≤ reduced_weight < 1`. Default `0`: the vertex no longer counts.
    pub reduced_weight: f64,
    /// The most solves, the first included: `1` only measures (the mismatch is reported and
    /// nothing is solved again). At least 1. Default `4`.
    pub max_rounds: usize,
    /// The most vertices that may be reweighted in total, as a fraction of the vertices of the
    /// mesh, between 0 and 1. If a round would exceed it, the loop stops without applying that
    /// round ([`FeedbackStop::TooManyVertices`]). Default `0.01`: a mesh of fewer than 100
    /// vertices is never reweighted.
    pub max_reweighted_fraction: f64,
    /// Whether the loop first looks at the surface **without cones** (only the prescribed
    /// ones), flattened, before the first solve ([`Feedback::observe_probe`]): a cone that the
    /// solver puts at a needle tip absorbs its curvature and hides the failure of the
    /// linearization there. Only on meshes whose every component has a boundary. Default
    /// `true`.
    pub cone_free_probe: bool,
}

impl Default for FeedbackOptions {
    fn default() -> Self {
        FeedbackOptions {
            threshold: 20.0,
            reduced_weight: 0.0,
            max_rounds: 4,
            max_reweighted_fraction: 0.01,
            cone_free_probe: true,
        }
    }
}

impl FeedbackOptions {
    /// Checks every option.
    ///
    /// # Errors
    /// [`ConeError::InvalidOptions`] naming the option, its value and the requirement.
    pub fn validate(&self) -> Result<(), ConeError> {
        let bad = |reason: String| Err(ConeError::InvalidOptions { reason });
        if !(self.threshold.is_finite() && self.threshold > 0.0) {
            return bad(format!(
                "feedback threshold = {}: must be positive and finite",
                self.threshold
            ));
        }
        if !(self.reduced_weight >= 0.0 && self.reduced_weight < 1.0) {
            return bad(format!(
                "feedback reduced_weight = {}: must be at least 0 and below 1",
                self.reduced_weight
            ));
        }
        if self.max_rounds == 0 {
            return bad("feedback max_rounds = 0: must be at least 1".to_owned());
        }
        if !(0.0..=1.0).contains(&self.max_reweighted_fraction) {
            return bad(format!(
                "feedback max_reweighted_fraction = {}: must be between 0 and 1",
                self.max_reweighted_fraction
            ));
        }
        Ok(())
    }
}

/// The comparison of the solver's log scale factor with a realized one ([`scale_mismatch`]).
#[derive(Debug, Clone, PartialEq)]
pub struct ScaleMismatch {
    /// `δ_i` per vertex; NaN where no face at the vertex was measured.
    pub mismatch: Vec<f64>,
    /// `R_i`, the area-weighted mean of the realized log area factor of the faces at `i`; NaN
    /// where none was measured.
    pub realized: Vec<f64>,
    /// The offset `c` of every connected component (numbered as
    /// `TriMesh::face_components`); NaN for a component without a measured vertex.
    pub offsets: Vec<f64>,
    /// The largest `|δ_i|` over the measured vertices (0 if there is none).
    pub max_abs: f64,
    /// The mass-weighted root mean square of `δ` over the measured vertices (0 if there is
    /// none).
    pub rms: f64,
    /// Number of vertices without a measured face.
    pub n_unmeasured: usize,
}

/// The weighted quantile: with the values in ascending order (equal values in the order of
/// their index), the first value at which the weights summed so far reach `q` times the sum
/// of all weights. `q = ½` is the lower weighted median. An empty list gives 0. The values must
/// not be NaN. This is the rule of `ddg-variational-cuts` (`weighted_quantile`) under the same
/// id; `lean/proofs` characterizes it (`IsWeightedQuantile`) and proves it shift-equivariant
/// (`weightedQuantile_shift`).
///
/// [formula: cut-weighted-quantile]
pub(crate) fn weighted_quantile(values: &[f64], weights: &[f64], q: f64) -> f64 {
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

/// The mismatch between the solver's log scale factor `u` and the one a flattening realizes.
///
/// `realized_face_log_scale[t]` is the realized log area factor `r_t = ½ ln(|B_t| / A_t)` of
/// face `t` (`A_t` its area on the surface, `B_t` in the plane; the absolute value, so that a
/// flipped face still says how much it was scaled), or NaN (any non-finite value) where it was
/// not measured. Then, for every vertex `i`,
///
/// ```text
/// R_i = Σ_{t ∋ i} A_t r_t / Σ_{t ∋ i} A_t        (measured faces with A_t > 0 only)
/// δ_i = (u_i − R_i) − c_k,   c_k = weighted median of (u_j − R_j) over the measured vertices j
///                                   of the component k of i, weights the lumped masses m_j
/// ```
///
/// `δ` does not change if a constant is added to `u` (the free constant on a closed mesh) or
/// to `r` (a uniform scale of the plane) (`lean/proofs`, `cone-feedback-mismatch`), and it is
/// zero at every vertex where the flattening realizes `u` up to such constants. The vertex
/// value `u_i` is compared with the mean of the faces around it, not with the mean of the
/// solver's own face values `(u_a + u_b + u_c)/3`: the latter would spread a tip's huge `u`
/// over its neighbours. The price is a bias where `u` has a kink: at a cone of angle `µ`
/// `u_i` lies below the average of its neighbours by `(µ − Ω_i) / Σ_j w_ij` (the cotan weights),
/// about `1.7` for `µ = 0.95 · 2π` on a regular mesh, and the realized face values are
/// averages of that; the default threshold has room for it (crate README).
///
/// [formula: cone-feedback-mismatch]
///
/// # Errors
/// [`ConeError::Geometry`] if the geometry is not that of the mesh, [`ConeError::WrongLength`]
/// if `log_scale` does not have one entry per vertex, [`ConeError::InvalidArgument`] if
/// `realized_face_log_scale` does not have one entry per face or `log_scale` is not finite.
pub fn scale_mismatch(
    mesh: &TriMesh,
    geometry: &Geometry,
    log_scale: &[f64],
    realized_face_log_scale: &[f64],
) -> Result<ScaleMismatch, ConeError> {
    geometry.check_mesh(mesh).map_err(ConeError::Geometry)?;
    let nv = mesh.n_vertices();
    let nf = mesh.n_faces();
    if log_scale.len() != nv {
        return Err(ConeError::WrongLength {
            what: "log_scale",
            expected: nv,
            got: log_scale.len(),
        });
    }
    if realized_face_log_scale.len() != nf {
        return Err(ConeError::InvalidArgument {
            reason: format!(
                "realized_face_log_scale has {} entries, but the mesh has {nf} faces",
                realized_face_log_scale.len()
            ),
        });
    }
    if let Some(v) = log_scale.iter().position(|x| !x.is_finite()) {
        return Err(ConeError::InvalidArgument {
            reason: format!("log_scale[{v}] = {} is not finite", log_scale[v]),
        });
    }
    let areas = geometry.face_areas();
    let mut area_sum = vec![0.0; nv];
    let mut weighted_sum = vec![0.0; nv];
    for (f, face) in mesh.faces().iter().enumerate() {
        let (a, r) = (areas[f], realized_face_log_scale[f]);
        if !(a > 0.0 && r.is_finite()) {
            continue;
        }
        for &v in face {
            area_sum[v as usize] += a;
            weighted_sum[v as usize] += a * r;
        }
    }
    let realized: Vec<f64> = (0..nv)
        .map(|v| {
            if area_sum[v] > 0.0 {
                weighted_sum[v] / area_sum[v]
            } else {
                f64::NAN
            }
        })
        .collect();
    let raw: Vec<f64> = (0..nv).map(|v| log_scale[v] - realized[v]).collect();
    let mass = lumped_mass(mesh, geometry).map_err(ConeError::Geometry)?;
    let components = mesh.vertex_components();
    let nc = mesh.n_components();
    let mut values: Vec<Vec<f64>> = vec![Vec::new(); nc];
    let mut weights: Vec<Vec<f64>> = vec![Vec::new(); nc];
    for v in (0..nv).filter(|&v| raw[v].is_finite()) {
        let c = components[v] as usize;
        values[c].push(raw[v]);
        weights[c].push(mass[v]);
    }
    let offsets: Vec<f64> = (0..nc)
        .map(|c| {
            if values[c].is_empty() {
                f64::NAN
            } else {
                weighted_quantile(&values[c], &weights[c], 0.5)
            }
        })
        .collect();
    let mismatch: Vec<f64> = (0..nv)
        .map(|v| raw[v] - offsets[components[v] as usize])
        .collect();
    let (mut max_abs, mut sum_m, mut sum_md2, mut n_unmeasured) = (0.0f64, 0.0, 0.0, 0);
    for v in 0..nv {
        let d = mismatch[v];
        if d.is_finite() {
            max_abs = max_abs.max(d.abs());
            sum_m += mass[v];
            sum_md2 += mass[v] * d * d;
        } else {
            n_unmeasured += 1;
        }
    }
    let rms = if sum_m > 0.0 {
        (sum_md2 / sum_m).sqrt()
    } else {
        0.0
    };
    Ok(ScaleMismatch {
        mismatch,
        realized,
        offsets,
        max_abs,
        rms,
        n_unmeasured,
    })
}

/// Why the feedback loop ended.
#[derive(Debug, Clone, PartialEq)]
pub enum FeedbackStop {
    /// The last round found no vertex to reweight: the solution is consistent with its
    /// flattening at the threshold.
    Stable,
    /// The last round found vertices to reweight, but [`FeedbackOptions::max_rounds`] solves
    /// were done; they are listed in [`FeedbackReport::pending`] and were not applied.
    RoundLimit,
    /// The last round would have reweighted more vertices in total than
    /// [`FeedbackOptions::max_reweighted_fraction`] allows; nothing of that round was applied.
    TooManyVertices {
        /// How many vertices would have been reweighted in total.
        wanted: usize,
        /// How many are allowed.
        allowed: usize,
    },
    /// With `reduced_weight = 0` the last round would have left a connected component without
    /// any positive distortion weight; nothing of that round was applied.
    NoWeightLeft {
        /// The component (numbered as `TriMesh::face_components`).
        component: u32,
    },
    /// The solve or the flattening after the last reweighting failed; the reweighting of that
    /// round was undone and the previous solution kept.
    RoundFailed {
        /// The round that failed (0 is the first solve).
        round: usize,
        /// What the caller reported.
        message: String,
    },
    /// The caller ended the loop ([`Feedback::finish`]) although it asked for another solve;
    /// the reweighting of the last round was undone.
    NotContinued,
}

/// One round of the loop: what the comparison of one solve with its flattening gave.
#[derive(Debug, Clone, PartialEq)]
pub struct FeedbackRound {
    /// Whether this round is the cone-free probe ([`Feedback::observe_probe`]), not a solve.
    pub probe: bool,
    /// The largest `|δ_i|` over the measured vertices.
    pub max_abs_mismatch: f64,
    /// The mass-weighted root mean square of `δ` over the measured vertices.
    pub rms_mismatch: f64,
    /// The largest `|δ_i|` over the vertices the rule looks at (interior, positive weight, not
    /// reweighted before this round): the margin to the threshold on a mesh where nothing is
    /// reweighted. 0 if there is none.
    pub max_abs_eligible: f64,
    /// The offsets `c` per connected component ([`ScaleMismatch::offsets`]).
    pub offsets: Vec<f64>,
    /// Vertices without a measured face.
    pub n_unmeasured: usize,
    /// Vertices the rule picked in this round (`|δ| > threshold`, eligible).
    pub n_over_threshold: usize,
    /// The vertices reweighted after this round, ascending; empty if the loop ended here.
    pub reweighted: Vec<u32>,
    /// Seconds of the solve of this round, as the caller reported them.
    pub seconds_solve: f64,
    /// Seconds of the flattening and its measurement, as the caller reported them.
    pub seconds_realize: f64,
    /// Seconds of the comparison and the rule.
    pub seconds_mismatch: f64,
}

/// What the feedback loop did.
#[derive(Debug, Clone, PartialEq)]
pub struct FeedbackReport {
    /// One entry per solve, the first included, after the cone-free probe if it ran.
    pub rounds: Vec<FeedbackRound>,
    /// Every reweighted vertex, ascending.
    pub reweighted: Vec<u32>,
    /// Vertices the last round picked but that were not reweighted because the loop ended
    /// ([`FeedbackStop::RoundLimit`], [`FeedbackStop::TooManyVertices`],
    /// [`FeedbackStop::NoWeightLeft`], [`FeedbackStop::RoundFailed`],
    /// [`FeedbackStop::NotContinued`]), ascending.
    pub pending: Vec<u32>,
    /// The distortion weights of the last kept solve, one per vertex.
    pub weights: Vec<f64>,
    /// Why the loop ended.
    pub stop: FeedbackStop,
    /// Seconds from [`Feedback::new`] to [`Feedback::finish`], the caller's solves and
    /// flattenings included.
    pub seconds: f64,
    /// Why the cone-free probe did not run or failed, if it was asked for (set by the drivers,
    /// [`solve_with_feedback`] and `ddg-uv`; `None` from [`Feedback::finish`]).
    pub probe_note: Option<String>,
}

/// The feedback rule as a state machine. The caller solves with [`Feedback::weights`], flattens,
/// and passes both to [`Feedback::observe`], which says whether to solve again.
///
/// ```text
/// let mut loop_ = Feedback::new(mesh, initial_weights, &options)?;
/// let mut kept = solve_and_flatten(loop_.weights());
/// while loop_.observe(mesh, geometry, &kept.u, &kept.realized, t_solve, t_flatten)? {
///     match solve_and_flatten(loop_.weights()) {
///         Ok(next) => kept = next,
///         Err(e) => { loop_.round_failed(e.to_string()); break; }
///     }
/// }
/// let report = loop_.finish();   // report.weights are those of `kept`
/// ```
#[derive(Debug, Clone)]
pub struct Feedback {
    options: FeedbackOptions,
    weights: Vec<f64>,
    eligible: Vec<bool>,
    is_reweighted: Vec<bool>,
    components: Vec<u32>,
    allowed: usize,
    rounds: Vec<FeedbackRound>,
    /// `(vertex, weight before)` of the last applied round, to undo it.
    last_applied: Vec<(u32, f64)>,
    pending: Vec<u32>,
    stop: Option<FeedbackStop>,
    start: Instant,
}

impl Feedback {
    /// Sets the loop up for `mesh`, with the distortion weights of the first solve
    /// (`None`: 1 everywhere).
    ///
    /// # Errors
    /// [`ConeError::InvalidOptions`] for invalid options, [`ConeError::WrongLength`] and
    /// [`ConeError::InvalidDistortionWeight`] for invalid weights.
    pub fn new(
        mesh: &TriMesh,
        initial_weights: Option<&[f64]>,
        options: &FeedbackOptions,
    ) -> Result<Feedback, ConeError> {
        options.validate()?;
        let nv = mesh.n_vertices();
        let weights = match initial_weights {
            None => vec![1.0; nv],
            Some(w) => {
                if w.len() != nv {
                    return Err(ConeError::WrongLength {
                        what: "distortion_weights",
                        expected: nv,
                        got: w.len(),
                    });
                }
                if let Some(v) = w.iter().position(|x| !(x.is_finite() && *x >= 0.0)) {
                    return Err(ConeError::InvalidDistortionWeight {
                        vertex: v as u32,
                        value: w[v],
                    });
                }
                w.to_vec()
            }
        };
        let eligible = (0..nv)
            .map(|v| !mesh.is_boundary_vertex(v) && weights[v] > 0.0)
            .collect();
        Ok(Feedback {
            options: *options,
            eligible,
            is_reweighted: vec![false; nv],
            components: mesh.vertex_components(),
            allowed: (options.max_reweighted_fraction * nv as f64).floor() as usize,
            weights,
            rounds: Vec::new(),
            last_applied: Vec::new(),
            pending: Vec::new(),
            stop: None,
            start: Instant::now(),
        })
    }

    /// The distortion weights for the next solve.
    pub fn weights(&self) -> &[f64] {
        &self.weights
    }

    /// Number of rounds observed so far.
    pub fn n_rounds(&self) -> usize {
        self.rounds.len()
    }

    /// Whether the loop has ended.
    pub fn is_finished(&self) -> bool {
        self.stop.is_some()
    }

    /// Observes one round: `log_scale` is the `u` the solve with [`Feedback::weights`] returned
    /// for the cones that were flattened (`Solution::log_scale`), `realized_face_log_scale` the
    /// log area factor of every face of their flattening (see [`scale_mismatch`]). The seconds
    /// are only reported. Returns `true` if the weights were lowered and the caller should
    /// solve and flatten again, `false` if the loop has ended.
    ///
    /// The rule is a design choice, not a formula (`lean/proofs/CITED.md`).
    ///
    /// [formula: cone-feedback-rule]
    ///
    /// # Errors
    /// Those of [`scale_mismatch`]; [`ConeError::InvalidArgument`] if the loop has already
    /// ended or `mesh` has another number of vertices than the one it was set up for.
    pub fn observe(
        &mut self,
        mesh: &TriMesh,
        geometry: &Geometry,
        log_scale: &[f64],
        realized_face_log_scale: &[f64],
        seconds_solve: f64,
        seconds_realize: f64,
    ) -> Result<bool, ConeError> {
        self.step(
            mesh,
            geometry,
            log_scale,
            realized_face_log_scale,
            [seconds_solve, seconds_realize],
            false,
        )
    }

    /// Observes the **cone-free probe**, before the first solve: `log_scale` is the scale factor
    /// of the surface without cones (only prescribed ones; `evaluate_cones`, `u = 0` on the
    /// boundary), `realized_face_log_scale` that of its flattening. A cone at the tip of a needle
    /// absorbs its curvature and hides the failure of the linearization there; without cones
    /// the tip's `|δ|` is of order `h/ρ`. The rule picks as in [`Feedback::observe`]; the picks
    /// are applied if the limits allow it (otherwise nothing is applied), and the loop never
    /// ends here: the caller solves next, with [`Feedback::weights`]. Returns whether weights
    /// were lowered (then a solve made with the previous weights is out of date). The probe
    /// does not count as a solve for [`FeedbackOptions::max_rounds`]. Only on meshes whose
    /// every component has a boundary (a closed surface cannot be flattened without cones).
    ///
    /// [formula: cone-feedback-rule]
    ///
    /// # Errors
    /// As [`Feedback::observe`], and [`ConeError::InvalidArgument`] if a round was observed
    /// already.
    pub fn observe_probe(
        &mut self,
        mesh: &TriMesh,
        geometry: &Geometry,
        log_scale: &[f64],
        realized_face_log_scale: &[f64],
        seconds_evaluate: f64,
        seconds_realize: f64,
    ) -> Result<bool, ConeError> {
        if !self.rounds.is_empty() {
            return Err(ConeError::InvalidArgument {
                reason: "the cone-free probe must come before the first round".into(),
            });
        }
        self.step(
            mesh,
            geometry,
            log_scale,
            realized_face_log_scale,
            [seconds_evaluate, seconds_realize],
            true,
        )?;
        Ok(!self.last_applied.is_empty())
    }

    fn step(
        &mut self,
        mesh: &TriMesh,
        geometry: &Geometry,
        log_scale: &[f64],
        realized_face_log_scale: &[f64],
        seconds: [f64; 2],
        probe: bool,
    ) -> Result<bool, ConeError> {
        if self.stop.is_some() {
            return Err(ConeError::InvalidArgument {
                reason: "the feedback loop has ended; it cannot observe another round".into(),
            });
        }
        let nv = self.weights.len();
        if mesh.n_vertices() != nv {
            return Err(ConeError::InvalidArgument {
                reason: format!(
                    "the mesh has {} vertices, but the feedback loop was set up for {nv}",
                    mesh.n_vertices()
                ),
            });
        }
        let start = Instant::now();
        let m = scale_mismatch(mesh, geometry, log_scale, realized_face_log_scale)?;
        let looked_at = |v: usize| self.eligible[v] && !self.is_reweighted[v];
        let max_abs_eligible = (0..nv)
            .filter(|&v| looked_at(v) && m.mismatch[v].is_finite())
            .map(|v| m.mismatch[v].abs())
            .fold(0.0f64, f64::max);
        let candidates: Vec<u32> = (0..nv)
            .filter(|&v| looked_at(v) && m.mismatch[v].abs() > self.options.threshold)
            .map(|v| v as u32)
            .collect();
        let n_over_threshold = candidates.len();
        let mut round = FeedbackRound {
            probe,
            max_abs_mismatch: m.max_abs,
            rms_mismatch: m.rms,
            max_abs_eligible,
            offsets: m.offsets,
            n_unmeasured: m.n_unmeasured,
            n_over_threshold,
            reweighted: Vec::new(),
            seconds_solve: seconds[0],
            seconds_realize: seconds[1],
            seconds_mismatch: 0.0,
        };
        self.last_applied.clear();
        let stop = if candidates.is_empty() {
            Some(FeedbackStop::Stable)
        } else if !probe && self.n_solves() + 1 >= self.options.max_rounds {
            Some(FeedbackStop::RoundLimit)
        } else if self.n_reweighted() + candidates.len() > self.allowed {
            Some(FeedbackStop::TooManyVertices {
                wanted: self.n_reweighted() + candidates.len(),
                allowed: self.allowed,
            })
        } else if self.options.reduced_weight == 0.0 {
            self.component_left_without_weight(&candidates)
                .map(|component| FeedbackStop::NoWeightLeft { component })
        } else {
            None
        };
        match stop {
            // The probe never ends the loop: a limit only means that nothing is applied.
            Some(_) if probe => {}
            Some(stop) => {
                if stop != FeedbackStop::Stable {
                    self.pending = candidates;
                }
                self.stop = Some(stop);
            }
            None => {
                for &v in &candidates {
                    let v = v as usize;
                    self.last_applied.push((v as u32, self.weights[v]));
                    self.weights[v] *= self.options.reduced_weight;
                    self.is_reweighted[v] = true;
                }
                round.reweighted = candidates;
            }
        }
        round.seconds_mismatch = start.elapsed().as_secs_f64();
        self.rounds.push(round);
        Ok(self.stop.is_none())
    }

    /// Number of solves observed (rounds other than the probe).
    fn n_solves(&self) -> usize {
        self.rounds.iter().filter(|r| !r.probe).count()
    }

    /// The solve or the flattening with the lowered weights failed: the last reweighting is
    /// undone (the caller keeps the previous solution) and the loop ends with
    /// [`FeedbackStop::RoundFailed`]. Does nothing if the loop has already ended.
    pub fn round_failed(&mut self, message: impl Into<String>) {
        if self.stop.is_none() {
            self.undo_last();
            self.stop = Some(FeedbackStop::RoundFailed {
                round: self.rounds.len(),
                message: message.into(),
            });
        }
    }

    /// Ends the loop and returns the report. If the loop asked for another solve that was not
    /// observed, the last reweighting is undone and the loop ends with
    /// [`FeedbackStop::NotContinued`].
    pub fn finish(mut self) -> FeedbackReport {
        if self.stop.is_none() {
            self.undo_last();
            self.stop = Some(FeedbackStop::NotContinued);
        }
        let mut reweighted: Vec<u32> = (0..self.weights.len())
            .filter(|&v| self.is_reweighted[v])
            .map(|v| v as u32)
            .collect();
        reweighted.sort_unstable();
        self.pending.sort_unstable();
        FeedbackReport {
            rounds: self.rounds,
            reweighted,
            pending: self.pending,
            weights: self.weights,
            stop: self.stop.unwrap_or(FeedbackStop::NotContinued),
            seconds: self.start.elapsed().as_secs_f64(),
            probe_note: None,
        }
    }

    fn n_reweighted(&self) -> usize {
        self.is_reweighted.iter().filter(|&&r| r).count()
    }

    fn undo_last(&mut self) {
        let undone: Vec<u32> = self.last_applied.iter().map(|&(v, _)| v).collect();
        for &(v, w) in &self.last_applied {
            self.weights[v as usize] = w;
            self.is_reweighted[v as usize] = false;
        }
        self.last_applied.clear();
        if let Some(last) = self.rounds.last_mut() {
            last.reweighted.clear();
        }
        self.pending = undone;
    }

    /// The first component in which zeroing the weights of `candidates` leaves no positive
    /// weight.
    fn component_left_without_weight(&self, candidates: &[u32]) -> Option<u32> {
        let nc = self
            .components
            .iter()
            .map(|&c| c as usize + 1)
            .max()
            .unwrap_or(0);
        let mut positive = vec![0usize; nc];
        for (v, &w) in self.weights.iter().enumerate() {
            if w > 0.0 {
                positive[self.components[v] as usize] += 1;
            }
        }
        for &v in candidates {
            positive[self.components[v as usize] as usize] -= 1;
        }
        candidates
            .iter()
            .map(|&v| self.components[v as usize])
            .find(|&c| positive[c as usize] == 0)
    }
}

/// The result of [`solve_with_feedback`].
#[derive(Debug, Clone)]
pub struct FeedbackSolution<R> {
    /// The solution of the last kept solve (certified for [`FeedbackReport::weights`]).
    pub solution: Solution,
    /// What the closure returned for it besides the face values (e.g. the UV map).
    pub realized: R,
    /// What the loop did.
    pub feedback: FeedbackReport,
}

/// The error of [`solve_with_feedback`]: invalid input, or a failure of the first solve or the
/// first flattening. Failures in later rounds are not errors: the loop keeps the previous
/// solution and reports [`FeedbackStop::RoundFailed`].
#[derive(Debug, Clone, PartialEq)]
pub enum FeedbackError<E> {
    /// The problem, the options or the first solve.
    Cones(ConeError),
    /// The closure failed on the first solution.
    Realize(E),
}

impl<E: fmt::Display> fmt::Display for FeedbackError<E> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            FeedbackError::Cones(e) => write!(f, "{e}"),
            FeedbackError::Realize(e) => {
                write!(f, "the flattening of the first solution failed: {e}")
            }
        }
    }
}

impl<E: fmt::Debug + fmt::Display> std::error::Error for FeedbackError<E> {}

/// Whether every connected component of the mesh has a boundary vertex (then the surface
/// without cones can be flattened, and the cone-free probe applies).
pub fn every_component_has_boundary(mesh: &TriMesh) -> bool {
    let components = mesh.vertex_components();
    let mut has = vec![false; mesh.n_components()];
    for (v, &c) in components.iter().enumerate() {
        if mesh.is_boundary_vertex(v) {
            has[c as usize] = true;
        }
    }
    has.iter().all(|&h| h)
}

/// Solves `problem` with the feedback loop. `realize(cones, log_scale)` flattens the surface
/// with these cones and scale factors and returns the realized log area factor of every face
/// (see [`scale_mismatch`]) plus anything the caller wants back.
///
/// With [`FeedbackOptions::cone_free_probe`] (and a mesh whose every component has a boundary)
/// the loop first flattens the surface without cones (the prescribed ones only, `u` from
/// [`evaluate_cones`](crate::evaluate_cones), zero on the boundary) and applies the rule to it
/// ([`Feedback::observe_probe`]). Then: [`solve_with_boundary`], `realize` with the solution's
/// `cones` and `log_scale`, the rule, until it is stable or a limit ends it. The distortion
/// weights of `problem` are those the loop starts from.
///
/// # Errors
/// [`FeedbackError::Cones`] for invalid options or problems and errors of the first solve,
/// [`FeedbackError::Realize`] if `realize` fails on the first solution. A failure of the probe
/// is not an error: the loop goes on without it ([`FeedbackReport::probe_note`]).
pub fn solve_with_feedback<R, E: fmt::Display>(
    problem: &ConeProblem<'_>,
    boundary: BoundaryTreatment,
    options: &SolveOptions<'_>,
    feedback: &FeedbackOptions,
    mut realize: impl FnMut(&[Cone], &[f64]) -> Result<(Vec<f64>, R), E>,
) -> Result<FeedbackSolution<R>, FeedbackError<E>> {
    let mut state = Feedback::new(
        problem.mesh,
        problem.distortion_weights.as_deref(),
        feedback,
    )
    .map_err(FeedbackError::Cones)?;
    let mut probe_note = None;
    if feedback.cone_free_probe {
        if every_component_has_boundary(problem.mesh) {
            let start = Instant::now();
            match evaluate_cones(problem, &problem.prescribed_cones, Parallelism::Sequential) {
                Err(e) => probe_note = Some(format!("the cone-free state failed: {e}")),
                Ok(free) => {
                    let seconds_evaluate = start.elapsed().as_secs_f64();
                    let start = Instant::now();
                    match realize(&problem.prescribed_cones, &free.log_scale) {
                        Err(e) => {
                            probe_note = Some(format!("the cone-free flattening failed: {e}"));
                        }
                        Ok((faces, _)) => {
                            let seconds_realize = start.elapsed().as_secs_f64();
                            state
                                .observe_probe(
                                    problem.mesh,
                                    problem.geometry,
                                    &free.log_scale,
                                    &faces,
                                    seconds_evaluate,
                                    seconds_realize,
                                )
                                .map_err(FeedbackError::Cones)?;
                        }
                    }
                }
            }
        } else {
            probe_note = Some("a component has no boundary".to_owned());
        }
    }
    let weighted = |state: &Feedback| ConeProblem {
        distortion_weights: Some(state.weights().to_vec()),
        ..problem.clone()
    };
    let start = Instant::now();
    let mut solution =
        solve_with_boundary(&weighted(&state), boundary, options).map_err(FeedbackError::Cones)?;
    let mut seconds_solve = start.elapsed().as_secs_f64();
    let start = Instant::now();
    let (mut faces, mut realized) =
        realize(&solution.cones, &solution.log_scale).map_err(FeedbackError::Realize)?;
    let mut seconds_realize = start.elapsed().as_secs_f64();
    loop {
        let again = state
            .observe(
                problem.mesh,
                problem.geometry,
                &solution.log_scale,
                &faces,
                seconds_solve,
                seconds_realize,
            )
            .map_err(FeedbackError::Cones)?;
        if !again {
            break;
        }
        let start = Instant::now();
        let next = match solve_with_boundary(&weighted(&state), boundary, options) {
            Ok(next) => next,
            Err(e) => {
                state.round_failed(format!("the solve failed: {e}"));
                break;
            }
        };
        seconds_solve = start.elapsed().as_secs_f64();
        let start = Instant::now();
        match realize(&next.cones, &next.log_scale) {
            Ok((f, r)) => {
                faces = f;
                realized = r;
                solution = next;
            }
            Err(e) => {
                state.round_failed(format!("the flattening failed: {e}"));
                break;
            }
        }
        seconds_realize = start.elapsed().as_secs_f64();
    }
    let mut report = state.finish();
    report.probe_note = probe_note;
    Ok(FeedbackSolution {
        solution,
        realized,
        feedback: report,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use ddg_core::generators::{bumpy_grid_disk, icosphere};

    fn grid(n: usize) -> (TriMesh, Geometry, Vec<[f64; 3]>) {
        let g = bumpy_grid_disk(n, n, 0.0).unwrap();
        let mesh = TriMesh::new(g.positions.len(), g.faces.clone()).unwrap();
        let geometry = Geometry::from_positions(&mesh, &g.positions).unwrap();
        (mesh, geometry, g.positions)
    }

    /// The model face values `(u_a + u_b + u_c)/3`.
    fn face_means(mesh: &TriMesh, u: &[f64]) -> Vec<f64> {
        mesh.faces()
            .iter()
            .map(|f| f.iter().map(|&v| u[v as usize]).sum::<f64>() / 3.0)
            .collect()
    }

    #[test]
    fn weighted_quantile_is_the_lower_weighted_median() {
        assert_eq!(weighted_quantile(&[], &[], 0.5), 0.0);
        assert_eq!(
            weighted_quantile(&[3.0, 1.0, 2.0], &[1.0, 1.0, 1.0], 0.5),
            2.0
        );
        // Half the weight exactly at the first value: the lower median.
        assert_eq!(weighted_quantile(&[1.0, 2.0], &[1.0, 1.0], 0.5), 1.0);
        // A heavy value wins.
        assert_eq!(
            weighted_quantile(&[5.0, 1.0, 2.0], &[10.0, 1.0, 1.0], 0.5),
            5.0
        );
        // An outlier of tiny weight does not move it.
        assert_eq!(
            weighted_quantile(&[0.0, 0.1, -0.1, -1e4], &[1.0, 1.0, 1.0, 1e-6], 0.5),
            0.0
        );
        // Shifting the values shifts the median.
        let values = [0.3, -2.0, 7.5, 1.25, 1.25];
        let weights = [0.5, 2.0, 0.1, 1.0, 0.7];
        let shifted: Vec<f64> = values.iter().map(|x| x + 3.5).collect();
        assert_eq!(
            weighted_quantile(&shifted, &weights, 0.5),
            weighted_quantile(&values, &weights, 0.5) + 3.5
        );
    }

    #[test]
    fn a_realized_scale_equal_to_the_solvers_has_no_mismatch() {
        // u constant on every component: the face values equal it, and δ = 0.
        let (mesh, geometry, _) = grid(8);
        let u = vec![0.7; mesh.n_vertices()];
        let realized = vec![-1.3; mesh.n_faces()];
        let m = scale_mismatch(&mesh, &geometry, &u, &realized).unwrap();
        assert_eq!(m.n_unmeasured, 0);
        assert!(m.max_abs < 1e-12 && m.rms < 1e-12, "{m:?}");
        assert!((m.offsets[0] - 2.0).abs() < 1e-12);
    }

    #[test]
    fn the_mismatch_does_not_depend_on_constants() {
        let (mesh, geometry, positions) = grid(9);
        let u: Vec<f64> = positions
            .iter()
            .map(|p| (3.0 * p[0]).sin() + p[1] * p[1])
            .collect();
        let mut realized = face_means(&mesh, &u);
        realized[7] += 0.4; // some disagreement
        realized[30] = f64::NAN; // and a face that was not measured
        let a = scale_mismatch(&mesh, &geometry, &u, &realized).unwrap();
        let u2: Vec<f64> = u.iter().map(|x| x - 11.0).collect();
        let r2: Vec<f64> = realized.iter().map(|x| x + 2.5).collect();
        let b = scale_mismatch(&mesh, &geometry, &u2, &r2).unwrap();
        for (x, y) in a.mismatch.iter().zip(&b.mismatch) {
            assert!((x - y).abs() < 1e-9, "{x} {y}");
        }
        assert!(a.max_abs > 0.0);
    }

    #[test]
    fn unmeasured_vertices_are_nan_and_counted() {
        let (mesh, geometry, _) = grid(4);
        let u = vec![0.0; mesh.n_vertices()];
        let realized = vec![f64::NAN; mesh.n_faces()];
        let m = scale_mismatch(&mesh, &geometry, &u, &realized).unwrap();
        assert_eq!(m.n_unmeasured, mesh.n_vertices());
        assert!(m.mismatch.iter().all(|x| x.is_nan()));
        assert_eq!((m.max_abs, m.rms), (0.0, 0.0));
        assert!(m.offsets[0].is_nan());
    }

    #[test]
    fn wrong_lengths_and_bad_values_are_errors() {
        let (mesh, geometry, _) = grid(4);
        let nv = mesh.n_vertices();
        let nf = mesh.n_faces();
        assert!(matches!(
            scale_mismatch(&mesh, &geometry, &vec![0.0; nv - 1], &vec![0.0; nf]),
            Err(ConeError::WrongLength { .. })
        ));
        assert!(matches!(
            scale_mismatch(&mesh, &geometry, &vec![0.0; nv], &vec![0.0; nf + 1]),
            Err(ConeError::InvalidArgument { .. })
        ));
        let mut u = vec![0.0; nv];
        u[3] = f64::INFINITY;
        let error = scale_mismatch(&mesh, &geometry, &u, &vec![0.0; nf]).unwrap_err();
        assert!(error.to_string().contains("log_scale[3]"), "{error}");
    }

    #[test]
    fn invalid_options_are_errors() {
        let (mesh, _, _) = grid(4);
        let d = FeedbackOptions::default();
        for bad in [
            FeedbackOptions {
                threshold: 0.0,
                ..d
            },
            FeedbackOptions {
                threshold: -1.0,
                ..d
            },
            FeedbackOptions {
                threshold: f64::NAN,
                ..d
            },
            FeedbackOptions {
                threshold: f64::INFINITY,
                ..d
            },
            FeedbackOptions {
                reduced_weight: 1.0,
                ..d
            },
            FeedbackOptions {
                reduced_weight: -0.1,
                ..d
            },
            FeedbackOptions {
                reduced_weight: f64::NAN,
                ..d
            },
            FeedbackOptions { max_rounds: 0, ..d },
            FeedbackOptions {
                max_reweighted_fraction: 1.5,
                ..d
            },
            FeedbackOptions {
                max_reweighted_fraction: f64::NAN,
                ..d
            },
        ] {
            let error = bad.validate().unwrap_err();
            assert!(matches!(error, ConeError::InvalidOptions { .. }), "{bad:?}");
            assert!(Feedback::new(&mesh, None, &bad).is_err(), "{bad:?}");
        }
        assert!(d.validate().is_ok());
        let nv = mesh.n_vertices();
        assert!(matches!(
            Feedback::new(&mesh, Some(vec![1.0; nv + 1].as_slice()), &d),
            Err(ConeError::WrongLength { .. })
        ));
        let mut w = vec![1.0; nv];
        w[2] = f64::NAN;
        assert!(matches!(
            Feedback::new(&mesh, Some(w.as_slice()), &d),
            Err(ConeError::InvalidDistortionWeight { vertex: 2, .. })
        ));
    }

    /// One vertex whose solver value is far from the faces around it: it is reweighted, once,
    /// and the next round with an agreeing flattening is stable.
    #[test]
    fn the_rule_reweights_a_far_vertex_once() {
        let (mesh, geometry, _) = grid(20); // 441 vertices: up to 4 may be reweighted
        let nv = mesh.n_vertices();
        let tip = (10 * 21 + 10) as usize;
        let mut u = vec![0.0; nv];
        u[tip] = -1.0e4;
        let realized = vec![0.0; mesh.n_faces()];
        let mut fb = Feedback::new(&mesh, None, &FeedbackOptions::default()).unwrap();
        assert!(
            fb.observe(&mesh, &geometry, &u, &realized, 0.0, 0.0)
                .unwrap()
        );
        assert_eq!(fb.weights()[tip], 0.0);
        assert_eq!(fb.weights().iter().filter(|&&w| w == 0.0).count(), 1);
        // The same mismatch again: the tip is not picked a second time, nothing else is.
        assert!(
            !fb.observe(&mesh, &geometry, &u, &realized, 0.0, 0.0)
                .unwrap()
        );
        let report = fb.finish();
        assert_eq!(report.reweighted, vec![tip as u32]);
        assert_eq!(report.stop, FeedbackStop::Stable);
        assert_eq!(report.rounds.len(), 2);
        assert_eq!(report.rounds[0].reweighted, vec![tip as u32]);
        assert!(report.rounds[0].max_abs_mismatch > 9.0e3);
        assert!(report.rounds[1].max_abs_eligible < 1.0);
    }

    #[test]
    fn the_limits_end_the_loop_without_applying_the_round() {
        let (mesh, geometry, _) = grid(20);
        let nv = mesh.n_vertices();
        let tips: Vec<usize> = (0..5).map(|k| (4 + 3 * k) * 21 + 10).collect();
        let mut u = vec![0.0; nv];
        for &t in &tips {
            u[t] = 50.0;
        }
        let realized = vec![0.0; mesh.n_faces()];
        // Five far vertices, four allowed.
        let mut fb = Feedback::new(&mesh, None, &FeedbackOptions::default()).unwrap();
        assert!(
            !fb.observe(&mesh, &geometry, &u, &realized, 0.0, 0.0)
                .unwrap()
        );
        let report = fb.finish();
        assert_eq!(
            report.stop,
            FeedbackStop::TooManyVertices {
                wanted: 5,
                allowed: 4
            }
        );
        assert!(report.reweighted.is_empty() && report.pending.len() == 5);
        assert!(report.weights.iter().all(|&w| w == 1.0));
        // One round only: measure, do not solve again.
        let one = FeedbackOptions {
            max_rounds: 1,
            ..FeedbackOptions::default()
        };
        let mut fb = Feedback::new(&mesh, None, &one).unwrap();
        u[tips[4]] = 0.0;
        assert!(
            !fb.observe(&mesh, &geometry, &u, &realized, 0.0, 0.0)
                .unwrap()
        );
        let report = fb.finish();
        assert_eq!(report.stop, FeedbackStop::RoundLimit);
        assert_eq!(report.pending.len(), 4);
        assert!(report.weights.iter().all(|&w| w == 1.0));
        // A failed round is undone.
        let mut fb = Feedback::new(&mesh, None, &FeedbackOptions::default()).unwrap();
        assert!(
            fb.observe(&mesh, &geometry, &u, &realized, 0.0, 0.0)
                .unwrap()
        );
        fb.round_failed("no flattening");
        let report = fb.finish();
        assert!(matches!(
            report.stop,
            FeedbackStop::RoundFailed { round: 1, .. }
        ));
        assert!(report.reweighted.is_empty());
        assert_eq!(report.pending.len(), 4);
        assert!(report.weights.iter().all(|&w| w == 1.0));
    }

    #[test]
    fn a_reduced_weight_factor_and_initial_weights() {
        let (mesh, geometry, _) = grid(20);
        let nv = mesh.n_vertices();
        let (a, b) = (5 * 21 + 5, 15 * 21 + 15);
        let mut u = vec![0.0; nv];
        u[a] = 30.0;
        u[b] = -30.0;
        let mut initial = vec![2.0; nv];
        initial[b] = 0.0; // already excluded: not eligible
        let options = FeedbackOptions {
            reduced_weight: 0.25,
            ..FeedbackOptions::default()
        };
        let mut fb = Feedback::new(&mesh, Some(initial.as_slice()), &options).unwrap();
        let realized = vec![0.0; mesh.n_faces()];
        assert!(
            fb.observe(&mesh, &geometry, &u, &realized, 0.0, 0.0)
                .unwrap()
        );
        assert_eq!(fb.weights()[a], 0.5);
        assert_eq!(fb.weights()[b], 0.0);
        assert!(
            !fb.observe(&mesh, &geometry, &u, &realized, 0.0, 0.0)
                .unwrap()
        );
        assert_eq!(fb.finish().reweighted, vec![a as u32]);
    }

    #[test]
    fn boundary_vertices_are_never_reweighted() {
        let (mesh, geometry, _) = grid(20);
        let nv = mesh.n_vertices();
        let corner = 0usize;
        assert!(mesh.is_boundary_vertex(corner));
        let mut u = vec![0.0; nv];
        u[corner] = 100.0;
        let mut fb = Feedback::new(&mesh, None, &FeedbackOptions::default()).unwrap();
        let realized = vec![0.0; mesh.n_faces()];
        assert!(
            !fb.observe(&mesh, &geometry, &u, &realized, 0.0, 0.0)
                .unwrap()
        );
        let report = fb.finish();
        assert_eq!(report.stop, FeedbackStop::Stable);
        assert!(report.rounds[0].max_abs_mismatch > 50.0);
        assert_eq!(report.rounds[0].max_abs_eligible, 0.0);
    }

    /// A closed mesh where every weight would go: refused, not a solve without weights.
    #[test]
    fn a_component_keeps_some_weight() {
        let ico = icosphere(0).unwrap(); // 12 vertices
        let mesh = TriMesh::new(ico.positions.len(), ico.faces.clone()).unwrap();
        let geometry = Geometry::from_positions(&mesh, &ico.positions).unwrap();
        let options = FeedbackOptions {
            max_reweighted_fraction: 1.0,
            ..FeedbackOptions::default()
        };
        let mut initial = vec![0.0; 12];
        initial[0] = 1.0;
        let mut fb = Feedback::new(&mesh, Some(initial.as_slice()), &options).unwrap();
        let mut u = vec![0.0; 12];
        u[0] = 100.0;
        let realized = vec![0.0; mesh.n_faces()];
        assert!(
            !fb.observe(&mesh, &geometry, &u, &realized, 0.0, 0.0)
                .unwrap()
        );
        assert_eq!(
            fb.finish().stop,
            FeedbackStop::NoWeightLeft { component: 0 }
        );
    }

    /// The probe applies its picks but never ends the loop, does not count as a solve, and must
    /// come first; picks beyond the limit are not applied.
    #[test]
    fn the_cone_free_probe() {
        let (mesh, geometry, _) = grid(20);
        let nv = mesh.n_vertices();
        let tip = 10 * 21 + 10;
        let realized = vec![0.0; mesh.n_faces()];
        let mut u = vec![0.0; nv];
        u[tip] = -1.0e3;
        let one = FeedbackOptions {
            max_rounds: 1,
            ..FeedbackOptions::default()
        };
        let mut fb = Feedback::new(&mesh, None, &one).unwrap();
        assert!(
            fb.observe_probe(&mesh, &geometry, &u, &realized, 0.0, 0.0)
                .unwrap()
        );
        assert!(!fb.is_finished());
        assert_eq!(fb.weights()[tip], 0.0);
        // The first solve agrees: stable, with the one solve that max_rounds = 1 allows.
        let flat = vec![0.0; nv];
        assert!(
            !fb.observe(&mesh, &geometry, &flat, &realized, 0.0, 0.0)
                .unwrap()
        );
        let report = fb.finish();
        assert_eq!(report.stop, FeedbackStop::Stable);
        assert_eq!(report.reweighted, vec![tip as u32]);
        assert!(report.rounds[0].probe && !report.rounds[1].probe);
        // Nothing found: nothing applied, the loop goes on.
        let mut fb = Feedback::new(&mesh, None, &FeedbackOptions::default()).unwrap();
        assert!(
            !fb.observe_probe(&mesh, &geometry, &flat, &realized, 0.0, 0.0)
                .unwrap()
        );
        assert!(!fb.is_finished());
        // Too many: nothing applied, the loop goes on.
        let mut many = vec![0.0; nv];
        for k in 0..5 {
            many[(4 + 3 * k) * 21 + 10] = 100.0;
        }
        let mut fb = Feedback::new(&mesh, None, &FeedbackOptions::default()).unwrap();
        assert!(
            !fb.observe_probe(&mesh, &geometry, &many, &realized, 0.0, 0.0)
                .unwrap()
        );
        assert!(!fb.is_finished() && fb.weights().iter().all(|&w| w == 1.0));
        // Not after a round.
        assert!(
            fb.observe(&mesh, &geometry, &flat, &realized, 0.0, 0.0)
                .is_ok()
        );
        let mut fb = Feedback::new(&mesh, None, &FeedbackOptions::default()).unwrap();
        fb.observe(&mesh, &geometry, &many, &realized, 0.0, 0.0)
            .unwrap();
        assert!(matches!(
            fb.observe_probe(&mesh, &geometry, &u, &realized, 0.0, 0.0),
            Err(ConeError::InvalidArgument { .. })
        ));
    }

    #[test]
    fn observing_after_the_end_is_an_error() {
        let (mesh, geometry, _) = grid(4);
        let mut fb = Feedback::new(&mesh, None, &FeedbackOptions::default()).unwrap();
        let u = vec![0.0; mesh.n_vertices()];
        let r = vec![0.0; mesh.n_faces()];
        assert!(!fb.observe(&mesh, &geometry, &u, &r, 0.0, 0.0).unwrap());
        assert!(matches!(
            fb.observe(&mesh, &geometry, &u, &r, 0.0, 0.0),
            Err(ConeError::InvalidArgument { .. })
        ));
    }
}
