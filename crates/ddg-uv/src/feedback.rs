//! Stage 3c (optional): the feedback loop of the cone solver.
//!
//! The cone solver measures distortion with the **linearized** log scale factor `u`. Where the
//! linearization is poor, typically at the tip of a thin spike, it optimizes a distortion that
//! the flattening does not have and spends its cones there (`ddg-optimal-cones`, README,
//! "Needle-like spikes"). With [`UnwrapOptions::cone_feedback`] set, the pipeline compares, on
//! every piece, the solver's `u` with the log area factor that the flattening of stage 3b
//! realizes, lowers the distortion weight of the vertices where they are far apart, and solves
//! and flattens that piece again, until nothing changes. The rule is
//! [`ddg_optimal_cones::Feedback`]; this module drives it with the pipeline's own solve
//! (policies, merge guard, slit) and flattening.
//!
//! On a piece with a boundary the loop first looks at the piece without cones (the cone-free
//! probe: the state the spike check evaluated, flattened), because a cone that the solver puts at
//! a needle tip hides the needle afterwards. Where nothing is reweighted the loop costs that probe
//! and the comparison only (no extra solve): the first round compares the solve and the
//! flattening that stages 3a and 3b made anyway.
//! Vertices that the spike check ([`crate::spikes`]) already excluded start with weight zero
//! and are not looked at.
//!
//! The crate README ("The cone feedback loop") and the cone crate's README ("The feedback
//! loop") hold the measurements.

use std::time::Instant;

use ddg_bff::ConeFlattenOptions;
use ddg_bff::quality::{log_scale_factor, signed_double_area};
use ddg_core::Geometry;
use ddg_core::io::Cone;
use ddg_optimal_cones::{
    ConeProblem, Feedback, FeedbackReport, FeedbackStop, Parallelism, evaluate_cones,
    every_component_has_boundary,
};

use crate::cones::{ConeStage, PieceCones, piece_solve, piece_weights};
use crate::error::{UvError, internal};
use crate::flatten::{FlatStage, flatten_options, flatten_piece};
use crate::options::UnwrapOptions;
use crate::prepare::{Piece, Prepared};

/// What the feedback loop did on one chart.
#[derive(Debug, Clone, PartialEq)]
pub struct ChartFeedback {
    /// The vertices whose distortion weight was lowered, as vertices of the caller, ascending.
    pub reweighted: Vec<u32>,
    /// The vertices the last round picked but that were not reweighted because the loop ended
    /// (see [`FeedbackStop`]), as vertices of the caller, ascending.
    pub pending: Vec<u32>,
    /// The report of the loop, in the vertex numbering of the piece's mesh (the piece with its
    /// slit, if it has one): rounds, the mismatch per round, the final weights, why it ended.
    pub report: FeedbackReport,
}

/// The realized log area factor `½ ln(|B_t| / A_t)` of every face of a flattening
/// (`corner_uvs[3 t + k]` is corner `k` of face `t`), or NaN where the face or its image has no
/// area or the value is not finite. The absolute value of the image's area is used, so that a
/// flipped face still says how much it was scaled.
///
/// [formula: face-log-scale-factor]
pub fn realized_face_log_scale(geometry: &Geometry, corner_uvs: &[[f64; 2]]) -> Vec<f64> {
    geometry
        .face_areas()
        .iter()
        .enumerate()
        .map(|(f, &a)| {
            let Some(q) = corner_uvs.get(3 * f..3 * f + 3) else {
                return f64::NAN;
            };
            let b = 0.5 * signed_double_area(q[0], q[1], q[2]).abs();
            let r = log_scale_factor(a, b);
            if a > 0.0 && b > 0.0 && r.is_finite() {
                r
            } else {
                f64::NAN
            }
        })
        .collect()
}

/// Runs the feedback loop on every piece that was solved, if [`UnwrapOptions::cone_feedback`]
/// is set (otherwise nothing happens). `cones` and `flat` must be the results of
/// [`crate::solve_cones`] and [`crate::flatten()`] for `prepared`; the pieces that the loop solves
/// again are replaced in both. See [`crate::unwrap`].
///
/// A solve or flattening that fails after a reweighting is not an error: the piece keeps its
/// previous cones and flattening, and the report and a warning say so.
///
/// # Errors
/// [`UvError::InvalidOption`] for invalid options, [`UvError::Internal`] if the stages do not
/// belong to `prepared`.
pub fn cone_feedback(
    prepared: &Prepared,
    cones: &mut ConeStage,
    flat: &mut FlatStage,
    options: &UnwrapOptions,
) -> Result<(), UvError> {
    let Some(feedback_options) = options.cone_feedback else {
        return Ok(());
    };
    options.validate()?;
    let n = prepared.pieces.len();
    if cones.pieces.len() != n || flat.pieces.len() != n || cones.piece_warnings.len() != n {
        return Err(internal(
            "the stages passed to the feedback loop were not computed from this input",
        ));
    }
    let start = Instant::now();
    let flattener = flatten_options(options);
    for (k, piece) in prepared.pieces.iter().enumerate() {
        let (Some(first_log_scale), Some(first_solve)) = (
            cones.pieces[k].log_scale.clone(),
            cones.pieces[k].solve.as_ref(),
        ) else {
            continue;
        };
        let mut seconds_solve = first_solve.seconds;
        // The mesh that was solved and flattened; a slit is cut the same way in every round.
        let slit = cones.pieces[k].slit_piece.clone();
        let on = slit.as_ref().unwrap_or(piece);
        if (0..on.mesh.n_vertices()).all(|v| on.mesh.is_boundary_vertex(v)) {
            continue;
        }
        // The weights of the first solve: zero at the vertices the spike check excluded.
        let spike_weights: Vec<(u32, f64)> = cones.pieces[k]
            .spikes
            .as_ref()
            .map(|s| s.excluded.iter().map(|&v| (v, 0.0)).collect())
            .unwrap_or_default();
        let initial = piece_weights(on, &spike_weights);
        let mut state = Feedback::new(&on.mesh, initial.as_deref(), &feedback_options)
            .map_err(|e| internal(format!("the feedback loop could not be set up: {e}")))?;
        let mut log_scale = first_log_scale;
        let mut realized = realized_face_log_scale(&on.geometry, &flat.pieces[k].corner_uvs);
        let mut seconds_flatten = flat.pieces[k].seconds;
        // The cone-free probe: the piece without cones (the prescribed ones only), flattened.
        let mut probe_note = None;
        let mut out_of_date = false;
        if feedback_options.cone_free_probe {
            if every_component_has_boundary(&on.mesh) {
                match probe(prepared, piece, on, &slit, &cones.pieces[k], &flattener) {
                    Ok((u0, faces, seconds_evaluate, seconds_probe_flatten)) => {
                        out_of_date = state
                            .observe_probe(
                                &on.mesh,
                                &on.geometry,
                                &u0,
                                &faces,
                                seconds_evaluate,
                                seconds_probe_flatten,
                            )
                            .map_err(|e| internal(format!("the feedback loop failed: {e}")))?;
                    }
                    Err(message) => probe_note = Some(message),
                }
            } else {
                probe_note = Some("the piece is closed".to_owned());
            }
        }
        loop {
            if out_of_date {
                // The weights, per vertex of the caller, where they are not 1. Interior vertices
                // of `on` are distinct vertices of the caller; the two copies of a vertex on a
                // slit are boundary vertices, which the loop never reweights, so they carry the
                // same weight.
                let mut weights: Vec<(u32, f64)> = (0..on.mesh.n_vertices())
                    .filter(|&v| state.weights()[v] != 1.0)
                    .map(|v| (on.vertex_input[v], state.weights()[v]))
                    .collect();
                weights.sort_unstable_by_key(|p| p.0);
                weights.dedup_by_key(|p| p.0);
                let solve_start = Instant::now();
                let mut solve_warnings = Vec::new();
                let solved =
                    match piece_solve(prepared, piece, options, &weights, &mut solve_warnings) {
                        Ok(solved) => solved,
                        Err(e) => {
                            state.round_failed(format!("the cone solve failed: {e}"));
                            break;
                        }
                    };
                let same_mesh = solved.slit_piece.as_ref().map(|p| p.mesh.n_vertices())
                    == slit.as_ref().map(|p| p.mesh.n_vertices());
                let Some(next_log_scale) = solved.log_scale.clone().filter(|_| same_mesh) else {
                    state.round_failed("the piece was not solved on the same mesh again");
                    break;
                };
                seconds_solve = solve_start.elapsed().as_secs_f64();
                let mut flatten_warnings = Vec::new();
                let flattened =
                    match flatten_piece(piece, &solved, &flattener, &mut flatten_warnings) {
                        Ok(flattened) => flattened,
                        Err(e) => {
                            state.round_failed(format!("the flattening failed: {e}"));
                            break;
                        }
                    };
                seconds_flatten = flattened.seconds;
                realized = realized_face_log_scale(&on.geometry, &flattened.corner_uvs);
                log_scale = next_log_scale;
                let spikes = cones.pieces[k].spikes.take();
                let cone_free_log_scale = cones.pieces[k].cone_free_log_scale.take();
                cones.pieces[k] = PieceCones {
                    spikes,
                    cone_free_log_scale,
                    ..solved
                };
                cones.piece_warnings[k].solve = solve_warnings;
                flat.pieces[k] = flattened;
                flat.piece_warnings[k] = flatten_warnings;
            }
            let again = state
                .observe(
                    &on.mesh,
                    &on.geometry,
                    &log_scale,
                    &realized,
                    seconds_solve,
                    seconds_flatten,
                )
                .map_err(|e| internal(format!("the feedback loop failed: {e}")))?;
            if !again {
                break;
            }
            out_of_date = true;
        }
        let mut report = state.finish();
        report.probe_note = probe_note;
        let caller = |vertices: &[u32]| {
            let mut out: Vec<u32> = vertices
                .iter()
                .map(|&v| on.vertex_input[v as usize])
                .collect();
            out.sort_unstable();
            out
        };
        let chart = ChartFeedback {
            reweighted: caller(&report.reweighted),
            pending: caller(&report.pending),
            report,
        };
        if let Some(warning) = warning(piece.first_polygon, on.mesh.n_vertices(), &chart) {
            cones.feedback_warnings.push(warning);
        }
        cones.pieces[k].feedback = Some(chart);
    }
    cones.seconds_feedback = start.elapsed().as_secs_f64();
    Ok(())
}

/// The cone-free probe of one piece: the scale factor without cones (the prescribed cones of
/// the piece only, `u = 0` on the boundary) and the realized log area factor of its
/// flattening, with the seconds of both. An error message if either fails.
fn probe(
    prepared: &Prepared,
    piece: &Piece,
    on: &Piece,
    slit: &Option<Piece>,
    solved: &PieceCones,
    flattener: &ConeFlattenOptions,
) -> Result<(Vec<f64>, Vec<f64>, f64, f64), String> {
    let start = Instant::now();
    let n = on.mesh.n_vertices();
    let prescribed: Vec<Cone> = (0..n)
        .filter(|&v| !on.mesh.is_boundary_vertex(v))
        .filter(|&v| prepared.steering_role[on.vertex_input[v] as usize] == 2)
        .map(|v| Cone {
            vertex: v as u32,
            angle: prepared.prescribed_angle[on.vertex_input[v] as usize],
        })
        .collect();
    // The spike check evaluated the same state on the piece without slit; reuse it.
    let u0 = match (&solved.cone_free_log_scale, slit) {
        (Some(u0), None) => u0.clone(),
        _ => {
            let problem = ConeProblem::new(&on.mesh, &on.geometry, 0.0);
            evaluate_cones(&problem, &prescribed, Parallelism::Sequential)
                .map_err(|e| format!("the cone-free state failed: {e}"))?
                .log_scale
        }
    };
    let seconds_evaluate = start.elapsed().as_secs_f64();
    let cone_free = PieceCones {
        slit_piece: slit.clone(),
        slit_edges: Vec::new(),
        cones: prescribed,
        log_scale: Some(u0.clone()),
        boundary_optimized: false,
        closed: solved.closed,
        positive_only: solved.positive_only,
        solve: None,
        spikes: None,
        feedback: None,
        cone_free_log_scale: None,
    };
    let mut ignored = Vec::new();
    let flattened = flatten_piece(piece, &cone_free, flattener, &mut ignored)
        .map_err(|e| format!("the cone-free flattening failed: {e}"))?;
    let faces = realized_face_log_scale(&on.geometry, &flattened.corner_uvs);
    Ok((u0, faces, seconds_evaluate, flattened.seconds))
}

/// The warning for a chart where the loop changed something or could not finish.
fn warning(first_polygon: u32, n_vertices: usize, chart: &ChartFeedback) -> Option<String> {
    const LISTED: usize = 16;
    let list = |vertices: &[u32]| {
        let mut text: Vec<String> = vertices.iter().take(LISTED).map(u32::to_string).collect();
        if vertices.len() > LISTED {
            text.push(format!("and {} more", vertices.len() - LISTED));
        }
        text.join(", ")
    };
    let report = &chart.report;
    let largest = report.rounds.first().map_or(0.0, |r| r.max_abs_eligible);
    let mut parts = Vec::new();
    if !chart.reweighted.is_empty() {
        parts.push(format!(
            "the cone feedback lowered the distortion weight of {} of the {n_vertices} vertices \
             of the piece of the surface that contains face {first_polygon} (vertices {}): there \
             the cone solver's linearized scale factor was off by up to {largest:.3e} (log \
             units) from the one the flattening realized, the signature of a needle-like \
             spike; {} solves",
            chart.reweighted.len(),
            list(&chart.reweighted),
            report.rounds.len()
        ));
    }
    let pending = list(&chart.pending);
    match &report.stop {
        FeedbackStop::Stable | FeedbackStop::NotContinued => {}
        FeedbackStop::RoundLimit => parts.push(format!(
            "the cone feedback on the piece that contains face {first_polygon} stopped at its \
             round limit with {} vertices still far from the flattening (vertices {pending})",
            chart.pending.len()
        )),
        FeedbackStop::TooManyVertices { wanted, allowed } => parts.push(format!(
            "the cone feedback on the piece that contains face {first_polygon} found {wanted} \
             vertices where the linearized scale factor is far from the flattening, more than \
             the {allowed} it may reweight; it left them alone (vertices {pending})"
        )),
        FeedbackStop::NoWeightLeft { .. } => parts.push(format!(
            "the cone feedback on the piece that contains face {first_polygon} would have left \
             it without distortion weight; it left vertices {pending} alone"
        )),
        FeedbackStop::RoundFailed { message, .. } => parts.push(format!(
            "the cone feedback on the piece that contains face {first_polygon} kept its \
             previous cones: {message}"
        )),
    }
    (!parts.is_empty()).then(|| parts.join("; "))
}
