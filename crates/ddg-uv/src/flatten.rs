//! Stage 3b: cutting every piece through its cones along polygon edges, and flattening it.

use std::time::Instant;

use ddg_bff::quality::signed_double_area;
use ddg_bff::{BffError, ConeFlattenInput, ConeFlattenOptions, CutSource, flatten_with_cones};

use crate::cones::ConeStage;
use crate::error::{UvError, internal};
use crate::measure::Accumulator;
use crate::options::UnwrapOptions;
use crate::prepare::Prepared;
use crate::report::FlattenReport;

/// The flattening of one piece.
#[derive(Debug, Clone)]
pub struct PieceFlat {
    /// The position in the plane of every triangle corner of the piece: entry `3 f + k` is
    /// corner `k` of face `f`. Unpacked: the chart is where the flattener put it, at the
    /// scale `e^u`.
    pub corner_uvs: Vec<[f64; 2]>,
    /// The edges the planner cut (through cones, around handles), as vertex pairs of the
    /// caller, `[min, max]`, ascending. Seams and slits are not in this list.
    pub cut_edges: Vec<[u32; 2]>,
    /// What the flattener reported.
    pub report: FlattenReport,
    /// Overlapped and covered area of the chart in the plane, if estimated.
    pub(crate) coverage: Option<(f64, f64)>,
}

/// The result of stage 3b.
#[derive(Debug, Clone)]
pub struct FlatStage {
    /// One entry per piece.
    pub pieces: Vec<PieceFlat>,
    pub(crate) warnings: Vec<String>,
    pub(crate) seconds: f64,
}

/// Runs stage 3b on every piece. See [`crate::unwrap`].
///
/// # Errors
/// [`UvError::Flatten`] with the flattener's message, [`UvError::Internal`] if a piece does
/// not come out as exactly one chart.
pub fn flatten(
    prepared: &Prepared,
    cones: &ConeStage,
    options: &UnwrapOptions,
) -> Result<FlatStage, UvError> {
    let start = Instant::now();
    options.validate()?;
    if cones.pieces.len() != prepared.pieces.len() {
        return Err(internal(format!(
            "{} pieces but cones for {}",
            prepared.pieces.len(),
            cones.pieces.len()
        )));
    }
    let flatten_options = ConeFlattenOptions {
        overlap_samples_per_triangle: options.overlap_samples_per_triangle,
        ..ConeFlattenOptions::default()
    };
    let mut warnings = Vec::new();
    let mut pieces = Vec::with_capacity(prepared.pieces.len());
    for (piece, piece_cones) in prepared.pieces.iter().zip(&cones.pieces) {
        let on = piece_cones.slit_piece.as_ref().unwrap_or(piece);
        let run = |scale_factors: Option<&[f64]>| {
            flatten_with_cones(
                &ConeFlattenInput {
                    mesh: &on.mesh,
                    geometry: &on.geometry,
                    cones: &piece_cones.cones,
                    scale_factors,
                    cut: CutSource::Plan {
                        may_cut: Some(&on.is_polygon_edge),
                        must_cut: &[],
                    },
                    anchors: &[],
                },
                &flatten_options,
            )
        };
        let flat = match run(piece_cones.log_scale.as_deref()) {
            Err(BffError::ScaleFactorsNotFlat {
                vertex, curvature, ..
            }) if piece_cones.log_scale.is_some() => {
                // The solver's scale factors leave curvature inside the chart: let the
                // flattener solve for its own.
                warnings.push(format!(
                    "the cone solver's scale factors leave a curvature of {curvature:.3e} at \
                     vertex {} away from the cones; the flattener solved for the scale factors \
                     again",
                    on.vertex_input[vertex as usize]
                ));
                run(None)
            }
            other => other,
        }
        .map_err(|e| on.flatten_error(&e))?;
        if flat.n_charts() != 1 || flat.charts.len() != 1 {
            return Err(internal(format!(
                "the piece that contains face {} was flattened into {} charts instead of one",
                piece.first_polygon,
                flat.n_charts()
            )));
        }
        let chart = &flat.charts[0];
        let lengths = on.geometry.edge_lengths();
        let mut cut_edges: Vec<[u32; 2]> = Vec::with_capacity(flat.cut_edges().len());
        let mut cut_length = 0.0;
        for &e in flat.cut_edges() {
            let [a, b] = on.mesh.edge_vertices(e as usize);
            let (a, b) = (on.vertex_input[a], on.vertex_input[b]);
            cut_edges.push([a.min(b), a.max(b)]);
            cut_length += lengths[e as usize];
        }
        cut_edges.sort_unstable();
        let mut accumulator = Accumulator::default();
        for (f, &area) in on.geometry.face_areas().iter().enumerate() {
            let q = &flat.corner_uvs[3 * f..3 * f + 3];
            accumulator.add(area, 0.5 * signed_double_area(q[0], q[1], q[2]));
        }
        let report = FlattenReport {
            n_cut_edges: cut_edges.len(),
            cut_length,
            flipped_triangles: chart.distortion.flipped_triangles,
            degenerate_triangles: chart.distortion.degenerate_triangles,
            wedges_over_2pi: chart.overlap.wedges_over_2pi,
            boundary_crossings: chart.overlap.boundary_crossings,
            overlapped_fraction: chart.overlap.overlapped_fraction(),
            quasi_conformal_max: chart.distortion.quasi_conformal.max,
            quasi_conformal_mean: chart.distortion.quasi_conformal.mean,
            closure_length_change: chart.closure.max_length_change,
            seam_length_mismatch: flat.seams.max_length_mismatch,
            surface_area: chart.distortion.surface_area,
            plane_area: chart.distortion.plane_area,
            area_distortion: accumulator.finish(),
        };
        pieces.push(PieceFlat {
            coverage: chart
                .overlap
                .coverage
                .map(|c| (c.overlapped_area, c.covered_area)),
            corner_uvs: flat.corner_uvs,
            cut_edges,
            report,
        });
    }
    Ok(FlatStage {
        pieces,
        warnings,
        seconds: start.elapsed().as_secs_f64(),
    })
}
