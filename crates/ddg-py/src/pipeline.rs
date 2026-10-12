//! Where the binding calls the unwrap pipeline.
//!
//! [`run_pipeline`] is the **only** function that knows which pipeline is linked in. Today it
//! calls the stand-in [`crate::placeholder::planar_projection`]. To wire in the real pipeline:
//!
//! 1. add `ddg-uv` to `[dependencies]` in `Cargo.toml`;
//! 2. replace the body of [`run_pipeline`] by a call into `ddg-uv` that builds its input from
//!    the accessors of [`UnwrapInput`], and converts its result and error into
//!    [`UnwrapOutput`] and [`UnwrapError`];
//! 3. set [`PIPELINE_NAME`] and [`IS_PLACEHOLDER`];
//! 4. delete `src/placeholder.rs` and its `mod` line.
//!
//! Everything around it stays: [`unwrap`] validates the options, checks the pipeline's output
//! against the input, measures the distortion of the result and assembles the report.

use std::time::Instant;

use crate::contract::{Report, UnwrapInput, UnwrapOptions, UnwrapOutput, check_output};
use crate::error::UnwrapError;
use crate::measure::measure_distortion;

/// Name of the linked pipeline: `report["pipeline"]` and `ddg_tools.PIPELINE`.
pub const PIPELINE_NAME: &str = "placeholder-planar-projection";

/// Whether the linked pipeline is the stand-in: `report["placeholder"]` and
/// `ddg_tools.IS_PLACEHOLDER`.
pub const IS_PLACEHOLDER: bool = true;

/// Runs the linked pipeline. **This is the seam**: see the [module documentation](self).
///
/// A pipeline must
///
/// * return one UV per polygon corner of `input`, in corner order, all finite;
/// * return the chart of every polygon;
/// * return the cones it used and the polygon edges it cut along (including honoured seams);
/// * report unusable meshes as [`UnwrapError::Mesh`] with a message in the caller's vertex and
///   polygon indices, and numerical failure as [`UnwrapError::Solve`];
/// * never panic and never modify the input.
pub fn run_pipeline(
    input: &UnwrapInput,
    options: &UnwrapOptions,
) -> Result<UnwrapOutput, UnwrapError> {
    // ---------------------------------------------------------------------------------------
    // STAND-IN. Replace this line by the call into `ddg-uv`.
    crate::placeholder::planar_projection(input, options)
    // ---------------------------------------------------------------------------------------
}

/// The complete Rust side of `ddg_tools.unwrap`: options check, pipeline, output check,
/// distortion measurement, report.
///
/// The returned report contains, in this order: `pipeline`, `placeholder`; the sizes of the
/// input (`n_vertices`, `n_polygons`, `n_corners`, `n_triangles`); the steering as received
/// (`n_seam_edges`, `n_pinned`, `n_prescribed`, `n_forbidden`, `lambda`,
/// `positive_cones_only`); the result counts (`n_charts`, `n_cones`, `n_cut_edges`); whatever
/// the pipeline reported; the measured distortion (see [`crate::measure`]); and the timings
/// `seconds_pipeline` and `seconds_measure`.
///
/// # Errors
/// [`UnwrapError::Input`] for bad options, and whatever the pipeline returns.
pub fn unwrap(input: &UnwrapInput, options: &UnwrapOptions) -> Result<UnwrapOutput, UnwrapError> {
    options.validate()?;

    let start = Instant::now();
    let mut output = run_pipeline(input, options)?;
    let seconds_pipeline = start.elapsed().as_secs_f64();
    check_output(input, &output)?;

    let start = Instant::now();
    let distortion = measure_distortion(input, &output.uvs);
    let seconds_measure = start.elapsed().as_secs_f64();

    let mut charts = output.polygon_charts.clone();
    charts.sort_unstable();
    charts.dedup();

    let mut report = Report::new();
    report.set("pipeline", PIPELINE_NAME);
    report.set("placeholder", IS_PLACEHOLDER);
    report.set("n_vertices", input.n_vertices());
    report.set("n_polygons", input.n_polygons());
    report.set("n_corners", input.n_corners());
    report.set("n_triangles", input.n_triangles());
    report.set("n_seam_edges", input.seam_edges().len());
    report.set("n_pinned", input.pinned_vertices().len());
    report.set("n_prescribed", input.prescribed_cones().len());
    report.set("n_forbidden", input.forbidden_vertices().len());
    report.set("lambda", options.lambda);
    report.set("positive_cones_only", options.positive_cones_only);
    report.set("n_charts", charts.len());
    report.set("n_cones", output.cones.len());
    report.set("n_cut_edges", output.cut_edges.len());
    report.extend(std::mem::take(&mut output.report));
    distortion.write_to(&mut report);
    report.set("seconds_pipeline", seconds_pipeline);
    report.set("seconds_measure", seconds_measure);
    output.report = report;
    Ok(output)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::contract::test_meshes::*;
    use crate::contract::{Cone, ReportValue};

    #[test]
    fn report_has_the_documented_keys_in_order() {
        let mut b = cube();
        b.seam_edges = vec![4, 5];
        b.pinned_vertices = vec![3];
        b.forbidden_vertices = vec![0, 1];
        let input = UnwrapInput::from_buffers(b).unwrap();
        let out = unwrap(&input, &UnwrapOptions::default()).unwrap();
        let keys: Vec<&str> = out
            .report
            .entries()
            .iter()
            .map(|(k, _)| k.as_str())
            .collect();
        assert_eq!(
            keys,
            [
                "pipeline",
                "placeholder",
                "n_vertices",
                "n_polygons",
                "n_corners",
                "n_triangles",
                "n_seam_edges",
                "n_pinned",
                "n_prescribed",
                "n_forbidden",
                "lambda",
                "positive_cones_only",
                "n_charts",
                "n_cones",
                "n_cut_edges",
                "n_unreferenced_vertices",
                "n_seam_edges_ignored",
                "n_fold_edges",
                "n_degenerate_triangles",
                "n_collapsed_uv_triangles",
                "n_flipped_uv_triangles",
                "measured_triangles",
                "measured_angle_error_max",
                "measured_angle_error_mean",
                "measured_log_scale_std",
                "measured_log_scale_min",
                "measured_log_scale_max",
                "seconds_pipeline",
                "seconds_measure",
            ]
        );
        let get = |k: &str| out.report.get(k).cloned().unwrap();
        assert_eq!(get("pipeline"), ReportValue::Text(PIPELINE_NAME.to_owned()));
        assert_eq!(get("placeholder"), ReportValue::Bool(true));
        assert_eq!(get("n_vertices"), ReportValue::Int(8));
        assert_eq!(get("n_corners"), ReportValue::Int(24));
        assert_eq!(get("n_seam_edges"), ReportValue::Int(1));
        assert_eq!(get("n_pinned"), ReportValue::Int(1));
        assert_eq!(get("n_forbidden"), ReportValue::Int(2));
        assert_eq!(get("lambda"), ReportValue::Float(0.5));
        assert_eq!(get("n_charts"), ReportValue::Int(1));
        assert_eq!(get("n_cones"), ReportValue::Int(1));
        // Seam (4, 5) plus the four fold edges of the bottom face.
        assert_eq!(get("n_cut_edges"), ReportValue::Int(5));
        assert_eq!(
            out.cones,
            vec![Cone {
                vertex: 3,
                angle: 0.0
            }]
        );
        // The projection of a cube along z collapses the four side faces.
        assert_eq!(get("n_collapsed_uv_triangles"), ReportValue::Int(8));
        assert_eq!(get("n_flipped_uv_triangles"), ReportValue::Int(2));
        assert_eq!(get("measured_triangles"), ReportValue::Int(4));
        assert_eq!(out.report.warnings().len(), 1);
    }

    #[test]
    fn bad_options_and_bad_meshes_are_errors() {
        let input = UnwrapInput::from_buffers(cube()).unwrap();
        let options = UnwrapOptions {
            lambda: f64::INFINITY,
            ..UnwrapOptions::default()
        };
        assert!(matches!(
            unwrap(&input, &options),
            Err(UnwrapError::Input(_))
        ));
        let input = UnwrapInput::from_buffers(three_fins()).unwrap();
        assert_eq!(
            unwrap(&input, &UnwrapOptions::default())
                .unwrap_err()
                .to_string(),
            "edge (0, 1) is shared by 3 faces: 0, 1, 2"
        );
    }

    #[test]
    fn results_are_deterministic() {
        let mut b = cube();
        b.seam_edges = vec![0, 1, 2, 3, 4, 5];
        let input = UnwrapInput::from_buffers(b).unwrap();
        let strip_times = |mut out: UnwrapOutput| {
            out.report.set("seconds_pipeline", 0.0);
            out.report.set("seconds_measure", 0.0);
            out
        };
        let a = strip_times(unwrap(&input, &UnwrapOptions::default()).unwrap());
        let b = strip_times(unwrap(&input, &UnwrapOptions::default()).unwrap());
        assert_eq!(a, b);
    }
}
