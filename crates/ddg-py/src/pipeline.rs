//! Where the binding calls the unwrap pipeline.
//!
//! [`run_pipeline`] is the **only** function that knows which pipeline is linked in: it calls
//! [`ddg_uv::unwrap`] (seams, optimal cones, cut, boundary first flattening) and converts its
//! result and errors into [`UnwrapOutput`] and [`UnwrapError`].
//!
//! Everything around it stays independent of the pipeline: [`unwrap`] validates the options,
//! checks the pipeline's output against the input, measures the distortion of the result and
//! assembles the report.

use std::time::Instant;

use crate::contract::{Cone, Report, UnwrapInput, UnwrapOptions, UnwrapOutput, check_output};
use crate::error::{InputError, UnwrapError};
use crate::measure::measure_distortion;

/// Name of the linked pipeline: `report["pipeline"]` and `ddg_tools.PIPELINE`.
pub const PIPELINE_NAME: &str = "ddg-uv";

/// Whether the linked pipeline is the stand-in: `report["placeholder"]` and
/// `ddg_tools.IS_PLACEHOLDER`.
pub const IS_PLACEHOLDER: bool = false;

/// Runs the linked pipeline, [`ddg_uv::unwrap`]. **This is the seam**: see the [module documentation](self).
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
    let prescribed: Vec<ddg_uv::Cone> = input
        .prescribed_cones()
        .iter()
        .map(|c| ddg_uv::Cone {
            vertex: c.vertex,
            angle: c.angle,
        })
        .collect();
    let mesh = ddg_uv::MeshInput {
        positions: input.positions(),
        corner_vertices: input.corner_vertices(),
        polygon_offsets: input.polygon_offsets(),
        triangles: input.triangles(),
        seam_edges: input.seam_edges(),
        pinned_vertices: input.pinned_vertices(),
        prescribed_cones: &prescribed,
        forbidden_vertices: input.forbidden_vertices(),
    };
    let uv_options = ddg_uv::UnwrapOptions {
        lambda: options.lambda,
        cone_signs: if options.positive_cones_only {
            ddg_uv::ConeSigns::PositiveWhereFeasible
        } else {
            ddg_uv::ConeSigns::Any
        },
        boundary: if options.optimize_boundary {
            ddg_uv::BoundaryValues::Optimized
        } else {
            ddg_uv::BoundaryValues::Zero
        },
        // The contract of `ddg_tools.unwrap` is "not packed": `ddg_tools.pack_charts` packs.
        packing: None,
        split_overlaps: options.split_overlaps.then(|| ddg_uv::SplitOptions {
            max_splits: options.split_max,
            min_overlap_fraction: options.split_threshold,
            ..ddg_uv::SplitOptions::default()
        }),
        cone_feedback: options.cone_feedback.then(|| ddg_uv::FeedbackOptions {
            threshold: options.feedback_threshold,
            max_rounds: options.feedback_max_rounds,
            ..ddg_uv::FeedbackOptions::default()
        }),
        ..ddg_uv::UnwrapOptions::default()
    };
    let map = ddg_uv::unwrap(&mesh, &uv_options).map_err(|e| match e {
        ddg_uv::UvError::InvalidOption {
            name,
            value,
            requirement,
        } => UnwrapError::Input(InputError::BadOption {
            name,
            value,
            requirement,
        }),
        e => match e.kind() {
            ddg_uv::ErrorKind::Mesh => UnwrapError::Mesh(e.to_string()),
            _ => UnwrapError::Solve(e.to_string()),
        },
    })?;
    let mut report = Report::new();
    for (key, value) in map.report.summary() {
        match value {
            ddg_uv::Value::Bool(v) => report.set(key, v),
            ddg_uv::Value::Int(v) => report.set(key, v),
            ddg_uv::Value::Float(v) => report.set(key, v),
            ddg_uv::Value::Text(v) => report.set(key, v),
        }
    }
    for warning in &map.report.warnings {
        report.warn(warning.clone());
    }
    Ok(UnwrapOutput {
        uvs: map.uvs,
        polygon_charts: map.polygon_charts,
        cones: map
            .cones
            .iter()
            .map(|c| Cone {
                vertex: c.vertex,
                angle: c.angle,
            })
            .collect(),
        cut_edges: map.cut_edges,
        report,
    })
}

/// The complete Rust side of `ddg_tools.unwrap`: options check, pipeline, output check,
/// distortion measurement, report.
///
/// The returned report contains, in this order: `pipeline`, `placeholder`; the sizes of the
/// input (`n_vertices`, `n_polygons`, `n_corners`, `n_triangles`); the steering as received
/// (`n_seam_edges`, `n_pinned`, `n_prescribed`, `n_forbidden`, `lambda`,
/// `positive_cones_only`, `optimize_boundary`); the result counts (`n_charts`, `n_cones`, `n_cut_edges`); whatever
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
    report.set("optimize_boundary", options.optimize_boundary);
    report.set("cone_feedback", options.cone_feedback);
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
    use crate::contract::ReportValue;
    use crate::contract::test_meshes::*;

    /// Keys the binding itself writes, in order: before and after whatever the pipeline reports.
    const LEADING: [&str; 17] = [
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
        "optimize_boundary",
        "cone_feedback",
        "n_charts",
        "n_cones",
        "n_cut_edges",
    ];
    const TRAILING: [&str; 11] = [
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
    ];

    #[test]
    fn report_frames_the_pipeline_report_with_the_documented_keys() {
        let input = UnwrapInput::from_buffers(cube()).unwrap();
        let out = unwrap(&input, &UnwrapOptions::default()).unwrap();
        let keys: Vec<&str> = out
            .report
            .entries()
            .iter()
            .map(|(k, _)| k.as_str())
            .collect();
        assert_eq!(keys[..LEADING.len()], LEADING);
        assert_eq!(keys[keys.len() - TRAILING.len()..], TRAILING);
        // No key is written twice (the pipeline's keys must not collide with the binding's).
        let mut sorted = keys.clone();
        sorted.sort_unstable();
        sorted.dedup();
        assert_eq!(sorted.len(), keys.len());

        let get = |k: &str| out.report.get(k).cloned().unwrap();
        assert_eq!(get("pipeline"), ReportValue::Text(PIPELINE_NAME.to_owned()));
        assert_eq!(get("placeholder"), ReportValue::Bool(false));
        assert_eq!(get("n_vertices"), ReportValue::Int(8));
        assert_eq!(get("n_corners"), ReportValue::Int(24));
        assert_eq!(get("lambda"), ReportValue::Float(0.1));
        assert_eq!(get("positive_cones_only"), ReportValue::Bool(true));
        assert_eq!(get("n_charts"), ReportValue::Int(1));
    }

    /// A cube is developable: its eight corners are cones of `π/2`, and cutting through them
    /// unfolds it without any distortion.
    #[test]
    fn a_cube_unfolds_isometrically() {
        let input = UnwrapInput::from_buffers(cube()).unwrap();
        let out = unwrap(&input, &UnwrapOptions::default()).unwrap();
        assert_eq!(out.cones.len(), 8);
        for cone in &out.cones {
            assert!(
                (cone.angle - std::f64::consts::FRAC_PI_2).abs() < 1e-9,
                "{cone:?}"
            );
        }
        // A spanning tree of the eight corners.
        assert_eq!(out.cut_edges.len(), 7);
        let float = |k: &str| match out.report.get(k) {
            Some(ReportValue::Float(v)) => *v,
            other => panic!("{k}: {other:?}"),
        };
        assert!(float("measured_angle_error_max") < 1e-9);
        assert!(float("measured_log_scale_std") < 1e-9);
        assert_eq!(
            out.report.get("n_flipped_uv_triangles"),
            Some(&ReportValue::Int(0))
        );
    }

    #[test]
    fn steering_reaches_the_pipeline() {
        // Seams around the face `z = 0` turn the cube into an open box plus its lid: two charts,
        // the seams among the cuts.
        let mut b = cube();
        b.seam_edges = vec![0, 2, 2, 3, 3, 1, 1, 0];
        let input = UnwrapInput::from_buffers(b).unwrap();
        let seams: Vec<[u32; 2]> = input.seam_edges().to_vec();
        let out = unwrap(&input, &UnwrapOptions::default()).unwrap();
        for seam in &seams {
            let mut edge = *seam;
            edge.sort_unstable();
            assert!(out.cut_edges.contains(&edge), "{edge:?} not cut");
        }
        let mut charts = out.polygon_charts.clone();
        charts.sort_unstable();
        charts.dedup();
        assert_eq!(charts.len(), 2);

        // A forbidden vertex gets no cone.
        let mut b = cube();
        b.forbidden_vertices = vec![0];
        let input = UnwrapInput::from_buffers(b).unwrap();
        let out = unwrap(&input, &UnwrapOptions::default()).unwrap();
        assert!(out.cones.iter().all(|c| c.vertex != 0));
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
        let error = unwrap(&input, &UnwrapOptions::default()).unwrap_err();
        assert!(matches!(error, UnwrapError::Mesh(_)));
        assert_eq!(
            error.to_string(),
            "edge (0, 1) is shared by 3 faces: 0, 1, 2"
        );
    }

    #[test]
    fn the_cone_feedback_reaches_the_pipeline() {
        let input = UnwrapInput::from_buffers(cube()).unwrap();
        assert!(UnwrapOptions::default().cone_feedback);
        let off_options = UnwrapOptions {
            cone_feedback: false,
            ..UnwrapOptions::default()
        };
        let off = unwrap(&input, &off_options).unwrap();
        assert_eq!(
            off.report.get("cone_feedback"),
            Some(&ReportValue::Bool(false))
        );
        assert!(off.report.get("cone_feedback_rounds").is_none());
        let on = unwrap(&input, &UnwrapOptions::default()).unwrap();
        assert_eq!(
            on.report.get("cone_feedback"),
            Some(&ReportValue::Bool(true))
        );
        // A cube has nothing to reweight: one round, the same cones.
        assert_eq!(
            on.report.get("cone_feedback_rounds"),
            Some(&ReportValue::Int(1))
        );
        assert_eq!(
            on.report.get("n_cone_feedback_reweighted"),
            Some(&ReportValue::Int(0))
        );
        assert_eq!(on.cones, off.cones);
        let bad = UnwrapOptions {
            cone_feedback: true,
            feedback_max_rounds: 0,
            ..UnwrapOptions::default()
        };
        assert!(matches!(unwrap(&input, &bad), Err(UnwrapError::Input(_))));
    }

    #[test]
    fn results_are_deterministic() {
        let mut b = cube();
        b.seam_edges = vec![0, 2, 2, 3];
        let input = UnwrapInput::from_buffers(b).unwrap();
        let strip_times = |mut out: UnwrapOutput| {
            let keys: Vec<String> = out
                .report
                .entries()
                .iter()
                .map(|(k, _)| k.clone())
                .filter(|k| k.starts_with("seconds_"))
                .collect();
            for key in keys {
                out.report.set(&key, 0.0);
            }
            out
        };
        let a = strip_times(unwrap(&input, &UnwrapOptions::default()).unwrap());
        let b = strip_times(unwrap(&input, &UnwrapOptions::default()).unwrap());
        assert_eq!(a, b);
    }
}
