//! Reads OBJ files, flips each to its intrinsic Delaunay triangulation, traces the intrinsic
//! edges over the input mesh, and prints mesh quality before and after with timings.
//!
//! ```text
//! cargo run --release -p ddg-intrinsic --example delaunay_report -- [options] a.obj b.obj ...
//!   --reference        use the Delaunay test of the authors' implementation (cotan weight, 1e-6)
//!   --angle-tol X      tolerance of the angle-sum test in radians (default 1e-10)
//!   --trace-tol X      tolerance of the edge traces, relative to the edge length (default 1e-6)
//!   --repeat N         repeat each timing N times and print the fastest (default 3)
//!   --failures N       print the first N failed edge traces of each mesh (default 3)
//! ```
//!
//! Polygons are split into a fan. For every mesh one block is printed: the flip report, a
//! quality table (before / after) and the edge-trace report.

use std::time::Instant;

use ddg_core::io::read_obj;
use ddg_core::{Geometry, TriMesh};
use ddg_intrinsic::{
    DelaunayCriterion, DelaunayOptions, EdgeTraceOptions, IntrinsicError, IntrinsicTriangulation,
    Quality,
};

struct Config {
    criterion: DelaunayCriterion,
    trace_tolerance: f64,
    repeat: usize,
    failures: usize,
    paths: Vec<String>,
}

fn parse_args() -> Result<Config, String> {
    let mut config = Config {
        criterion: DelaunayCriterion::default(),
        trace_tolerance: EdgeTraceOptions::default().tolerance,
        repeat: 3,
        failures: 3,
        paths: Vec::new(),
    };
    let mut args = std::env::args().skip(1);
    while let Some(arg) = args.next() {
        let mut value = |name: &str| {
            args.next()
                .ok_or_else(|| format!("{name} needs a value"))
                .and_then(|v| v.parse::<f64>().map_err(|e| format!("{name}: {e}")))
        };
        match arg.as_str() {
            "--reference" => config.criterion = DelaunayCriterion::REFERENCE,
            "--angle-tol" => {
                config.criterion = DelaunayCriterion::AngleSum {
                    tolerance: value("--angle-tol")?,
                }
            }
            "--trace-tol" => config.trace_tolerance = value("--trace-tol")?,
            "--repeat" => config.repeat = (value("--repeat")? as usize).max(1),
            "--failures" => config.failures = value("--failures")? as usize,
            other if other.starts_with("--") => return Err(format!("unknown option {other}")),
            _ => config.paths.push(arg),
        }
    }
    if config.paths.is_empty() {
        return Err("no OBJ files given".to_string());
    }
    Ok(config)
}

fn load(path: &str) -> Result<(TriMesh, Geometry, usize), String> {
    let obj = read_obj(path).map_err(|e| e.to_string())?;
    let mut faces = Vec::new();
    let mut polygons = 0;
    for f in 0..obj.n_faces() {
        let corners = &obj.corner_verts[obj.face_corners(f)];
        if corners.len() != 3 {
            polygons += 1;
        }
        for k in 1..corners.len().saturating_sub(1) {
            faces.push([corners[0], corners[k], corners[k + 1]]);
        }
    }
    let mesh = TriMesh::new(obj.positions.len(), faces).map_err(|e| e.to_string())?;
    let geometry = Geometry::from_positions(&mesh, &obj.positions).map_err(|e| e.to_string())?;
    Ok((mesh, geometry, polygons))
}

fn quality_row(label: &str, q: &Quality) -> String {
    let deg = f64::to_degrees;
    format!(
        "  {label:<7} {:>9.4} {:>8.2} {:>8.2} {:>8.2} {:>8.2} {:>7.2} {:>7.2} {:>7.2} {:>9.3} {:>11.3e}",
        deg(q.min_angle),
        deg(q.largest_angle_quantiles[0]),
        deg(q.largest_angle_quantiles[1]),
        deg(q.largest_angle_quantiles[2]),
        deg(q.largest_angle_quantiles[3]),
        100.0 * q.largest_angle_above[0],
        100.0 * q.largest_angle_above[1],
        100.0 * q.largest_angle_above[2],
        100.0 * q.negative_weight_fraction(),
        q.min_cotan_weight,
    )
}

/// The fastest of `repeat` runs of `f`, in seconds, and the last result.
fn timed<T>(repeat: usize, mut f: impl FnMut() -> T) -> (f64, T) {
    let mut best = f64::INFINITY;
    let mut out = None;
    for _ in 0..repeat {
        let start = Instant::now();
        let value = f();
        best = best.min(start.elapsed().as_secs_f64());
        out = Some(value);
    }
    (best, out.expect("repeat is at least one"))
}

fn report(path: &str, config: &Config) -> Result<(), String> {
    let (mesh, geometry, polygons) = load(path)?;
    let name = std::path::Path::new(path)
        .file_stem()
        .and_then(|s| s.to_str())
        .unwrap_or(path);
    println!(
        "== {name}: {} vertices, {} edges, {} faces, {} boundary loops{}",
        mesh.n_vertices(),
        mesh.n_edges(),
        mesh.n_faces(),
        mesh.boundary_loops().len(),
        if polygons > 0 {
            format!(" ({polygons} polygons split into fans)")
        } else {
            String::new()
        }
    );
    let (t_init, input) = timed(config.repeat, || {
        IntrinsicTriangulation::new(&mesh, &geometry)
    });
    let input = input.map_err(|e| e.to_string())?;
    // Ties (cocircular quadrilaterals) have a weight of the order of the rounding error; they
    // are not counted as negative.
    let weight_threshold = 1e-12;
    let before = input.quality(weight_threshold);

    let options = DelaunayOptions {
        criterion: config.criterion,
        ..DelaunayOptions::default()
    };
    let (t_flip, flipped) = timed(config.repeat, || {
        let mut triangulation = input.clone();
        let report = triangulation.flip_to_delaunay(&options);
        (triangulation, report)
    });
    // The clone is part of the timed closure; measure it alone and subtract.
    let (t_clone, _) = timed(config.repeat, || input.clone());
    let t_flip = (t_flip - t_clone).max(0.0);
    let (triangulation, flip_report) = flipped;
    let flip_report = flip_report.map_err(|e| e.to_string())?;
    let after = triangulation.quality(weight_threshold);
    let irregular = triangulation.irregularities();
    let consistency = triangulation.check().map_err(|e| e.to_string())?;

    println!(
        "  flips {} ({:?}), tested {}, refused {}, not Delaunay {}, ties {}, largest angle excess {:.2e} rad",
        flip_report.n_flips,
        flip_report.status,
        flip_report.n_tested,
        flip_report.n_refused,
        flip_report.n_not_delaunay,
        flip_report.n_ties,
        flip_report.max_angle_excess,
    );
    println!(
        "  after: {} self-edges, {} edges in multi-edges, {} faces with a repeated vertex, {} vertices of degree one",
        irregular.self_edges.len(),
        irregular.multi_edges.len(),
        irregular.faces_with_repeated_vertex.len(),
        irregular.degree_one_vertices.len(),
    );
    println!(
        "  consistency: area {:.1e} (relative), angle sums {:.1e} rad, signposts {:.1e} rad, Gauss-Bonnet {:.1e}",
        consistency.relative_area_error,
        consistency.max_angle_sum_error,
        consistency.max_signpost_error,
        consistency.gauss_bonnet_residual,
    );
    println!(
        "  {:<7} {:>9} {:>8} {:>8} {:>8} {:>8} {:>7} {:>7} {:>7} {:>9} {:>11}",
        "",
        "min ang",
        "big 50%",
        "big 90%",
        "big 99%",
        "big max",
        "%>90",
        "%>120",
        "%>150",
        "%neg cot",
        "min weight"
    );
    println!("{}", quality_row("before", &before));
    println!("{}", quality_row("after", &after));

    let trace_options = EdgeTraceOptions {
        tolerance: config.trace_tolerance,
        ..EdgeTraceOptions::default()
    };
    let (t_trace, traces) = timed(config.repeat, || triangulation.trace_edges(&trace_options));
    let n_edges = triangulation.n_edges();
    let n_on_input = (0..n_edges)
        .filter(|&e| triangulation.edge_on_input(e).is_some())
        .count();
    let n_to_trace = n_edges - n_on_input;
    let (n_failed, n_crossings) = match &traces {
        Ok(traces) => {
            let r = &traces.report;
            println!(
                "  traces: {} edges traced, {} along input edges, 0 failed; {} crossings (at most {} per edge), largest mismatch {:.2e} of the edge length, largest clamp {:.1e}",
                r.n_traced,
                r.n_on_input,
                r.n_crossings,
                r.max_crossings,
                r.max_relative_mismatch,
                r.max_clamp
            );
            (0, r.n_crossings)
        }
        Err(IntrinsicError::EdgeTraces { failures, .. }) => {
            println!(
                "  traces: {} of {} traced edges FAILED the end-point check ({:.3} %); tolerance {:.0e} of the edge length",
                failures.len(),
                n_to_trace,
                100.0 * failures.len() as f64 / n_to_trace.max(1) as f64,
                config.trace_tolerance
            );
            for failure in failures.iter().take(config.failures) {
                println!("    {failure}");
            }
            (failures.len(), 0)
        }
        Err(e) => return Err(e.to_string()),
    };
    // The rate of the flip alone: every edge that can be flipped is flipped there and back.
    let (t_pure, n_pure) = timed(config.repeat, || {
        let mut scratch = triangulation.clone();
        let mut flips = 0usize;
        for e in 0..scratch.n_edges() {
            if scratch.flip_edge(e).is_ok() {
                flips += 1 + usize::from(scratch.flip_edge(e).is_ok());
            }
        }
        flips
    });
    let t_pure = (t_pure - t_clone).max(1e-12);
    println!(
        "  time: init {:.2} ms, flip to Delaunay {:.2} ms ({:.2} M edge tests/s), trace {:.2} ms ({:.2} M crossings/s); flips alone: {} in {:.2} ms ({:.2} M flips/s)",
        1e3 * t_init,
        1e3 * t_flip,
        1e-6 * flip_report.n_tested as f64 / t_flip.max(1e-12),
        1e3 * t_trace,
        1e-6 * n_crossings as f64 / t_trace.max(1e-12),
        n_pure,
        1e3 * t_pure,
        1e-6 * n_pure as f64 / t_pure,
    );
    // One line for tables: name V F flips min-angle before/after, max angle before/after,
    // % negative weights before/after, traced, failed, times.
    println!(
        "  row | {name} | {} | {} | {} | {:.3} | {:.3} | {:.2} | {:.2} | {:.2} | {:.2} | {} | {} | {} | {} | {:.2} | {:.2} | {:.2}",
        mesh.n_vertices(),
        mesh.n_faces(),
        flip_report.n_flips,
        before.min_angle.to_degrees(),
        after.min_angle.to_degrees(),
        before.max_angle.to_degrees(),
        after.max_angle.to_degrees(),
        100.0 * before.negative_weight_fraction(),
        100.0 * after.negative_weight_fraction(),
        irregular.self_edges.len(),
        irregular.multi_edges.len(),
        n_to_trace,
        n_failed,
        1e3 * t_init,
        1e3 * t_flip,
        1e3 * t_trace,
    );
    Ok(())
}

fn main() {
    let config = match parse_args() {
        Ok(config) => config,
        Err(message) => {
            eprintln!(
                "{message}\nusage: delaunay_report [--reference] [--angle-tol X] [--trace-tol X] [--repeat N] [--failures N] mesh.obj ..."
            );
            std::process::exit(2);
        }
    };
    let mut failed = false;
    for path in &config.paths {
        if let Err(message) = report(path, &config) {
            println!("== {path}: {message}");
            failed = true;
        }
    }
    if failed {
        std::process::exit(1);
    }
}
