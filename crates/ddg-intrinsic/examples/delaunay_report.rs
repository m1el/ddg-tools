//! Reads OBJ files and, for each, flips to the intrinsic Delaunay triangulation, traces the
//! intrinsic edges over the input mesh, refines to an angle bound, builds the overlay and the
//! simplicial export, and prints mesh quality before and after with timings.
//!
//! ```text
//! cargo run --release -p ddg-intrinsic --example delaunay_report -- [options] a.obj b.obj ...
//!   --reference        use the Delaunay test of the authors' implementation (cotan weight, 1e-6)
//!   --angle-tol X      tolerance of the angle-sum test in radians (default 1e-10)
//!   --trace-tol X      tolerance of the edge traces, relative to the edge length (default 1e-6)
//!   --repeat N         repeat each timing N times and print the fastest (default 3)
//!   --failures N       print the first N failed edge traces of each mesh (default 3)
//!   --refine DEG       refine to this angle bound in degrees (may be given several times;
//!                      default 25 and 30; `--refine 0` alone: no refinement)
//! ```
//!
//! Polygons are split into a fan. For every mesh one block is printed: the flip report, a
//! quality table (before / after flipping / after each refinement), the edge-trace report,
//! and for each angle bound the refinement report, the overlay and the export. Lines that
//! start with `row |` and `refine |` are for tables. Timings are single-thread: run one
//! process at a time when they matter.

use std::time::Instant;

use ddg_core::io::read_obj;
use ddg_core::{Geometry, TriMesh};
use ddg_intrinsic::{
    BelowReason, DelaunayCriterion, DelaunayOptions, EdgeTraceOptions, ExportOptions,
    IntrinsicError, IntrinsicTriangulation, Quality, RefineOptions,
};

struct Config {
    criterion: DelaunayCriterion,
    trace_tolerance: f64,
    repeat: usize,
    failures: usize,
    refine: Vec<f64>,
    paths: Vec<String>,
}

fn parse_args() -> Result<Config, String> {
    let mut config = Config {
        criterion: DelaunayCriterion::default(),
        trace_tolerance: EdgeTraceOptions::default().tolerance,
        repeat: 3,
        failures: 3,
        refine: Vec::new(),
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
            "--refine" => config.refine.push(value("--refine")?),
            other if other.starts_with("--") => return Err(format!("unknown option {other}")),
            _ => config.paths.push(arg),
        }
    }
    if config.paths.is_empty() {
        return Err("no OBJ files given".to_string());
    }
    if config.refine.is_empty() {
        config.refine = vec![25.0, 30.0];
    }
    config.refine.retain(|&degrees| degrees > 0.0);
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
    for &degrees in &config.refine {
        refine_report(name, &mesh, &geometry, degrees, config)?;
    }
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

/// Refinement to `degrees`, then the traces, the overlay and the export of the result.
fn refine_report(
    name: &str,
    mesh: &TriMesh,
    geometry: &Geometry,
    degrees: f64,
    config: &Config,
) -> Result<(), String> {
    let text = |e: IntrinsicError| e.to_string();
    let options = RefineOptions {
        min_angle: degrees.to_radians(),
        criterion: config.criterion,
        ..RefineOptions::default()
    };
    let fresh = IntrinsicTriangulation::new(mesh, geometry).map_err(text)?;
    let (t_clone, _) = timed(config.repeat, || fresh.clone());
    let (t_refine, result) = timed(config.repeat, || {
        let mut tri = fresh.clone();
        let report = tri.refine(&options);
        (tri, report)
    });
    let t_refine = (t_refine - t_clone).max(1e-12);
    let (mut tri, report) = result;
    let report = report.map_err(text)?;
    let quality = tri.quality(1e-12);
    println!("{}", quality_row(&format!("ref {degrees:.0}"), &quality));
    // What is left below the bound, by reason.
    let mut reasons = [0usize; 4];
    for b in &report.below_bound {
        reasons[match b.reason {
            BelowReason::AngleSum { .. } => 0,
            BelowReason::NeedleVertex { .. } => 1,
            BelowReason::NotInsertable { .. } => 2,
            BelowReason::InsertionLimit => 3,
        }] += 1;
    }
    println!(
        "  refine {degrees:.0}: {:?}; {} inserted ({} in faces, {} on interior edges, {} on boundary edges), {} flips, {} places moved onto input edges; smallest angle {:.3} ({:.3} where the bound can be reached), largest {:.2}; below the bound: {} corners at vertices with a smaller angle sum, {} next to a needle vertex, {} not insertable, {} left by the limit; {:.2} ms",
        report.status,
        report.n_inserted,
        report.n_face_insertions,
        report.n_edge_splits,
        report.n_boundary_splits,
        report.n_flips,
        report.n_snapped,
        report.min_angle.to_degrees(),
        report.min_angle_reachable.to_degrees(),
        report.max_angle.to_degrees(),
        reasons[0],
        reasons[1],
        reasons[2],
        reasons[3],
        1e3 * t_refine,
    );
    let trace_options = EdgeTraceOptions {
        tolerance: config.trace_tolerance,
        ..EdgeTraceOptions::default()
    };
    let (t_trace, traces) = timed(config.repeat, || tri.trace_edges(&trace_options));
    let traces = traces.map_err(text)?;
    let (t_overlay, overlay) = timed(config.repeat, || tri.overlay(&traces));
    let overlay = overlay.map_err(text)?;
    println!(
        "    traces: {} edges, {} crossings, largest mismatch {:.2e}, {:.2} ms; overlay: {} vertices, {} faces, {:.2} ms",
        traces.report.n_traced,
        traces.report.n_crossings,
        traces.report.max_relative_mismatch,
        1e3 * t_trace,
        overlay.vertices.len(),
        overlay.n_faces(),
        1e3 * t_overlay,
    );
    // The export, by splits and flips only, and with refinement between the rounds.
    let irregular = tri.irregularities();
    let mut refined_between = tri.clone();
    let (t_export, plain) = timed(1, || tri.export_simplicial(&ExportOptions::default()));
    let plain = plain.map_err(text)?;
    let (t_export_refined, with) = timed(1, || {
        refined_between.export_simplicial(&ExportOptions {
            refine: Some(options),
            ..ExportOptions::default()
        })
    });
    let with = with.map_err(text)?;
    println!(
        "    export: {} self-edges, {} edges in groups; by splits and flips: {} splits, {} flips, angles {:.3} to {:.2}, {:.2} ms; with refinement between the rounds: {} rounds, {} splits, {} inserted, angles {:.3} to {:.2}, {}, {:.2} ms; fingerprint {:016x}",
        irregular.self_edges.len(),
        irregular.multi_edges.len(),
        plain.n_splits,
        plain.n_flips,
        plain.min_angle.to_degrees(),
        plain.max_angle.to_degrees(),
        1e3 * t_export,
        with.n_refine_rounds,
        with.n_splits,
        with.n_refined,
        with.min_angle.to_degrees(),
        with.max_angle.to_degrees(),
        if with.refinement.is_some() {
            "ended refined"
        } else {
            "not ended refined"
        },
        1e3 * t_export_refined,
        with.fingerprint,
    );
    // One line for tables: name bound V inserted min reachable max below-by-reason time,
    // overlay faces, export splits and angles (plain, then with refinement).
    println!(
        "  refine | {name} | {degrees:.0} | {} | {} | {:.3} | {:.3} | {:.2} | {} | {} | {} | {} | {:.2} | {} | {} | {:.3} | {:.2} | {} | {:.3} | {:.2}",
        mesh.n_vertices(),
        report.n_inserted,
        report.min_angle.to_degrees(),
        report.min_angle_reachable.to_degrees(),
        report.max_angle.to_degrees(),
        reasons[0],
        reasons[1],
        reasons[2],
        reasons[3],
        1e3 * t_refine,
        overlay.n_faces(),
        plain.n_splits,
        plain.min_angle.to_degrees(),
        plain.max_angle.to_degrees(),
        with.n_splits,
        with.min_angle.to_degrees(),
        with.max_angle.to_degrees(),
    );
    Ok(())
}

fn main() {
    let config = match parse_args() {
        Ok(config) => config,
        Err(message) => {
            eprintln!(
                "{message}\nusage: delaunay_report [--reference] [--angle-tol X] [--trace-tol X] [--repeat N] [--failures N] [--refine DEG]... mesh.obj ..."
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
