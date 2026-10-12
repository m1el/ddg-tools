//! Verification tables on the optimal-cones dataset (real output of the reference BFF).
//!
//! ```text
//! cargo run --release -p ddg-bff --example dataset_report            # all meshes
//! cargo run --release -p ddg-bff --example dataset_report -- bunny   # some meshes
//! ```
//!
//! Prints, per mesh: agreement with the stored UV map when flattening with the reference's
//! own cut, agreement of the flat metric when flattening with our cut, overlap and flipped
//! triangles, and timings. Everything runs on one thread.

#[path = "../tests/common/mod.rs"]
mod common;

use std::time::Instant;

use common::{
    DatasetMesh, compare_layouts, count_flipped, dataset_names, length_ratio_quantiles, load,
    total_length,
};
use ddg_bff::{
    ConeFlattenInput, ConeFlattenOptions, ConeFlattening, CutSource, GaussBonnetPolicy,
    REFERENCE_DIAGONAL_SHIFT, REFERENCE_MOLLIFY_DELTA, ScaleFactorOptions, flatten_with_cones,
    mollify_edge_lengths,
};
use ddg_core::Geometry;

fn options(shift: f64) -> ConeFlattenOptions {
    ConeFlattenOptions {
        scale_factors: ScaleFactorOptions {
            gauss_bonnet: GaussBonnetPolicy::SpreadOverVertices,
            diagonal_shift: shift,
            ..ScaleFactorOptions::default()
        },
        ..ConeFlattenOptions::default()
    }
}

fn flatten(
    data: &DatasetMesh,
    geometry: &Geometry,
    cut: CutSource<'_>,
    shift: f64,
) -> Result<ConeFlattening, ddg_bff::BffError> {
    flatten_cones(data, geometry, &data.cones, cut, shift)
}

fn flatten_cones(
    data: &DatasetMesh,
    geometry: &Geometry,
    cones: &[ddg_bff::Cone],
    cut: CutSource<'_>,
    shift: f64,
) -> Result<ConeFlattening, ddg_bff::BffError> {
    flatten_with_cones(
        &ConeFlattenInput {
            mesh: &data.mesh,
            geometry,
            cones,
            scale_factors: None,
            cut,
            anchors: &[],
        },
        &options(shift),
    )
}

struct Row {
    name: String,
    same_cut: String,
    own_cut: String,
    overlap: String,
    timing: String,
}

fn report(name: &str) -> Row {
    let start = Instant::now();
    let data = load(name);
    let load_seconds = start.elapsed().as_secs_f64();
    let n_faces = data.mesh.n_faces();
    let closed = data.mesh.boundary_loops().is_empty();
    let cone_sum: f64 = data.cones.iter().map(|c| c.angle).sum();
    let discrepancy = if closed {
        format!(
            "{:9.1e}",
            cone_sum - 2.0 * std::f64::consts::PI * data.mesh.euler_characteristic() as f64
        )
    } else {
        "        -".to_string()
    };

    // 1. Same cut, same map. Where the stored UVs do not show the whole cut (cones of
    //    nearly zero angle), the planner completes it.
    let (completed_cut, n_completed) = data.completed_reference_cut();
    let reference_cut = CutSource::Edges(&completed_cut);
    let same_cut = match flatten(&data, &data.geometry, reference_cut, 0.0) {
        Err(e) => format!("FAILED: {e}"),
        Ok(exact) => {
            let a = compare_layouts(
                &exact.face_chart,
                exact.n_charts(),
                &exact.corner_uvs,
                &data.reference_uvs,
            );
            let shifted = flatten(&data, &data.geometry, reference_cut, REFERENCE_DIAGONAL_SHIFT)
                .expect("the shifted solve works where the exact one does");
            let b = compare_layouts(
                &shifted.face_chart,
                shifted.n_charts(),
                &shifted.corner_uvs,
                &data.reference_uvs,
            );
            // The reference mollifies its edge lengths; this matters only if a triangle is
            // nearly degenerate.
            let mollified = mollify_edge_lengths(
                &data.mesh,
                data.geometry.edge_lengths(),
                REFERENCE_MOLLIFY_DELTA,
            )
            .unwrap();
            let c = if mollified.epsilon > 0.0 {
                let geometry =
                    Geometry::from_edge_lengths(&data.mesh, mollified.edge_lengths).unwrap();
                let m = flatten(&data, &geometry, reference_cut, REFERENCE_DIAGONAL_SHIFT).unwrap();
                let c = compare_layouts(
                    &m.face_chart,
                    m.n_charts(),
                    &m.corner_uvs,
                    &data.reference_uvs,
                );
                format!("{:8.1e} {:8.1e}", c.rms, c.max)
            } else {
                "       =        =".to_string()
            };
            // How much does rounding the input to six digits, as the files do, change the
            // result? Worst of three random perturbations of that size.
            let mut sensitivity = (0.0f64, 0.0f64);
            for seed in 1..=3 {
                let (geometry, cones) = common::perturbed_like_file_rounding(&data, seed);
                let p = flatten_cones(&data, &geometry, &cones, reference_cut, REFERENCE_DIAGONAL_SHIFT)
                    .expect("perturbed input flattens");
                let d = compare_layouts(
                    &p.face_chart,
                    p.n_charts(),
                    &p.corner_uvs,
                    &shifted.corner_uvs,
                );
                sensitivity = (sensitivity.0.max(d.rms), sensitivity.1.max(d.max));
            }
            format!(
                "{:6} {:6}{:+3} {:2} {} | {:8.1e} {:8.1e} | {:8.1e} {:8.1e} | {} | {:8.1e} {:8.1e} | {:8.1e} {:8.1e}",
                data.cones.len(),
                data.reference_cut.len(),
                n_completed,
                exact.n_charts(),
                discrepancy,
                a.rms,
                a.max,
                b.rms,
                b.max,
                c,
                sensitivity.0,
                sensitivity.1,
                shifted.seams.max_length_mismatch,
                shifted
                    .charts
                    .iter()
                    .map(|c| c.closure.max_length_change)
                    .fold(0.0, f64::max),
            )
        }
    };

    // 2. Own cut, same metric; 6. timing (best of three).
    let own_cut = CutSource::Plan {
        may_cut: None,
        must_cut: &[],
    };
    let mut best: Option<ConeFlattening> = None;
    for _ in 0..3 {
        let result = flatten(&data, &data.geometry, own_cut, 0.0);
        match result {
            Err(e) => {
                return Row {
                    name: name.to_string(),
                    same_cut,
                    own_cut: format!("FAILED: {e}"),
                    overlap: String::new(),
                    timing: String::new(),
                };
            }
            Ok(r) => {
                if best.as_ref().is_none_or(|b| r.timings.total < b.timings.total) {
                    best = Some(r);
                }
            }
        }
    }
    let own = best.expect("three runs");
    let q = length_ratio_quantiles(&own.corner_uvs, &data.reference_uvs);
    let plan = own.plan.as_ref().expect("planned");
    let own_length = total_length(&data.geometry, &plan.cut_edges);
    let reference_length = total_length(&data.geometry, &data.reference_cut);
    let own_cut_text = format!(
        "{:6} {:6} {:7.3} | {:7.4} {:7.4} {:7.4} {:7.4} | {:8.1e} {:8.1e} {:8.1e} | {:6.4} {:7.4}",
        data.reference_cut.len(),
        plan.cut_edges.len(),
        own_length / reference_length,
        q[0],
        q[1],
        q[3],
        q[4],
        own.seams.max_length_mismatch,
        own.charts.iter().map(|c| c.closure.gap_before).fold(0.0, f64::max),
        own.charts
            .iter()
            .map(|c| c.closure.max_length_change)
            .fold(0.0, f64::max),
        own.charts[0].distortion.quasi_conformal.mean,
        own.charts[0].distortion.quasi_conformal.max.min(9999.0),
    );

    // Overlap and flips: the stored map, ours with the reference's cut, ours with our cut.
    let with_reference_cut = flatten(&data, &data.geometry, reference_cut, 0.0).ok();
    let describe = |r: &ConeFlattening| {
        let negative: Vec<_> = r.cones.iter().filter(|c| c.angle < 0.0).collect();
        format!(
            "{:5} {:5} {:5} {:5} {:5}",
            r.charts.iter().map(|c| c.distortion.flipped_triangles).sum::<usize>(),
            r.charts.iter().map(|c| c.overlap.wedges_over_2pi).sum::<usize>(),
            r.charts.iter().map(|c| c.overlap.boundary_crossings).sum::<usize>(),
            negative.iter().filter(|c| c.n_wedges == 1).count(),
            negative
                .iter()
                .filter(|c| c.max_wedge_angle > 2.0 * std::f64::consts::PI)
                .count(),
        )
    };
    let negative = data.cones.iter().filter(|c| c.angle < 0.0).count();
    let overlap = format!(
        "{:6} {:6} | {:5} | {} | {}",
        data.cones.len(),
        negative,
        count_flipped(&data.reference_uvs),
        with_reference_cut
            .as_ref()
            .map_or("FAILED".to_string(), describe),
        describe(&own),
    );

    let t = own.timings;
    let timing = format!(
        "{:7} {:6} | {:8.1} | {:7.1} {:7.1} {:7.1} {:7.1} {:7.1} | {:8.1}",
        n_faces,
        data.cones.len(),
        1e3 * t.total,
        1e3 * t.scale_factors,
        1e3 * t.cut,
        1e3 * t.curves,
        1e3 * t.extension,
        1e3 * t.report,
        1e3 * load_seconds,
    );
    Row {
        name: name.to_string(),
        same_cut,
        own_cut: own_cut_text,
        overlap,
        timing,
    }
}

fn main() {
    let Some(mut names) = dataset_names() else {
        eprintln!("dataset not found at {}", common::dataset_dir().display());
        std::process::exit(1);
    };
    let wanted: Vec<String> = std::env::args().skip(1).collect();
    if !wanted.is_empty() {
        names.retain(|n| wanted.contains(n));
    }
    let rows: Vec<Row> = names.iter().map(|n| report(n)).collect();

    println!("## 1. Same cut, same map");
    println!("Flattened with the cut recovered from the stored UV map; per chart the best similarity is");
    println!("fitted; residuals (rms, max) are relative to the diagonal of the chart's bounding box in the");
    println!("stored map. 'input rounding': change of our result when positions and cone angles are");
    println!("perturbed by half a unit in the sixth digit (worst of 3 samples), i.e. what the six-digit");
    println!("files can cost.");
    println!(
        "{:<20} {:>6} {:>9} {:>2} {:>9} | {:>17} | {:>17} | {:>17} | {:>17} | {:>17}",
        "mesh", "cones", "cutE+add", "ch", "GB discr",
        "exact: rms, max", "shift 1e-8", "shift+mollify", "input rounding", "seam mism, max adj"
    );
    for r in &rows {
        println!("{:<20} {}", r.name, r.same_cut);
    }

    println!("\n## 2. Own cut, same metric");
    println!("Flattened with the planner's cut (exact solves). Ratio of UV edge lengths ours / stored, per");
    println!("triangle edge, divided by its median: minimum, 1 %, 99 %, maximum.");
    println!(
        "{:<20} {:>6} {:>6} {:>7} | {:>7} {:>7} {:>7} {:>7} | {:>8} {:>8} {:>8} | {:>6} {:>7}",
        "mesh", "refE", "ownE", "len/ref", "min", "p1", "p99", "max",
        "seam", "gap", "max adj", "Qmean", "Qmax"
    );
    for r in &rows {
        println!("{:<20} {}", r.name, r.own_cut);
    }

    println!("\n## 3. Flipped triangles and overlap");
    println!("flip = triangles with negative orientation; w>2pi = boundary vertices with more than a full");
    println!("turn around them; cross = pairs of crossing boundary edges; leaf- = negative cones at the end");
    println!("of a cut; ovl- = negative cones with a wedge over 2 pi.");
    println!(
        "{:<20} {:>6} {:>6} | {:>5} | {:>29} | {:>29}",
        "mesh", "cones", "neg", "ref", "reference cut", "own cut"
    );
    println!(
        "{:<20} {:>6} {:>6} | {:>5} | {:>5} {:>5} {:>5} {:>5} {:>5} | {:>5} {:>5} {:>5} {:>5} {:>5}",
        "", "", "", "flip", "flip", "w>2pi", "cross", "leaf-", "ovl-", "flip", "w>2pi", "cross", "leaf-", "ovl-"
    );
    for r in &rows {
        println!("{:<20} {}", r.name, r.overlap);
    }

    println!("\n## 4. Timing (milliseconds, one thread, best of 3)");
    println!("flatten_with_cones with the planner's cut: scale factors, cut, curves, extension, report.");
    println!(
        "{:<20} {:>7} {:>6} | {:>8} | {:>7} {:>7} {:>7} {:>7} {:>7} | {:>8}",
        "mesh", "faces", "cones", "total", "scale", "cut", "curves", "extend", "report", "load"
    );
    for r in &rows {
        println!("{:<20} {}", r.name, r.timing);
    }
}
