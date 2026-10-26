//! The feedback loop of the cone solver on the dataset: what it costs, and how far the solver's
//! linearized scale factor is from the flattening's on every mesh.
//!
//! ```text
//! cargo run --release -p ddg-uv --example cone_feedback [name ...]
//! ```
//!
//! For every mesh of `data/optimal-cones/optimal-cones-data/data` (or `DDG_DATASET`; the
//! arguments restrict to exact names) and every `λ` of `DDG_LAMBDAS` (default `0.1,0.5`), with
//! `u = 0` on the boundary and no angle bounds (the problem of the crate README's tables):
//!
//! * `plain`: one solve and one Boundary First Flattening with its cones and scale factors,
//!   their seconds, the cones, and `D*` of the flattening (the dataset's measure);
//! * `loop`: `solve_with_feedback` with the default `FeedbackOptions` (or `DDG_FB_THRESHOLD`,
//!   `DDG_FB_ROUNDS`): rounds, reweighted vertices, why it stopped, the seconds of every round
//!   (solve, flattening, comparison), the total, the total over the plain solve, cones, `D*`;
//! * `mismatch`: the distribution of `|δ|` over the interior vertices in the first round (the
//!   plain solve): largest, 99.9 % and 99 % quantiles, and how many vertices exceed 2, 3, 4, 6,
//!   10 and 100: what a threshold would pick.
//!
//! The solver runs single-threaded. For timing numbers run one process at a time on a quiet
//! machine; for the mismatch tables the meshes may run in parallel.

use std::path::PathBuf;
use std::time::Instant;

use ddg_bff::{
    ConeFlattenInput, ConeFlattenOptions, CutSource, GaussBonnetPolicy, ScaleFactorOptions,
    flatten_with_cones,
};
use ddg_core::io::{Cone, read_obj};
use ddg_core::{Geometry, TriMesh};
use ddg_optimal_cones::{
    BoundaryTreatment, ConeProblem, FeedbackOptions, SolveOptions, scale_mismatch,
    solve_with_boundary, solve_with_feedback,
};
use ddg_uv::realized_face_log_scale;

fn dataset_dir() -> PathBuf {
    std::env::var_os("DDG_DATASET").map_or_else(
        || {
            PathBuf::from(env!("CARGO_MANIFEST_DIR"))
                .join("../../data/optimal-cones/optimal-cones-data/data")
        },
        PathBuf::from,
    )
}

fn env_f64(name: &str, default: f64) -> f64 {
    std::env::var(name)
        .ok()
        .and_then(|s| s.parse().ok())
        .unwrap_or(default)
}

/// The flattening of the solver's cones with its scale factors: corner UVs, seconds.
fn flatten(
    mesh: &TriMesh,
    geometry: &Geometry,
    cones: &[Cone],
    u: &[f64],
) -> Result<(Vec<[f64; 2]>, f64), String> {
    let start = Instant::now();
    let flat = flatten_with_cones(
        &ConeFlattenInput {
            mesh,
            geometry,
            cones,
            scale_factors: Some(u),
            cut: CutSource::Plan {
                may_cut: None,
                must_cut: &[],
            },
            anchors: &[],
        },
        &ConeFlattenOptions {
            scale_factors: ScaleFactorOptions {
                gauss_bonnet: GaussBonnetPolicy::Reject,
                ..ScaleFactorOptions::default()
            },
            overlap_samples_per_triangle: 0.0,
            ..ConeFlattenOptions::default()
        },
    )
    .map_err(|e| e.to_string())?;
    Ok((flat.corner_uvs, start.elapsed().as_secs_f64()))
}

/// The dataset's measure `D*`: `Σ a_t (u_t − ū)²` over the triangles whose image has the
/// majority orientation, `a_t` the area share, `u_t = ½ ln(B_t/A_t)`.
fn dstar(geometry: &Geometry, uvs: &[[f64; 2]]) -> f64 {
    let areas = geometry.face_areas();
    let signed: Vec<f64> = (0..areas.len())
        .map(|f| {
            let (p, q, r) = (uvs[3 * f], uvs[3 * f + 1], uvs[3 * f + 2]);
            0.5 * ((q[0] - p[0]) * (r[1] - p[1]) - (q[1] - p[1]) * (r[0] - p[0]))
        })
        .collect();
    let positive = signed.iter().filter(|b| **b > 0.0).count();
    let sign = if 2 * positive >= signed.len() {
        1.0
    } else {
        -1.0
    };
    let total: f64 = areas.iter().sum();
    let (mut sw, mut swu, mut swuu) = (0.0, 0.0, 0.0);
    for (f, &a) in areas.iter().enumerate() {
        let b = sign * signed[f];
        if !(a > 0.0 && b > 0.0) {
            continue;
        }
        let u = 0.5 * (b / a).ln();
        if !u.is_finite() {
            continue;
        }
        let w = a / total;
        sw += w;
        swu += w * u;
        swuu += w * u * u;
    }
    if sw > 0.0 {
        (swuu - swu * swu / sw).max(0.0)
    } else {
        f64::NAN
    }
}

/// The quantile `q` of the sorted values (nearest rank).
fn quantile(sorted: &[f64], q: f64) -> f64 {
    if sorted.is_empty() {
        return f64::NAN;
    }
    let k = ((q * sorted.len() as f64).ceil() as usize).clamp(1, sorted.len()) - 1;
    sorted[k]
}

fn run(name: &str, lambdas: &[f64], feedback: &FeedbackOptions) -> Result<(), String> {
    let obj = read_obj(dataset_dir().join(format!("{name}.obj"))).map_err(|e| e.to_string())?;
    let faces = obj.triangles().ok_or("not a triangle mesh")?;
    let mesh = TriMesh::new(obj.positions.len(), faces).map_err(|e| e.to_string())?;
    let geometry = Geometry::from_positions(&mesh, &obj.positions).map_err(|e| e.to_string())?;
    let bt = BoundaryTreatment::Fixed;
    let options = SolveOptions::default();
    for &lambda in lambdas {
        let tag = format!("{name} λ={lambda}");
        let problem = ConeProblem::new(&mesh, &geometry, lambda);
        // Plain: one solve, one flattening.
        let start = Instant::now();
        let plain = solve_with_boundary(&problem, bt, &options).map_err(|e| e.to_string())?;
        let plain_solve = start.elapsed().as_secs_f64();
        let (plain_uvs, plain_flatten) =
            match flatten(&mesh, &geometry, &plain.cones, &plain.log_scale) {
                Ok(x) => x,
                Err(e) => {
                    println!("{tag} | plain flattening FAILED: {e}");
                    continue;
                }
            };
        println!(
            "{tag} | plain | {} vertices | {} cones | solve {plain_solve:.3}s | flatten {plain_flatten:.3}s | D* {:.5e}",
            mesh.n_vertices(),
            plain.cones.len(),
            dstar(&geometry, &plain_uvs)
        );
        // The mismatch of the plain solve.
        let faces = realized_face_log_scale(&geometry, &plain_uvs);
        let m = scale_mismatch(&mesh, &geometry, &plain.log_scale, &faces)
            .map_err(|e| e.to_string())?;
        let mut interior: Vec<f64> = (0..mesh.n_vertices())
            .filter(|&v| !mesh.is_boundary_vertex(v) && m.mismatch[v].is_finite())
            .map(|v| m.mismatch[v].abs())
            .collect();
        interior.sort_by(f64::total_cmp);
        let above = |t: f64| interior.iter().filter(|&&x| x > t).count();
        println!(
            "{tag} | mismatch | max {:.3e} | q99.9 {:.3e} | q99 {:.3e} | rms {:.3e} | >2: {} | >3: {} | >4: {} | >6: {} | >10: {} | >100: {} | unmeasured {}",
            interior.last().copied().unwrap_or(0.0),
            quantile(&interior, 0.999),
            quantile(&interior, 0.99),
            m.rms,
            above(2.0),
            above(3.0),
            above(4.0),
            above(6.0),
            above(10.0),
            above(100.0),
            m.n_unmeasured
        );
        // The loop.
        let start = Instant::now();
        let result = solve_with_feedback(
            &problem,
            bt,
            &options,
            feedback,
            |cones: &[Cone], u: &[f64]| {
                let (uvs, _) = flatten(&mesh, &geometry, cones, u)?;
                Ok::<_, String>((realized_face_log_scale(&geometry, &uvs), uvs))
            },
        );
        let total = start.elapsed().as_secs_f64();
        match result {
            Err(e) => println!("{tag} | loop FAILED: {e}"),
            Ok(r) => {
                let report = &r.feedback;
                let rounds: Vec<String> = report
                    .rounds
                    .iter()
                    .map(|x| {
                        format!(
                            "{}solve {:.3}s flatten {:.3}s compare {:.4}s max|δ| {:.2e} over {}",
                            if x.probe { "probe: " } else { "" },
                            x.seconds_solve,
                            x.seconds_realize,
                            x.seconds_mismatch,
                            x.max_abs_eligible,
                            x.n_over_threshold
                        )
                    })
                    .collect();
                println!(
                    "{tag} | loop | {} rounds | {} reweighted | {:?} | total {total:.3}s = {:.2} x plain solve, {:.2} x plain solve + flatten | {} cones | D* {:.5e} | {}",
                    report.rounds.len(),
                    report.reweighted.len(),
                    report.stop,
                    total / plain_solve,
                    total / (plain_solve + plain_flatten),
                    r.solution.cones.len(),
                    dstar(&geometry, &r.realized),
                    rounds.join(" ; ")
                );
            }
        }
    }
    Ok(())
}

fn main() {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let lambdas: Vec<f64> = std::env::var("DDG_LAMBDAS")
        .unwrap_or_else(|_| "0.1,0.5".to_owned())
        .split(',')
        .filter_map(|s| s.trim().parse().ok())
        .collect();
    let d = FeedbackOptions::default();
    let feedback = FeedbackOptions {
        threshold: env_f64("DDG_FB_THRESHOLD", d.threshold),
        max_rounds: env_f64("DDG_FB_ROUNDS", d.max_rounds as f64) as usize,
        max_reweighted_fraction: env_f64("DDG_FB_FRACTION", d.max_reweighted_fraction),
        reduced_weight: env_f64("DDG_FB_REDUCED", d.reduced_weight),
        cone_free_probe: env_f64("DDG_FB_PROBE", 1.0) != 0.0,
    };
    let mut names: Vec<String> = std::fs::read_dir(dataset_dir())
        .map(|rd| {
            rd.filter_map(|e| e.ok())
                .filter_map(|e| e.file_name().into_string().ok())
                .filter_map(|n| n.strip_suffix(".obj").map(str::to_owned))
                .collect()
        })
        .unwrap_or_default();
    names.sort();
    if !args.is_empty() {
        names.retain(|n| args.contains(n));
    }
    if names.is_empty() {
        eprintln!("no dataset meshes found in {}", dataset_dir().display());
        return;
    }
    for name in &names {
        if let Err(e) = run(name, &lambdas, &feedback) {
            println!("{name} | FAILED: {e}");
        }
    }
}
