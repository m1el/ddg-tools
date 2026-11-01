//! Distortion (Fig. 18) and area error (Fig. 19) of ICE against QEM at equal vertex counts
//! (paper Sec. 8.3, 8.4; SPEC §10.3, §10.4, §10 "Results of record").
//!
//! ```text
//! cargo run --release -p ddg-intrinsic-simplification --example evaluate -- MESH.obj \
//!     [--vertices 200,500 | --fractions 0.1,0.02] [--methods ice,qem-closest,qem-successive] \
//!     [--sign paper|cetm] [--threads 8] [--dump DIR]
//! ```
//!
//! One line per (target, method): `ice` is `simplify` with the defaults (SPEC §5; with
//! `--sign cetm` the Eq. 10 sign of S25), `qem-closest` the QEM baseline with the closest-point
//! projection, `qem-successive` with the successive self-parameterization of Liu et al. 2021
//! (closed meshes only). `--dump DIR` writes the per-face distortions (`aniso`, `log σ₁σ₂`,
//! `σ₁/σ₂`) to `DIR/<mesh>_<n>_<method>.tsv`.

use std::time::Instant;

use ddg_core::io::read_obj;
use ddg_core::sparse::ThreadPool;
use ddg_intrinsic_simplification::applications::CoarseApplications;
use ddg_intrinsic_simplification::evaluation::{
    DistortionOptions, DistortionReport, ExtrinsicMap, coarse_area, distortion,
    extrinsic_applications, qem_to_target,
};
use ddg_intrinsic_simplification::{
    CoarseMesh, CoarseMeshOptions, ProjectiveSign, SimplifyOptions, Target, simplify,
};

fn die(msg: String) -> ! {
    eprintln!("{msg}");
    std::process::exit(1);
}

fn arg(args: &[String], name: &str) -> Option<String> {
    args.iter()
        .position(|a| a == name)
        .and_then(|i| args.get(i + 1))
        .cloned()
}

fn area(positions: &[[f64; 3]], faces: &[[u32; 3]]) -> f64 {
    faces
        .iter()
        .map(|f| {
            let p = f.map(|v| positions[v as usize]);
            let a = [p[1][0] - p[0][0], p[1][1] - p[0][1], p[1][2] - p[0][2]];
            let b = [p[2][0] - p[0][0], p[2][1] - p[0][1], p[2][2] - p[0][2]];
            let c = [
                a[1] * b[2] - a[2] * b[1],
                a[2] * b[0] - a[0] * b[2],
                a[0] * b[1] - a[1] * b[0],
            ];
            0.5 * (c[0] * c[0] + c[1] * c[1] + c[2] * c[2]).sqrt()
        })
        .sum()
}

fn quantile(v: &[f64], q: f64) -> f64 {
    let mut x: Vec<f64> = v.iter().copied().filter(|x| x.is_finite()).collect();
    x.sort_by(|a, b| a.total_cmp(b));
    if x.is_empty() {
        f64::NAN
    } else {
        x[((x.len() - 1) as f64 * q) as usize]
    }
}

#[allow(clippy::too_many_arguments)]
fn line(
    mesh: &str,
    n_fine: usize,
    n_coarse: usize,
    method: &str,
    t: f64,
    area_err: f64,
    r: &DistortionReport,
    t_dist: f64,
) {
    let aniso: Vec<f64> = r.per_face.iter().map(|x| x[0]).collect();
    println!(
        "{mesh}\t{n_fine}\t{n_coarse}\t{method}\t{t:.3}\t{:.4}\t{:.4}\t{:.4}\t{:.4}\t{:.4}\t{:.4}\t{:.4}\t{:.4}\t{}\t{:.2e}\t{}/{}/{}/{}/{}\t{:.1}",
        100.0 * area_err,
        r.aniso_mean,
        quantile(&aniso, 0.5),
        quantile(&aniso, 0.99),
        r.aniso_max,
        r.area_mean_abs,
        r.area_max_abs,
        r.area_mean_signed,
        r.n_invalid_faces,
        r.invalid_area_fraction,
        r.n_same_face,
        r.n_unfolded,
        r.n_via_vertex,
        r.n_fallback,
        r.n_unresolved,
        t_dist
    );
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let Some(path) = args.get(1) else {
        die("usage: evaluate MESH.obj [--vertices a,b | --fractions a,b] [--methods ...] [--sign paper|cetm] [--threads n] [--dump DIR]".into());
    };
    let obj = read_obj(path).unwrap_or_else(|e| die(format!("{path}: {e}")));
    let mut faces = Vec::new();
    for f in 0..obj.n_faces() {
        let c = &obj.corner_verts[obj.face_corners(f)];
        for k in 1..c.len().saturating_sub(1) {
            faces.push([c[0], c[k], c[k + 1]]);
        }
    }
    let positions = obj.positions;
    let n = positions.len();
    let name = std::path::Path::new(path)
        .file_stem()
        .map(|s| s.to_string_lossy().to_string())
        .unwrap_or_default();
    let sign = match arg(&args, "--sign").as_deref() {
        None | Some("paper") => ProjectiveSign::Paper,
        Some("cetm") => ProjectiveSign::Cetm,
        Some(x) => die(format!("--sign: paper or cetm, not {x}")),
    };
    let sign_name = if sign == ProjectiveSign::Paper {
        "paper"
    } else {
        "cetm"
    };
    let threads: usize = arg(&args, "--threads")
        .and_then(|x| x.parse().ok())
        .unwrap_or(8);
    let dump = arg(&args, "--dump");
    let methods: Vec<String> = arg(&args, "--methods")
        .unwrap_or_else(|| "ice,qem-closest,qem-successive".to_string())
        .split(',')
        .map(String::from)
        .collect();
    let targets: Vec<usize> = if let Some(v) = arg(&args, "--vertices") {
        v.split(',').filter_map(|x| x.parse().ok()).collect()
    } else {
        arg(&args, "--fractions")
            .unwrap_or_else(|| "0.1,0.02".to_string())
            .split(',')
            .filter_map(|x| x.parse::<f64>().ok())
            .map(|f| ((n as f64) * f).round() as usize)
            .collect()
    };
    let pool = ThreadPool::new(threads).unwrap_or_else(|e| die(format!("{e}")));
    let fine_area = area(&positions, &faces);
    println!(
        "# {path}: {n} vertices, {} faces, fine area {fine_area:.6e}, sign {sign_name}",
        faces.len()
    );
    println!(
        "# mesh\tn_fine\tn_coarse\tmethod\tt_coarsen_s\tarea_err_%\taniso_mean\taniso_median\taniso_p99\taniso_max\tlogarea_abs_mean\tlogarea_abs_max\tlogarea_signed_mean\tinvalid_faces\tinvalid_area\tedges same/unfolded/via/fallback/unresolved\tt_distortion_s"
    );
    let dopts = DistortionOptions::default();
    for &target in &targets {
        for m in &methods {
            let t0 = Instant::now();
            let app: Result<(CoarseApplications, usize), String> = match m.as_str() {
                "ice" => (|| {
                    let mut mesh = CoarseMesh::from_positions(
                        &positions,
                        &faces,
                        &[],
                        CoarseMeshOptions {
                            projective_sign: sign,
                            ..CoarseMeshOptions::default()
                        },
                    )
                    .map_err(|e| e.to_string())?;
                    let opts = SimplifyOptions {
                        target: Target::Vertices(target),
                        ..SimplifyOptions::default()
                    };
                    simplify(&mut mesh, &opts).map_err(|e| e.to_string())?;
                    let app = CoarseApplications::new(&mesh).map_err(|e| e.to_string())?;
                    let nc = app.output().n_vertices();
                    Ok((app, nc))
                })(),
                "qem-closest" | "qem-successive" => (|| {
                    let successive = m == "qem-successive";
                    let (q, relaxed) = qem_to_target(&positions, &faces, target, successive)
                        .map_err(|e| e.to_string())?;
                    if relaxed > 0 {
                        eprintln!("# {m}: relaxation rung {relaxed} needed to reach {target}");
                    }
                    if let Some(st) = q.successive_stats {
                        eprintln!("# {m}: successive map {st:?}");
                    }
                    if !q.report.reached_target {
                        eprintln!(
                            "# {m}: stopped at {} vertices (target {target}): {:?}",
                            q.report.n_output_vertices, q.report
                        );
                    }
                    let map = if successive {
                        ExtrinsicMap::Successive
                    } else {
                        ExtrinsicMap::ClosestPoint
                    };
                    let app = extrinsic_applications(&q, n, map).map_err(|e| e.to_string())?;
                    Ok((app, q.report.n_output_vertices))
                })(),
                other => Err(format!("unknown method {other}")),
            };
            let (app, nc) = match app {
                Ok(x) => x,
                Err(e) => {
                    println!("# {name}\t{target}\t{m}\tFAILED: {e}");
                    continue;
                }
            };
            let t_coarse = t0.elapsed().as_secs_f64();
            let ca = coarse_area(&app).unwrap_or(f64::NAN);
            let t1 = Instant::now();
            match distortion(&positions, &faces, &app, &pool, &dopts) {
                Ok(r) => {
                    let label = if m == "ice" {
                        format!("ice-{sign_name}")
                    } else {
                        m.clone()
                    };
                    line(
                        &name,
                        n,
                        nc,
                        &label,
                        t_coarse,
                        (ca - fine_area).abs() / fine_area,
                        &r,
                        t1.elapsed().as_secs_f64(),
                    );
                    if let Some(dir) = &dump {
                        let _ = std::fs::create_dir_all(dir);
                        let mut s = String::from("aniso\tlogarea\tsigma_ratio\n");
                        for x in &r.per_face {
                            s.push_str(&format!("{}\t{}\t{}\n", x[0], x[1], x[2]));
                        }
                        let file = format!("{dir}/{name}_{nc}_{label}.tsv");
                        if let Err(e) = std::fs::write(&file, s) {
                            eprintln!("{file}: {e}");
                        }
                    }
                }
                Err(e) => println!("# {name}\t{target}\t{m}\tdistortion FAILED: {e}"),
            }
        }
    }
}
