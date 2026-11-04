//! WP14 density study (SPEC.md §16.4; paper Fig. 7 and Fig. 29).
//!
//! ```text
//! density <mesh.obj> <out-dir> [--counts 1024,8192,65536] [--grid 128] [--mc 0|res] [--plane 512]
//!                                [--threads T] [--seed S]
//! ```
//!
//! The mesh is centered (vertex mean) and scaled to `[-1, 1]³`, checked closed and outward, and
//! sampled uniformly by area at each count. For each count: the default PAT field (32-NN,
//! `PerQuery` C = 64, `SignRule::Reference`), its error on the `grid³` grid over `[-1, 1]³`
//! against the exact mesh SDF (overall, narrow band `|φ_true| < 0.05`, interior `φ_true < −0.05`,
//! exterior `φ_true > 0.05`), the torus size statistics in units of the point's `σ_i` (the paper's
//! hypothesis for Fig. 29: dense neighborhoods give compact tori), optionally the zero level set by
//! marching cubes (`--mc`, OBJ written to the output directory), and `plane × plane` SDF values on
//! the mid-plane `z = 0` (`slices.npy`, truth first) for the figures. Writes `density.json`.

#[path = "eval_util/util.rs"]
mod util;

use ddg_points_as_tori::generate::rng::Rng;
use ddg_points_as_tori::generate::sampling::{center_and_scale, sample_mesh_uniform};
use ddg_points_as_tori::mc::{Band, McOptions, MeshCheck, Resolution, marching_cubes};
use ddg_points_as_tori::truth::{MeshTruth, TruthOptions, check_mesh};
use ddg_points_as_tori::{BuildOptions, EvalOut, Model, NormalPolicy, PointCloud, PointsAsTori};
use std::path::PathBuf;
use std::time::Instant;
use util::*;

fn main() {
    if let Err(e) = run() {
        eprintln!("error: {e}");
        std::process::exit(1);
    }
}

fn opt<T: std::str::FromStr>(a: &[String], k: &str, d: T) -> R<T> {
    match a.iter().position(|x| x == k) {
        None => Ok(d),
        Some(i) => a.get(i + 1).ok_or(format!("{k} needs a value"))?.parse().map_err(|_| format!("{k}: invalid value")),
    }
}

fn run() -> R<()> {
    let a: Vec<String> = std::env::args().skip(1).collect();
    if a.len() < 2 {
        return Err("usage: density <mesh.obj> <out-dir> [--counts ..] [--grid 128] [--mc 0] [--plane 512] [--threads T] [--seed S]".into());
    }
    let mesh_path = PathBuf::from(&a[0]);
    let out = PathBuf::from(&a[1]);
    std::fs::create_dir_all(&out).map_err(err)?;
    let counts: Vec<usize> = opt::<String>(&a, "--counts", "1024,8192,65536".into())?.split(',').map(|s| s.parse().map_err(err)).collect::<R<_>>()?;
    let g: usize = opt(&a, "--grid", 128)?;
    let mc: usize = opt(&a, "--mc", 0)?;
    let plane: usize = opt(&a, "--plane", 512)?;
    let t: usize = opt(&a, "--threads", 0)?;
    let th = if t == 0 { None } else { Some(t) };
    let seed: u64 = opt(&a, "--seed", 7)?;
    let pool = rayon::ThreadPoolBuilder::new().num_threads(t).build().map_err(err)?;

    let obj = ddg_core::io::read_obj(&mesh_path).map_err(err)?;
    let f = obj.triangles().ok_or("the mesh has non-triangle faces")?;
    let v = center_and_scale(&obj.positions);
    let rep = check_mesh(&v, &f).map_err(err)?;
    if !rep.is_closed() || !(rep.signed_volume > 0.0) {
        return Err(format!("the mesh is not closed and outward: {rep:?}"));
    }
    let t0 = Instant::now();
    let (truth, _) = MeshTruth::new(&v, &f, &TruthOptions { threads: th, ..Default::default() }).map_err(err)?;
    let xs = grid(g);
    let mut tv = vec![0.0; xs.len()];
    truth.evaluate(&xs, &mut tv, None).map_err(err)?;
    let ax = axis(plane);
    let pxs: Vec<[f64; 3]> = ax.iter().flat_map(|&y| ax.iter().map(move |&x| [x, y, 0.0])).collect();
    let mut ptv = vec![0.0; pxs.len()];
    truth.evaluate(&pxs, &mut ptv, None).map_err(err)?;
    let mut sl = ptv.clone();
    eprintln!("truth {} faces, {:.1}s", f.len(), t0.elapsed().as_secs_f64());
    let model = Model::reference().map_err(err)?;
    let mut rows = Vec::new();
    for &n in &counts {
        let t1 = Instant::now();
        let mut rng = Rng::new(seed ^ n as u64);
        let (p, nr) = sample_mesh_uniform(&v, &f, n, &mut rng).map_err(err)?;
        let (cloud, _) = PointCloud::new(p, nr, NormalPolicy::Normalize).map_err(err)?;
        let (field, br) = PointsAsTori::build(&cloud, &model, &BuildOptions::default().with_threads(th)).map_err(err)?;
        let build_s = t1.elapsed().as_secs_f64();
        let mut phi = vec![0.0; xs.len()];
        let er = field.evaluate(&xs, EvalOut { distance: Some(&mut phi), ..Default::default() }).map_err(err)?;
        let m = metrics(g, &phi, &tv);
        let (mut si, mut ni, mut se, mut ne) = (0.0, 0usize, 0.0, 0usize);
        for (x, y) in phi.iter().zip(&tv) {
            if *y < -BAND {
                si += (x - y).abs();
                ni += 1;
            } else if *y > BAND {
                se += (x - y).abs();
                ne += 1;
            }
        }
        let mut pp = vec![0.0; pxs.len()];
        field.evaluate(&pxs, EvalOut { distance: Some(&mut pp), ..Default::default() }).map_err(err)?;
        sl.extend_from_slice(&pp);
        // torus size in units of σ_i
        let tori = field.tori();
        let mut big: Vec<f64> = tori.tori.iter().zip(&tori.sigma).map(|(t, s)| t.major / s).collect();
        let mut small: Vec<f64> = tori.tori.iter().zip(&tori.sigma).map(|(t, s)| t.minor.abs() / s).collect();
        let spindle = tori.tori.iter().filter(|t| t.major < t.minor.abs()).count();
        let neg = tori.tori.iter().filter(|t| t.minor < 0.0).count();
        let mut sig = tori.sigma.clone();
        let mut j = Json::new();
        j.int("n", n as i64)
            .raw("m", &metrics_json(&m))
            .num("interior_mae", si / ni.max(1) as f64)
            .num("exterior_mae", se / ne.max(1) as f64)
            .num("build_s", build_s)
            .num("pat_build_s", br.seconds)
            .num("grid_eval_s", er.seconds)
            .num("sigma_median", quantile(&mut sig, 0.5))
            .num("major_over_sigma_median", quantile(&mut big, 0.5))
            .num("major_over_sigma_p10", quantile(&mut big, 0.1))
            .num("minor_over_sigma_median", quantile(&mut small, 0.5))
            .num("minor_over_sigma_p10", quantile(&mut small, 0.1))
            .num("spindle_frac", spindle as f64 / n as f64)
            .num("negative_minor_frac", neg as f64 / n as f64)
            .int("saddles", br.fit.saddle as i64);
        if mc > 0 {
            let (mesh, mr) = pool
                .install(|| {
                    marching_cubes(
                        &field,
                        &McOptions { bbox: Some([[-1.1; 3], [1.1; 3]]), resolution: Resolution::Cells([mc; 3]), band: Band::Dense, threads: th, ..Default::default() },
                    )
                })
                .map_err(err)?;
            if let MeshCheck::Manifold { components, euler_characteristic, .. } = mr.mesh {
                j.int("mc_components", components as i64).int("mc_chi", euler_characteristic);
            }
            let tri = mesh.to_trimesh().ok();
            let mut s = String::new();
            for p in &mesh.positions {
                s.push_str(&format!("v {} {} {}\n", p[0], p[1], p[2]));
            }
            for t in &mesh.faces {
                s.push_str(&format!("f {} {} {}\n", t[0] + 1, t[1] + 1, t[2] + 1));
            }
            std::fs::write(out.join(format!("mc_{n}.obj")), s).map_err(err)?;
            j.int("mc_triangles", mesh.faces.len() as i64).bool("mc_trimesh_ok", tri.is_some());
        }
        eprintln!("n = {n}: MAE {:.5} interior {:.5} exterior {:.5} band {:.5} (build {:.1}s)", m.mae, si / ni.max(1) as f64, se / ne.max(1) as f64, m.band_mae, build_s);
        rows.push(j.finish());
    }
    write_npy_f32(&out.join("slices.npy"), &[counts.len() + 1, plane, plane], &sl)?;
    let res = Json::new()
        .str("mesh", &mesh_path.display().to_string())
        .int("faces", f.len() as i64)
        .int("grid", g as i64)
        .int("plane", plane as i64)
        .str("threads", &format!("{th:?}"))
        .raw("rows", &jarr(&rows))
        .finish();
    std::fs::write(out.join("density.json"), res).map_err(err)?;
    Ok(())
}
