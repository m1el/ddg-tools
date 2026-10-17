//! Re-optimized angles at the published cone positions of the dataset, and the time of the
//! direct solve against the generic solve of the same problem.
//!
//! `cargo run --release -p ddg-optimal-cones --example reoptimize [-- <mesh>...]`: one line per
//! mesh, sequentially and on one thread (the times are single-thread times; run nothing else
//! meanwhile). `--repeat N` takes the fastest of N runs.
use std::f64::consts::PI;
use std::path::PathBuf;
use std::time::Instant;

use ddg_core::io::{read_cones, read_obj};
use ddg_core::{Geometry, TriMesh};
use ddg_optimal_cones::*;

fn best<T>(repeat: usize, mut f: impl FnMut() -> T) -> (T, f64) {
    let mut out = None;
    let mut fastest = f64::INFINITY;
    for _ in 0..repeat {
        let start = Instant::now();
        let x = f();
        fastest = fastest.min(start.elapsed().as_secs_f64());
        out = Some(x);
    }
    (out.unwrap(), fastest)
}

fn main() {
    let mut args: Vec<String> = std::env::args().skip(1).collect();
    let mut repeat = 1;
    if let Some(i) = args.iter().position(|a| a == "--repeat") {
        repeat = args[i + 1].parse().unwrap();
        args.drain(i..i + 2);
    }
    let dir = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../data/optimal-cones/optimal-cones-data/data");
    let mut names: Vec<String> = if args.is_empty() {
        std::fs::read_dir(&dir)
            .unwrap()
            .filter_map(|e| e.ok())
            .filter_map(|e| {
                e.file_name()
                    .to_str()
                    .and_then(|n| n.strip_suffix(".obj"))
                    .map(str::to_string)
            })
            .collect()
    } else {
        args
    };
    names.sort();
    println!(
        "| file (cones) | D, published angles | D, re-optimized | ratio | largest angle change | \
         D, published, boundary optimized | D, re-optimized, boundary optimized | ratio | \
         reoptimize_cones (s) | method | solve, costs 0/∞ (s) |"
    );
    println!("|---|---:|---:|---:|---:|---:|---:|---:|---:|---|---:|");
    for name in &names {
        let obj = read_obj(dir.join(format!("{name}.obj"))).unwrap();
        let mesh = TriMesh::new(obj.positions.len(), obj.triangles().unwrap()).unwrap();
        let geometry = Geometry::from_positions(&mesh, &obj.positions).unwrap();
        let published =
            read_cones(dir.join(format!("{name}-cones.txt")), mesh.n_vertices()).unwrap();
        let problem = ConeProblem::new(&mesh, &geometry, 0.5);
        let vertices: Vec<u32> = published.iter().map(|c| c.vertex).collect();
        let closed = mesh.n_boundary_edges() == 0;
        let mut angles = published.clone();
        if closed {
            let defect: f64 = angles.iter().map(|c| c.angle).sum::<f64>()
                - 2.0 * PI * mesh.euler_characteristic() as f64;
            for c in &mut angles {
                c.angle -= defect / published.len() as f64;
            }
        }
        let evaluated = evaluate_cones(&problem, &angles, Parallelism::Sequential).unwrap();
        let (r, seconds) = best(repeat, || {
            reoptimize_cones(&problem, &vertices, &ReoptimizeOptions::default()).unwrap()
        });
        assert_eq!(r.report.status, Status::Converged, "{name}");
        let reopt = reoptimization_problem(&problem, &vertices).unwrap();
        let (s, seconds_solve) = best(repeat, || solve(&reopt, &SolveOptions::default()).unwrap());
        assert_eq!(s.report.status, Status::Converged, "{name}");
        let max_change = r
            .cones
            .iter()
            .zip(&published)
            .map(|(x, y)| (x.angle - y.angle).abs())
            .fold(0.0, f64::max);
        let free = if closed {
            "– | – | –".to_string()
        } else {
            let fixed = ConeProblem {
                cone_costs: Some(vec![f64::INFINITY; mesh.n_vertices()]),
                prescribed_cones: published.clone(),
                ..problem.clone()
            };
            let p = solve_with_boundary(&fixed, BoundaryTreatment::Free, &SolveOptions::default())
                .unwrap();
            let o = ReoptimizeOptions {
                boundary: BoundaryTreatment::Free,
                ..ReoptimizeOptions::default()
            };
            let q = reoptimize_cones(&problem, &vertices, &o).unwrap();
            format!(
                "{:.4e} | {:.4e} | {:.4}",
                p.distortion_optimal,
                q.distortion,
                q.distortion / p.distortion_optimal
            )
        };
        println!(
            "| {name}{} ({}) | {:.4e} | {:.4e} | {:.4} | {:.1e} | {free} | {:.3} | {:?} | {:.3} |",
            if closed { ", closed" } else { "" },
            published.len(),
            evaluated.distortion,
            r.distortion,
            r.distortion / evaluated.distortion,
            max_change,
            seconds,
            r.report.method,
            seconds_solve
        );
    }
}
