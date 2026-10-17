//! Adaptive against uniform refinement (paper Sec. 5.4, Fig. 19 and 20).
//!
//! ```text
//! cargo run --release -p ddg-intrinsic-apps --example amr_study -- <mesh> <problem> [options]
//!   <mesh>     a path below the workspace's data/ directory, or `bumpy_disk`, `sheared_disk`,
//!              `jittered_disk`
//!   <problem>  green | heat | disk
//!   --source V          the source vertex (green, heat; default 0)
//!   --heat-time F       t = F · (mean input edge length)² (default 1)
//!   --reference N       vertices of the uniform part of the reference (default 400000)
//!   --reference-adaptive N  vertices of the adaptive part of the reference (default N/4)
//!   --max-vertices N    budget of the adaptive run (default 30000)
//!   --squares           Dörfler marking on squared indicators instead of the supplement's sums
//! ```
//!
//! All computation is on one thread. The reference is an adaptive refinement with several times
//! the vertices of any adaptive run, refined further with a bound on the circumradius that gives
//! a uniform refinement with many more vertices than any uniform run: finer than both, near a
//! point source and everywhere else. Its own accuracy is estimated by a second reference of
//! half the size. Errors are relative L² errors on the surface (the reference's quadrature). For
//! every step of the adaptive run the uniform refinement reaching the same error is found by
//! interpolating (log–log) between uniform runs, and the vertex counts and times compared.

use std::path::PathBuf;
use std::time::Instant;

use ddg_core::generators;
use ddg_core::io::read_obj;
use ddg_core::{Geometry, TriMesh};
use ddg_intrinsic::{InsertionLimit, IntrinsicTriangulation, RefineOptions};
use ddg_intrinsic_apps::amr::*;

fn load(name: &str) -> (TriMesh, Geometry, Vec<[f64; 3]>) {
    let (positions, faces) = match name {
        "bumpy_disk" | "sheared_disk" => {
            let mut g = generators::bumpy_grid_disk(12, 12, if name == "bumpy_disk" { 0.3 } else { 0.0 })
                .unwrap();
            for p in &mut g.positions {
                p[0] += 1.6 * p[1];
            }
            (g.positions, g.faces)
        }
        "jittered_disk" => {
            // A flat square grid, jittered and sheared (the harmonic map of the input flips).
            let n = 10;
            let mut g = generators::grid_disk(n, n).unwrap();
            let mut rng = generators::SplitMix64::new(5);
            let h = 1.0 / n as f64;
            for p in &mut g.positions {
                let inside = p[0] > 1e-9 && p[0] < 1.0 - 1e-9 && p[1] > 1e-9 && p[1] < 1.0 - 1e-9;
                if inside {
                    p[0] += rng.uniform(-0.45, 0.45) * h;
                    p[1] += rng.uniform(-0.45, 0.45) * h;
                }
                p[0] += 2.5 * p[1];
            }
            (g.positions, g.faces)
        }
        _ => {
            let path = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../data").join(name);
            let obj = read_obj(&path).unwrap_or_else(|e| panic!("{}: {e}", path.display()));
            let mut faces = Vec::new();
            for f in 0..obj.n_faces() {
                let c = &obj.corner_verts[obj.face_corners(f)];
                for k in 1..c.len() - 1 {
                    faces.push([c[0], c[k], c[k + 1]]);
                }
            }
            (obj.positions, faces)
        }
    };
    let mesh = TriMesh::new(positions.len(), faces).unwrap();
    let geometry = Geometry::from_positions(&mesh, &positions).unwrap();
    (mesh, geometry, positions)
}

struct Run {
    vertices: usize,
    seconds: f64,
    error: f64,
}

fn uniform(
    mesh: &TriMesh,
    geometry: &Geometry,
    problem: &Problem,
    target_vertices: usize,
) -> (IntrinsicTriangulation, Solution, f64) {
    let t = Instant::now();
    let mut tri = IntrinsicTriangulation::new(mesh, geometry).unwrap();
    let r = (geometry.total_area() / (2.6 * target_vertices as f64)).sqrt();
    tri.refine(&RefineOptions {
        max_circumradius: Some(r),
        max_insertions: InsertionLimit::Unlimited,
        ..Default::default()
    })
    .unwrap();
    let sol = solve(&tri, problem).unwrap();
    (tri, sol, t.elapsed().as_secs_f64())
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    if args.len() < 3 {
        eprintln!("usage: amr_study <mesh> <green|heat|disk> [options]");
        std::process::exit(2);
    }
    let mut source = 0u32;
    let mut heat_factor = 1.0;
    let mut reference_vertices = 400_000usize;
    let mut reference_adaptive: Option<usize> = None;
    let mut max_vertices = 30_000usize;
    let mut bulk_sum = BulkSum::Indicators;
    let mut i = 3;
    while i < args.len() {
        match args[i].as_str() {
            "--source" => source = args[i + 1].parse().unwrap(),
            "--heat-time" => heat_factor = args[i + 1].parse().unwrap(),
            "--reference" => reference_vertices = args[i + 1].parse().unwrap(),
            "--reference-adaptive" => reference_adaptive = Some(args[i + 1].parse().unwrap()),
            "--max-vertices" => max_vertices = args[i + 1].parse().unwrap(),
            "--squares" => {
                bulk_sum = BulkSum::Squares;
                i += 1;
                continue;
            }
            other => panic!("unknown option {other}"),
        }
        i += 2;
    }
    let (mesh, geometry, _) = load(&args[1]);
    let mean_edge = geometry.edge_lengths().iter().sum::<f64>() / geometry.n_edges() as f64;
    let problem = match args[2].as_str() {
        "green" => Problem::GreenFunction { source },
        "heat" => Problem::HeatKernel {
            source,
            time: heat_factor * mean_edge * mean_edge,
        },
        "disk" => Problem::HarmonicDisk,
        other => panic!("unknown problem {other}"),
    };
    let exclude: Vec<u32> = match problem {
        Problem::HarmonicDisk => vec![],
        _ => vec![source],
    };
    println!(
        "mesh {} ({} vertices, {} faces), problem {:?}, marking {:?}",
        args[1],
        mesh.n_vertices(),
        mesh.n_faces(),
        problem,
        bulk_sum
    );
    let input_tri = IntrinsicTriangulation::new(&mesh, &geometry).unwrap();

    // The reference, and a second one of half the size to estimate its accuracy.
    let build_reference = |uniform_vertices: usize, adaptive_vertices: usize| {
        let t = Instant::now();
        let mut tri = IntrinsicTriangulation::new(&mesh, &geometry).unwrap();
        adaptive_solve(
            &mut tri,
            &problem,
            &AmrOptions {
                bulk_sum,
                max_vertices: adaptive_vertices,
                max_iterations: 1000,
                ..Default::default()
            },
        )
        .unwrap();
        let r = (geometry.total_area() / (2.6 * uniform_vertices as f64)).sqrt();
        tri.refine(&RefineOptions {
            max_circumradius: Some(r),
            max_insertions: InsertionLimit::Unlimited,
            ..Default::default()
        })
        .unwrap();
        let sol = solve(&tri, &problem).unwrap();
        (tri, sol, t.elapsed().as_secs_f64())
    };
    let adaptive_reference = reference_adaptive.unwrap_or(reference_vertices / 4);
    let (ref_tri, reference, ref_seconds) = build_reference(reference_vertices, adaptive_reference);
    let (half_tri, half, _) = build_reference(reference_vertices / 2, adaptive_reference / 2);
    println!(
        "reference: {} vertices, {:.2} s; the half-size reference ({} vertices) differs by {:.3e}",
        ref_tri.n_vertices(),
        ref_seconds,
        half_tri.n_vertices(),
        l2_error(&ref_tri, &reference, &half_tri, &half).unwrap()
    );

    // The input mesh and its intrinsic Delaunay triangulation.
    for (what, delaunay) in [("input", false), ("iDT", true)] {
        let t = Instant::now();
        let mut tri = IntrinsicTriangulation::new(&mesh, &geometry).unwrap();
        if delaunay {
            tri.flip_to_delaunay(&Default::default()).unwrap();
        }
        let sol = solve(&tri, &problem).unwrap();
        let s = t.elapsed().as_secs_f64();
        let err = l2_error(&ref_tri, &reference, &tri, &sol).unwrap();
        let verr = input_vertex_error(&input_tri, &sol, &reference, &exclude).unwrap();
        let est = estimate(&tri, &problem, &sol).unwrap();
        print!(
            "{what}: {} vertices, {:.3} s, L2 error {:.4e}, error at input vertices {:.4e}, estimate {:.4e}",
            tri.n_vertices(),
            s,
            err,
            verr,
            est.total
        );
        if matches!(problem, Problem::HarmonicDisk) {
            let q = map_quality(&tri, &sol);
            print!(
                ", flipped {}, conformal mean {:.4} q99 {:.3} max {:.3}",
                q.n_flipped, q.mean_conformal, q.conformal_q99, q.max_conformal
            );
        }
        println!();
    }

    // Uniform runs.
    let mut uniform_runs: Vec<Run> = Vec::new();
    let mut target = (mesh.n_vertices() as f64 * 1.5) as usize;
    while target * 2 <= reference_vertices {
        let (tri, sol, seconds) = uniform(&mesh, &geometry, &problem, target);
        let error = l2_error(&ref_tri, &reference, &tri, &sol).unwrap();
        let verr = input_vertex_error(&input_tri, &sol, &reference, &exclude).unwrap();
        let est = estimate(&tri, &problem, &sol).unwrap();
        println!(
            "uniform: {:>7} vertices {:>8.3} s L2 error {:.4e} vertex error {:.4e} estimate {:.4e}",
            tri.n_vertices(),
            seconds,
            error,
            verr,
            est.total
        );
        if matches!(problem, Problem::HarmonicDisk) {
            let q = map_quality(&tri, &sol);
            println!(
                "         flipped {}, conformal mean {:.4} q99 {:.3} max {:.3}",
                q.n_flipped, q.mean_conformal, q.conformal_q99, q.max_conformal
            );
        }
        uniform_runs.push(Run {
            vertices: tri.n_vertices(),
            seconds,
            error,
        });
        target = (target as f64 * 1.6) as usize;
    }

    // The adaptive run, step by step.
    let options = AmrOptions {
        bulk_sum,
        ..Default::default()
    };
    let mut t = Instant::now();
    let mut tri = IntrinsicTriangulation::new(&mesh, &geometry).unwrap();
    tri.refine(&options.refine).unwrap();
    let mut amr_runs: Vec<Run> = Vec::new();
    let mut last_quality = None;
    loop {
        let sol = solve(&tri, &problem).unwrap();
        let est = estimate(&tri, &problem, &sol).unwrap();
        let seconds = t.elapsed().as_secs_f64();
        let paused = Instant::now();
        let error = l2_error(&ref_tri, &reference, &tri, &sol).unwrap();
        let verr = input_vertex_error(&input_tri, &sol, &reference, &exclude).unwrap();
        println!(
            "adaptive: {:>7} vertices {:>8.3} s L2 error {:.4e} vertex error {:.4e} estimate {:.4e} oscillation {:.3e}",
            tri.n_vertices(),
            seconds,
            error,
            verr,
            est.total,
            est.total_oscillation
        );
        if matches!(problem, Problem::HarmonicDisk) {
            let q = map_quality(&tri, &sol);
            println!(
                "         flipped {}, conformal mean {:.4} q99 {:.3} max {:.3}",
                q.n_flipped, q.mean_conformal, q.conformal_q99, q.max_conformal
            );
            last_quality = Some(q);
        }
        amr_runs.push(Run {
            vertices: tri.n_vertices(),
            seconds,
            error,
        });
        let evaluation = paused.elapsed();
        t += evaluation;
        let marked = mark(&est, options.bulk, options.bulk_sum);
        if tri.n_vertices() + marked.len() > max_vertices {
            break;
        }
        refine_marked(&mut tri, &options.refine, &marked).unwrap();
    }
    let _ = last_quality;

    // Equal accuracy.
    println!();
    println!("| adaptive V | adaptive s | error | uniform V (same error) | uniform s | V ratio | time ratio |");
    println!("|---:|---:|---:|---:|---:|---:|---:|");
    for a in &amr_runs {
        let mut found = None;
        for w in uniform_runs.windows(2) {
            let (p, q) = (&w[0], &w[1]);
            if p.error >= a.error && a.error >= q.error {
                let s = (a.error.ln() - p.error.ln()) / (q.error.ln() - p.error.ln());
                let v = ((p.vertices as f64).ln() * (1.0 - s) + (q.vertices as f64).ln() * s).exp();
                let time = (p.seconds.ln() * (1.0 - s) + q.seconds.ln() * s).exp();
                found = Some((v, time));
            }
        }
        match found {
            Some((v, time)) => println!(
                "| {} | {:.3} | {:.3e} | {:.0} | {:.3} | {:.1} | {:.1} |",
                a.vertices,
                a.seconds,
                a.error,
                v,
                time,
                v / a.vertices as f64,
                time / a.seconds
            ),
            None => println!(
                "| {} | {:.3} | {:.3e} | (outside the uniform runs) | | | |",
                a.vertices, a.seconds, a.error
            ),
        }
    }
}
