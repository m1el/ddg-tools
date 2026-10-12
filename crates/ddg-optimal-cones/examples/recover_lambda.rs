//! Recovers the parameter `λ` behind each cone file of the optimal-cones dataset.
//!
//! ```text
//! cargo run --release -p ddg-optimal-cones --example recover_lambda -- [--data DIR] [names...]
//! ```
//!
//! If a cone file is an exact optimum of the canonical problem (unit weights, no bounds), the
//! adjoint state `φ` computed from its cones satisfies `φ_i = λ·sign(µ_i)` at every cone and
//! `|φ_i| ≤ λ` elsewhere; on a closed mesh the same holds for `φ + c` with one constant `c`.
//! The example prints, per mesh, the `λ` (and `c`) that fit the cones best, the spread of
//! `sign(µ_i)·(φ_i + c)` over the cones, and the largest violation of `|φ + c| ≤ λ` elsewhere.

use std::path::PathBuf;

use ddg_core::io::{read_cones, read_obj};
use ddg_core::{Geometry, TriMesh};
use ddg_optimal_cones::{ConeProblem, Parallelism, evaluate_cones};

fn main() {
    let mut dir = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../data/optimal-cones/optimal-cones-data/data");
    let mut names: Vec<String> = Vec::new();
    let mut args = std::env::args().skip(1);
    while let Some(arg) = args.next() {
        match arg.as_str() {
            "--data" => dir = PathBuf::from(args.next().expect("--data needs a value")),
            name => names.push(name.to_string()),
        }
    }
    if !dir.is_dir() {
        eprintln!("dataset directory {} not found", dir.display());
        std::process::exit(2);
    }
    if names.is_empty() {
        names = std::fs::read_dir(&dir)
            .expect("readable dataset directory")
            .filter_map(|e| e.ok())
            .filter_map(|e| {
                e.file_name()
                    .to_str()
                    .and_then(|n| n.strip_suffix(".obj"))
                    .map(str::to_string)
            })
            .collect();
        names.sort();
    }
    println!(
        "{:<20} {:>7} {:>6} {:>6} {:>5} {:>5} {:>11} {:>11} {:>10} {:>10} {:>10} {:>10} {:>9}",
        "mesh",
        "V",
        "closed",
        "cones",
        "pos",
        "neg",
        "lambda",
        "shift c",
        "min s*phi",
        "max s*phi",
        "spread/lam",
        "excess/lam",
        "GB defect"
    );
    for name in &names {
        let obj = read_obj(dir.join(format!("{name}.obj"))).unwrap_or_else(|e| panic!("{e}"));
        let faces = obj.triangles().expect("triangle mesh");
        let mesh = TriMesh::new(obj.positions.len(), faces).unwrap_or_else(|e| panic!("{e}"));
        let geometry =
            Geometry::from_positions(&mesh, &obj.positions).unwrap_or_else(|e| panic!("{e}"));
        let cones = read_cones(dir.join(format!("{name}-cones.txt")), mesh.n_vertices())
            .unwrap_or_else(|e| panic!("{e}"));
        let problem = ConeProblem::new(&mesh, &geometry, 0.0);
        let eval = evaluate_cones(&problem, &cones, Parallelism::Sequential)
            .unwrap_or_else(|e| panic!("{name}: {e}"));
        let closed = mesh.boundary_loops().is_empty();
        let phi = &eval.adjoint;
        let pos: Vec<f64> = cones
            .iter()
            .filter(|c| c.angle > 0.0)
            .map(|c| phi[c.vertex as usize])
            .collect();
        let neg: Vec<f64> = cones
            .iter()
            .filter(|c| c.angle < 0.0)
            .map(|c| phi[c.vertex as usize])
            .collect();
        let mean = |x: &[f64]| x.iter().sum::<f64>() / x.len() as f64;
        // Least squares fit of  φ_i + c = λ s_i  over the cones.
        let (lambda, shift) = if !closed {
            let all: Vec<f64> = cones
                .iter()
                .map(|c| c.angle.signum() * phi[c.vertex as usize])
                .collect();
            (mean(&all), 0.0)
        } else if !pos.is_empty() && !neg.is_empty() {
            (
                0.5 * (mean(&pos) - mean(&neg)),
                -0.5 * (mean(&pos) + mean(&neg)),
            )
        } else {
            // One sign only: λ and c are not separately identifiable from the cones. Take the
            // smallest λ for which |φ + c| ≤ λ can hold everywhere: the cones sit at the
            // maximum (resp. minimum) of φ, so λ = (max φ − min φ) / 2.
            let (lo, hi) = phi
                .iter()
                .fold((f64::INFINITY, f64::NEG_INFINITY), |(lo, hi), &p| {
                    (lo.min(p), hi.max(p))
                });
            if neg.is_empty() {
                let top = mean(&pos);
                (0.5 * (top - lo), -0.5 * (top + lo))
            } else {
                let bottom = mean(&neg);
                (0.5 * (hi - bottom), -0.5 * (hi + bottom))
            }
        };
        let signed: Vec<f64> = cones
            .iter()
            .map(|c| c.angle.signum() * (phi[c.vertex as usize] + shift))
            .collect();
        let min = signed.iter().cloned().fold(f64::INFINITY, f64::min);
        let max = signed.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
        let mut is_cone = vec![false; mesh.n_vertices()];
        for c in &cones {
            is_cone[c.vertex as usize] = true;
        }
        let excess = (0..mesh.n_vertices())
            .filter(|&v| !is_cone[v] && !mesh.is_boundary_vertex(v))
            .map(|v| (phi[v] + shift).abs() - lambda)
            .fold(f64::NEG_INFINITY, f64::max);
        println!(
            "{:<20} {:>7} {:>6} {:>6} {:>5} {:>5} {:>11.6} {:>11.3e} {:>10.6} {:>10.6} {:>10.2e} {:>10.2e} {:>9.1e}",
            name,
            mesh.n_vertices(),
            closed,
            cones.len(),
            pos.len(),
            neg.len(),
            lambda,
            shift,
            min,
            max,
            (max - min) / lambda,
            excess / lambda,
            eval.gauss_bonnet_defect,
        );
    }
}
