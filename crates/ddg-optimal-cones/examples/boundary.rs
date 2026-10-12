//! Experiments with the boundary treatments on the published dataset.
//!
//! ```text
//! cargo run --release -p ddg-optimal-cones --example boundary -- table [lambda]
//! cargo run --release -p ddg-optimal-cones --example boundary -- spikes
//! ```
//!
//! * `table`: for every mesh with boundary, the optimum at `λ` (default 0.5) with `u = 0` on the
//!   boundary and with optimized boundary values (and the penalized and convex variants):
//!   distortion, objective, number of cones.
//! * `spikes`: the weights hypothesis for `camelshotnoise` (see the README, "Needle-like
//!   spikes"): the angles at the published positions re-optimized under various distortion
//!   weights, with `u = 0` and with optimized boundary values, against the published angles.

use std::path::PathBuf;

use ddg_core::io::{read_cones, read_obj};
use ddg_core::operators::{angle_defects, lumped_mass_normalized};
use ddg_core::{Geometry, TriMesh};
use ddg_optimal_cones::{
    BoundaryTreatment, Cone, ConeProblem, FeatureSizeConvention, Solution, SolveOptions, Status,
    feature_size_weights, solve_with_boundary,
};

struct Data {
    mesh: TriMesh,
    geometry: Geometry,
    cones: Vec<Cone>,
}

fn dir() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../data/optimal-cones/optimal-cones-data/data")
}

fn load(name: &str) -> Data {
    let obj = read_obj(dir().join(format!("{name}.obj"))).unwrap();
    let mesh = TriMesh::new(obj.positions.len(), obj.triangles().unwrap()).unwrap();
    let geometry = Geometry::from_positions(&mesh, &obj.positions).unwrap();
    let cones = read_cones(dir().join(format!("{name}-cones.txt")), mesh.n_vertices()).unwrap();
    Data {
        mesh,
        geometry,
        cones,
    }
}

fn run(p: &ConeProblem<'_>, b: BoundaryTreatment) -> Solution {
    solve_with_boundary(p, b, &SolveOptions::default()).unwrap()
}

fn table(lambda: f64) {
    println!(
        "lambda = {lambda}; distortion = ½ Σ wE m u² of the optimum, cones = number of cones of the optimum"
    );
    println!(
        "{:<20} {:>7} | {:>10} {:>5} | {:>10} {:>5} {:>7} | {:>10} {:>5} | {:>10} {:>5}",
        "mesh",
        "V",
        "D (u=0)",
        "cones",
        "D (free)",
        "cones",
        "ratio",
        "D (pen.)",
        "cones",
        "D (conv.)",
        "cones"
    );
    let mut names: Vec<String> = std::fs::read_dir(dir())
        .unwrap()
        .filter_map(|e| {
            e.ok()?
                .file_name()
                .to_str()?
                .strip_suffix(".obj")
                .map(str::to_string)
        })
        .collect();
    names.sort();
    for name in names {
        let d = load(&name);
        if d.mesh.boundary_loops().is_empty() {
            continue;
        }
        let p = ConeProblem::new(&d.mesh, &d.geometry, lambda);
        let count = |s: &Solution| {
            s.measure
                .iter()
                .enumerate()
                .filter(|(v, m)| **m != 0.0 && !d.mesh.is_boundary_vertex(*v))
                .count()
        };
        let mut cells = Vec::new();
        for b in [
            BoundaryTreatment::Fixed,
            BoundaryTreatment::Free,
            BoundaryTreatment::Penalized,
            BoundaryTreatment::Convex,
        ] {
            let s = run(&p, b);
            let status = if s.report.status == Status::Converged {
                ""
            } else {
                "!"
            };
            cells.push((s.distortion_optimal, count(&s), status, s.objective));
        }
        println!(
            "{:<20} {:>7} | {:>10.5} {:>5} {:>9.5} | {:>10.5} {:>5} {:>9.5} {:>7.3} | {:>10.5} {:>5} | {:>10.5} {:>5} {}",
            name,
            d.mesh.n_vertices(),
            cells[0].0,
            cells[0].1,
            cells[0].3,
            cells[1].0,
            cells[1].1,
            cells[1].3,
            cells[1].0 / cells[0].0,
            cells[2].0,
            cells[2].1,
            cells[3].0,
            cells[3].1,
            cells.iter().map(|c| c.2).collect::<String>()
        );
    }
}

/// Re-optimizes the angles at `cones` (free sites there, forbidden elsewhere) with the distortion
/// weights `w`; returns the largest difference from the file.
fn refit(d: &Data, w: Option<Vec<f64>>, b: BoundaryTreatment) -> (f64, Solution) {
    let mut costs = vec![f64::INFINITY; d.mesh.n_vertices()];
    for c in &d.cones {
        costs[c.vertex as usize] = 0.0;
    }
    let p = ConeProblem {
        cone_costs: Some(costs),
        distortion_weights: w,
        ..ConeProblem::new(&d.mesh, &d.geometry, 1.0)
    };
    let s = run(&p, b);
    let worst = d
        .cones
        .iter()
        .map(|c| (s.measure[c.vertex as usize] - c.angle).abs())
        .fold(0.0, f64::max);
    if s.report.status != Status::Converged {
        println!("    (not converged: {:?})", s.report.message);
    }
    (worst, s)
}

fn spikes() {
    let only = std::env::args().nth(2);
    for name in ["camel", "camelshotnoise"] {
        if only.as_deref().is_some_and(|o| o != name) {
            continue;
        }
        let d = load(name);
        let nv = d.mesh.n_vertices();
        let omega = angle_defects(&d.mesh, &d.geometry).unwrap();
        let mass = lumped_mass_normalized(&d.mesh, &d.geometry).unwrap();
        let mut dens: Vec<f64> = (0..nv)
            .filter(|&v| !d.mesh.is_boundary_vertex(v))
            .map(|v| omega[v].abs() / mass[v])
            .collect();
        dens.sort_by(f64::total_cmp);
        let q = |f: f64| dens[((dens.len() - 1) as f64 * f) as usize];
        let tips: Vec<usize> = (0..nv).filter(|&v| omega[v] > 6.0).collect();
        println!(
            "== {name}: {nv} vertices, {} published cones, {} needle tips (defect > 6)",
            d.cones.len(),
            tips.len()
        );
        println!(
            "   |Omega|/M over interior vertices: min {:.3e}, median {:.3e}, 99% {:.3e}, max {:.3e}",
            dens[0],
            q(0.5),
            q(0.99),
            dens[dens.len() - 1]
        );
        let mut settings: Vec<(String, Option<Vec<f64>>)> = vec![("uniform".into(), None)];
        if !tips.is_empty() {
            let mut w = vec![1.0; nv];
            for &t in &tips {
                w[t] = 0.0;
            }
            settings.push(("zero weight at the tips".into(), Some(w)));
        }
        if !tips.is_empty() {
            let mut w = vec![1.0; nv];
            for &t in &tips {
                w[t] = 0.0;
                for n in d.mesh.vertex_neighbors(t) {
                    w[n] = 0.0;
                }
            }
            settings.push(("zero weight at the tips and their 1-rings".into(), Some(w)));
        }
        let median = q(0.5);
        let mut eps_list: Vec<f64> = [1e-2, 1.0, 1e2, 1e4].iter().map(|f| f * median).collect();
        eps_list.push(dens[dens.len() - 1]);
        for &eps in &eps_list {
            for (conv, label) in [
                (FeatureSizeConvention::AsPrinted, "printed"),
                (FeatureSizeConvention::AsDescribed, "described"),
            ] {
                let f = feature_size_weights(&d.mesh, &d.geometry, eps, conv).unwrap();
                settings.push((
                    format!("feature size {label}, eps = {eps:.2e}"),
                    Some(f.distortion_weights),
                ));
            }
        }
        println!(
            "   {:<48} {:>12} {:>12}",
            "distortion weights", "u = 0", "optimized"
        );
        let mut best = (f64::INFINITY, String::new());
        for (label, w) in &settings {
            let (fixed, _) = refit(&d, w.clone(), BoundaryTreatment::Fixed);
            let (free, _) = refit(&d, w.clone(), BoundaryTreatment::Free);
            println!("   {:<48} {:>12.3e} {:>12.3e}", label, fixed, free);
            if free < best.0 {
                best = (free, label.clone());
            }
        }
        println!(
            "   closest with optimized boundary values: {} ({:.3e})",
            best.1, best.0
        );
        // Penalized solve with the best weights (feature-size ones carry wR = 1/wE too).
        if name == "camelshotnoise" {
            let published: std::collections::HashSet<u32> =
                d.cones.iter().map(|c| c.vertex).collect();
            let ring: std::collections::HashSet<usize> = tips
                .iter()
                .flat_map(|&t| std::iter::once(t).chain(d.mesh.vertex_neighbors(t)))
                .collect();
            // candidates: the closest setting, and zero weight at tips
            let mut cands: Vec<(String, ConeProblem<'_>)> = Vec::new();
            for (label, w) in &settings {
                if *label == best.1 || label.starts_with("zero") {
                    let mut p = ConeProblem::new(&d.mesh, &d.geometry, 0.5);
                    p.distortion_weights = w.clone();
                    if label.starts_with("feature size") {
                        let eps: f64 = label.rsplit("= ").next().unwrap().parse().unwrap();
                        let conv = if label.contains("printed") {
                            FeatureSizeConvention::AsPrinted
                        } else {
                            FeatureSizeConvention::AsDescribed
                        };
                        p.cone_costs = Some(
                            feature_size_weights(&d.mesh, &d.geometry, eps, conv)
                                .unwrap()
                                .cone_costs,
                        );
                    }
                    cands.push((label.clone(), p));
                }
            }
            for (label, base) in cands {
                for b in [BoundaryTreatment::Fixed, BoundaryTreatment::Free] {
                    println!("   penalized solves, {label}, boundary {b:?}");
                    println!(
                        "   {:>8} {:>6} {:>10} {:>9} {:>11}",
                        "lambda", "cones", "published", "on tips", "elsewhere"
                    );
                    for lambda in [
                        2.0, 1.0, 0.5, 0.2, 0.1, 0.05, 0.02, 0.01, 0.005, 0.002, 0.001,
                    ] {
                        let p = ConeProblem {
                            lambda,
                            ..base.clone()
                        };
                        let s = run(&p, b);
                        let on: Vec<u32> = s.cones.iter().map(|c| c.vertex).collect();
                        let hit = on.iter().filter(|v| published.contains(v)).count();
                        let tip = on.iter().filter(|&&v| ring.contains(&(v as usize))).count();
                        println!(
                            "   {:>8} {:>6} {:>10} {:>9} {:>11} {}",
                            lambda,
                            on.len(),
                            format!("{hit}/{}", published.len()),
                            tip,
                            on.len() - hit - tip,
                            if s.report.status == Status::Converged {
                                ""
                            } else {
                                "NOT CONVERGED"
                            }
                        );
                    }
                }
            }
        }
    }
}

fn main() {
    let mut args = std::env::args().skip(1);
    match args.next().as_deref() {
        Some("table") => table(args.next().map_or(0.5, |a| a.parse().unwrap())),
        Some("spikes") => spikes(),
        _ => eprintln!("usage: boundary table [lambda] | spikes"),
    }
}
