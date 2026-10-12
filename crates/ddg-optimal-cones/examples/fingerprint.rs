//! Prints a hash of every bit of the results of default-mode solves (44 solves: synthetic meshes and
//! dataset meshes, with free sites, bounds, weights and merging), to compare two commits:
//! `cargo run --release -p ddg-optimal-cones --example fingerprint > a.log`, then `diff`.
//! Used to check that adding the boundary treatments left `solve` bit for bit unchanged.
use ddg_core::generators::{
    SplitMix64, annulus, bumpy_grid_disk, icosphere, open_cylinder, perturb, torus,
};
use ddg_core::io::{read_cones, read_obj};
use ddg_core::{Geometry, TriMesh};
use ddg_optimal_cones::*;

fn fnv(h: &mut u64, x: f64) {
    for b in x.to_bits().to_le_bytes() {
        *h ^= b as u64;
        *h = h.wrapping_mul(0x100000001b3);
    }
}

fn fingerprint(s: &Solution) -> u64 {
    let mut h = 0xcbf29ce484222325u64;
    for v in s
        .measure
        .iter()
        .chain(&s.log_scale_optimal)
        .chain(&s.adjoint)
        .chain(&s.log_scale)
    {
        fnv(&mut h, *v);
    }
    for v in [
        s.objective,
        s.distortion_optimal,
        s.cone_cost,
        s.distortion,
        s.total_cone_angle,
    ] {
        fnv(&mut h, v);
    }
    let c = &s.report.certificate;
    for v in [
        c.state_residual,
        c.adjoint_residual,
        c.cone_residual,
        c.bound_violation,
        c.dual_infeasibility,
        c.gauss_bonnet_residual,
        c.primal_objective,
        c.dual_objective,
        c.duality_gap,
        c.complementarity_gap,
        c.gap_uncertainty,
    ] {
        fnv(&mut h, v);
    }
    for cone in &s.cones {
        fnv(&mut h, cone.vertex as f64);
        fnv(&mut h, cone.angle);
    }
    fnv(&mut h, s.report.stages.len() as f64);
    fnv(&mut h, s.report.regularized_iterations as f64);
    fnv(&mut h, s.report.exact_iterations as f64);
    h
}

fn run(name: &str, mesh: &TriMesh, geometry: &Geometry) {
    let nv = mesh.n_vertices();
    for (label, lambda, merge, tweak) in [
        ("plain", 0.05, false, 0),
        ("merge", 0.01, true, 0),
        ("free-site", 0.05, false, 1),
        ("bounded", 0.02, false, 2),
    ] {
        let mut p = ConeProblem::new(mesh, geometry, lambda);
        match tweak {
            1 => {
                let mut c = vec![1.0; nv];
                c[nv / 2] = 0.0;
                for v in 0..nv / 5 {
                    c[v] = f64::INFINITY;
                }
                c[nv / 2] = 0.0;
                p.cone_costs = Some(c);
            }
            2 => {
                p.lower_bounds = Some(vec![0.0; nv]);
                p.upper_bounds = Some(vec![1.5; nv]);
                p.distortion_weights = Some((0..nv).map(|v| 0.5 + (v % 3) as f64).collect());
            }
            _ => {}
        }
        let options = SolveOptions {
            merge_cone_clusters: merge,
            ..SolveOptions::default()
        };
        match solve(&p, &options) {
            Ok(s) => println!(
                "{name:<12} {label:<10} {:?} {:016x} obj {:.15e}",
                s.report.status,
                fingerprint(&s),
                s.objective
            ),
            Err(e) => println!("{name:<12} {label:<10} error {e}"),
        }
    }
}

fn main() {
    let mut meshes: Vec<(String, Vec<[f64; 3]>, Vec<[u32; 3]>)> = Vec::new();
    let g = bumpy_grid_disk(14, 11, 0.3).unwrap();
    meshes.push(("disk".into(), g.positions.clone(), g.faces.clone()));
    let g = annulus(14, 4, 0.4, 1.0).unwrap();
    meshes.push(("annulus".into(), g.positions.clone(), g.faces.clone()));
    let g = open_cylinder(10, 6, 0.7, 1.5).unwrap();
    meshes.push(("cylinder".into(), g.positions.clone(), g.faces.clone()));
    let g = icosphere(2).unwrap();
    let mut pos = g.positions.clone();
    perturb(&mut pos, [0.05; 3], 7);
    meshes.push(("sphere".into(), pos, g.faces.clone()));
    let g = torus(9, 6, 1.0, 0.4).unwrap();
    let mut pos = g.positions.clone();
    perturb(&mut pos, [0.03; 3], 3);
    meshes.push(("torus".into(), pos, g.faces.clone()));
    let _ = SplitMix64::new(1);
    for (name, pos, faces) in &meshes {
        let mesh = TriMesh::new(pos.len(), faces.clone()).unwrap();
        let geometry = Geometry::from_positions(&mesh, pos).unwrap();
        run(name, &mesh, &geometry);
    }
    let dir = std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../data/optimal-cones/optimal-cones-data/data");
    for name in ["bump", "camel", "cow1", "hand", "axes", "airplane"] {
        let Ok(obj) = read_obj(dir.join(format!("{name}.obj"))) else {
            continue;
        };
        let faces = obj.triangles().unwrap();
        let mesh = TriMesh::new(obj.positions.len(), faces).unwrap();
        let geometry = Geometry::from_positions(&mesh, &obj.positions).unwrap();
        let _ = read_cones(dir.join(format!("{name}-cones.txt")), mesh.n_vertices());
        run(name, &mesh, &geometry);
    }
}
