//! The measurements behind the spike check (`ddg_uv::spikes`): how the distortion `m_i u_i²` of
//! the surface without cones is spread over its vertices.
//!
//! ```text
//! cargo run --release -p ddg-uv --example spikes [name filter | prim]
//! ```
//!
//! One row per dataset mesh (from `data/optimal-cones/optimal-cones-data/data`, or
//! `DDG_DATASET`) and per generated primitive: the cone-free distortion `E0 = ½ Σ m u²`, the
//! share of it carried by the 1, 8 and 32 largest vertices, the number of vertices needed for
//! 90 % and 99 %, the area share of the 8 largest, and the time of the cone-free evaluation
//! (one Laplacian factorization). It uses the solver's public `evaluate_cones` only.
use std::time::Instant;

use ddg_core::io::read_obj;
use ddg_core::operators::lumped_mass_normalized;
use ddg_optimal_cones::{ConeProblem, Parallelism, evaluate_cones};
use ddg_uv::primitives::{capsule, cone, cube, cylinder, grid_cube, torus, uv_sphere};
use ddg_uv::{PolygonMesh, UnwrapOptions, prepare};

fn probe(name: &str, mesh: &PolygonMesh) {
    let options = UnwrapOptions::default();
    let Ok(prepared) = prepare(&mesh.input(), &options) else {
        println!("{name}: prepare failed");
        return;
    };
    for (pi, piece) in prepared.pieces.iter().enumerate() {
        let n = piece.mesh.n_vertices();
        let t = Instant::now();
        let problem = ConeProblem::new(&piece.mesh, &piece.geometry, 0.1);
        let Ok(ev) = evaluate_cones(&problem, &[], Parallelism::Sequential) else {
            println!("{name}: eval failed");
            continue;
        };
        let secs = t.elapsed().as_secs_f64();
        let mass = lumped_mass_normalized(&piece.mesh, &piece.geometry).unwrap();
        let mut c: Vec<(f64, usize)> = (0..n)
            .map(|i| (mass[i] * ev.log_scale[i].powi(2), i))
            .collect();
        let total: f64 = c.iter().map(|x| x.0).sum();
        c.sort_by(|a, b| b.0.partial_cmp(&a.0).unwrap());
        let share = |k: usize| c.iter().take(k).map(|x| x.0).sum::<f64>() / total;
        let mshare = |k: usize| c.iter().take(k).map(|x| mass[x.1]).sum::<f64>();
        // smallest k reaching 90 %, 99 %
        let reach = |f: f64| {
            let mut acc = 0.0;
            for (k, x) in c.iter().enumerate() {
                acc += x.0;
                if acc >= f * total {
                    return k + 1;
                }
            }
            n
        };
        println!(
            "{name} p{pi} n={n} E0={:.4e} top1={:.4} top8={:.4} top32={:.4} k90={} k99={} mass(top8)={:.2e} t={:.3}s",
            ev.distortion,
            share(1),
            share(8),
            share(32),
            reach(0.9),
            reach(0.99),
            mshare(8),
            secs
        );
    }
}

fn main() {
    let filter = std::env::args().nth(1).unwrap_or_default();
    let dir = std::env::var("DDG_DATASET")
        .unwrap_or("../../data/optimal-cones/optimal-cones-data/data".into());
    let mut names: Vec<String> = std::fs::read_dir(&dir)
        .unwrap()
        .filter_map(|e| e.ok())
        .filter_map(|e| {
            e.file_name()
                .into_string()
                .ok()?
                .strip_suffix(".obj")
                .map(str::to_owned)
        })
        .collect();
    names.sort();
    for name in names.iter().filter(|n| n.contains(&filter)) {
        let obj = read_obj(format!("{dir}/{name}.obj")).unwrap();
        let Some(faces) = obj.triangles() else {
            continue;
        };
        let mesh = PolygonMesh::from_triangles(obj.positions.clone(), &faces);
        probe(name, &mesh);
    }
    if filter.is_empty() || filter == "prim" {
        let ico = ddg_core::generators::icosphere(2).unwrap();
        let list: Vec<(&str, PolygonMesh)> = vec![
            ("uv_sphere", uv_sphere(32, 16, 1.0).unwrap()),
            (
                "icosphere",
                PolygonMesh::from_triangles(ico.positions.clone(), &ico.faces),
            ),
            ("cube", cube(2.0)),
            ("grid_cube8", grid_cube(8, 2.0)),
            ("cylinder", cylinder(32, 1, 1.0, 2.0).unwrap()),
            ("cone", cone(32, 1.0, 2.0).unwrap()),
            ("capsule", capsule(32, 8, 0.5, 1.0).unwrap()),
            ("torus", torus(48, 12, 1.0, 0.25).unwrap()),
        ];
        for (n, m) in &list {
            probe(n, m);
        }
    }
}
