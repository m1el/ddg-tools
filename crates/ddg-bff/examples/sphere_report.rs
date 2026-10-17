//! The table of the README section "Mapping to the sphere": `map_to_sphere` on generated
//! spheres and on the closed meshes of `data/` (one thread, best of three).
//!
//! ```sh
//! cargo run --release -p ddg-bff --example sphere_report               # all meshes
//! cargo run --release -p ddg-bff --example sphere_report -- a.obj b.obj  # these files
//! ```
//!
//! Meshes that are not closed of genus 0 are listed with the error.

use std::path::{Path, PathBuf};

use ddg_bff::{SphereMap, SphereOptions, map_to_sphere};
use ddg_core::generators::{icosphere, perturb};
use ddg_core::io::read_obj;
use ddg_core::{Geometry, TriMesh};

fn root() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../..")
}

fn load(path: &Path) -> Result<(TriMesh, Vec<[f64; 3]>), String> {
    let obj = read_obj(path).map_err(|e| e.to_string())?;
    let faces = obj.triangles().ok_or("not a triangle mesh")?;
    let mesh = TriMesh::new(obj.positions.len(), faces).map_err(|e| e.to_string())?;
    Ok((mesh, obj.positions))
}

fn best_of_three(mesh: &TriMesh, geometry: &Geometry) -> Result<SphereMap, String> {
    let mut best: Option<SphereMap> = None;
    for _ in 0..3 {
        let map =
            map_to_sphere(mesh, geometry, &SphereOptions::default()).map_err(|e| e.to_string())?;
        if best
            .as_ref()
            .is_none_or(|b| map.report.timings.total < b.report.timings.total)
        {
            best = Some(map);
        }
    }
    Ok(best.expect("three runs"))
}

fn row(name: &str, mesh: &TriMesh, positions: &[[f64; 3]]) {
    let geometry = match Geometry::from_positions(mesh, positions) {
        Ok(g) => g,
        Err(e) => {
            println!("| {name} | {} | – | {e} |", mesh.n_vertices());
            return;
        }
    };
    match best_of_three(mesh, &geometry) {
        Err(e) => println!("| {name} | {} | – | {e} |", mesh.n_vertices()),
        Ok(map) => {
            let r = &map.report;
            let d = &r.distortion.faces;
            let c = r.centering.expect("centering is on by default");
            let t = &r.timings;
            println!(
                "| {name} | {} | {} | {:.4} / {:.3} | {:.3} | {:.2} / {:.2} | {} | {:.9} | {} ({:.0e}) | {:.0} | {:.0} / {:.0} / {:.0} / {:.0} / {:.0} |",
                mesh.n_vertices(),
                r.pole,
                d.quasi_conformal.mean,
                d.quasi_conformal.max,
                d.log_scale_deviation,
                d.log_scale.min,
                d.log_scale.max,
                d.flipped_triangles + d.degenerate_triangles,
                r.distortion.covering_degree,
                c.iterations,
                c.final_center,
                1e3 * t.total,
                1e3 * t.puncture,
                1e3 * t.disk,
                1e3 * t.projection,
                1e3 * t.centering,
                1e3 * t.report,
            );
        }
    }
}

fn main() {
    println!(
        "| mesh | vertices | pole | QC mean / max | area (rms u) | u min / max | flipped | degree | centering its (|µ|) | ms | puncture / disk / projection / centering / report |"
    );
    println!("|---|---:|---:|---|---:|---|---:|---:|---|---:|---|");
    let args: Vec<String> = std::env::args().skip(1).collect();
    if !args.is_empty() {
        for a in &args {
            let path = PathBuf::from(a);
            match load(&path) {
                Ok((mesh, positions)) => row(a, &mesh, &positions),
                Err(e) => println!("| {a} | – | – | {e} |"),
            }
        }
        return;
    }
    for level in [3, 5] {
        let g = icosphere(level).unwrap();
        row(
            &format!("icosphere {level}"),
            &g.topology().unwrap(),
            &g.positions,
        );
    }
    {
        let mut g = icosphere(4).unwrap();
        // A squashed, noisy sphere: far from round, still genus 0.
        for p in &mut g.positions {
            p[0] *= 2.0;
            p[2] *= 0.5;
        }
        perturb(&mut g.positions, [0.01; 3], 7);
        row(
            "ellipsoid 2:1:0.5, noisy",
            &g.topology().unwrap(),
            &g.positions,
        );
    }
    let data = root().join("data");
    let files = [
        "models/spot/spot/spot_triangulated.obj",
        "models/blub/blub/blub_triangulated.obj",
        "models/bob/bob/bob_tri.obj",
        "optimal-cones/optimal-cones-data/data/bunny.obj",
        "optimal-cones/optimal-cones-data/data/axes.obj",
        "optimal-cones/optimal-cones-data/data/chair.obj",
        "optimal-cones/optimal-cones-data/data/brain.obj",
    ];
    for f in files {
        let path = data.join(f);
        if !path.exists() {
            continue;
        }
        let name = path.file_stem().unwrap().to_string_lossy().into_owned();
        match load(&path) {
            Ok((mesh, positions)) => row(&name, &mesh, &positions),
            Err(e) => println!("| {name} | – | – | {e} |"),
        }
    }
}
