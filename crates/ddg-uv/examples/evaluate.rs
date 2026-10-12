//! The measurements behind the defaults of `ddg-uv`, as tables.
//!
//! ```text
//! cargo run --release -p ddg-uv --example evaluate -- closed    # closed surfaces without seams
//! cargo run --release -p ddg-uv --example evaluate -- lambda    # cones and distortion against λ
//! cargo run --release -p ddg-uv --example evaluate -- signs     # negative cones
//! cargo run --release -p ddg-uv --example evaluate -- merge     # cluster merging and its guard
//! cargo run --release -p ddg-uv --example evaluate -- dataset   # the 23 dataset meshes
//! ```
//!
//! A further argument restricts the meshes to those whose name contains it. The dataset is
//! read from `data/optimal-cones/optimal-cones-data/data` (or `DDG_DATASET`); without it the
//! dataset rows are left out.

use std::f64::consts::{PI, TAU};
use std::path::PathBuf;

use ddg_core::generators::icosphere;
use ddg_core::io::read_obj;
use ddg_uv::primitives::{capsule, cone, cube, cylinder, grid_cube, holed_slab, torus, uv_sphere};
use ddg_uv::{
    ClosedSurfaces, ConeSigns, MergeClusters, MergeOutcome, PolygonMesh, SlitPlacement,
    UnwrapOptions, area_distortion, unwrap,
};

fn dataset_dir() -> PathBuf {
    std::env::var_os("DDG_DATASET").map_or_else(
        || {
            PathBuf::from(env!("CARGO_MANIFEST_DIR"))
                .join("../../data/optimal-cones/optimal-cones-data/data")
        },
        PathBuf::from,
    )
}

/// A dataset mesh and the area distortion of its published UV map (if it has one).
fn load(name: &str) -> Option<(PolygonMesh, Option<f64>)> {
    let obj = read_obj(dataset_dir().join(format!("{name}.obj"))).ok()?;
    let faces = obj.triangles()?;
    let mesh = PolygonMesh::from_triangles(obj.positions.clone(), &faces);
    let published = obj.corner_uvs().map(|uvs| {
        area_distortion(&mesh.positions, &mesh.corner_vertices, &mesh.triangles, &uvs).scale_free
    });
    Some((mesh, published))
}

fn dataset_names() -> Vec<String> {
    let Ok(entries) = std::fs::read_dir(dataset_dir()) else {
        return Vec::new();
    };
    let mut names: Vec<String> = entries
        .filter_map(|e| e.ok())
        .filter_map(|e| {
            let name = e.file_name().into_string().ok()?;
            name.strip_suffix(".obj").map(str::to_owned)
        })
        .collect();
    names.sort();
    names
}

fn primitives() -> Vec<(String, PolygonMesh)> {
    let ico = icosphere(2).unwrap();
    vec![
        ("uv_sphere 32x16".into(), uv_sphere(32, 16, 1.0).unwrap()),
        ("icosphere 320".into(), PolygonMesh::from_triangles(ico.positions.clone(), &ico.faces)),
        ("cube".into(), cube(2.0)),
        ("grid_cube 8".into(), grid_cube(8, 2.0)),
        ("cylinder 32".into(), cylinder(32, 1, 1.0, 2.0).unwrap()),
        ("cone 32".into(), cone(32, 1.0, 2.0).unwrap()),
        ("capsule 32x8".into(), capsule(32, 8, 0.5, 1.0).unwrap()),
        ("torus 48x12".into(), torus(48, 12, 1.0, 0.25).unwrap()),
    ]
}

/// One row: what a run with these options gives.
fn row(mesh: &PolygonMesh, options: &UnwrapOptions) -> String {
    match unwrap(&mesh.input(), options) {
        Err(e) => {
            let mut message = e.to_string();
            message.truncate(110);
            format!("FAILED: {message}")
        }
        Ok(map) => {
            let r = &map.report;
            let solve = r.charts[0].solve.as_ref();
            format!(
                "{:>5} | {:>3} | {:>7.3} | {:>6.3} | {:>8.4} | {:>8.4} | {:>6.2} | {:>4} | {:>5} | {:>5.3} | {:>4} | {:>7.3} | {}",
                r.cones.n_cones,
                r.cones.n_negative,
                r.cones.total_angle,
                r.cones.max_angle,
                r.area_distortion.scale_free,
                solve.map_or(f64::NAN, |s| s.distortion),
                100.0 * r.overlapped_fraction.unwrap_or(f64::NAN),
                r.flipped_triangles,
                r.n_cut_edges,
                r.packing.map_or(f64::NAN, |p| p.fill),
                solve.map_or(0, |s| s.n_cones_optimal),
                r.timings.total,
                match solve.map(|s| s.merge) {
                    Some(MergeOutcome::Merged) => "merged",
                    Some(MergeOutcome::Rejected) => "merge rejected",
                    Some(MergeOutcome::NothingToMerge) => "-",
                    Some(MergeOutcome::NotRequested) => "off",
                    None => "no solve",
                },
            )
        }
    }
}

const HEADER: &str = "cones | neg |  Σ|µ|   | max µ  | distort. | solver D | ovl %  | flip | cuts  | fill  | raw  | seconds | merge";

fn closed(filter: &str) {
    let slit = |placement, relative_length| ClosedSurfaces::Slit {
        placement,
        relative_length,
    };
    let sites = |count, max_angle, forbid_elsewhere| ClosedSurfaces::FreeSites {
        count,
        max_angle,
        forbid_elsewhere,
    };
    let big = 0.95 * TAU;
    let policies: Vec<(&str, ClosedSurfaces, f64)> = vec![
        ("(a) raw optimum, merged      λ 0.5", ClosedSurfaces::RawOptimum, 0.5),
        ("(b) 1 free site, unbounded   λ 0.5", sites(1, f64::INFINITY, false), 0.5),
        ("(b) 1 free site, unbounded   λ 0.05", sites(1, f64::INFINITY, false), 0.05),
        ("(b) 1 free site, µ<=0.95·2π  λ 0.5", sites(1, big, false), 0.5),
        ("(b) 4 sites only, µ<=0.95·2π λ 0.5", sites(4, big, true), 0.5),
        ("(b) 8 sites only, µ<=π       λ 0.5", sites(8, PI, true), 0.5),
        ("(b) 8 sites only, µ<=0.95·2π λ 0.5", sites(8, big, true), 0.5),
        ("(b) 16 sites only, µ<=π      λ 0.5", sites(16, PI, true), 0.5),
        ("(c) slit at max |Ω|, 2 edges λ 0.5", slit(SlitPlacement::AtLargestCurvature, 0.0), 0.5),
        ("(c) slit at max |Ω|, 2 edges λ 0.1", slit(SlitPlacement::AtLargestCurvature, 0.0), 0.1),
        ("(c) slit at max |Ω|, 2 edges λ 0.02", slit(SlitPlacement::AtLargestCurvature, 0.0), 0.02),
        ("(c) slit mid-diameter, 0.25  λ 0.5", slit(SlitPlacement::DiameterMiddle, 0.25), 0.5),
        ("(c) slit mid-diameter, 0.5   λ 0.5", slit(SlitPlacement::DiameterMiddle, 0.5), 0.5),
        ("(c) slit mid-diameter, 0.5   λ 0.1", slit(SlitPlacement::DiameterMiddle, 0.5), 0.1),
        ("(c) slit mid-diameter, 0.5   λ 0.02", slit(SlitPlacement::DiameterMiddle, 0.5), 0.02),
        ("(c) slit mid-diameter, 1.0   λ 0.5", slit(SlitPlacement::DiameterMiddle, 1.0), 0.5),
        ("(c) slit mid-diameter, 1.0   λ 0.1", slit(SlitPlacement::DiameterMiddle, 1.0), 0.1),
        ("(c) slit mid-diameter, 1.0   λ 0.02", slit(SlitPlacement::DiameterMiddle, 1.0), 0.02),
    ];
    let mut meshes = primitives();
    for name in ["axes", "bunny", "chair", "brain"] {
        if !filter.is_empty()
            && name.contains(filter)
            && let Some((mesh, _)) = load(name)
        {
            meshes.push((name.to_string(), mesh));
        } else if filter.is_empty()
            && name != "brain"
            && let Some((mesh, _)) = load(name)
        {
            meshes.push((name.to_string(), mesh));
        }
    }
    for (name, mesh) in meshes.iter().filter(|(n, _)| n.contains(filter)) {
        println!("\n## {name} ({} polygons)", mesh.n_polygons());
        println!("{:<36} | {HEADER}", "policy");
        for (label, policy, lambda) in &policies {
            let options = UnwrapOptions {
                closed_surfaces: *policy,
                lambda: *lambda,
                merge: MergeClusters::Guarded {
                    max_cluster_size: usize::MAX,
                    max_factor: 2.0,
                    max_increase: 0.01,
                },
                ..UnwrapOptions::default()
            };
            println!("{label:<36} | {}", row(mesh, &options));
        }
    }
}

fn lambda(filter: &str) {
    let lambdas = [2.0, 1.0, 0.5, 0.2, 0.1, 0.05, 0.02, 0.01, 0.005, 0.002];
    let mut meshes = primitives();
    meshes.push(("holed_slab genus 2".into(), holed_slab(2, 4, 1.0)));
    for name in dataset_names() {
        if name != "brain"
            && name != "cow00001"
            && let Some((mesh, _)) = load(&name)
        {
            meshes.push((name, mesh));
        }
    }
    for (name, mesh) in meshes.iter().filter(|(n, _)| n.contains(filter)) {
        println!("\n## {name} ({} polygons)", mesh.n_polygons());
        println!("{:<8} | {HEADER}", "λ");
        for &lambda in &lambdas {
            let options = UnwrapOptions {
                lambda,
                ..UnwrapOptions::default()
            };
            println!("{lambda:<8} | {}", row(mesh, &options));
        }
    }
}

fn signs(filter: &str) {
    let mut meshes = vec![
        ("torus 48x12".to_string(), torus(48, 12, 1.0, 0.25).unwrap()),
        ("torus 48x12 thick".to_string(), torus(48, 12, 1.0, 0.6).unwrap()),
        ("holed_slab genus 2".to_string(), holed_slab(2, 4, 1.0)),
        ("holed_slab genus 3".to_string(), holed_slab(3, 4, 1.0)),
    ];
    for name in dataset_names() {
        if name != "brain"
            && name != "cow00001"
            && let Some((mesh, _)) = load(&name)
        {
            meshes.push((name, mesh));
        }
    }
    for (name, mesh) in meshes.iter().filter(|(n, _)| n.contains(filter)) {
        println!("\n## {name} ({} polygons)", mesh.n_polygons());
        println!("{:<22} | {HEADER}", "signs, λ");
        for lambda in [0.5, 0.1, 0.02] {
            for (label, cone_signs) in [
                ("positive", ConeSigns::PositiveWhereFeasible),
                ("any", ConeSigns::Any),
            ] {
                let options = UnwrapOptions {
                    lambda,
                    cone_signs,
                    ..UnwrapOptions::default()
                };
                println!("{:<22} | {}", format!("{label}, {lambda}"), row(mesh, &options));
            }
        }
    }
}

fn merge(filter: &str) {
    let mut meshes = Vec::new();
    for name in dataset_names() {
        if name != "brain"
            && name != "cow00001"
            && let Some((mesh, _)) = load(&name)
        {
            meshes.push((name, mesh));
        }
    }
    println!(
        "mesh | λ | signs | raw cones | clusters | merged cones | solver D before | after | ratio | guard"
    );
    for (name, mesh) in meshes.iter().filter(|(n, _)| n.contains(filter)) {
        for lambda in [0.5, 0.02] {
            for (label, cone_signs) in [
                ("positive", ConeSigns::PositiveWhereFeasible),
                ("any", ConeSigns::Any),
            ] {
                let options = UnwrapOptions {
                    lambda,
                    cone_signs,
                    closed_surfaces: ClosedSurfaces::RawOptimum,
                    merge: MergeClusters::Guarded {
                        max_cluster_size: usize::MAX,
                        max_factor: 2.0,
                        max_increase: 0.01,
                    },
                    max_cone_angle: TAU,
                    ..UnwrapOptions::default()
                };
                let line = match ddg_uv::prepare(&mesh.input(), &options)
                    .and_then(|prepared| ddg_uv::solve_cones(&prepared, &options))
                {
                    Err(e) => format!("FAILED: {e}"),
                    Ok(stage) => match stage.pieces[0].solve.as_ref() {
                        None => "no solve".to_string(),
                        Some(s) => format!(
                            "{} | {} | {} | {:.5} | {:.5} | {:.3} | {:?}",
                            s.n_cones_optimal,
                            s.n_clusters,
                            stage.pieces[0].cones.len(),
                            s.distortion_optimal,
                            s.distortion_merged,
                            s.distortion_merged / s.distortion_optimal,
                            s.merge
                        ),
                    },
                };
                println!("{name} | {lambda} | {label} | {line}");
            }
        }
    }
}

fn dataset(filter: &str) {
    let names = dataset_names();
    if names.is_empty() {
        println!("dataset not found at {}", dataset_dir().display());
        return;
    }
    println!(
        "mesh | triangles | cones | Σ|µ| | distortion | published | overlap % | flipped | charts | \
         validate | seams | cones s | flatten | pack | verify | total | warnings"
    );
    for name in names.iter().filter(|n| n.contains(filter)) {
        let Some((mesh, published)) = load(name) else {
            println!("{name} | not a triangle mesh");
            continue;
        };
        let options = UnwrapOptions::default();
        match unwrap(&mesh.input(), &options) {
            Err(e) => println!("{name} | {} | FAILED: {e}", mesh.n_polygons()),
            Ok(map) => {
                let r = &map.report;
                let t = &r.timings;
                println!(
                    "{name} | {} | {} | {:.3} | {:.5} | {} | {:.3} | {} | {} | {:.3} | {:.3} | \
                     {:.3} | {:.3} | {:.3} | {:.3} | {:.3} | {}",
                    mesh.n_polygons(),
                    r.cones.n_cones,
                    r.cones.total_angle,
                    r.area_distortion.scale_free,
                    published.map_or("-".to_string(), |d| format!("{d:.5}")),
                    100.0 * r.overlapped_fraction.unwrap_or(f64::NAN),
                    r.flipped_triangles,
                    r.charts.len(),
                    t.validate,
                    t.seams,
                    t.cones,
                    t.flatten,
                    t.pack,
                    t.verify,
                    t.total,
                    r.warnings.len(),
                );
                for warning in &r.warnings {
                    println!("    warning: {warning}");
                }
            }
        }
    }
}

fn main() {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let filter = args.get(1).map_or("", String::as_str);
    match args.first().map(String::as_str) {
        Some("closed") => closed(filter),
        Some("lambda") => lambda(filter),
        Some("signs") => signs(filter),
        Some("merge") => merge(filter),
        Some("dataset") => dataset(filter),
        _ => eprintln!("usage: evaluate closed|lambda|signs|merge|dataset [name filter]"),
    }
}
