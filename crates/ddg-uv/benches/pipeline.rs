//! End-to-end benchmarks: from an in-memory polygon mesh to packed UVs, and the stages apart.
//!
//! ```text
//! cargo bench -p ddg-uv                                   # everything but the largest mesh
//! cargo bench -p ddg-uv -- unwrap/                        # end to end only
//! cargo bench -p ddg-uv -- stages/                        # per-stage breakdown
//! DDG_BENCH_LARGE=1 cargo bench -p ddg-uv -- brain        # the 306k-triangle mesh (minutes)
//! cargo bench -p ddg-uv -- --save-baseline before         # record ...
//! cargo bench -p ddg-uv -- --baseline before              # ... and compare after a change
//! ```
//!
//! The generated meshes need no files. The dataset meshes are read from
//! `data/optimal-cones/optimal-cones-data/data` (or `DDG_DATASET`) and skipped with a notice
//! when they are absent. All runs use `UnwrapOptions::default()` and one thread.

use std::path::PathBuf;
use std::time::Duration;

use criterion::{BatchSize, Criterion, criterion_group, criterion_main};
use ddg_core::generators::bumpy_grid_disk;
use ddg_core::io::read_obj;
use ddg_uv::primitives::{grid_cube, torus, uv_sphere};
use ddg_uv::{PolygonMesh, UnwrapOptions, assemble, flatten, prepare, solve_cones, unwrap};

fn dataset(name: &str) -> Option<PolygonMesh> {
    let dir = std::env::var_os("DDG_DATASET").map_or_else(
        || {
            PathBuf::from(env!("CARGO_MANIFEST_DIR"))
                .join("../../data/optimal-cones/optimal-cones-data/data")
        },
        PathBuf::from,
    );
    let Ok(obj) = read_obj(dir.join(format!("{name}.obj"))) else {
        eprintln!("dataset mesh {name} not found in {}: skipped", dir.display());
        return None;
    };
    Some(PolygonMesh::from_triangles(
        obj.positions.clone(),
        &obj.triangles()?,
    ))
}

/// The meshes: name, mesh, and whether it runs only when `DDG_BENCH_LARGE` is set.
fn meshes() -> Vec<(&'static str, PolygonMesh, bool)> {
    let disk = bumpy_grid_disk(64, 64, 0.3).expect("valid parameters");
    let mut meshes = vec![
        // 8 192 triangles with boundary: cones from the penalized problem.
        (
            "bumpy_disk_64",
            PolygonMesh::from_triangles(disk.positions.clone(), &disk.faces),
            false,
        ),
        // 2 048 polygons (quads and fans), closed and round: the free-site fallback, two solves.
        ("uv_sphere_64x32", uv_sphere(64, 32, 1.0).expect("valid"), false),
        // 1 536 quads, closed with corners: the raw optimum.
        ("grid_cube_16", grid_cube(16, 2.0), false),
        // 2 304 quads, genus 1: no cone solve, handle cuts.
        ("torus_96x24", torus(96, 24, 1.0, 0.25).expect("valid"), false),
    ];
    for (name, large) in [
        ("hand", false),   // 6 660 triangles, boundary
        ("camel", false),  // 22 704 triangles, boundary
        ("bunny", false),  // 28 576 triangles, closed
        ("brain", true),   // 305 782 triangles, closed
    ] {
        if let Some(mesh) = dataset(name) {
            meshes.push((name, mesh, large));
        }
    }
    meshes
}

fn large_enabled() -> bool {
    std::env::var_os("DDG_BENCH_LARGE").is_some()
}

fn end_to_end(c: &mut Criterion) {
    let options = UnwrapOptions::default();
    let mut group = c.benchmark_group("unwrap");
    group.sample_size(10).warm_up_time(Duration::from_secs(1));
    for (name, mesh, large) in meshes() {
        if large && !large_enabled() {
            eprintln!("unwrap/{name}: skipped (set DDG_BENCH_LARGE=1)");
            continue;
        }
        let input = mesh.input();
        group.bench_function(name, |b| {
            b.iter(|| unwrap(&input, &options).expect("the benchmark meshes unwrap"))
        });
    }
    group.finish();
}

fn stages(c: &mut Criterion) {
    let options = UnwrapOptions::default();
    let mut group = c.benchmark_group("stages");
    group.sample_size(10).warm_up_time(Duration::from_secs(1));
    for (name, mesh, large) in meshes() {
        if !matches!(name, "bumpy_disk_64" | "uv_sphere_64x32" | "camel" | "brain") {
            continue;
        }
        if large && !large_enabled() {
            continue;
        }
        let input = mesh.input();
        let prepared = prepare(&input, &options).expect("prepare");
        let cones = solve_cones(&prepared, &options).expect("cones");
        let flat = flatten(&prepared, &cones, &options).expect("flatten");
        group.bench_function(format!("{name}/1_prepare"), |b| {
            b.iter(|| prepare(&input, &options).expect("prepare"))
        });
        group.bench_function(format!("{name}/2_cones"), |b| {
            b.iter(|| solve_cones(&prepared, &options).expect("cones"))
        });
        group.bench_function(format!("{name}/3_flatten"), |b| {
            b.iter(|| flatten(&prepared, &cones, &options).expect("flatten"))
        });
        group.bench_function(format!("{name}/4_verify_pack"), |b| {
            b.iter_batched(
                || (),
                |()| assemble(&input, &prepared, &cones, &flat, &options).expect("assemble"),
                BatchSize::LargeInput,
            )
        });
    }
    group.finish();
}

criterion_group!(benches, end_to_end, stages);
criterion_main!(benches);
