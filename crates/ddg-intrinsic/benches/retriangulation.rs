//! Benchmarks of the retriangulations of paper Sec. 4: Delaunay flipping (iDT), Delaunay
//! refinement to 25° (iDR) and the optimal Delaunay triangulation (iODT, default options), and
//! one sweep of vertex moves.
//!
//! ```text
//! cargo bench -p ddg-intrinsic --bench retriangulation                 # everything available
//! cargo bench -p ddg-intrinsic --bench retriangulation -- iodt         # one operation
//! ```
//!
//! The generated meshes need no files; spot and camelbad are read from `data/` of the
//! workspace (`DDG_MODELS`, `DDG_DATASET` override the directories) and skipped with a notice if
//! absent. Single thread; run nothing else at the same time when the numbers matter.

use std::path::PathBuf;
use std::time::Duration;

use criterion::{BatchSize, Criterion, criterion_group, criterion_main};
use ddg_core::generators::{SplitMix64, grid_disk, icosphere};
use ddg_core::io::read_obj;
use ddg_core::{Geometry, TriMesh};
use ddg_intrinsic::{
    DelaunayOptions, IntrinsicTriangulation, OptimalDelaunayOptions, RefineOptions,
};

struct Input {
    name: String,
    mesh: TriMesh,
    geometry: Geometry,
}

fn input(name: &str, positions: Vec<[f64; 3]>, faces: Vec<[u32; 3]>) -> Input {
    let mesh = TriMesh::new(positions.len(), faces).expect("a valid mesh");
    let geometry = Geometry::from_positions(&mesh, &positions).expect("a valid geometry");
    Input {
        name: name.to_string(),
        mesh,
        geometry,
    }
}

fn from_file(name: &str, path: PathBuf) -> Option<Input> {
    if !path.is_file() {
        eprintln!("SKIPPED: {} not found", path.display());
        return None;
    }
    let obj = read_obj(&path).ok()?;
    let mut faces = Vec::new();
    for f in 0..obj.n_faces() {
        let c = &obj.corner_verts[obj.face_corners(f)];
        for k in 1..c.len().saturating_sub(1) {
            faces.push([c[0], c[k], c[k + 1]]);
        }
    }
    Some(input(name, obj.positions, faces))
}

fn inputs() -> Vec<Input> {
    let root = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../data");
    let models = std::env::var("DDG_MODELS")
        .map(PathBuf::from)
        .unwrap_or_else(|_| root.join("models/spot/spot"));
    let dataset = std::env::var("DDG_DATASET")
        .map(PathBuf::from)
        .unwrap_or_else(|_| root.join("optimal-cones/optimal-cones-data/data"));
    let mut out = Vec::new();
    let g = icosphere(4).expect("icosphere");
    let mut positions = g.positions;
    ddg_core::generators::perturb(&mut positions, [0.004; 3], 7);
    out.push(input("noisy icosphere 4", positions, g.faces));
    let g = grid_disk(48, 48).expect("grid");
    let mut positions = g.positions;
    let mut rng = SplitMix64::new(3);
    for p in &mut positions {
        p[0] += 1.0 * p[1] + 1e-3 * rng.uniform(-1.0, 1.0);
    }
    out.push(input("sheared grid 48", positions, g.faces));
    out.extend(from_file("spot", models.join("spot_triangulated.obj")));
    out.extend(from_file("camelbad", dataset.join("camelbad.obj")));
    out
}

fn bench(c: &mut Criterion) {
    for inp in inputs() {
        let fresh = IntrinsicTriangulation::new(&inp.mesh, &inp.geometry).expect("valid");
        let mut group = c.benchmark_group(format!("retriangulation/{}", inp.name));
        group
            .sample_size(10)
            .measurement_time(Duration::from_secs(10));
        group.bench_function("idt", |b| {
            b.iter_batched(
                || fresh.clone(),
                |mut t| t.flip_to_delaunay(&DelaunayOptions::default()).unwrap(),
                BatchSize::LargeInput,
            )
        });
        group.bench_function("idr25", |b| {
            b.iter_batched(
                || fresh.clone(),
                |mut t| t.refine(&RefineOptions::default()).unwrap(),
                BatchSize::LargeInput,
            )
        });
        group.bench_function("iodt", |b| {
            b.iter_batched(
                || fresh.clone(),
                |mut t| t.optimal_delaunay(&OptimalDelaunayOptions::default()).unwrap(),
                BatchSize::LargeInput,
            )
        });
        // One sweep of moves of every inserted vertex towards its target, after iODT.
        let mut converged = fresh.clone();
        converged
            .optimal_delaunay(&OptimalDelaunayOptions::default())
            .unwrap();
        group.bench_function("move sweep", |b| {
            b.iter_batched(
                || converged.clone(),
                |mut t| {
                    let mut moved = 0usize;
                    for v in t.n_input_vertices()..t.n_vertices() {
                        if t.is_boundary_vertex(v) {
                            continue;
                        }
                        if let Ok(goal) = t.odt_target(v, true) {
                            moved += usize::from(t.move_vertex(v, goal).is_ok());
                        }
                    }
                    moved
                },
                BatchSize::LargeInput,
            )
        });
        group.finish();
    }
}

criterion_group!(benches, bench);
criterion_main!(benches);
