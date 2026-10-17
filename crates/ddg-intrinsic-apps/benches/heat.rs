//! Benchmarks of the heat method (paper Sec. 5.3): factoring the two systems, and one distance
//! query with the factorizations reused, on the input mesh, its intrinsic Delaunay
//! triangulation and its Delaunay refinement to 25°.
//!
//! ```text
//! cargo bench -p ddg-intrinsic-apps --bench heat                 # everything available
//! cargo bench -p ddg-intrinsic-apps --bench heat -- spot         # one mesh
//! ```
//!
//! An icosphere needs no file; spot and camelbad are read from `data/` of the workspace
//! (`DDG_MODELS`, `DDG_DATASET` override the directories) and skipped with a notice if absent.
//! Single thread; run nothing else at the same time when the numbers matter.

use std::path::PathBuf;
use std::time::Duration;

use criterion::{Criterion, criterion_group, criterion_main};
use ddg_core::generators::icosphere;
use ddg_core::io::read_obj;
use ddg_core::{Geometry, TriMesh};
use ddg_intrinsic::{DelaunayOptions, IntrinsicTriangulation, RefineOptions};
use ddg_intrinsic_apps::heat::{HeatOptions, HeatSolver, HeatSource};

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
        .unwrap_or(root.join("models"));
    let dataset = std::env::var("DDG_DATASET")
        .map(PathBuf::from)
        .unwrap_or(root.join("optimal-cones/optimal-cones-data/data"));
    let sphere = icosphere(4).expect("icosphere");
    let mut list = vec![input("icosphere4", sphere.positions, sphere.faces)];
    list.extend(from_file("spot", models.join("spot/spot/spot_triangulated.obj")));
    list.extend(from_file("camelbad", dataset.join("camelbad.obj")));
    list
}

fn bench(c: &mut Criterion) {
    for input in inputs() {
        let fresh = IntrinsicTriangulation::new(&input.mesh, &input.geometry).expect("valid");
        let mut idt = fresh.clone();
        idt.flip_to_delaunay(&DelaunayOptions::default())
            .expect("flips");
        let mut idr = idt.clone();
        idr.refine(&RefineOptions::default()).expect("refinement");
        for (mode, tri) in [("input", &fresh), ("idt", &idt), ("idr25", &idr)] {
            let options = HeatOptions::default();
            c.bench_function(&format!("{}/{mode}/factor", input.name), |b| {
                b.iter(|| HeatSolver::new(tri, &options).expect("factors"))
            });
            let solver = HeatSolver::new(tri, &options).expect("factors");
            c.bench_function(&format!("{}/{mode}/distance", input.name), |b| {
                b.iter(|| solver.distance(&[HeatSource::Vertex(0)]).expect("solves"))
            });
        }
    }
}

criterion_group! {
    name = benches;
    config = Criterion::default().sample_size(10).measurement_time(Duration::from_secs(3));
    targets = bench
}
criterion_main!(benches);
