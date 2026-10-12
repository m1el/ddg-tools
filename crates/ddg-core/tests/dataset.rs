//! Loads and validates every mesh of the optimal-cones dataset and prints its topology.
//!
//! The dataset is not part of the repository. When `data/` is absent the test prints a notice
//! and passes without checking anything.
//!
//! Run with `cargo test -p ddg-core --test dataset -- --nocapture` to see the table.

use std::path::PathBuf;

use ddg_core::geometry::Geometry;
use ddg_core::io::{read_cones, read_obj};
use ddg_core::operators::gauss_bonnet_residual;
use ddg_core::topology::TriMesh;

fn dataset_dir() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../data/optimal-cones/optimal-cones-data/data")
}

#[test]
fn all_dataset_meshes_load_and_validate() {
    let dir = dataset_dir();
    if !dir.is_dir() {
        eprintln!("SKIPPED: dataset directory {} not found; nothing was checked", dir.display());
        return;
    }
    let mut names: Vec<String> = std::fs::read_dir(&dir)
        .unwrap()
        .filter_map(|e| e.ok())
        .filter_map(|e| e.file_name().to_str().and_then(|n| n.strip_suffix(".obj")).map(str::to_string))
        .collect();
    names.sort();
    assert_eq!(names.len(), 23, "expected the 23 meshes of the dataset, found {names:?}");

    println!(
        "{:<20} {:>8} {:>8} {:>8} {:>6} {:>6} {:>5} {:>5} {:>6} {:>10}",
        "mesh", "V", "E", "F", "loops", "bdryV", "comp", "genus", "cones", "GB resid"
    );
    for name in &names {
        let obj = read_obj(dir.join(format!("{name}.obj"))).unwrap_or_else(|e| panic!("{e}"));
        let faces = obj.triangles().unwrap_or_else(|| panic!("{name}: not a triangle mesh"));
        let mesh = TriMesh::new(obj.positions.len(), faces).unwrap_or_else(|e| panic!("{name}: {e}"));
        let geometry = Geometry::from_positions(&mesh, &obj.positions).unwrap_or_else(|e| panic!("{name}: {e}"));
        let residual = gauss_bonnet_residual(&mesh, &geometry).unwrap();
        let cones =
            read_cones(dir.join(format!("{name}-cones.txt")), mesh.n_vertices()).unwrap_or_else(|e| panic!("{e}"));
        println!(
            "{:<20} {:>8} {:>8} {:>8} {:>6} {:>6} {:>5} {:>5} {:>6} {:>10.1e}",
            name,
            mesh.n_vertices(),
            mesh.n_edges(),
            mesh.n_faces(),
            mesh.boundary_loops().len(),
            mesh.n_boundary_edges(),
            mesh.n_components(),
            mesh.genus(),
            cones.len(),
            residual
        );
        assert!(residual.abs() < 1e-8, "{name}: Gauss-Bonnet residual {residual}");
        assert!(obj.corner_uvs().is_some(), "{name}: every corner of the dataset meshes has a UV");
        assert!(!cones.is_empty(), "{name}: no cones");
        // Cones of the dataset sit at interior vertices.
        for cone in &cones {
            assert!((cone.vertex as usize) < mesh.n_vertices());
        }
    }
}
