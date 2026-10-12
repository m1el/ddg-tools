//! Loads and validates every mesh of the optimal-cones dataset and prints its topology.
//!
//! The dataset is not part of the repository. When `data/` is absent the test prints a notice
//! and passes without checking anything.
//!
//! Run with `cargo test -p ddg-core --test dataset -- --nocapture` to see the table.

use std::path::PathBuf;

use ddg_core::geometry::Geometry;
use ddg_core::io::{read_cones, read_obj};
use ddg_core::operators::{cotan_laplacian, gauss_bonnet_residual, lumped_mass_normalized};
use ddg_core::sparse::{
    Method, Ordering, Parallelism, RefineOptions, SparseMatrix, Symbolic, ThreadPool,
    fill_reducing_ordering, relative_residual,
};
use ddg_core::topology::TriMesh;

fn dataset_dir() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../data/optimal-cones/optimal-cones-data/data")
}

#[test]
fn all_dataset_meshes_load_and_validate() {
    let dir = dataset_dir();
    if !dir.is_dir() {
        eprintln!(
            "SKIPPED: dataset directory {} not found; nothing was checked",
            dir.display()
        );
        return;
    }
    let mut names: Vec<String> = std::fs::read_dir(&dir)
        .unwrap()
        .filter_map(|e| e.ok())
        .filter_map(|e| {
            e.file_name()
                .to_str()
                .and_then(|n| n.strip_suffix(".obj"))
                .map(str::to_string)
        })
        .collect();
    names.sort();
    assert_eq!(
        names.len(),
        23,
        "expected the 23 meshes of the dataset, found {names:?}"
    );

    println!(
        "{:<20} {:>8} {:>8} {:>8} {:>6} {:>6} {:>5} {:>5} {:>6} {:>10}",
        "mesh", "V", "E", "F", "loops", "bdryV", "comp", "genus", "cones", "GB resid"
    );
    for name in &names {
        let obj = read_obj(dir.join(format!("{name}.obj"))).unwrap_or_else(|e| panic!("{e}"));
        let faces = obj
            .triangles()
            .unwrap_or_else(|| panic!("{name}: not a triangle mesh"));
        let mesh =
            TriMesh::new(obj.positions.len(), faces).unwrap_or_else(|e| panic!("{name}: {e}"));
        let geometry = Geometry::from_positions(&mesh, &obj.positions)
            .unwrap_or_else(|e| panic!("{name}: {e}"));
        let residual = gauss_bonnet_residual(&mesh, &geometry).unwrap();
        let cones = read_cones(dir.join(format!("{name}-cones.txt")), mesh.n_vertices())
            .unwrap_or_else(|e| panic!("{e}"));
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
        assert!(
            residual.abs() < 1e-8,
            "{name}: Gauss-Bonnet residual {residual}"
        );
        assert!(
            obj.corner_uvs().is_some(),
            "{name}: every corner of the dataset meshes has a UV"
        );
        assert!(!cones.is_empty(), "{name}: no cones");
        // Cones of the dataset sit at interior vertices.
        for cone in &cones {
            assert!((cone.vertex as usize) < mesh.n_vertices());
        }
    }
}

/// The linear systems the later crates will solve, on real meshes: the pinned Laplacian
/// (flattening) and the saddle-point matrix of the cone solver's Newton step.
#[test]
fn sparse_solvers_work_on_dataset_meshes() {
    let dir = dataset_dir();
    if !dir.is_dir() {
        eprintln!(
            "SKIPPED: dataset directory {} not found; nothing was checked",
            dir.display()
        );
        return;
    }
    let pool = ThreadPool::new(2).unwrap();
    // One closed mesh and one with boundary.
    for name in ["bunny", "face"] {
        let obj = read_obj(dir.join(format!("{name}.obj"))).unwrap();
        let mesh = TriMesh::new(obj.positions.len(), obj.triangles().unwrap()).unwrap();
        let geometry = Geometry::from_positions(&mesh, &obj.positions).unwrap();
        let laplacian = cotan_laplacian(&mesh, &geometry).unwrap();
        let mass = lumped_mass_normalized(&mesh, &geometry).unwrap();
        let n = mesh.n_vertices();

        // Interior block L_II (boundary pinned), or all but one vertex on a closed mesh.
        let free: Vec<u32> = if mesh.boundary_loops().is_empty() {
            (1..n as u32).collect()
        } else {
            mesh.interior_vertices()
        };
        let a = laplacian.submatrix(&free, &free).unwrap();
        let x_true: Vec<f64> = (0..free.len()).map(|i| (0.01 * i as f64).sin()).collect();
        let mut b = vec![0.0; free.len()];
        a.mul_vec(&x_true, &mut b).unwrap();
        let symbolic = Symbolic::analyze(&a, Method::Llt).unwrap();
        for par in [Parallelism::Sequential, Parallelism::Pool(&pool)] {
            let factorization = symbolic.factor(&a, par).unwrap();
            let mut x = b.clone();
            factorization.solve_in_place(&mut x, par).unwrap();
            let residual = relative_residual(&a, &x, &b).unwrap();
            assert!(residual < 1e-12, "{name}: LLt residual {residual}");
        }
        println!("{name}: LLt of the pinned Laplacian: {:?}", symbolic.info());

        // K = [[-M, L], [L, D]] with a sparse non-negative D, by LDLt in the paired order.
        let mut triplets = Vec::new();
        for i in 0..n {
            triplets.push((i as u32, i as u32, -mass[i]));
            triplets.push((
                (n + i) as u32,
                (n + i) as u32,
                if i % 50 == 0 { 1.0 } else { 0.0 },
            ));
        }
        for (i, j, v) in laplacian.triplets() {
            triplets.push((n as u32 + i, j, v));
            triplets.push((i, n as u32 + j, v));
        }
        let k = SparseMatrix::from_triplets(2 * n, 2 * n, &triplets).unwrap();
        let paired: Vec<u32> = fill_reducing_ordering(&laplacian)
            .unwrap()
            .iter()
            .flat_map(|&i| [i, n as u32 + i])
            .collect();
        let symbolic =
            Symbolic::analyze_with_ordering(&k, Method::Ldlt, Ordering::Custom(&paired)).unwrap();
        let factorization = symbolic.factor(&k, Parallelism::Sequential).unwrap();
        let x_true: Vec<f64> = (0..2 * n).map(|i| (0.02 * i as f64).cos()).collect();
        let mut b = vec![0.0; 2 * n];
        k.mul_vec(&x_true, &mut b).unwrap();
        let mut x = vec![0.0; 2 * n];
        let report = factorization
            .solve_refined(
                &k,
                &b,
                &mut x,
                RefineOptions {
                    max_steps: 5,
                    tolerance: 1e-12,
                },
                Parallelism::Sequential,
            )
            .unwrap();
        println!(
            "{name}: LDLt of K in paired order: {:?}, {report:?}",
            symbolic.info()
        );
        assert!(report.converged, "{name}: {report:?}");
    }
}
