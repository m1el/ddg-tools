//! Finite elements on intrinsic triangulations (paper Sec. 5.2).

mod common;

use common::{Shape, jittered_grid, sheared_parallelogram, sheared_torus, vertex_positions};
use ddg_core::SparseMatrix;
use ddg_core::generators::{icosphere, torus};
use ddg_intrinsic::{
    DelaunayOptions, IntrinsicTriangulation, OptimalDelaunayOptions, RefineOptions,
};
use ddg_intrinsic_apps::fem::{
    EigenOptions, FemError, MassMatrix, PinnedLaplacian, condition_number, galerkin_mass_matrix,
    lumped_mass, regular_subdivision, stiffness_matrix,
};

fn refined(shape: &Shape) -> Vec<(&'static str, IntrinsicTriangulation)> {
    let mut idt = shape.triangulation();
    idt.flip_to_delaunay(&DelaunayOptions::default()).unwrap();
    let mut idr = idt.clone();
    idr.refine(&RefineOptions::default()).unwrap();
    let mut iodt = idt.clone();
    iodt.optimal_delaunay(&OptimalDelaunayOptions::default())
        .unwrap();
    vec![
        ("input", shape.triangulation()),
        ("idt", idt),
        ("idr", idr),
        ("iodt", iodt),
    ]
}

#[test]
fn the_stiffness_matrix_is_exact_on_linear_functions_of_a_flat_triangulation() {
    for shape in [jittered_grid(10, 0.35, 2), sheared_parallelogram(8, 3.0)] {
        for (mode, tri) in refined(&shape) {
            let l = stiffness_matrix(&tri);
            let p = vertex_positions(&shape, &tri);
            for (a, b) in [(1.0, 0.0), (0.3, -1.7)] {
                let f: Vec<f64> = p.iter().map(|q| a * q[0] + b * q[1] + 0.25).collect();
                let mut lf = vec![0.0; f.len()];
                l.mul_vec(&f, &mut lf).unwrap();
                for v in 0..tri.n_vertices() {
                    if !tri.is_boundary_vertex(v) {
                        assert!(lf[v].abs() < 1e-11, "{} {mode}: vertex {v}: {}", shape.name, lf[v]);
                    }
                }
                // The Dirichlet energy of a linear function is |∇f|² times the area.
                let energy: f64 = f.iter().zip(&lf).map(|(x, y)| x * y).sum();
                let area = tri.total_area();
                assert!((energy - (a * a + b * b) * area).abs() < 1e-10 * area, "{mode}");
            }
        }
    }
}

#[test]
fn mass_matrices() {
    let shape = sheared_torus(20, 8, 1.0, 0.4, 3);
    for (mode, tri) in refined(&shape) {
        let m = galerkin_mass_matrix(&tri);
        assert!(m.symmetry_defect().unwrap() == 0.0, "{mode}");
        let lumped = lumped_mass(&tri);
        let mut rows = vec![0.0; tri.n_vertices()];
        m.mul_vec(&vec![1.0; tri.n_vertices()], &mut rows).unwrap();
        for v in 0..tri.n_vertices() {
            assert!((rows[v] - lumped[v]).abs() < 1e-14, "{mode}: vertex {v}");
        }
        // ∫ f² for f = 1 at one vertex: its diagonal entry; positive definiteness by Cholesky.
        assert!(
            ddg_core::sparse::Symbolic::analyze(&m, ddg_core::sparse::Method::Llt)
                .and_then(|s| s.factor(&m, ddg_core::sparse::Parallelism::Sequential))
                .is_ok()
        );
    }
}

#[test]
fn regular_subdivision_keeps_the_surface() {
    let shape = sheared_torus(12, 6, 1.0, 0.4, 2);
    let (mesh, geometry) = regular_subdivision(&shape.mesh, &shape.geometry, 2).unwrap();
    assert_eq!(mesh.n_faces(), 16 * shape.mesh.n_faces());
    let once = shape.mesh.n_vertices() + shape.mesh.n_edges();
    assert_eq!(mesh.n_vertices(), once + 2 * shape.mesh.n_edges() + 3 * shape.mesh.n_faces());
    assert!((geometry.total_area() - shape.geometry.total_area()).abs() < 1e-12);
    // Every small triangle is similar to its parent: the same angles.
    let (one, g1) = regular_subdivision(&shape.mesh, &shape.geometry, 1).unwrap();
    for f in 0..shape.mesh.n_faces() {
        let mut parent = shape.geometry.corner_angles()[f];
        parent.sort_by(f64::total_cmp);
        for c in 0..4 {
            let mut child = g1.corner_angles()[4 * f + c];
            child.sort_by(f64::total_cmp);
            for k in 0..3 {
                assert!((parent[k] - child[k]).abs() < 1e-12);
            }
        }
    }
    assert_eq!(one.n_vertices(), once);
    // The midpoints are where they should be in space: the lengths agree with the positions.
    let mut positions = shape.positions.clone();
    for e in 0..shape.mesh.n_edges() {
        let [a, b] = shape.mesh.edge_vertices(e);
        positions.push([0, 1, 2].map(|k| 0.5 * (shape.positions[a][k] + shape.positions[b][k])));
    }
    let from_positions = ddg_core::Geometry::from_positions(&one, &positions).unwrap();
    for e in 0..one.n_edges() {
        let (x, y) = (g1.edge_lengths()[e], from_positions.edge_lengths()[e]);
        // Chords of a curved surface: the subdivided surface is the input polyhedron, so the
        // half-edges are exact and the edges between midpoints are exact in each flat face.
        assert!((x - y).abs() < 1e-12, "edge {e}: {x} {y}");
    }
}

/// The eigenvalues of `L x = λ M x` for small dense matrices, by the Jacobi method on
/// `M^{-1/2} L M^{-1/2}` (diagonal `M`).
fn dense_eigenvalues(l: &SparseMatrix, m: &[f64]) -> Vec<f64> {
    let n = l.nrows();
    let mut a = vec![0.0; n * n];
    for (i, j, x) in l.triplets() {
        let (i, j) = (i as usize, j as usize);
        a[i * n + j] += x / (m[i] * m[j]).sqrt();
    }
    for _ in 0..100 {
        let mut off = 0.0;
        for p in 0..n {
            for r in p + 1..n {
                off += a[p * n + r] * a[p * n + r];
                if a[p * n + r].abs() < 1e-300 {
                    continue;
                }
                let theta = (a[r * n + r] - a[p * n + p]) / (2.0 * a[p * n + r]);
                let t = theta.signum() / (theta.abs() + (theta * theta + 1.0).sqrt());
                let c = 1.0 / (t * t + 1.0).sqrt();
                let s = t * c;
                for i in 0..n {
                    let (x, y) = (a[i * n + p], a[i * n + r]);
                    a[i * n + p] = c * x - s * y;
                    a[i * n + r] = s * x + c * y;
                }
                for i in 0..n {
                    let (x, y) = (a[p * n + i], a[r * n + i]);
                    a[p * n + i] = c * x - s * y;
                    a[r * n + i] = s * x + c * y;
                }
            }
        }
        if off < 1e-28 {
            break;
        }
    }
    let mut e: Vec<f64> = (0..n).map(|i| a[i * n + i]).collect();
    e.sort_by(f64::total_cmp);
    e
}

#[test]
fn condition_numbers_agree_with_dense_eigenvalues() {
    let shapes = [
        Shape::generated("icosphere 1", icosphere(1).unwrap()),
        jittered_grid(5, 0.3, 4),
        sheared_torus(9, 5, 1.0, 0.4, 2),
    ];
    for shape in &shapes {
        let tri = shape.triangulation();
        let l = stiffness_matrix(&tri);
        let n = tri.n_vertices();
        for (name, mass, diag) in [
            ("identity", MassMatrix::Identity, vec![1.0; n]),
            ("lumped", MassMatrix::Lumped(lumped_mass(&tri)), lumped_mass(&tri)),
        ] {
            let report = condition_number(&l, &mass, &EigenOptions::default()).unwrap();
            let e = dense_eigenvalues(&l, &diag);
            assert!(e[0].abs() < 1e-10 * e[n - 1]);
            let (max, min) = (e[n - 1], e[1]);
            assert!(report.lambda_max.converged && report.lambda_min.converged, "{report:?}");
            assert!((report.lambda_max.value - max).abs() < 1e-8 * max, "{} {name}: {report:?} {max}", shape.name);
            assert!((report.lambda_min.value - min).abs() < 1e-8 * min, "{} {name}: {report:?} {min}", shape.name);
            assert_eq!(report.n_components, 1);
        }
        // Galerkin mass: bounded by the lumped one (below).
        let galerkin = condition_number(
            &l,
            &MassMatrix::Galerkin(galerkin_mass_matrix(&tri)),
            &EigenOptions::default(),
        )
        .unwrap();
        let lumped =
            condition_number(&l, &MassMatrix::Lumped(lumped_mass(&tri)), &EigenOptions::default())
                .unwrap();
        // ¼ M_L ≤ M_G ≤ M_L (per triangle the Galerkin mass has the eigenvalues A/12, A/12,
        // A/3 and the lumped one A/3, A/3, A/3), so by the min-max principle
        // λ_k(L, M_L) ≤ λ_k(L, M_G) ≤ 4 λ_k(L, M_L) for every k.
        for (g, l) in [
            (galerkin.lambda_max.value, lumped.lambda_max.value),
            (galerkin.lambda_min.value, lumped.lambda_min.value),
        ] {
            assert!(g >= l * (1.0 - 1e-9) && g <= 4.0 * l, "{} {g} {l}", shape.name);
        }
    }
}

#[test]
fn intrinsic_triangulations_have_better_condition_numbers() {
    // Paper Fig. 17: a coarse mesh of bad triangles.
    let shape = sheared_torus(40, 10, 1.0, 0.4, 5);
    let mut kappa = Vec::new();
    for (mode, tri) in refined(&shape) {
        let l = stiffness_matrix(&tri);
        let report = condition_number(
            &l,
            &MassMatrix::Galerkin(galerkin_mass_matrix(&tri)),
            &EigenOptions::default(),
        )
        .unwrap();
        kappa.push((mode, tri.n_vertices(), report.condition_number));
    }
    println!("{kappa:?}");
    let (_, _, input) = kappa[0];
    let (_, _, idt) = kappa[1];
    assert!(input > 2.0 * idt, "{kappa:?}");
}

#[test]
fn the_pinned_laplacian_solves_exactly_and_bad_input_is_an_error() {
    let shape = Shape::generated("torus", torus(10, 5, 1.0, 0.3).unwrap());
    let tri = shape.triangulation();
    let l = stiffness_matrix(&tri);
    let solver = PinnedLaplacian::new(&l).unwrap();
    let n = tri.n_vertices();
    let mut b: Vec<f64> = (0..n).map(|i| ((i * 7919) % 13) as f64 - 6.0).collect();
    let mean = b.iter().sum::<f64>() / n as f64;
    b.iter_mut().for_each(|x| *x -= mean);
    assert!(solver.inconsistency(&b) < 1e-12);
    let mut x = b.clone();
    solver.solve_in_place(&mut x).unwrap();
    let mut lx = vec![0.0; n];
    l.mul_vec(&x, &mut lx).unwrap();
    for v in 0..n {
        assert!((lx[v] - b[v]).abs() < 1e-10, "vertex {v}");
    }
    assert!(matches!(
        solver.solve_in_place(&mut [0.0; 3]),
        Err(FemError::Invalid(_))
    ));
    assert!(matches!(
        condition_number(&l, &MassMatrix::Lumped(vec![1.0; 3]), &EigenOptions::default()),
        Err(FemError::Invalid(_))
    ));
    let mut bad = vec![1.0; n];
    bad[4] = 0.0;
    assert!(matches!(
        condition_number(&l, &MassMatrix::Lumped(bad), &EigenOptions::default()),
        Err(FemError::Invalid(_))
    ));
}
