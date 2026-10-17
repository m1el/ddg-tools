//! Tests of the connection Laplacians, harmonic vector fields, the vector maximum principle and
//! smoothest direction fields (SPEC §5.5).

mod common;

use std::f64::consts::{PI, TAU};

use common::*;
use ddg_core::sparse::{Method, Parallelism, Symbolic};
use ddg_core::SparseMatrix;
use ddg_intrinsic::{DelaunayOptions, IntrinsicTriangulation, RefineOptions};
use ddg_intrinsic_apps::vector::*;

fn delaunay(shape: &Shape) -> IntrinsicTriangulation {
    let mut tri = shape.triangulation();
    tri.flip_to_delaunay(&DelaunayOptions::default()).unwrap();
    tri
}

fn wrap_pi(a: f64) -> f64 {
    let x = (a + PI).rem_euclid(TAU) - PI;
    if x <= -PI { x + TAU } else { x }
}

/// The global angle of the frame of every vertex of a flat triangulation whose vertices are
/// the input vertices (no insertions): the angle, in the plane, of the frame direction 0.
fn frame_angles_in_plane(shape: &Shape, tri: &IntrinsicTriangulation) -> Vec<f64> {
    let dirs = halfedge_directions(tri);
    (0..tri.n_vertices())
        .map(|v| {
            let h = tri.outgoing_halfedge(v);
            let (a, b) = (shape.positions[v], shape.positions[tri.to_vertex(h)]);
            (b[1] - a[1]).atan2(b[0] - a[0]) - dirs[h]
        })
        .collect()
}

/// Solves the scalar Laplace equation with the cotan Laplacian of `tri` and boundary values.
fn scalar_harmonic(tri: &IntrinsicTriangulation, boundary: &[f64]) -> Vec<f64> {
    let l = tri.cotan_laplacian();
    let interior: Vec<u32> = (0..tri.n_vertices() as u32)
        .filter(|&v| !tri.is_boundary_vertex(v as usize))
        .collect();
    let bnd: Vec<u32> = (0..tri.n_vertices() as u32)
        .filter(|&v| tri.is_boundary_vertex(v as usize))
        .collect();
    let l_ii = l.submatrix(&interior, &interior).unwrap();
    let l_ib = l.submatrix(&interior, &bnd).unwrap();
    let xb: Vec<f64> = bnd.iter().map(|&v| boundary[v as usize]).collect();
    let mut rhs = vec![0.0; interior.len()];
    l_ib.mul_vec(&xb, &mut rhs).unwrap();
    for x in &mut rhs {
        *x = -*x;
    }
    let f = Symbolic::analyze(&l_ii, Method::Llt)
        .unwrap()
        .factor(&l_ii, Parallelism::Sequential)
        .unwrap();
    f.solve_in_place(&mut rhs, Parallelism::Sequential).unwrap();
    let mut out = boundary.to_vec();
    for (k, &v) in interior.iter().enumerate() {
        out[v as usize] = rhs[k];
    }
    out
}

/// On a flat domain the connection is trivial: transport is the change of frame, every face has
/// holonomy zero, and the harmonic vector field is the scalar harmonic function in each
/// Cartesian component (SPEC N41: boundary frames are not normalized).
#[test]
fn on_a_flat_mesh_the_vector_laplacian_is_the_scalar_one_componentwise() {
    let shape = sheared_grid(10, 0.8, 0.3, 7);
    for (what, tri) in [("input", shape.triangulation()), ("iDT", delaunay(&shape))] {
        let beta = frame_angles_in_plane(&shape, &tri);
        for h in 0..tri.n_halfedges() {
            let (i, j) = (tri.from_vertex(h), tri.to_vertex(h));
            let d = wrap_pi(transport_angle(&tri, h) - (beta[i] - beta[j]));
            assert!(d.abs() < 1e-12, "{what}: halfedge {h}: transport off by {d}");
        }
        for f in 0..tri.n_faces() {
            assert!(face_holonomy(&tri, f).abs() < 1e-12);
        }
        // Boundary vectors G(p) = (x² − y², 2xy + y), given in the frames.
        let g = |p: [f64; 3]| C64::new(p[0] * p[0] - p[1] * p[1], 2.0 * p[0] * p[1] + p[1]);
        let given: Vec<(u32, C64)> = (0..tri.n_vertices())
            .filter(|&v| tri.is_boundary_vertex(v))
            .map(|v| (v as u32, g(shape.positions[v]) * C64::from_angle(-beta[v])))
            .collect();
        for laplacian in [ConnectionLaplacian::FiniteDifference, ConnectionLaplacian::FiniteElement] {
            let options = HarmonicFieldOptions {
                n_sym: 1,
                laplacian,
                boundary: BoundaryValues::AtInputVertices(given.clone()),
            };
            let (field, report) = harmonic_field(&tri, &options).unwrap();
            assert!(report.residual < 1e-12, "{what} {laplacian:?}: residual {}", report.residual);
            let bx: Vec<f64> = shape.positions.iter().map(|&p| g(p).re).collect();
            let by: Vec<f64> = shape.positions.iter().map(|&p| g(p).im).collect();
            let (sx, sy) = (scalar_harmonic(&tri, &bx), scalar_harmonic(&tri, &by));
            for v in 0..tri.n_vertices() {
                let global = field.values[v] * C64::from_angle(beta[v]);
                assert!(
                    (global.re - sx[v]).abs() < 1e-10 && (global.im - sy[v]).abs() < 1e-10,
                    "{what} {laplacian:?} vertex {v}: {global:?} vs ({}, {})",
                    sx[v],
                    sy[v]
                );
            }
        }
    }
}

/// The 3D output of a flat field is the field in the plane.
#[test]
fn input_vectors_in_3d_are_the_field_in_the_plane() {
    let shape = sheared_grid(6, 0.5, 0.2, 3);
    let tri = delaunay(&shape);
    let beta = frame_angles_in_plane(&shape, &tri);
    let values: Vec<C64> = (0..tri.n_vertices())
        .map(|v| C64::from_angle(0.3 * v as f64).scale(1.0 + 0.01 * v as f64))
        .collect();
    for n_sym in [1u32, 2, 4] {
        let field = TangentField {
            n_sym,
            values: values.clone(),
        };
        let out = field.input_vectors_3d(&tri, &shape.positions).unwrap();
        for v in 0..tri.n_vertices() {
            let u = values[v];
            let angle = u.arg() / n_sym as f64 + beta[v];
            let expected = [u.abs() * angle.cos(), u.abs() * angle.sin(), 0.0];
            for k in 0..3 {
                assert!((out[v][k] - expected[k]).abs() < 1e-12, "n {n_sym} vertex {v}");
            }
        }
    }
    let short = TangentField { n_sym: 1, values: vec![C64::default(); 3] };
    assert!(short.input_vectors_3d(&tri, &shape.positions).is_err());
}

fn symmetric_defect(a: &SparseMatrix) -> f64 {
    let scale = a.values().iter().fold(0.0f64, |m, v| m.max(v.abs()));
    a.symmetry_defect().unwrap() / scale
}

fn is_positive_definite(a: &SparseMatrix) -> bool {
    Symbolic::analyze(a, Method::Llt)
        .unwrap()
        .factor(a, Parallelism::Sequential)
        .is_ok()
}

fn shifted(a: &SparseMatrix, m: &SparseMatrix, eps: f64) -> SparseMatrix {
    let mut t: Vec<(u32, u32, f64)> = a.triplets().collect();
    t.extend(m.triplets().map(|(i, j, v)| (i, j, eps * v)));
    SparseMatrix::from_triplets(a.nrows(), a.ncols(), &t).unwrap()
}

/// Both connection Laplacians are Hermitian; the finite element one is positive semidefinite on
/// every triangulation, the finite-difference one on the intrinsic Delaunay triangulation.
#[test]
fn connection_laplacians_are_hermitian_and_positive_semidefinite() {
    let mut shapes = vec![jittered_icosphere(2, 0.12, 5), bumpy_sheared_grid(8, 0.2, 1.5)];
    shapes.extend(dataset("camelbad"));
    for shape in &shapes {
        for (what, tri) in [("input", shape.triangulation()), ("iDT", delaunay(shape))] {
            for n_sym in [1u32, 2, 4] {
                for kind in [ConnectionLaplacian::FiniteDifference, ConnectionLaplacian::FiniteElement] {
                    let m = connection_matrices(&tri, n_sym, kind).unwrap();
                    assert!(symmetric_defect(&m.stiffness) < 1e-14, "{} {what} {kind:?}", shape.name);
                    assert!(symmetric_defect(&m.mass) < 1e-14);
                    let tr = m.stiffness.diagonal().iter().sum::<f64>()
                        / m.mass.diagonal().iter().sum::<f64>();
                    let pd = is_positive_definite(&shifted(&m.stiffness, &m.mass, 1e-9 * tr));
                    let must = kind == ConnectionLaplacian::FiniteElement || what == "iDT";
                    if must {
                        assert!(pd, "{} {what} n {n_sym} {kind:?}: A + εM is not definite", shape.name);
                    } else {
                        eprintln!("{} {what} n {n_sym} {kind:?}: A + εM definite: {pd}", shape.name);
                    }
                    assert!(is_positive_definite(&m.mass));
                }
            }
        }
    }
}

/// The holonomies of the faces add up to `2πχ` on a closed surface, and to the interior
/// angle defects on a surface with boundary (boundary frames are not normalized).
#[test]
fn face_holonomies_add_up_to_the_curvature() {
    let mut shapes = vec![icosphere(2), jittered_icosphere(2, 0.12, 9), bumpy_sheared_grid(6, 0.3, 0.5)];
    shapes.extend(dataset("camelbad"));
    for shape in &shapes {
        let mut tri = shape.triangulation();
        tri.refine(&RefineOptions::default()).unwrap();
        let total: f64 = (0..tri.n_faces()).map(|f| face_holonomy(&tri, f)).sum();
        let expected: f64 = tri.angle_defects().iter().sum();
        assert!((total - expected).abs() < 1e-9, "{}: {total} vs {expected}", shape.name);
        if (0..tri.n_vertices()).all(|v| !tri.is_boundary_vertex(v)) {
            assert!((total - TAU * tri.euler_characteristic() as f64).abs() < 1e-9);
        }
    }
}

/// On a sphere the smoothest `n`-direction field has singularities whose indices add up to
/// `n χ = 2n` (units of `1/n`), and, on a round sphere, exactly the minimal number of them.
#[test]
fn smoothest_fields_on_spheres_satisfy_poincare_hopf() {
    for shape in [icosphere(3), jittered_icosphere(3, 0.03, 4)] {
        for (what, tri) in [("input", shape.triangulation()), ("iDT", delaunay(&shape))] {
            for n_sym in [1u32, 2, 4] {
                for laplacian in [ConnectionLaplacian::FiniteElement, ConnectionLaplacian::FiniteDifference] {
                    // On a round sphere the smallest eigenvalues of n-vectors form a cluster
                    // of 2n + 1 (spin-weighted harmonics of degree n), split only by the
                    // discretization: the block must be larger.
                    let options = DirectionFieldOptions {
                        n_sym,
                        laplacian,
                        block: 2 * n_sym as usize + 2,
                        ..Default::default()
                    };
                    let (field, report) = match smoothest_field(&tri, &options) {
                        Ok(x) => x,
                        Err(VectorError::NotPositiveDefinite { .. })
                            if laplacian == ConnectionLaplacian::FiniteDifference && what == "input" =>
                        {
                            continue;
                        }
                        Err(e) => panic!("{e}"),
                    };
                    assert_eq!(report.status, EigenStatus::Converged, "{} {what} n {n_sym}: {report:?}", shape.name);
                    let s = singularities(&tri, &field).unwrap();
                    assert!(s.max_rounding < 1e-6);
                    assert_eq!(s.index_sum, 2 * n_sym as i64, "{} {what} n {n_sym} {laplacian:?}", shape.name);
                    eprintln!(
                        "{} {what} n {n_sym} {laplacian:?}: λ {:.4} ({} iterations), {} singular faces",
                        shape.name, report.eigenvalue, report.iterations, s.faces.len()
                    );
                    // The smooth value: the rough Laplacian of n-vectors on the unit sphere has
                    // the smallest eigenvalue n (spin-weighted harmonics, l = n).
                    assert!(
                        (report.eigenvalue - n_sym as f64).abs() < 0.03 * n_sym as f64,
                        "{} {what} n {n_sym}: λ {}",
                        shape.name,
                        report.eigenvalue
                    );
                    if shape.name == "icosphere_3" {
                        assert!(s.faces.iter().all(|x| x.index > 0));
                        assert_eq!(s.faces.len(), 2 * n_sym as usize, "{what} n {n_sym}");
                    }
                }
            }
        }
    }
}

/// VHM Sec. 5.4: on the intrinsic Delaunay triangulation every vector of the harmonic field is
/// in the convex cone of its transported neighbours.
#[test]
fn the_vector_maximum_principle_holds_on_the_intrinsic_delaunay_triangulation() {
    let mut shapes = vec![
        bumpy_sheared_grid(10, 0.25, 1.5),
        sheared_grid(12, 1.5, 0.45, 2),
        sheared_grid(12, 4.0, 0.45, 2),
    ];
    shapes.extend(dataset("face"));
    for shape in &shapes {
        for n_sym in [1u32, 2] {
            for boundary in [BoundaryValues::InwardNormal, BoundaryValues::Tangent] {
                let options = HarmonicFieldOptions {
                    n_sym,
                    boundary: boundary.clone(),
                    ..Default::default()
                };
                let input = shape.triangulation();
                let (f_in, _) = harmonic_field(&input, &options).unwrap();
                let r_in = check_maximum_principle(&input, &f_in, 1e-9);
                let mut tri = delaunay(shape);
                let (f, report) = harmonic_field(&tri, &options).unwrap();
                assert_eq!(report.factorization, Factorized::Cholesky);
                let r = check_maximum_principle(&tri, &f, 1e-9);
                eprintln!(
                    "{} n {n_sym} {boundary:?}: input {} of {} outside (margin {:.2e}), iDT {} (margin {:.2e})",
                    shape.name, r_in.n_outside_cone, r_in.n_checked, r_in.min_margin, r.n_outside_cone, r.min_margin
                );
                assert_eq!(r.n_outside_cone, 0, "{} n {n_sym}", shape.name);
                if shape.name == "sheared_grid_12" && n_sym == 1 && r_in.n_checked == 121 {
                    // Strongly sheared (shear 4): the input violates the principle.
                    let shear4 = shapes.iter().position(|x| std::ptr::eq(x, shape)) == Some(2);
                    if shear4 {
                        assert!(r_in.n_outside_cone > 0, "{boundary:?}");
                    }
                }
                // And after refinement.
                tri.refine(&RefineOptions::default()).unwrap();
                let (f, _) = harmonic_field(&tri, &options).unwrap();
                assert_eq!(check_maximum_principle(&tri, &f, 1e-9).n_outside_cone, 0);
            }
        }
    }
}

/// Boundary values at the input vertices are interpolated to inserted boundary vertices by
/// transport along the boundary: a constant field on a flat domain stays constant.
#[test]
fn boundary_values_are_interpolated_along_the_boundary() {
    let shape = sheared_grid(5, 0.4, 0.2, 11);
    let input = shape.triangulation();
    let beta = frame_angles_in_plane(&shape, &input);
    let mut tri = shape.triangulation();
    tri.refine(&RefineOptions {
        max_circumradius: Some(0.05),
        ..Default::default()
    })
    .unwrap();
    assert!(tri.n_vertices() > 3 * input.n_vertices());
    let constant = C64::from_angle(0.7);
    let given: Vec<(u32, C64)> = (0..input.n_vertices())
        .filter(|&v| input.is_boundary_vertex(v))
        .map(|v| (v as u32, constant * C64::from_angle(-beta[v])))
        .collect();
    let options = HarmonicFieldOptions {
        boundary: BoundaryValues::AtInputVertices(given),
        ..Default::default()
    };
    let (field, _) = harmonic_field(&tri, &options).unwrap();
    // The field is parallel: zero energy, and at the input vertices the given constant.
    let m = connection_matrices(&tri, 1, ConnectionLaplacian::FiniteDifference).unwrap();
    let energy = m.energy(&field.values) / m.norm_squared(&field.values);
    assert!(energy < 1e-12, "{energy}");
    for (v, &b) in beta.iter().enumerate().take(input.n_vertices()) {
        let global = field.values[v] * C64::from_angle(b);
        assert!((global - constant).abs() < 1e-10);
    }
    // A missing vertex is an error naming it.
    let options = HarmonicFieldOptions {
        boundary: BoundaryValues::AtInputVertices(vec![(0, constant)]),
        ..Default::default()
    };
    assert!(matches!(harmonic_field(&tri, &options), Err(VectorError::InvalidInput { .. })));
}
