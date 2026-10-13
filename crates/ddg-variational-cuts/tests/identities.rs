//! SPEC §6 item 1 and the structural properties of an evaluation: the curvature identity, the
//! areas and Euler characteristics of the cut mesh, the consistency of the cut mesh that is
//! handed out with the one the evaluation used, and the properties of the flux recovery.

mod common;

use std::f64::consts::PI;

use common::{Shape, bumpy, sphere};
use ddg_core::generators::{annulus, icosphere, torus};
use ddg_core::operators::cotan_laplacian;
use ddg_core::sparse::{Method, Parallelism, Symbolic};
use ddg_core::{Geometry, TriMesh};
use ddg_variational_cuts::{
    CutProblem, CutSegment, CutVertex, EdgeKind, Energy, EvalOptions, Evaluation, Evaluator,
    FluxRecovery, NONE, face_curvature_sum_residual, face_curvatures,
};

fn shapes() -> Vec<Shape> {
    vec![
        sphere(2),
        Shape::from_generated("torus", torus(24, 12, 1.0, 0.4).unwrap()),
        bumpy(12, 0.2),
        Shape::from_generated("annulus", annulus(30, 5, 0.4, 1.0).unwrap()),
    ]
}

/// A level-set function that cuts every test shape (and meets the boundary where there is one).
fn tilted(p: [f64; 3]) -> f64 {
    0.31 * p[0] + 0.52 * p[1] + 0.4 * p[2] - 0.37 + 0.2 * p[0] * p[1]
}

fn evaluate(shape: &Shape, phi: &[f64], energy: Energy, options: &EvalOptions<'_>) -> Evaluation {
    shape
        .evaluator()
        .evaluate(phi, &CutProblem::new(energy, 0.5), options)
        .unwrap()
}

#[test]
fn face_curvatures_redistribute_the_angle_defects() {
    for shape in shapes() {
        let residual = face_curvature_sum_residual(&shape.mesh, &shape.geometry).unwrap();
        assert!(residual.abs() < 1e-11, "{}: {residual}", shape.name);
        let total: f64 = face_curvatures(&shape.mesh, &shape.geometry)
            .unwrap()
            .iter()
            .sum();
        if shape.mesh.n_boundary_edges() == 0 {
            let expect = 2.0 * PI * shape.mesh.euler_characteristic() as f64;
            assert!((total - expect).abs() < 1e-11, "{}: {total}", shape.name);
        }
    }
    // On a sphere every face has positive curvature (SPEC D1).
    let ball = sphere(2);
    assert!(
        face_curvatures(&ball.mesh, &ball.geometry)
            .unwrap()
            .iter()
            .all(|&omega| omega > 0.0)
    );
    // The density of the evaluator is per unit area of the unit-area surface: it integrates to
    // the same total.
    let evaluator = ball.evaluator();
    let scale = evaluator.scale();
    let total: f64 = evaluator
        .gaussian_curvature()
        .iter()
        .zip(ball.geometry.face_areas())
        .map(|(k, &a)| k * scale.area_to_unit_area(a))
        .sum();
    assert!((total - 4.0 * PI).abs() < 1e-11);
    assert!((evaluator.vertex_areas().iter().sum::<f64>() - 1.0).abs() < 1e-13);
}

#[test]
fn cut_mesh_areas_counts_and_euler_characteristics() {
    for shape in shapes() {
        let phi = shape.sample(tilted);
        let e = evaluate(&shape, &phi, Energy::Dirichlet, &EvalOptions::default());
        let (nv, nf) = (shape.mesh.n_vertices(), shape.mesh.n_faces());
        let (n_points, n_segments) = (e.cut.points.len(), e.cut.segments.len());
        assert!(n_points > 0, "{}", shape.name);
        assert_eq!(e.report.n_cut_points, n_points);
        assert_eq!(e.report.n_segments, n_segments);

        // Areas: the regions share the unit area; the triangles of every cut face sum to it.
        let total: f64 =
            e.energy.area.iter().sum::<f64>() + e.report.excluded_area.iter().sum::<f64>();
        assert!((total - 1.0).abs() < 1e-13, "{}: {total}", shape.name);
        let tri = &e.triangles;
        let mut per_parent = vec![0.0; nf];
        for (t, &parent) in tri.parent.iter().enumerate() {
            per_parent[parent as usize] += tri.area[t];
        }
        for (f, area) in per_parent.iter().enumerate() {
            let expect = e.scale.area_to_unit_area(shape.geometry.face_areas()[f]);
            assert!((area - expect).abs() < 1e-14 * expect.max(1e-3));
        }

        // Counts: two copies per cut point, two more triangles per cut face.
        assert_eq!(tri.n_vertices(), nv + 2 * n_points);
        assert_eq!(tri.n_faces(), nf + 2 * n_segments);
        // Every face has zero or one segment; an interior cut point has two segments, one on
        // a boundary edge has one.
        let mut per_point = vec![0; n_points];
        for segment in &e.cut.segments {
            for p in segment.points {
                per_point[p as usize] += 1;
            }
        }
        let mut on_boundary = 0;
        for (k, point) in e.cut.points.iter().enumerate() {
            assert_eq!(e.cut.edge_point[point.edge as usize], k as u32);
            assert!(point.t > 0.0 && point.t < 1.0);
            if shape.mesh.is_boundary_edge(point.edge as usize) {
                on_boundary += 1;
                assert_eq!(per_point[k], 1);
            } else {
                assert_eq!(per_point[k], 2);
            }
        }
        assert_eq!(
            e.cut.edge_point.iter().filter(|&&p| p != NONE).count(),
            n_points
        );

        // The cut mesh as a mesh: valid, with the regions separated.
        let cut_mesh = e.cut_mesh().unwrap();
        assert_eq!(cut_mesh.mesh.n_vertices(), nv + 2 * n_points);
        assert_eq!(cut_mesh.face_parent, tri.parent);
        // Cutting along a closed curve keeps the Euler characteristic, every arc that ends on
        // the boundary adds one.
        assert_eq!(on_boundary % 2, 0);
        assert_eq!(
            cut_mesh.mesh.euler_characteristic(),
            shape.mesh.euler_characteristic() + (on_boundary / 2) as i64,
            "{}",
            shape.name
        );
        assert!(cut_mesh.mesh.n_components() >= 2);
        // No edge of the cut mesh joins two regions, and faces have the region of their corners.
        for (f, face) in cut_mesh.mesh.faces().iter().enumerate() {
            for &v in face {
                assert_eq!(cut_mesh.vertex_region[v as usize], cut_mesh.face_region[f]);
            }
        }
        for v in 0..cut_mesh.mesh.n_vertices() {
            match cut_mesh.vertex_source[v] {
                CutVertex::Original(i) => {
                    assert_eq!(i as usize, v);
                    assert_eq!(cut_mesh.vertex_region[v], e.cut.vertex_region[v]);
                }
                CutVertex::Point(k) => {
                    assert!(tri.point_vertices(k as usize).contains(&v));
                    // Copies of cut points are boundary vertices of the cut mesh.
                    assert!(cut_mesh.mesh.is_boundary_vertex(v));
                }
            }
        }
        // The geometry handed out is in the caller's unit.
        let area = cut_mesh.geometry.total_area();
        let expect = shape.geometry.total_area();
        assert!((area / expect - 1.0).abs() < 1e-12, "{}", shape.name);
        // Positions of the cut mesh reproduce its edge lengths.
        let positions = cut_mesh
            .positions(&e.cut, &shape.mesh, &shape.positions)
            .unwrap();
        for (edge, length) in cut_mesh
            .mesh
            .edges()
            .iter()
            .zip(cut_mesh.geometry.edge_lengths())
        {
            let (p, q) = (positions[edge[0] as usize], positions[edge[1] as usize]);
            let d = (0..3).map(|k| (p[k] - q[k]).powi(2)).sum::<f64>().sqrt();
            assert!((d - length).abs() < 1e-12 * expect.sqrt(), "{}", shape.name);
        }
    }
}

#[test]
fn boundary_edges_dual_lengths_and_total_flux() {
    for shape in shapes() {
        let phi = shape.sample(tilted);
        let e = evaluate(&shape, &phi, Energy::Hencky, &EvalOptions::default());
        let tri = &e.triangles;
        let nv = tri.n_original;
        // Edge kinds against the mesh that is handed out.
        let cut_mesh = e.cut_mesh().unwrap();
        for t in 0..tri.n_faces() {
            for k in 0..3 {
                let boundary = cut_mesh.mesh.twin(3 * t + k).is_none();
                assert_eq!(boundary, tri.edge_kind[t][k] != EdgeKind::Interior);
                if tri.edge_kind[t][k] == EdgeKind::Cut {
                    assert!(tri.faces[t][k] as usize >= nv);
                    assert!(tri.faces[t][(k + 1) % 3] as usize >= nv);
                }
            }
        }
        // Dual lengths: positive exactly at the boundary vertices of the cut mesh; the two
        // copies of an interior cut point have the same one; they sum to the boundary length.
        for v in 0..tri.n_vertices() {
            assert_eq!(e.dual_length[v] > 0.0, cut_mesh.mesh.is_boundary_vertex(v));
            if e.dual_length[v] == 0.0 {
                assert_eq!((e.flux[v], e.normal_derivative[v]), (0.0, 0.0));
            }
        }
        let mut cut_length = 0.0;
        for (k, point) in e.cut.points.iter().enumerate() {
            let [v0, v1] = tri.point_vertices(k);
            if !shape.mesh.is_boundary_edge(point.edge as usize) {
                assert!((e.dual_length[v0] - e.dual_length[v1]).abs() < 1e-14);
                cut_length += e.dual_length[v0];
            }
        }
        if shape.mesh.n_boundary_edges() == 0 {
            assert!((cut_length - e.energy.length).abs() < 1e-13);
        }
        assert!((tri.segment_length.iter().sum::<f64>() - e.energy.length).abs() < 1e-14);

        // Gauss–Bonnet for u: the total flux through the boundary of a region is the total
        // source in it, Σ_B v_i = Σ_i g_i = Σ_t k_t A_t. The same for the adjoint with its
        // source r_i = −A_i b'(u_i) = −2 A_i u_i.
        let evaluator = shape.evaluator();
        let q = e.adjoint_flux.as_ref().unwrap();
        for region in 0..2 {
            let source: f64 = (0..tri.n_faces())
                .filter(|&t| tri.region[t] == region)
                .map(|t| evaluator.gaussian_curvature()[tri.parent[t] as usize] * tri.area[t])
                .sum();
            let in_region = |v: &usize| tri.vertex_region[*v] == region;
            let flux: f64 = (0..tri.n_vertices())
                .filter(in_region)
                .map(|v| e.flux[v])
                .sum();
            assert!(
                (flux - source).abs() < 1e-11,
                "{}: {flux} vs {source}",
                shape.name
            );
            let mut vertex_area = vec![0.0; tri.n_vertices()];
            for (t, face) in tri.faces.iter().enumerate() {
                for &v in face {
                    vertex_area[v as usize] += tri.area[t] / 3.0;
                }
            }
            let adjoint_source: f64 = (0..tri.n_vertices())
                .filter(in_region)
                .map(|v| -2.0 * vertex_area[v] * e.log_scale[v])
                .sum();
            let adjoint_flux: f64 = (0..tri.n_vertices()).filter(in_region).map(|v| q[v]).sum();
            assert!(
                (adjoint_flux + adjoint_source).abs() < 1e-11,
                "{}",
                shape.name
            );
        }
        assert!(e.report.yamabe_residual < 1e-10 && e.report.adjoint_residual.unwrap() < 1e-10);
    }
}

#[test]
fn the_cut_mesh_handed_out_reproduces_the_scale_factor() {
    // An independent path: the cotan Laplacian of ddg-core on the mesh and geometry that
    // `cut_mesh()` returns, its own lumped source, a solve on its interior vertices.
    for shape in shapes() {
        let phi = shape.sample(tilted);
        let e = evaluate(&shape, &phi, Energy::Dirichlet, &EvalOptions::default());
        let evaluator = shape.evaluator();
        let cut_mesh = e.cut_mesh().unwrap();
        let (mesh, geometry) = (&cut_mesh.mesh, &cut_mesh.geometry);
        let laplacian = cotan_laplacian(mesh, geometry).unwrap();
        let mut source = vec![0.0; mesh.n_vertices()];
        for (f, face) in mesh.faces().iter().enumerate() {
            let area = e.scale.area_to_unit_area(geometry.face_areas()[f]);
            let density = evaluator.gaussian_curvature()[cut_mesh.face_parent[f] as usize];
            for &v in face {
                source[v as usize] += density * area / 3.0;
            }
        }
        let interior = mesh.interior_vertices();
        let matrix = laplacian.submatrix(&interior, &interior).unwrap();
        let factorization = Symbolic::analyze(&matrix, Method::Llt)
            .unwrap()
            .factor(&matrix, Parallelism::Sequential)
            .unwrap();
        let mut u: Vec<f64> = interior.iter().map(|&v| -source[v as usize]).collect();
        factorization
            .solve_in_place(&mut u, Parallelism::Sequential)
            .unwrap();
        let scale = e.log_scale.iter().fold(0.0f64, |m, x| m.max(x.abs()));
        for (&v, value) in interior.iter().zip(&u) {
            assert!(
                (value - e.log_scale[v as usize]).abs() < 1e-9 * scale,
                "{}: vertex {v}: {value} vs {}",
                shape.name,
                e.log_scale[v as usize]
            );
        }
        // The unknowns of the evaluator are the interior vertices of the input mesh, and its
        // factorization is the one of this evaluation.
        assert_eq!(
            evaluator.vertex_of_unknown(),
            shape.mesh.interior_vertices()
        );
        assert!(evaluator.factorization().is_none());
    }
}

#[test]
fn the_factorization_of_the_last_evaluation_is_kept() {
    let shape = sphere(2);
    let mut evaluator = shape.evaluator();
    assert!(evaluator.factorization().is_none());
    let phi = shape.sample(tilted);
    let problem = CutProblem::new(Energy::Dirichlet, 0.0);
    let e = evaluator
        .evaluate(&phi, &problem, &EvalOptions::default())
        .unwrap();
    // Solve L_II x = −g_I again through the kept factorization: it reproduces u.
    let tri = &e.triangles;
    let mut source = vec![0.0; tri.n_vertices()];
    for (t, face) in tri.faces.iter().enumerate() {
        for &v in face {
            source[v as usize] +=
                evaluator.gaussian_curvature()[tri.parent[t] as usize] * tri.area[t] / 3.0;
        }
    }
    let mut x: Vec<f64> = evaluator
        .vertex_of_unknown()
        .iter()
        .map(|&v| -source[v as usize])
        .collect();
    let factorization = evaluator.factorization().unwrap();
    assert_eq!(factorization.dim(), evaluator.vertex_of_unknown().len());
    factorization
        .solve_in_place(&mut x, Parallelism::Sequential)
        .unwrap();
    for (&v, value) in evaluator.vertex_of_unknown().iter().zip(&x) {
        assert!(evaluator.unknown_of_vertex()[v as usize] != NONE);
        assert!((value - e.log_scale[v as usize]).abs() < 1e-13);
    }
    // A failed evaluation invalidates nothing that was not touched: bad input is rejected
    // before the matrix is assembled.
    assert!(
        evaluator
            .evaluate(&phi[1..], &problem, &EvalOptions::default())
            .is_err()
    );
    assert!(evaluator.factorization().is_some());
    // L applied to the extension of the speed is what the flow needs next (SPEC §3.8).
    let mut s = vec![0.0; tri.n_vertices()];
    for (k, &speed) in e.speed.iter().enumerate() {
        for v in tri.point_vertices(k) {
            s[v] = speed;
        }
    }
    let mut product = vec![0.0; tri.n_vertices()];
    tri.apply_laplacian(None, &s, &mut product).unwrap();
    assert!(product.iter().all(|x| x.is_finite()));
    // Zero column sums: the Laplacian annihilates constants on each side.
    assert!(
        product.iter().sum::<f64>().abs() < 1e-9 * product.iter().map(|x| x.abs()).sum::<f64>()
    );
}

#[test]
fn flux_recovery_conserves_and_stays_in_range() {
    for shape in shapes() {
        let phi = shape.sample(tilted);
        let nodal = evaluate(
            &shape,
            &phi,
            Energy::Hencky,
            &EvalOptions {
                flux_recovery: FluxRecovery::Nodal,
                ..EvalOptions::default()
            },
        );
        let smoothed = evaluate(&shape, &phi, Energy::Hencky, &EvalOptions::default());
        let zero = evaluate(
            &shape,
            &phi,
            Energy::Hencky,
            &EvalOptions {
                flux_recovery: FluxRecovery::Smoothed { length: 0.0 },
                ..EvalOptions::default()
            },
        );
        let tri = &smoothed.triangles;
        let nv = tri.n_original;
        // The integrated quantities do not depend on the recovery.
        assert_eq!(nodal.flux, smoothed.flux);
        assert_eq!(nodal.adjoint_flux, smoothed.adjoint_flux);
        // Length 0 is the nodal recovery, and the nodal recovery is v / l.
        assert_eq!(zero.normal_derivative, nodal.normal_derivative);
        assert_eq!(zero.speed, nodal.speed);
        for v in 0..tri.n_vertices() {
            if nodal.dual_length[v] > 0.0 {
                assert_eq!(
                    nodal.normal_derivative[v],
                    nodal.flux[v] / nodal.dual_length[v]
                );
            }
            // Off the cut the recovery changes nothing.
            if v < nv {
                assert_eq!(smoothed.normal_derivative[v], nodal.normal_derivative[v]);
            }
        }
        for (d, v) in [
            (&smoothed.normal_derivative, &smoothed.flux),
            (
                smoothed.adjoint_normal_derivative.as_ref().unwrap(),
                smoothed.adjoint_flux.as_ref().unwrap(),
            ),
        ] {
            for region in 0..2 {
                let copies = (nv..tri.n_vertices()).filter(|&c| tri.vertex_region[c] == region);
                // Σ l d = Σ v over the cut vertices of a region.
                let (mut lhs, mut rhs) = (0.0, 0.0);
                let (mut lo, mut hi) = (f64::INFINITY, f64::NEG_INFINITY);
                let (mut dlo, mut dhi) = (f64::INFINITY, f64::NEG_INFINITY);
                for c in copies {
                    lhs += smoothed.dual_length[c] * d[c];
                    rhs += v[c];
                    let ratio = v[c] / smoothed.dual_length[c];
                    (lo, hi) = (lo.min(ratio), hi.max(ratio));
                    (dlo, dhi) = (dlo.min(d[c]), dhi.max(d[c]));
                }
                assert!(
                    (lhs - rhs).abs() < 1e-11 * rhs.abs().max(1.0),
                    "{}",
                    shape.name
                );
                // A convex combination of the nodal values.
                let slack = 1e-9 * (hi - lo).abs().max(1.0);
                assert!(dlo >= lo - slack && dhi <= hi + slack, "{}", shape.name);
            }
        }
    }
}

#[test]
fn a_cut_through_vertices_does_not_blow_up() {
    // φ = 0 exactly on a ring of vertices: every crossing at such a vertex is moved to the
    // minimal offset, several of them around one vertex. The nodal quotient v/l is then wrong
    // by orders of magnitude; the default recovery is not (SPEC D12).
    let shape = sphere(3);
    let ring = shape
        .positions
        .iter()
        .map(|p| p[2])
        .min_by(|a, b| (a - 0.55).abs().total_cmp(&(b - 0.55).abs()))
        .unwrap();
    let phi = shape.sample(|p| ring - p[2]);
    assert!(phi.iter().filter(|&&x| x == 0.0).count() >= 2);
    let r = common::unit_area_radius();
    let psi0 = ring.acos();
    let exact = [(0.5 * psi0).tan() / r, 1.0 / ((0.5 * psi0).tan() * r)];
    let exact_speed = exact[0] * exact[0] - exact[1] * exact[1];
    let worst = |e: &Evaluation| -> f64 {
        e.speed
            .iter()
            .map(|s| (s / exact_speed - 1.0).abs())
            .fold(0.0, f64::max)
    };
    let default = evaluate(&shape, &phi, Energy::Dirichlet, &EvalOptions::default());
    assert!(default.report.n_clamped_crossings >= 4);
    assert_eq!(default.report.smallest_crossing_offset, 0.0);
    assert!(worst(&default) < 0.05, "{}", worst(&default));
    for k in 0..default.cut.points.len() {
        for (side, v) in default
            .triangles
            .point_vertices_by_region(k)
            .into_iter()
            .enumerate()
        {
            let error = (default.normal_derivative[v] / exact[side] - 1.0).abs();
            assert!(error < 0.03, "cut point {k}, side {side}: {error}");
        }
    }
    let paper = evaluate(&shape, &phi, Energy::Dirichlet, &EvalOptions::paper());
    assert!(worst(&paper) > 1e3, "{}", worst(&paper));
}

#[test]
fn closed_component_without_a_cut_is_reported_and_left_out() {
    // Two spheres; the cut touches only the first.
    let a = icosphere(2).unwrap();
    let n = a.positions.len();
    let mut positions = a.positions.clone();
    positions.extend(a.positions.iter().map(|p| [p[0] + 3.0, p[1], p[2]]));
    let mut faces = a.faces.clone();
    faces.extend(a.faces.iter().map(|f| f.map(|v| v + n as u32)));
    let shape = Shape::new("two spheres", positions, faces);
    let phi: Vec<f64> = shape
        .positions
        .iter()
        .enumerate()
        .map(|(v, p)| if v < n { 0.3 - p[2] } else { 1.0 })
        .collect();
    for energy in [Energy::Dirichlet, Energy::Hencky] {
        let e = evaluate(&shape, &phi, energy, &EvalOptions::default());
        assert_eq!(e.report.uncut_closed_components, vec![1]);
        // The second sphere lies in region 1 and has half of the area.
        assert!(
            e.report.excluded_area[0] == 0.0 && (e.report.excluded_area[1] - 0.5).abs() < 1e-13
        );
        assert!((e.energy.area.iter().sum::<f64>() - 0.5).abs() < 1e-13);
        assert!(e.log_scale[n..2 * n].iter().all(|&u| u == 0.0));
        if let Some(p) = &e.adjoint {
            assert!(p[n..2 * n].iter().all(|&x| x == 0.0));
        }
        // The first sphere is evaluated as if it were alone (it is half of the unit area, so
        // the dimensionless Dirichlet energy and u agree with the single sphere).
        let alone = Shape::from_generated("one sphere", a.clone());
        let single = evaluate(
            &alone,
            &phi[..n],
            Energy::Dirichlet,
            &EvalOptions::default(),
        );
        for v in 0..n {
            assert!((e.log_scale[v] - single.log_scale[v]).abs() < 1e-12);
        }
        assert!(e.report.yamabe_residual < 1e-10);
    }
    // No cut at all on a closed mesh: nothing to solve, nothing in the energy.
    let alone = Shape::from_generated("one sphere", a);
    let e = evaluate(
        &alone,
        &vec![1.0; n],
        Energy::Dirichlet,
        &EvalOptions::default(),
    );
    assert_eq!(e.report.uncut_closed_components, vec![0]);
    assert_eq!((e.cut.points.len(), e.speed.len()), (0, 0));
    assert_eq!(e.energy.total, 0.0);
    assert!(e.report.excluded_area[0] == 0.0 && (e.report.excluded_area[1] - 1.0).abs() < 1e-13);
    assert_eq!(e.energy.area, vec![0.0, 0.0]);
}

#[test]
fn no_cut_on_a_disk_is_the_plain_yamabe_problem() {
    let shape = bumpy(10, 0.25);
    let n = shape.positions.len();
    let e = evaluate(
        &shape,
        &vec![-1.0; n],
        Energy::Dirichlet,
        &EvalOptions::default(),
    );
    assert!(e.report.uncut_closed_components.is_empty());
    assert_eq!(e.cut.points.len(), 0);
    assert!((e.energy.area[0] - 1.0).abs() < 1e-13 && e.energy.area[1] == 0.0);
    assert!(e.energy.distortion[0] > 0.0 && e.energy.distortion[1] == 0.0);
    assert_eq!(e.energy.total, e.energy.distortion[0]);
    // u = 0 on the boundary of the mesh, and the flux through it balances the curvature.
    for v in 0..n {
        assert_eq!(shape.mesh.is_boundary_vertex(v), e.dual_length[v] > 0.0);
        if shape.mesh.is_boundary_vertex(v) {
            assert_eq!(e.log_scale[v], 0.0);
        }
    }
    let cut_mesh = e.cut_mesh().unwrap();
    assert_eq!(cut_mesh.mesh.faces(), shape.mesh.faces());
}

#[test]
fn one_cut_face_in_detail() {
    // An equilateral triangle with vertex 0 alone in region 0.
    let mesh = TriMesh::new(3, vec![[0, 1, 2]]).unwrap();
    let h = 3f64.sqrt() / 2.0;
    let positions = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.5, h, 0.0]];
    let geometry = Geometry::from_positions(&mesh, &positions).unwrap();
    let mut evaluator = Evaluator::new(&mesh, &geometry).unwrap();
    let problem = CutProblem::new(Energy::Dirichlet, 1.0);
    let options = EvalOptions::default();

    // phi = (-1, 1, 3): crossings at 1/2 of the edge [0, 1] and at 1/4 of [0, 2].
    let e = evaluator
        .evaluate(&[-1.0, 1.0, 3.0], &problem, &options)
        .unwrap();
    assert_eq!(e.cut.vertex_region, vec![0, 1, 1]);
    assert_eq!(e.cut.points.len(), 2);
    assert_eq!((e.cut.points[0].edge, e.cut.points[0].t), (0, 0.5));
    assert_eq!((e.cut.points[1].edge, e.cut.points[1].t), (1, 0.25));
    // a = 0, p on 0 -> 1 (point 0), q on 2 -> 0 (point 1); region 0 is left of p -> q.
    assert_eq!(
        e.cut.segments,
        vec![CutSegment {
            points: [0, 1],
            face: 0
        }]
    );
    // Copies: point k has the vertices 3 + 2k (side of the lower end of its edge, here vertex
    // 0) and 3 + 2k + 1. |pc|^2 = 3/4 < |bq|^2 = 13/16: the quad is split by p-c.
    let tri = &e.triangles;
    assert_eq!(tri.faces, vec![[0, 3, 5], [4, 1, 2], [4, 2, 6]]);
    assert_eq!(tri.region, vec![0, 1, 1]);
    assert_eq!(tri.vertex_region, vec![0, 1, 1, 0, 1, 0, 1]);
    assert_eq!(e.report.n_diagonal_near_ties, 0);
    let area = 3f64.sqrt() / 4.0;
    let expect = [0.125, 0.5, 0.375];
    // Areas on the unit-area surface: s_p s_q, 1 - s_p, s_p (1 - s_q).
    for (computed, exact) in tri.area.iter().zip(expect) {
        assert!((computed - exact).abs() < 1e-15);
    }
    // The segment p q: |pq|^2 = 1/4 + 1/16 - 1/8 = 3/16 in the caller's unit.
    let unit = e.scale.length_scale;
    assert!((unit * unit * area - 1.0).abs() < 1e-15);
    assert!((tri.segment_length[0] - unit * (3f64 / 16.0).sqrt()).abs() < 1e-15);
    assert!((e.energy.length - tri.segment_length[0]).abs() < 1e-15);
    assert_eq!(
        tri.edge_kind[0],
        [EdgeKind::Boundary, EdgeKind::Cut, EdgeKind::Boundary]
    );
    assert_eq!(
        tri.edge_kind[2],
        [EdgeKind::Interior, EdgeKind::Boundary, EdgeKind::Cut]
    );
    // Cotangents of the corner triangle (0, p, q) with sides 1/2, sqrt(3)/4, 1/4: a right
    // angle at q, 60 degrees at 0, 30 degrees at p.
    let cot = tri.cot[0];
    assert!((cot[0] - 1.0 / 3f64.sqrt()).abs() < 1e-14);
    assert!((cot[1] - 3f64.sqrt()).abs() < 1e-14);
    assert!(cot[2].abs() < 1e-14);
    // Dual lengths at the cut point on the boundary edge [0, 1], side of vertex 0: half of
    // the segment and half of the piece 0-p of the boundary edge (SPEC D9d).
    let expect = 0.5 * unit * ((3f64 / 16.0).sqrt() + 0.5);
    assert!((e.dual_length[3] - expect).abs() < 1e-15);

    // phi = (-1, 1, 1): the crossings are symmetric and the two diagonals are equal. The tie
    // goes to p-c, and it is counted.
    let e = evaluator
        .evaluate(&[-1.0, 1.0, 1.0], &problem, &options)
        .unwrap();
    assert_eq!(e.report.n_diagonal_near_ties, 1);
    assert_eq!(e.triangles.faces, vec![[0, 3, 5], [4, 1, 2], [4, 2, 6]]);
    // With vertex 1 alone in region 1 the frame is (a, b, c) = (1, 2, 0).
    let e = evaluator
        .evaluate(&[-1.0, 1.0, -3.0], &problem, &options)
        .unwrap();
    assert_eq!(e.cut.vertex_region, vec![0, 1, 0]);
    // Edges: [0, 1] -> point 0 at t = 1/2, [1, 2] -> point 1 at t = 1/4 from vertex 1.
    // p on 1 -> 2 is point 1, q on 0 -> 1 is point 0; region 0 is left of q -> p.
    assert_eq!(e.cut.segments[0].points, [0, 1]);
    assert_eq!(e.triangles.faces[0], [1, 5, 4]);
    assert_eq!(e.triangles.region, vec![1, 0, 0]);
}

/// Prints how long the stages take. Run with
/// `cargo test --release -p ddg-variational-cuts --test identities -- --ignored --nocapture`.
#[test]
#[ignore = "prints timings on icospheres of level 5 to 7"]
fn print_timings() {
    use std::time::Instant;
    println!(
        "{:>5} {:>8} {:>6} | {:>9} {:>11} {:>11} {:>11} {:>9}",
        "level", "vertices", "cut", "new (ms)", "Dirichlet", "Hencky", "scaled", "cut_mesh"
    );
    for level in 5..=7 {
        let shape = sphere(level);
        let phi = shape.sample(tilted);
        let start = Instant::now();
        let mut evaluator = shape.evaluator();
        let setup = start.elapsed().as_secs_f64() * 1e3;
        let weight =
            ddg_variational_cuts::ScaleWeight::PerVertex(shape.sample(|p| 1.0 + 0.3 * p[0]));
        let mut times = Vec::new();
        let mut last = None;
        for energy in [
            Energy::Dirichlet,
            Energy::Hencky,
            Energy::ScaledDirichlet(weight),
        ] {
            let problem = CutProblem::new(energy, 0.1);
            // The first evaluation allocates the factor; time the second.
            evaluator
                .evaluate(&phi, &problem, &EvalOptions::default())
                .unwrap();
            let start = Instant::now();
            let e = evaluator
                .evaluate(&phi, &problem, &EvalOptions::default())
                .unwrap();
            times.push(start.elapsed().as_secs_f64() * 1e3);
            last = Some(e);
        }
        let e = last.unwrap();
        let start = Instant::now();
        let cut_mesh = e.cut_mesh().unwrap();
        let mesh_time = start.elapsed().as_secs_f64() * 1e3;
        assert_eq!(cut_mesh.mesh.n_faces(), e.triangles.n_faces());
        println!(
            "{:>5} {:>8} {:>6} | {:>9.1} {:>11.1} {:>11.1} {:>11.1} {:>9.1}",
            level,
            shape.positions.len(),
            e.cut.points.len(),
            setup,
            times[0],
            times[1],
            times[2],
            mesh_time
        );
    }
}
