//! Tests of the applications on the coarse mesh (SPEC §10.1, WP6): operators, Poisson, geodesic
//! distance (exact and heat), all pairs, Voronoi, exponential map. Coarse meshes come from
//! `greedy_coarsening` (not the ICE order: numbers will change with `simplify`).

mod common;

use common::{Shape, load_data};
use ddg_core::generators::{grid_disk, icosphere};
use ddg_core::operators;
use ddg_core::sparse::ThreadPool;
use ddg_core::{Geometry, TriMesh};
use ddg_intrinsic_apps::exact::{ExactOptions, ExactSource, exact_distance};
use ddg_intrinsic_simplification::applications::{
    CoarseApplications, CoarseHeat, DistanceError, ExpStop, PoissonSolver, greedy_coarsening,
};
use ddg_intrinsic_simplification::{CoarseMesh, CoarsePoint};

fn shape_of(g: ddg_core::generators::GeneratedMesh, name: &str) -> Shape {
    Shape::new(name, g.positions, g.faces)
}

fn coarsen(shape: &Shape, fraction: f64) -> (CoarseMesh, CoarseApplications) {
    let mut mesh = shape.coarse(&[]);
    let target = ((shape.positions.len() as f64) * fraction).ceil() as usize;
    greedy_coarsening(&mut mesh, target).unwrap();
    mesh.certify(false).unwrap();
    let app = CoarseApplications::new(&mesh).unwrap();
    (mesh, app)
}

fn fine_geometry(shape: &Shape) -> (TriMesh, Geometry) {
    let mesh = TriMesh::new(shape.positions.len(), shape.faces.clone()).unwrap();
    let geo = Geometry::from_positions(&mesh, &shape.positions).unwrap();
    (mesh, geo)
}

/// The position in the input's space of a coarse point, through the input positions of the
/// corners (meaningful for a flat input, where the coarse metric is the planar one).
fn flat_position(shape: &Shape, app: &CoarseApplications, pt: &CoarsePoint) -> [f64; 3] {
    let out = app.output();
    match *pt {
        CoarsePoint::Vertex(v) => shape.positions[out.vertices[v as usize] as usize],
        CoarsePoint::Face { face, bary } => {
            let c = out.faces[face as usize];
            [0, 1, 2].map(|k| {
                (0..3)
                    .map(|i| bary[i] * shape.positions[out.vertices[c[i] as usize] as usize][k])
                    .sum()
            })
        }
    }
}

#[test]
fn operators_match_the_mesh_and_have_the_right_structure() {
    let shape = shape_of(icosphere(3).unwrap(), "ico3");
    let (mesh, app) = coarsen(&shape, 0.2);
    let l = app.laplacian();
    let l2 = mesh.cotan_laplacian();
    assert_eq!(l.nrows(), mesh.n_vertices());
    let mut worst: f64 = 0.0;
    for (i, j, v) in l2.triplets() {
        worst = worst.max((l.get(i as usize, j as usize) - v).abs());
    }
    assert!(worst < 1e-12, "{worst}");
    let ones = vec![1.0; l.nrows()];
    let mut y = vec![0.0; l.nrows()];
    l.mul_vec(&ones, &mut y).unwrap();
    assert!(y.iter().all(|v| v.abs() < 1e-10));
    assert!(l.symmetry_defect().unwrap() < 1e-12);
    let m = app.mass();
    let m2 = mesh.lumped_mass();
    assert!(m.iter().zip(&m2).all(|(a, b)| (a - b).abs() < 1e-14));
    let area: f64 = app.output().face_areas().unwrap().iter().sum();
    assert!((m.iter().sum::<f64>() - area).abs() < 1e-12 * area);
}

#[test]
fn poisson_neumann_and_dirichlet_on_a_sphere() {
    let shape = shape_of(icosphere(4).unwrap(), "ico4");
    let (fmesh, geo) = fine_geometry(&shape);
    let fl = operators::cotan_laplacian(&fmesh, &geo).unwrap();
    let fm = operators::lumped_mass(&fmesh, &geo).unwrap();
    let f: Vec<f64> = shape.positions.iter().map(|p| p[2] + p[0] * p[1]).collect();
    let fine = PoissonSolver::new(&fl, &fm, &[]).unwrap();
    let (xf, rep) = fine.solve(&f, &[]).unwrap();
    assert!(rep.residual_max < 1e-10, "{rep:?}");
    for fraction in [0.3, 0.1] {
        let (_mesh, app) = coarsen(&shape, fraction);
        let fbar = app.sample(&f).unwrap();
        let (_xbar, x, report) = app.poisson(&fbar, &[], &[]).unwrap();
        assert!(report.residual_max < 1e-9, "{report:?}");
        let err = DistanceError::of(&x, &xf).unwrap();
        eprintln!(
            "poisson sphere {:.0}% ({} verts): rel L2 {:.4}",
            100.0 * fraction,
            app.output().n_vertices(),
            err.rel_l2
        );
        assert!(err.rel_l2 < 0.15, "{err:?}");
        let b: Vec<f64> = fm.iter().zip(&f).map(|(m, f)| m * f).collect();
        let load = app.restrict_load(&b).unwrap();
        let solver = app.poisson_solver(&[]).unwrap();
        let (xg, _) = solver.solve_load(&load, &[]).unwrap();
        let err_g = DistanceError::of(&app.prolong(&xg).unwrap(), &xf).unwrap();
        eprintln!("  galerkin load: rel L2 {:.4}", err_g.rel_l2);
        assert!(err_g.rel_l2 < 0.3, "{err_g:?}");
    }
    // Dirichlet at a few coarse vertices (their fine ids are the fixed fine vertices).
    let (_mesh, app) = coarsen(&shape, 0.1);
    let fixed_coarse: Vec<u32> = (0..app.output().n_vertices() as u32).step_by(23).collect();
    let fixed_fine: Vec<u32> = fixed_coarse
        .iter()
        .map(|&c| app.output().vertices[c as usize])
        .collect();
    let values: Vec<f64> = fixed_fine.iter().map(|&v| f[v as usize]).collect();
    let fine = PoissonSolver::new(&fl, &fm, &fixed_fine).unwrap();
    let (xf, _) = fine.solve(&f, &values).unwrap();
    let fbar = app.sample(&f).unwrap();
    let (_, x, rep) = app.poisson(&fbar, &fixed_coarse, &values).unwrap();
    assert!(rep.residual_max < 1e-9);
    let err = DistanceError::of(&x, &xf).unwrap();
    eprintln!(
        "poisson dirichlet: rel L2 {:.4}, max {:.4}",
        err.rel_l2, err.max_over_max
    );
    assert!(err.rel_l2 < 0.3, "{err:?}");
}

#[test]
fn poisson_rejects_bad_input() {
    let shape = shape_of(icosphere(2).unwrap(), "ico2");
    let (_m, app) = coarsen(&shape, 0.5);
    assert!(app.poisson_solver(&[0, 0]).is_err());
    assert!(app.poisson_solver(&[u32::MAX]).is_err());
    let s = app.poisson_solver(&[0]).unwrap();
    assert!(s.solve(&[0.0; 3], &[0.0]).is_err());
    assert!(s.solve(&vec![0.0; app.output().n_vertices()], &[]).is_err());
}

fn exact_fine(shape: &Shape, source: usize) -> Vec<f64> {
    let (m, g) = fine_geometry(shape);
    exact_distance(&m, &g, &[source], &ExactOptions::default())
        .unwrap()
        .distances
}

#[test]
fn geodesic_distance_exact_and_heat_on_a_sphere() {
    let shape = shape_of(icosphere(4).unwrap(), "ico4");
    let src = 17;
    let reference = exact_fine(&shape, src);
    for fraction in [0.5, 0.1, 0.03] {
        let (_m, app) = coarsen(&shape, fraction);
        let (_, d) = app.geodesic_exact(src, &ExactOptions::default()).unwrap();
        let e = DistanceError::of(&d, &reference).unwrap();
        let (_, dh) = app.geodesic_heat(src, 1.0).unwrap();
        let eh = DistanceError::of(&dh, &reference).unwrap();
        eprintln!(
            "sphere {:.0}% ({} verts): exact rel L1 {:.4}, heat rel L1 {:.4}",
            100.0 * fraction,
            app.output().n_vertices(),
            e.rel_l1,
            eh.rel_l1
        );
        assert!(e.rel_l1 < 0.03, "exact {e:?}");
        assert!(eh.rel_l1 < 0.15, "heat {eh:?}");
    }
}

#[test]
fn geodesic_from_a_removed_vertex_uses_a_face_point_source() {
    let shape = shape_of(icosphere(4).unwrap(), "ico4");
    let (_m, app) = coarsen(&shape, 0.05);
    let removed = (0..shape.positions.len())
        .find(|&v| app.output().coarse_of_fine[v].is_none())
        .unwrap();
    assert!(matches!(
        app.output().points[removed],
        CoarsePoint::Face { .. }
    ));
    let reference = exact_fine(&shape, removed);
    let (dbar, d) = app
        .geodesic_exact(removed, &ExactOptions::default())
        .unwrap();
    assert!(dbar.iter().all(|x| x.is_finite() && *x >= 0.0));
    let e = DistanceError::of(&d, &reference).unwrap();
    eprintln!("face-point source: rel L1 {:.4}", e.rel_l1);
    assert!(e.rel_l1 < 0.05, "{e:?}");
    let (_, dh) = app.geodesic_heat(removed, 1.0).unwrap();
    let eh = DistanceError::of(&dh, &reference).unwrap();
    eprintln!("face-point source heat: rel L1 {:.4}", eh.rel_l1);
    assert!(eh.rel_l1 < 0.2, "{eh:?}");
}

#[test]
fn all_pairs_is_deterministic_and_accurate() {
    let shape = shape_of(icosphere(3).unwrap(), "ico3");
    let (_m, app) = coarsen(&shape, 0.25);
    let p1 = ThreadPool::new(1).unwrap();
    let p4 = ThreadPool::new(4).unwrap();
    let d1 = app.all_pairs(&p1, &ExactOptions::default()).unwrap();
    let d4 = app.all_pairs(&p4, &ExactOptions::default()).unwrap();
    assert_eq!(d1, d4, "all pairs depends on the thread count");
    let n = d1.n_rows;
    let mut asym: f64 = 0.0;
    for i in 0..n {
        assert_eq!(d1.get(i, i), 0.0);
        for j in 0..n {
            asym = asym.max((d1.get(i, j) - d1.get(j, i)).abs());
        }
    }
    assert!(asym < 1e-9, "{asym}");
    let dense = app.prolongate_dense(&d1, &p4).unwrap();
    let rows = [0usize, 5, 100, shape.positions.len() - 1];
    let part = app.prolongate_rows(&d1, &rows).unwrap();
    for (k, &r) in rows.iter().enumerate() {
        assert_eq!(part.row(k), dense.row(r));
    }
    let mut err_sum = 0.0;
    let mut ref_sum = 0.0;
    for &r in &rows {
        let reference = exact_fine(&shape, r);
        for (a, b) in dense.row(r).iter().zip(&reference) {
            err_sum += (a - b).abs();
            ref_sum += b;
        }
    }
    eprintln!(
        "all pairs {} coarse verts: rel L1 {:.4}",
        n,
        err_sum / ref_sum
    );
    assert!(err_sum / ref_sum < 0.05);
}

#[test]
fn voronoi_agrees_with_the_fine_diagram() {
    let shape = shape_of(icosphere(4).unwrap(), "ico4");
    let (_m, app) = coarsen(&shape, 0.1);
    let sources: Vec<usize> = (0..12).map(|k| 37 + 199 * k).collect();
    let pool = ThreadPool::new(4).unwrap();
    let v = app
        .voronoi(&sources, &pool, &ExactOptions::default())
        .unwrap();
    let fine: Vec<Vec<f64>> = sources.iter().map(|&s| exact_fine(&shape, s)).collect();
    let n = shape.positions.len();
    let mut wrong = 0;
    for i in 0..n {
        let mut best = 0;
        for k in 1..fine.len() {
            if fine[k][i] < fine[best][i] {
                best = k;
            }
        }
        if best as u32 != v.fine_labels[i] {
            wrong += 1;
        }
    }
    let frac = wrong as f64 / n as f64;
    eprintln!("voronoi: {:.2}% misclassified", 100.0 * frac);
    assert!(frac < 0.05, "{frac}");
}

#[test]
fn heat_distance_on_a_flat_grid_is_close_to_euclidean() {
    let shape = shape_of(grid_disk(24, 24).unwrap(), "grid");
    let mesh = shape.coarse(&[]);
    let app = CoarseApplications::new(&mesh).unwrap();
    let heat = CoarseHeat::new(&app, 1.0).unwrap();
    let c = 12 * 25 + 12;
    let d = heat.distance(&[ExactSource::Vertex(c)]).unwrap();
    let mut worst: f64 = 0.0;
    for (v, p) in shape.positions.iter().enumerate() {
        let e = ((p[0] - shape.positions[c][0]).powi(2) + (p[1] - shape.positions[c][1]).powi(2))
            .sqrt();
        worst = worst.max((d[v] - e).abs());
    }
    eprintln!("flat heat: max error {worst:.4} (domain size 1)");
    assert!(worst < 0.06, "{worst}");
    assert!(heat.distance(&[]).is_err());
    assert!(CoarseHeat::new(&app, 0.0).is_err());
}

#[test]
fn exp_map_on_a_flat_grid_is_the_straight_line() {
    let shape = shape_of(grid_disk(10, 10).unwrap(), "grid");
    let mesh = shape.coarse(&[]);
    let app = CoarseApplications::new(&mesh).unwrap();
    let out = app.output();
    let f = 85u32;
    let bary = [0.2, 0.3, 0.5];
    let zero = app.exp_map_face(f, bary, [0.0, 0.0]).unwrap();
    assert_eq!(zero.faces, vec![f]);
    let start = flat_position(&shape, &app, &CoarsePoint::Face { face: f, bary });
    // The face frame is a rotation of the plane: from the first edge.
    let c = out.faces[f as usize];
    let pa = shape.positions[out.vertices[c[0] as usize] as usize];
    let pb = shape.positions[out.vertices[c[1] as usize] as usize];
    let ang = (pb[1] - pa[1]).atan2(pb[0] - pa[0]);
    for (dx, dy) in [
        (0.15, 0.05),
        (-0.2, 0.1),
        (0.05, -0.12),
        (0.0, 0.3),
        (0.3, 0.3),
    ] {
        let r = app.exp_map_face(f, bary, [dx, dy]).unwrap();
        assert_eq!(r.stop, ExpStop::Completed);
        let (c, s) = (ang.cos(), ang.sin());
        let expected = [start[0] + c * dx - s * dy, start[1] + s * dx + c * dy];
        let got = flat_position(&shape, &app, &r.point);
        let err = ((got[0] - expected[0]).powi(2) + (got[1] - expected[1]).powi(2)).sqrt();
        assert!(err < 1e-10, "vector ({dx},{dy}): {got:?} vs {expected:?}");
        assert!((r.traced - (dx * dx + dy * dy).sqrt()).abs() < 1e-12);
        if let CoarsePoint::Face { face, bary } = r.point {
            let len = r.traced;
            let back = app
                .exp_map_face(face, bary, [-r.direction[0] * len, -r.direction[1] * len])
                .unwrap();
            let pbk = flat_position(&shape, &app, &back.point);
            let e = ((pbk[0] - start[0]).powi(2) + (pbk[1] - start[1]).powi(2)).sqrt();
            assert!(e < 1e-9, "round trip error {e}");
        }
    }
    let r = app.exp_map_face(f, bary, [5.0, 0.0]).unwrap();
    assert_eq!(r.stop, ExpStop::Boundary);
    assert!(r.traced < 2.0);
    assert!(app.exp_map_face(9999, bary, [0.1, 0.1]).is_err());
    assert!(app.exp_map_face(f, [0.5, 0.5, 0.5], [0.1, 0.1]).is_err());
    assert!(app.exp_map_face(f, bary, [f64::NAN, 0.0]).is_err());
    assert!(app.exp_map_vertex(9999, 0.0, 1.0).is_err());
    assert!(app.exp_map_vertex(60, 0.0, -1.0).is_err());
}

#[test]
fn exp_map_from_a_vertex_passes_straight_through_flat_vertices() {
    let shape = shape_of(grid_disk(10, 10).unwrap(), "grid");
    let mesh = shape.coarse(&[]);
    let app = CoarseApplications::new(&mesh).unwrap();
    let out = app.output();
    let v = (5 * 11 + 5) as u32;
    let pv = shape.positions[v as usize];
    let mut tested = 0;
    for k in 0..out.faces.len() {
        for slot in 0..3 {
            if out.faces[k][slot] != v {
                continue;
            }
            let w = out.faces[k][(slot + 1) % 3];
            let pw = shape.positions[w as usize];
            let len = out.face_lengths[k][slot];
            let angle = out.phi[3 * k + slot];
            let r = app.exp_map_vertex(v, angle, 2.0 * len).unwrap();
            let target = [2.0 * pw[0] - pv[0], 2.0 * pw[1] - pv[1]];
            let got = flat_position(&shape, &app, &r.point);
            let e = ((got[0] - target[0]).powi(2) + (got[1] - target[1]).powi(2)).sqrt();
            assert!(e < 1e-9, "edge {v}->{w}: got {got:?} want {target:?}");
            assert!(!r.vertices_hit.is_empty(), "the ray passes through {w}");
            tested += 1;
        }
    }
    assert!(tested >= 6);
}

#[test]
fn exp_map_on_a_coarse_delta_complex_is_reversible() {
    let shape = shape_of(icosphere(4).unwrap(), "ico4");
    let (_m, app) = coarsen(&shape, 0.02);
    let out = app.output();
    let mut checked = 0;
    for f in (0..out.n_faces() as u32).step_by(3) {
        for (bary, vec) in [
            ([0.3, 0.3, 0.4], [0.7, 0.4]),
            ([0.1, 0.2, 0.7], [-1.1, 0.5]),
        ] {
            let r = app.exp_map_face(f, bary, vec).unwrap();
            assert_eq!(r.stop, ExpStop::Completed);
            let (face, b) = match r.point {
                CoarsePoint::Face { face, bary } => (face, bary),
                CoarsePoint::Vertex(_) => continue,
            };
            let len = r.traced;
            let back = app
                .exp_map_face(face, b, [-r.direction[0] * len, -r.direction[1] * len])
                .unwrap();
            let p0 = app.face_position(f, bary).unwrap();
            if let CoarsePoint::Face { face: f2, bary: b2 } = back.point
                && f2 == f
            {
                let p = app.face_position(f2, b2).unwrap();
                let e = ((p[0] - p0[0]).powi(2) + (p[1] - p0[1]).powi(2)).sqrt();
                assert!(e < 1e-8, "face {f}: round trip error {e}");
                checked += 1;
            }
        }
    }
    assert!(checked > 5, "only {checked} round trips checked");
}

#[test]
fn spot_smoke_test() {
    // The full numbers are printed by `examples/applications.rs`.
    let Some(shape) = load_data("models/spot/spot/spot_triangulated.obj") else {
        return;
    };
    let (_m, app) = coarsen(&shape, 0.1);
    let reference = exact_fine(&shape, 100);
    let (_, d) = app.geodesic_exact(100, &ExactOptions::default()).unwrap();
    let e = DistanceError::of(&d, &reference).unwrap();
    eprintln!("spot 10%: exact rel L1 {:.4}", e.rel_l1);
    assert!(e.rel_l1 < 0.05, "{e:?}");
}
