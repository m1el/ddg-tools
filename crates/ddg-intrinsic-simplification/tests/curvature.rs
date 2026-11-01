//! Tests of the principal curvature estimator (SPEC §8, S22) and of the direction-to-angle
//! conversion.

use std::f64::consts::{PI, TAU};

use ddg_core::generators::{self, GeneratedMesh};
use ddg_core::{Geometry, TriMesh};
use ddg_intrinsic::IntrinsicTriangulation;
use ddg_intrinsic_simplification::curvature::{
    CurvatureError, CurvatureOptions, CurvatureStatus, direction_angle, direction_angles,
    principal_curvatures,
};

fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}
fn nearest(g: &GeneratedMesh, q: [f64; 3]) -> usize {
    let d2 = |p: &[f64; 3]| (0..3).map(|k| (p[k] - q[k]).powi(2)).sum::<f64>();
    (0..g.positions.len())
        .min_by(|&a, &b| d2(&g.positions[a]).total_cmp(&d2(&g.positions[b])))
        .unwrap()
}
fn mesh(g: &GeneratedMesh) -> TriMesh {
    g.topology().unwrap()
}

#[test]
fn sphere_curvature_and_orthogonal_directions() {
    let g = generators::icosphere(4).unwrap();
    let m = mesh(&g);
    let c = principal_curvatures(&m, &g.positions, &CurvatureOptions::default()).unwrap();
    assert_eq!(c.n_failed(), 0);
    for v in 0..m.n_vertices() {
        assert!((c.k_max[v] - 1.0).abs() < 0.02, "k_max {}", c.k_max[v]);
        assert!((c.k_min[v] - 1.0).abs() < 0.02, "k_min {}", c.k_min[v]);
        assert!(c.k_max[v] >= c.k_min[v]);
        assert!(dot(c.d_max[v], c.d_min[v]).abs() < 1e-9);
        let p = g.positions[v];
        assert!(dot(c.d_max[v], p).abs() < 0.05, "tangent");
        assert!((dot(c.d_max[v], c.d_max[v]) - 1.0).abs() < 1e-12);
    }
}

#[test]
fn sign_convention_can_be_flipped() {
    let g = generators::icosphere(3).unwrap();
    let m = mesh(&g);
    let opt = CurvatureOptions {
        convex_positive: false,
        ..Default::default()
    };
    let c = principal_curvatures(&m, &g.positions, &opt).unwrap();
    assert!((c.k_max[5] + 1.0).abs() < 0.05);
}

#[test]
fn cylinder_directions() {
    let r = 0.5;
    let g = generators::open_cylinder(96, 60, r, 4.0).unwrap();
    let m = mesh(&g);
    let opt = CurvatureOptions {
        scale: 2.0,
        ..Default::default()
    };
    let c = principal_curvatures(&m, &g.positions, &opt).unwrap();
    let mut worst = 0.0f64;
    let mut n = 0;
    for v in 0..m.n_vertices() {
        let p = g.positions[v];
        // Away from the two boundary loops (z in [-h/2, h/2] or [0, h]: use the middle half).
        let zs: Vec<f64> = g.positions.iter().map(|q| q[2]).collect();
        let (zmin, zmax) = (
            zs.iter().cloned().fold(f64::MAX, f64::min),
            zs.iter().cloned().fold(f64::MIN, f64::max),
        );
        if (p[2] - zmin) < 0.3 * (zmax - zmin) || (zmax - p[2]) < 0.3 * (zmax - zmin) {
            continue;
        }
        n += 1;
        assert_eq!(c.status[v], CurvatureStatus::Fitted);
        assert!(
            (c.k_max[v] - 1.0 / r).abs() < 0.03 / r,
            "k_max {} at {p:?}",
            c.k_max[v]
        );
        assert!(c.k_min[v].abs() < 0.03 / r, "k_min {}", c.k_min[v]);
        let rho = p[0].hypot(p[1]);
        let around = [-p[1] / rho, p[0] / rho, 0.0];
        let ang_max = dot(c.d_max[v], around).abs().min(1.0).acos();
        let ang_min = dot(c.d_min[v], [0.0, 0.0, 1.0]).abs().min(1.0).acos();
        worst = worst.max(ang_max).max(ang_min);
    }
    assert!(n > 500);
    assert!(worst.to_degrees() < 3.0, "worst {} deg", worst.to_degrees());
}

#[test]
fn saddle_at_origin() {
    // z = x^2 - y^2 on [-1,1]^2 from the unit-square grid: shift to the center.
    let mut g = generators::grid_disk(60, 60).unwrap();
    let h = 1.0 / 60.0;
    for p in &mut g.positions {
        p[0] = (p[0] - 0.5) * 2.0 * 0.5;
        p[1] = (p[1] - 0.5) * 2.0 * 0.5;
        // Domain [-0.5, 0.5]^2 (spacing h): z = x^2 - y^2.
        p[2] = p[0] * p[0] - p[1] * p[1];
    }
    let _ = h;
    let m = mesh(&g);
    let o = nearest(&g, [0.0; 3]);
    assert!(dot(g.positions[o], g.positions[o]) < 1e-12, "origin vertex exists");
    let c = principal_curvatures(&m, &g.positions, &CurvatureOptions::default()).unwrap();
    // Normal +z (counter-clockwise from +z): z = x^2 curves toward the normal along x ...
    // the surface bends toward +z along x, i.e. concave seen from outside: k = -2 along x,
    // convex k = +2 along y.
    assert!((c.k_max[o] - 2.0).abs() < 0.05, "k_max {}", c.k_max[o]);
    assert!((c.k_min[o] + 2.0).abs() < 0.05, "k_min {}", c.k_min[o]);
    assert!(c.d_max[o][1].abs() > 0.995, "d_max {:?}", c.d_max[o]);
    assert!(c.d_min[o][0].abs() > 0.995, "d_min {:?}", c.d_min[o]);
}

#[test]
fn torus_inner_and_outer() {
    let (big, small) = (2.0, 0.7);
    let g = generators::torus(128, 64, big, small).unwrap();
    let m = mesh(&g);
    let opt = CurvatureOptions {
        scale: 2.0,
        ..Default::default()
    };
    let c = principal_curvatures(&m, &g.positions, &opt).unwrap();
    // Outer equator: curvatures 1/small (meridian) and 1/(R + r) (parallel).
    let o = nearest(&g, [big + small, 0.0, 0.0]);
    assert!((c.k_max[o] - 1.0 / small).abs() < 0.05 / small, "outer {}", c.k_max[o]);
    assert!(
        (c.k_min[o] - 1.0 / (big + small)).abs() < 0.04,
        "outer min {}",
        c.k_min[o]
    );
    // Inner equator: saddle, 1/small and -1/(R - r).
    let i = nearest(&g, [big - small, 0.0, 0.0]);
    assert!((c.k_max[i] - 1.0 / small).abs() < 0.05 / small, "inner {}", c.k_max[i]);
    assert!(
        (c.k_min[i] + 1.0 / (big - small)).abs() < 0.08,
        "inner min {}",
        c.k_min[i]
    );
    // Top of the tube: z = r, curvatures 1/r and 0 at the top (cos phi = 0): k = cos/(R + r cos) = 0.
    let t = nearest(&g, [big, 0.0, small]);
    assert!(c.k_min[t].abs() < 0.06, "top min {}", c.k_min[t]);
    assert!((c.k_max[t] - 1.0 / small).abs() < 0.05 / small);
}

#[test]
fn convergence_under_refinement() {
    let mut errs = Vec::new();
    for level in 2..=4 {
        let g = generators::icosphere(level).unwrap();
        let m = mesh(&g);
        // A fixed number of edge lengths means a shrinking neighbourhood.
        let c = principal_curvatures(&m, &g.positions, &CurvatureOptions::default()).unwrap();
        let e = (0..m.n_vertices())
            .map(|v| (c.k_max[v] - 1.0).abs().max((c.k_min[v] - 1.0).abs()))
            .fold(0.0, f64::max);
        errs.push(e);
    }
    assert!(errs[1] < errs[0] && errs[2] < errs[1], "errors {errs:?}");
    assert!(errs[2] < 0.02);
}

#[test]
fn flat_region_is_not_degenerate() {
    let g = generators::grid_disk(10, 10).unwrap();
    let m = mesh(&g);
    let c = principal_curvatures(&m, &g.positions, &CurvatureOptions::default()).unwrap();
    for v in 0..m.n_vertices() {
        assert_eq!(c.status[v], CurvatureStatus::Fitted);
        assert!(c.k_max[v].abs() < 1e-9 && c.k_min[v].abs() < 1e-9);
        assert!((dot(c.d_max[v], c.d_max[v]) - 1.0).abs() < 1e-9);
        assert!(dot(c.d_max[v], c.d_min[v]).abs() < 1e-9);
        assert!(c.d_max[v][2].abs() < 1e-9);
    }
}

#[test]
fn deterministic_and_sampled_close_to_full() {
    let g = generators::icosphere(3).unwrap();
    let m = mesh(&g);
    let big = CurvatureOptions {
        scale: 6.0,
        ..Default::default()
    };
    let a = principal_curvatures(&m, &g.positions, &big).unwrap();
    let b = principal_curvatures(&m, &g.positions, &big).unwrap();
    assert_eq!(a, b);
    let full = principal_curvatures(
        &m,
        &g.positions,
        &CurvatureOptions {
            max_samples: None,
            ..big.clone()
        },
    )
    .unwrap();
    // A sphere is quadric-like at this scale: both are about 1 (slightly below: the large ball).
    for v in 0..m.n_vertices() {
        assert!((a.k_max[v] - full.k_max[v]).abs() < 0.1);
    }
}

#[test]
fn degenerate_inputs_do_not_panic() {
    // All vertices coincide: no normals.
    let g = generators::icosphere(1).unwrap();
    let m = mesh(&g);
    let pts = vec![[0.5, 0.5, 0.5]; m.n_vertices()];
    let c = principal_curvatures(&m, &pts, &CurvatureOptions::default()).unwrap();
    assert!(c.status.iter().all(|&s| s == CurvatureStatus::DegenerateNormal));
    assert!(c.k_max.iter().all(|k| *k == 0.0));
    assert!(c.d_max.iter().all(|d| d.iter().all(|t| t.is_finite())));
    // A single triangle: too few points for a fit.
    let tri = TriMesh::new(3, vec![[0, 1, 2]]).unwrap();
    let c = principal_curvatures(
        &tri,
        &[[0.0; 3], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]],
        &CurvatureOptions::default(),
    )
    .unwrap();
    assert!(c.status.iter().all(|&s| s == CurvatureStatus::Underdetermined));
    // Bad input.
    assert!(matches!(
        principal_curvatures(&m, &pts[..3], &CurvatureOptions::default()),
        Err(CurvatureError::PositionCount { .. })
    ));
    let mut nan = g.positions.clone();
    nan[2][1] = f64::NAN;
    assert!(matches!(
        principal_curvatures(&m, &nan, &CurvatureOptions::default()),
        Err(CurvatureError::NonFinitePosition { vertex: 2 })
    ));
    for bad in [0.0, -1.0, f64::NAN] {
        let o = CurvatureOptions {
            scale: bad,
            ..Default::default()
        };
        assert!(principal_curvatures(&m, &g.positions, &o).is_err());
    }
    let o = CurvatureOptions {
        max_samples: Some(3),
        ..Default::default()
    };
    assert!(principal_curvatures(&m, &g.positions, &o).is_err());
}

#[test]
fn boundary_vertices_are_fitted_one_sided() {
    let g = generators::open_cylinder(48, 12, 1.0, 2.0).unwrap();
    let m = mesh(&g);
    let c = principal_curvatures(&m, &g.positions, &CurvatureOptions::default()).unwrap();
    let mut nb = 0;
    for v in 0..m.n_vertices() {
        assert!(c.k_max[v].is_finite());
        if c.is_boundary[v] {
            nb += 1;
            // Cylinder is exactly representable by the one-sided quadric well enough.
            assert!((c.k_max[v] - 1.0).abs() < 0.1, "boundary k_max {}", c.k_max[v]);
        }
    }
    assert_eq!(nb, 2 * 48);
}

// ---- direction to angle ----

fn frames(g: &GeneratedMesh) -> (TriMesh, IntrinsicTriangulation) {
    let m = mesh(g);
    let geo = Geometry::from_positions(&m, &g.positions).unwrap();
    let t = IntrinsicTriangulation::new(&m, &geo).unwrap();
    (m, t)
}

#[test]
fn angle_of_edge_directions_is_the_signpost() {
    let g = generators::icosphere(2).unwrap();
    let (m, t) = frames(&g);
    for h in (0..m.n_halfedges()).step_by(7) {
        let v = m.from_vertex(h);
        let w = m.to_vertex(h);
        let d = [0, 1, 2].map(|k| g.positions[w][k] - g.positions[v][k]);
        let a = direction_angle(&m, &g.positions, &t, v, d).unwrap();
        let s = t.signpost(h);
        let diff = ((a - s + PI).rem_euclid(TAU) - PI).abs();
        assert!(diff < 1e-9, "h {h}: {a} vs {s}");
    }
}

#[test]
fn angle_on_flat_mesh_is_euclidean_and_boundary_in_range() {
    let g = generators::grid_disk(8, 8).unwrap();
    let (m, t) = frames(&g);
    // interior vertex: direction at angle phi from +x: the frame's zero is the
    // outgoing halfedge, so compare differences between two directions.
    let v = 4 * 9 + 4;
    let dir = |phi: f64| [phi.cos(), phi.sin(), 0.0];
    let a0 = direction_angle(&m, &g.positions, &t, v, dir(0.3)).unwrap();
    let a1 = direction_angle(&m, &g.positions, &t, v, dir(1.3)).unwrap();
    let diff = (a1 - a0).rem_euclid(TAU);
    assert!((diff - 1.0).abs() < 1e-9, "{diff}");
    // Boundary vertex: all directions, even out of the surface, give an angle in [0, pi].
    let b = 0usize;
    for k in 0..16 {
        let phi = k as f64 * TAU / 16.0;
        let a = direction_angle(&m, &g.positions, &t, b, dir(phi)).unwrap();
        assert!((0.0..=PI + 1e-12).contains(&a), "angle {a} for {phi}");
    }
    // Vector with a normal component projects to the same angle.
    let a = direction_angle(&m, &g.positions, &t, v, [dir(0.3)[0], dir(0.3)[1], 5.0]).unwrap();
    assert!((a - a0).abs() < 1e-9);
}

#[test]
fn angles_of_curvature_directions_and_errors() {
    let g = generators::icosphere(2).unwrap();
    let m = mesh(&g);
    let c = principal_curvatures(&m, &g.positions, &CurvatureOptions::default()).unwrap();
    let a = direction_angles(&m, &g.positions, &c.d_max).unwrap();
    assert!(a.iter().all(|x| (0.0..TAU).contains(x)));
    assert!(direction_angles(&m, &g.positions, &c.d_max[..3]).is_err());
    let mut zero = c.d_max.clone();
    zero[7] = [0.0; 3];
    assert!(matches!(
        direction_angles(&m, &g.positions, &zero),
        Err(CurvatureError::BadDirection { vertex: 7 })
    ));
}
