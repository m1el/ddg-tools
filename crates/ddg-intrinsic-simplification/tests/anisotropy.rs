#![allow(clippy::needless_range_loop)]
//! Tests of the anisotropic coarsening (paper Sec. 8.6.2, Fig. 27; SPEC §8).

mod common;

use ddg_core::generators;
use ddg_core::{Geometry, TriMesh};
use ddg_intrinsic_simplification::anisotropy::{
    AnisotropyError, AnisotropyOptions, PrincipalDirection, anisotropic_coarse_mesh,
    anisotropic_lengths, coarse_face_alignment, coarsen_by_area, edge_alignment, length_factor,
    max_valid_tau, principal_direction_field, violating_faces,
};
use ddg_intrinsic_simplification::curvature::{CurvatureOptions, direction_angles};
use ddg_intrinsic_simplification::{CoarseMesh, CoarseMeshOptions};

fn setup(positions: &[[f64; 3]], faces: &[[u32; 3]]) -> (TriMesh, Vec<f64>) {
    let mesh = TriMesh::new(positions.len(), faces.to_vec()).unwrap();
    let lengths = Geometry::from_positions(&mesh, positions)
        .unwrap()
        .edge_lengths()
        .to_vec();
    (mesh, lengths)
}

#[test]
fn tau_zero_is_bitwise_and_factor_is_in_range() {
    for g in [
        generators::icosphere(3).unwrap(),
        generators::torus(12, 6, 1.0, 0.4).unwrap(),
        generators::open_cylinder(24, 8, 0.5, 2.0).unwrap(),
    ] {
        let (mesh, lengths) = setup(&g.positions, &g.faces);
        let f = principal_direction_field(
            &g.positions,
            &g.faces,
            PrincipalDirection::Max,
            &CurvatureOptions::default(),
        )
        .unwrap();
        let (l0, r0) = anisotropic_lengths(
            &mesh,
            &lengths,
            &f.angles,
            &AnisotropyOptions::with_tau(0.0),
        )
        .unwrap();
        assert!(
            l0.iter()
                .zip(&lengths)
                .all(|(a, b)| a.to_bits() == b.to_bits())
        );
        assert!(!r0.clamped && r0.min_factor == 1.0 && r0.max_factor == 1.0);
        let a = edge_alignment(&mesh, &lengths, &f.angles).unwrap();
        assert!(a.iter().all(|&x| (0.0..=1.0).contains(&x)));
        for tau in [0.1, 0.3, 0.5] {
            let (l, r) = anisotropic_lengths(
                &mesh,
                &lengths,
                &f.angles,
                &AnisotropyOptions::with_tau(tau),
            )
            .unwrap();
            for e in 0..l.len() {
                let q = l[e] / lengths[e];
                assert!(q >= 1.0 - tau - 1e-14 && q <= 1.0 + 1e-14, "factor {q}");
            }
            assert!(r.min_factor >= 1.0 - tau - 1e-14 && r.max_factor <= 1.0 + 1e-14);
        }
    }
}

/// On a flat mesh the alignment of an edge is the analytic `cos²` of its angle with `u`.
#[test]
fn planar_grid_matches_the_formula() {
    let s = common::planar_grid(8, 0.4, 3);
    let (mesh, lengths) = setup(&s.positions, &s.faces);
    for dir in [[1.0, 0.0, 0.0], [0.6, 0.8, 0.0], [-0.6, -0.8, 0.0]] {
        let angles = direction_angles(&mesh, &s.positions, &vec![dir; s.positions.len()]).unwrap();
        let a = edge_alignment(&mesh, &lengths, &angles).unwrap();
        for e in 0..mesh.n_edges() {
            let [i, j] = mesh.edge_vertices(e);
            // At a boundary vertex a direction pointing out of the surface is snapped to the
            // nearest boundary edge (direction_angles), so only interior ends are analytic.
            if mesh.is_boundary_vertex(i) || mesh.is_boundary_vertex(j) {
                continue;
            }
            let d = [0, 1].map(|k| s.positions[j][k] - s.positions[i][k]);
            let l = (d[0] * d[0] + d[1] * d[1]).sqrt();
            let c = (d[0] * dir[0] + d[1] * dir[1]) / l;
            assert!(
                (a[e] - c * c).abs() < 1e-9,
                "edge {e}: {} vs {}",
                a[e],
                c * c
            );
        }
        // Parallel edges keep their length, orthogonal ones shrink by 1 - tau.
        assert_eq!(length_factor(0.4, 1.0), 1.0);
        assert!((length_factor(0.4, 0.0) - 0.6).abs() < 1e-15);
    }
    // u -> -u gives the same factors.
    let a1 = direction_angles(
        &mesh,
        &s.positions,
        &vec![[0.6, 0.8, 0.0]; s.positions.len()],
    )
    .unwrap();
    let a2 = direction_angles(
        &mesh,
        &s.positions,
        &vec![[-0.6, -0.8, 0.0]; s.positions.len()],
    )
    .unwrap();
    let (x, y) = (
        edge_alignment(&mesh, &lengths, &a1).unwrap(),
        edge_alignment(&mesh, &lengths, &a2).unwrap(),
    );
    for e in 0..mesh.n_edges() {
        let [i, j] = mesh.edge_vertices(e);
        if !mesh.is_boundary_vertex(i) && !mesh.is_boundary_vertex(j) {
            assert!((x[e] - y[e]).abs() < 1e-9, "edge {e}");
        }
    }
}

#[test]
fn invalid_tau_is_an_error_or_clamped() {
    let g = generators::icosphere(3).unwrap();
    let (mesh, lengths) = setup(&g.positions, &g.faces);
    let f = principal_direction_field(
        &g.positions,
        &g.faces,
        PrincipalDirection::Max,
        &CurvatureOptions::default(),
    )
    .unwrap();
    let a = edge_alignment(&mesh, &lengths, &f.angles).unwrap();
    let tau_max = max_valid_tau(&mesh, &lengths, &a);
    // Bisection on the validity agrees with the exact root.
    let (mut lo, mut hi) = (0.0f64, 1.0f64);
    if violating_faces(&mesh, &lengths, &a, hi).is_empty() {
        lo = 1.0;
    }
    while hi - lo > 1e-12 && lo < 1.0 {
        let mid = 0.5 * (lo + hi);
        if violating_faces(&mesh, &lengths, &a, mid).is_empty() {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    assert!(tau_max.is_finite() && tau_max < 1.0, "tau_max {tau_max}");
    assert!(
        (lo - tau_max).abs() < 1e-9,
        "bisection {lo} vs exact {tau_max}"
    );

    let too_strong = AnisotropyOptions::with_tau(1.0);
    match anisotropic_lengths(&mesh, &lengths, &f.angles, &too_strong) {
        Err(AnisotropyError::TriangleInequality {
            faces, tau_max: t, ..
        }) => {
            assert!(!faces.is_empty() && (t - tau_max).abs() < 1e-15);
            assert_eq!(faces, violating_faces(&mesh, &lengths, &a, 1.0));
        }
        other => panic!("expected the triangle inequality error, got {other:?}"),
    }
    let msg = anisotropic_lengths(&mesh, &lengths, &f.angles, &too_strong)
        .unwrap_err()
        .to_string();
    assert!(
        msg.contains("triangle inequality") && msg.contains("face"),
        "{msg}"
    );

    let clamp = AnisotropyOptions {
        tau: 1.0,
        clamp_tau: true,
        ..Default::default()
    };
    let (scaled, report) = anisotropic_lengths(&mesh, &lengths, &f.angles, &clamp).unwrap();
    assert!(report.clamped && report.tau_used < tau_max && report.tau_used >= 0.99999 * tau_max);
    assert!(report.n_violating_faces_requested > 0);
    // The result is accepted by the strict validation of the geometry.
    Geometry::from_edge_lengths(&mesh, scaled).unwrap();
    // A request below the limit is not changed.
    let ok = AnisotropyOptions {
        tau: 0.5 * tau_max,
        clamp_tau: true,
        ..Default::default()
    };
    let (_, r) = anisotropic_lengths(&mesh, &lengths, &f.angles, &ok).unwrap();
    assert!(!r.clamped && r.tau_used == 0.5 * tau_max);
}

#[test]
fn bad_input_is_an_error_not_a_panic() {
    let g = generators::icosphere(1).unwrap();
    let (mesh, lengths) = setup(&g.positions, &g.faces);
    let angles = vec![0.0; mesh.n_vertices()];
    for tau in [-0.1, 1.5, f64::NAN, f64::INFINITY] {
        assert!(matches!(
            anisotropic_lengths(&mesh, &lengths, &angles, &AnisotropyOptions::with_tau(tau)),
            Err(AnisotropyError::BadOption { .. })
        ));
    }
    let o = AnisotropyOptions::with_tau(0.2);
    assert!(matches!(
        anisotropic_lengths(&mesh, &lengths[1..], &angles, &o),
        Err(AnisotropyError::BadInput { .. })
    ));
    assert!(matches!(
        anisotropic_lengths(&mesh, &lengths, &angles[1..], &o),
        Err(AnisotropyError::BadInput { .. })
    ));
    let mut a = angles.clone();
    a[3] = f64::NAN;
    assert!(anisotropic_lengths(&mesh, &lengths, &a, &o).is_err());
    let mut l = lengths.clone();
    l[2] = -1.0;
    assert!(anisotropic_lengths(&mesh, &l, &angles, &o).is_err());
}

/// The coarse mesh of the scaled metric tracks the input points; with `tau = 0` it is the
/// plain coarse mesh.
#[test]
fn coarse_mesh_tracking_is_intact() {
    let g = generators::torus(20, 10, 1.0, 0.4).unwrap();
    let s = common::Shape::new("torus", g.positions.clone(), g.faces.clone());
    let (mesh, lengths) = setup(&g.positions, &g.faces);
    let f = principal_direction_field(
        &g.positions,
        &g.faces,
        PrincipalDirection::Max,
        &CurvatureOptions::default(),
    )
    .unwrap();
    let points = s.random_points(5);
    // tau = 0 : identical to the plain construction.
    let (c0, _) = anisotropic_coarse_mesh(
        &mesh,
        &lengths,
        &f.angles,
        &points,
        CoarseMeshOptions::default(),
        &AnisotropyOptions::with_tau(0.0),
    )
    .unwrap();
    let plain = CoarseMesh::from_lengths(
        &mesh,
        lengths.clone(),
        &points,
        CoarseMeshOptions::default(),
    )
    .unwrap();
    assert!(c0.bitwise_eq(&plain));

    let opts = AnisotropyOptions {
        tau: 0.5,
        clamp_tau: true,
        ..Default::default()
    };
    let (mut c, report) = anisotropic_coarse_mesh(
        &mesh,
        &lengths,
        &f.angles,
        &points,
        CoarseMeshOptions::default(),
        &opts,
    )
    .unwrap();
    c.certify(false).unwrap();
    // the metric being coarsened is the scaled one
    let (scaled, _) = anisotropic_lengths(&mesh, &lengths, &f.angles, &opts).unwrap();
    let area = Geometry::from_edge_lengths(&mesh, scaled)
        .unwrap()
        .total_area();
    assert!((c.total_area() - area).abs() < 1e-9 * area);
    let left = coarsen_by_area(&mut c, 40).unwrap();
    assert!(left < 60, "{left} vertices left");
    c.certify(false).unwrap();
    // every input point and fine vertex is somewhere on the coarse mesh
    let out = c.compact().unwrap();
    assert_eq!(out.points.len(), g.positions.len() + points.len());
    assert!(out.face_of_points().iter().all(|f| f.is_some()));
    assert!(!report.clamped);
}

/// Fig. 27 on a cylinder: with `u` the maximum principal direction (around the axis) the
/// scaled metric is short along the axis, so coarse triangles become long along the axis, i.e.
/// across `u` and along the minimum direction; and the other way round for `Min`.
#[test]
fn cylinder_triangles_elongate_across_the_field() {
    let g = generators::open_cylinder(96, 60, 0.5, 4.0).unwrap();
    let (mesh, lengths) = setup(&g.positions, &g.faces);
    let curv = ddg_intrinsic_simplification::curvature::principal_curvatures(
        &mesh,
        &g.positions,
        &CurvatureOptions::default(),
    )
    .unwrap();
    let measure = |which: PrincipalDirection, tau: f64| {
        let f =
            principal_direction_field(&g.positions, &g.faces, which, &CurvatureOptions::default())
                .unwrap();
        let (mut c, report) = anisotropic_coarse_mesh(
            &mesh,
            &lengths,
            &f.angles,
            &[],
            CoarseMeshOptions::default(),
            &AnisotropyOptions {
                tau,
                clamp_tau: true,
                ..Default::default()
            },
        )
        .unwrap();
        coarsen_by_area(&mut c, 300).unwrap();
        c.certify(false).unwrap();
        let out = c.compact().unwrap();
        let vs_u = coarse_face_alignment(&out, &g.positions, &f.directions, &curv.normal, 4);
        let other = match which {
            PrincipalDirection::Max => &curv.d_min,
            PrincipalDirection::Min => &curv.d_max,
        };
        let vs_other = coarse_face_alignment(&out, &g.positions, other, &curv.normal, 4);
        eprintln!(
            "{which:?} tau {tau} (used {:.3}): {} faces, angle to u {:.1} deg, to the other \
             direction {:.1} deg, aspect {:.2}",
            report.tau_used,
            vs_u.n_faces,
            vs_u.mean_angle_deg,
            vs_other.mean_angle_deg,
            vs_u.mean_aspect
        );
        (vs_u, vs_other)
    };
    let (base_u, _) = measure(PrincipalDirection::Max, 0.0);
    for which in [PrincipalDirection::Max, PrincipalDirection::Min] {
        let (u, other) = measure(which, 0.6);
        assert!(u.n_faces > 100);
        assert!(
            u.mean_angle_deg > 70.0,
            "{which:?}: across u: {}",
            u.mean_angle_deg
        );
        assert!(
            other.mean_angle_deg < 20.0,
            "{which:?}: along the other: {}",
            other.mean_angle_deg
        );
        assert!(u.mean_aspect > base_u.mean_aspect * 0.0 + 1.5);
    }
}
