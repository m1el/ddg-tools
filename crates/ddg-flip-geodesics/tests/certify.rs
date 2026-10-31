//! The geodesic certificate on any network state (SPEC §7 layer 1a).

mod common;

use std::f64::consts::PI;

use common::*;
use ddg_flip_geodesics::{
    CertifyOptions, FlipError, GeodesicNetwork, JointStatus, Side, certify::classify,
};

#[test]
fn straight_bent_and_pinned_joints() {
    let shape = grid(4, 4);
    let g = |i, j| gv(4, i, j);
    let options = CertifyOptions::default();

    let tri = shape.triangulation();
    let line = vpath(&tri, &[g(0, 2), g(1, 2), g(2, 2), g(3, 2), g(4, 2)], false);
    let gn = GeodesicNetwork::new(tri.clone(), &[line], &[], None).unwrap();
    let cert = gn.certify(&options).unwrap();
    assert!(cert.geodesic);
    assert_eq!(
        (
            cert.n_joints,
            cert.n_straight,
            cert.n_pinned,
            cert.n_blocked
        ),
        (3, 3, 0, 0)
    );
    assert!((cert.min_angle - PI).abs() < 1e-14);
    assert!(cert.failures.is_empty());

    let bent = vpath(&tri, &[g(1, 2), g(2, 2), g(2, 3)], false);
    let gn = GeodesicNetwork::new(tri.clone(), std::slice::from_ref(&bent), &[], None).unwrap();
    let cert = gn.certify(&options).unwrap();
    assert!(!cert.geodesic);
    assert_eq!(cert.failures.len(), 1);
    let f = cert.failures[0];
    assert_eq!((f.vertex, f.side), (g(2, 2), Side::Left));
    assert!((f.angle - PI / 2.0).abs() < 1e-14 && f.flexible);
    assert_eq!(f.joint, gn.network().path_joints(0)[0]);
    assert!((cert.min_angle - PI / 2.0).abs() < 1e-14);

    // Pinned there: exempt.
    let gn = GeodesicNetwork::new(tri.clone(), &[bent], &[g(2, 2)], None).unwrap();
    let cert = gn.certify(&options).unwrap();
    assert!(cert.geodesic);
    assert_eq!((cert.n_pinned, cert.min_angle), (1, f64::INFINITY));

    // A loop around one cell bends at every corner.
    let square = vpath(&tri, &[g(1, 1), g(2, 1), g(2, 2), g(1, 2)], true);
    let gn = GeodesicNetwork::new(tri, &[square], &[], None).unwrap();
    let cert = gn.certify(&options).unwrap();
    assert_eq!((cert.n_joints, cert.failures.len()), (4, 4));
    assert!(
        cert.failures
            .iter()
            .all(|f| f.side == Side::Left && f.flexible)
    );
}

/// A wedge blocked by a path end, on a flat hexagon (centre 0, ring 1..=6 counter-clockwise).
#[test]
fn blocked_by_an_end() {
    let shape = hexagon();
    let tri = shape.triangulation();
    let options = CertifyOptions::default();
    // A path that starts at the centre and comes back through it: the centre is an end, so
    // pinned, and the joint is exempt as such.
    let path = vpath(&tri, &[0, 1, 2, 0, 6, 5], false);
    let gn = GeodesicNetwork::new(tri.clone(), std::slice::from_ref(&path), &[], None).unwrap();
    let j = gn.joint_geometry(gn.network().path_joints(0)[2]).unwrap();
    assert_eq!(j.vertex(), 0);
    assert_eq!(
        classify(&j, gn.network(), 1e-4),
        (JointStatus::Pinned, None)
    );

    // Path 0 turns at the centre, 2 → 0 → 6, with a wedge of 120° on its left that holds
    // edge 0–1; path 1 runs out and back along that edge from vertex 1, so its two segments
    // are a path end: path 0 is blocked there (path 1 itself is not straight: α = 0).
    let p0 = vpath(&tri, &[2, 0, 6], false);
    let p1 = vpath(&tri, &[1, 0, 1], false);
    let gn = GeodesicNetwork::new(tri.clone(), &[p0, p1], &[], None).unwrap();
    let j = gn.joint_geometry(gn.network().path_joints(0)[0]).unwrap();
    assert!((j.left.angle - 2.0 * PI / 3.0).abs() < 1e-14);
    assert!((j.right.angle - 4.0 * PI / 3.0).abs() < 1e-14);
    assert!(j.is_blocked_by_end(gn.network(), Side::Left));
    assert!(!j.is_flexible(Side::Left));
    assert!(j.is_flexible(Side::Right) && !j.is_blocked_by_end(gn.network(), Side::Right));
    assert_eq!(
        classify(&j, gn.network(), 1e-4),
        (JointStatus::BlockedByEnd, None)
    );
    let cert = gn.certify(&options).unwrap();
    assert!(!cert.geodesic);
    assert_eq!((cert.n_joints, cert.n_pinned, cert.n_blocked), (2, 0, 1));
    assert_eq!(cert.failures.len(), 1);
    let f = cert.failures[0];
    assert_eq!((f.joint.path, f.vertex, f.angle), (1, 0, 0.0));
    assert!(f.flexible);

    // Without the first segment in the wedge (start at 1 instead): not blocked.
    let path = vpath(&tri, &[1, 2, 0, 6, 5], false);
    let gn = GeodesicNetwork::new(tri, &[path], &[2, 6], None).unwrap();
    let cert = gn.certify(&options).unwrap();
    assert!(!cert.geodesic);
    assert_eq!(cert.failures.len(), 1);
    assert_eq!(
        (cert.failures[0].vertex, cert.failures[0].side),
        (0, Side::Left)
    );
}

#[test]
fn saddle_and_cone() {
    let options = CertifyOptions::default();
    let saddle = saddle();
    let tri = saddle.triangulation();
    for ring in [[1, 0, 4], [1, 0, 3], [2, 0, 6]] {
        let gn =
            GeodesicNetwork::new(tri.clone(), &[vpath(&tri, &ring, false)], &[], None).unwrap();
        assert!(gn.certify(&options).unwrap().geodesic, "{ring:?}");
    }
    // Θ < π: both sides bend; the failure names the smaller one.
    let cone = pyramid(5.0);
    let tri = cone.triangulation();
    let gn =
        GeodesicNetwork::new(tri.clone(), &[vpath(&tri, &[0, 4, 1], false)], &[], None).unwrap();
    let cert = gn.certify(&options).unwrap();
    assert_eq!(cert.failures.len(), 1);
    assert_eq!(
        (cert.failures[0].vertex, cert.failures[0].side),
        (4, Side::Right)
    );
}

#[test]
fn options_and_determinism() {
    let shape = grid(3, 3);
    let tri = shape.triangulation();
    let g = |i, j| gv(3, i, j);
    let path = vpath(&tri, &[g(0, 1), g(1, 1), g(1, 2)], false);
    let gn = GeodesicNetwork::new(tri, &[path], &[], None).unwrap();
    for bad in [-1.0, f64::NAN, f64::INFINITY] {
        let err = gn.certify(&CertifyOptions { tolerance: bad }).unwrap_err();
        assert!(matches!(err, FlipError::InvalidOptions { .. }), "{err}");
        assert!(err.to_string().contains("tolerance"));
    }
    // The tolerance decides: π/2 passes with a tolerance above π/2.
    assert!(
        gn.certify(&CertifyOptions { tolerance: 1.6 })
            .unwrap()
            .geodesic
    );
    assert_eq!(
        gn.certify(&CertifyOptions::default()).unwrap(),
        gn.clone().certify(&CertifyOptions::default()).unwrap()
    );
}
