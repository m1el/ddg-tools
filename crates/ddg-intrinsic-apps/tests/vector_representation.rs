//! The tangent fields on integer coordinates and on signposts (Integer Coordinates I26): both
//! representations of the same input, flipped to Delaunay and refined by the same rules, give
//! the same triangulation, and on it the same frames to rounding, so the transport, the
//! connection matrices, the smoothest fields, their singularities (exactly), the harmonic fields
//! and the maximum principle agree. The differences are printed (`--nocapture`).

mod common;

use std::f64::consts::{PI, TAU};

use common::*;
use ddg_integer_coords::IntegerTriangulation;
use ddg_intrinsic::{
    DelaunayOptions, IntrinsicTriangulation, RefineOptions, RefineRules, Triangulation,
};
use ddg_intrinsic_apps::vector::*;

fn shapes() -> Vec<Shape> {
    let mut shapes = vec![
        jittered_icosphere(2, 0.12, 5),
        jittered_icosphere(3, 0.03, 4),
        bumpy_sheared_grid(10, 0.25, 1.5),
        sheared_grid(12, 1.5, 0.45, 2),
        sheared_torus(20, 8, 1.0, 0.4, 3),
        polar_disk(6, 0.3),
    ];
    shapes.extend(dataset("face"));
    shapes.extend(dataset("camelbad"));
    shapes.extend(load_data("models/spot/spot/spot_triangulated.obj"));
    shapes
}

/// Both representations of `shape`, Delaunay-flipped and (if `refine`) refined by Ruppert's
/// rules, or `None` where a rounding-dependent decision made different triangulations.
fn pair(shape: &Shape, refine: bool) -> Option<(IntrinsicTriangulation, IntegerTriangulation)> {
    let mut s = shape.triangulation();
    let mut i = IntegerTriangulation::new(&shape.mesh, &shape.geometry).unwrap();
    s.flip_to_delaunay(&DelaunayOptions::default()).unwrap();
    i.flip_to_delaunay(&DelaunayOptions::default()).unwrap();
    if refine {
        let options = RefineOptions {
            max_insertions: ddg_intrinsic::InsertionLimit::Vertices(4000),
            ..RefineOptions::default()
        };
        s.refine_with_rules(&options, RefineRules::Ruppert).unwrap();
        i.refine_with_rules(&options, RefineRules::Ruppert).unwrap();
    }
    (s.complex() == i.complex()).then_some((s, i))
}

fn angle_difference(a: f64, b: f64) -> f64 {
    let d = (a - b).rem_euclid(TAU);
    d.min(TAU - d)
}

fn max_diff(a: &[f64], b: &[f64]) -> f64 {
    assert_eq!(a.len(), b.len());
    a.iter()
        .zip(b)
        .map(|(x, y)| (x - y).abs())
        .fold(0.0, f64::max)
}

/// `max |a − e^{iψ} b| / max |a|` with the phase `ψ` that aligns `b` with `a` (eigenvectors are
/// defined up to a unit factor).
fn aligned_difference(a: &[C64], b: &[C64]) -> f64 {
    let mut dot = C64::default();
    for (x, y) in a.iter().zip(b) {
        dot = dot + *x * y.conj();
    }
    let phase = C64::from_angle(dot.arg());
    let scale = a.iter().fold(0.0f64, |m, x| m.max(x.abs()));
    a.iter()
        .zip(b)
        .map(|(x, y)| (*x - phase * *y).abs())
        .fold(0.0, f64::max)
        / scale
}

#[test]
fn frames_transport_and_matrices_agree() {
    let mut n_compared = 0;
    for shape in shapes() {
        for refine in [false, true] {
            let Some((s, i)) = pair(&shape, refine) else {
                eprintln!(
                    "{} refine={refine}: different triangulations, skipped",
                    shape.name
                );
                continue;
            };
            n_compared += 1;
            let (ds, di) = (
                halfedge_directions_of(&s).unwrap(),
                halfedge_directions_of(&i).unwrap(),
            );
            assert_eq!(ds, halfedge_directions(&s));
            let dd = ds
                .iter()
                .zip(&di)
                .map(|(&a, &b)| angle_difference(a, b))
                .fold(0.0, f64::max);
            let (ts, ti) = (transport_angles(&s).unwrap(), transport_angles(&i).unwrap());
            let dt = ts
                .iter()
                .zip(&ti)
                .map(|(&a, &b)| angle_difference(a, b))
                .fold(0.0, f64::max);
            let mut dm = 0.0f64;
            for n_sym in [1u32, 2, 4] {
                for kind in [
                    ConnectionLaplacian::FiniteDifference,
                    ConnectionLaplacian::FiniteElement,
                ] {
                    let (ms, mi) = (
                        connection_matrices(&s, n_sym, kind).unwrap(),
                        connection_matrices(&i, n_sym, kind).unwrap(),
                    );
                    assert_eq!(ms.stiffness.col_ptr(), mi.stiffness.col_ptr());
                    assert_eq!(ms.stiffness.row_idx(), mi.stiffness.row_idx());
                    let scale = ms
                        .stiffness
                        .values()
                        .iter()
                        .fold(0.0f64, |m, v| m.max(v.abs()));
                    dm = dm.max(max_diff(ms.stiffness.values(), mi.stiffness.values()) / scale);
                    assert!(max_diff(ms.mass.values(), mi.mass.values()) < 1e-14);
                }
            }
            eprintln!(
                "{} refine={refine} ({} vertices): directions {dd:.2e}, transport {dt:.2e}, \
                 stiffness (relative) {dm:.2e}",
                shape.name,
                s.n_vertices()
            );
            assert!(dd < 1e-8 && dt < 1e-8 && dm < 1e-8, "{}", shape.name);
        }
    }
    assert!(n_compared >= 10, "{n_compared}");
}

#[test]
fn smoothest_fields_and_singularities_agree() {
    let mut n_compared = 0;
    for shape in shapes() {
        if shape.name.starts_with("polar_disk") || shape.name.starts_with("sheared_grid") {
            continue;
        }
        for refine in [false, true] {
            let Some((s, i)) = pair(&shape, refine) else {
                continue;
            };
            // Large meshes: vectors and cross fields with finite elements only (time).
            let large = s.n_vertices() > 5000;
            for n_sym in [1u32, 2, 4] {
                for laplacian in [
                    ConnectionLaplacian::FiniteElement,
                    ConnectionLaplacian::FiniteDifference,
                ] {
                    if large && (n_sym == 2 || laplacian == ConnectionLaplacian::FiniteDifference) {
                        continue;
                    }
                    let options = DirectionFieldOptions {
                        n_sym,
                        laplacian,
                        ..Default::default()
                    };
                    let (fs, rs) = smoothest_field(&s, &options).unwrap();
                    let (fi, ri) = smoothest_field(&i, &options).unwrap();
                    let du = aligned_difference(&fs.values, &fi.values);
                    let dl =
                        (rs.eigenvalue - ri.eigenvalue).abs() / rs.eigenvalue.abs().max(1e-300);
                    let (ss, si) = (
                        singularities(&s, &fs).unwrap(),
                        singularities(&i, &fi).unwrap(),
                    );
                    eprintln!(
                        "{} refine={refine} n {n_sym} {laplacian:?}: λ {:.6} ({:?}, {} it) vs {:.6} \
                         ({:?}, {} it), rel {dl:.1e}; field {du:.1e}; {} vs {} singular faces, sum {} {}",
                        shape.name,
                        rs.eigenvalue,
                        rs.status,
                        rs.iterations,
                        ri.eigenvalue,
                        ri.status,
                        ri.iterations,
                        ss.faces.len(),
                        si.faces.len(),
                        ss.index_sum,
                        si.index_sum
                    );
                    assert!(dl < 1e-8, "{}: eigenvalues differ by {dl}", shape.name);
                    assert_eq!(ss.index_sum, si.index_sum);
                    if du < 1e-6 {
                        // The same eigenvector: the same singular faces and indices exactly.
                        assert_eq!(ss.faces, si.faces, "{} n {n_sym}", shape.name);
                    } else {
                        // A cluster of nearly equal eigenvalues (the iteration's own tolerance
                        // decides which vector of the cluster comes out): count only.
                        eprintln!("    eigenvector not determined to 1e-6 (cluster): {du:.1e}");
                    }
                    // The same field on both gives the same singularities, whatever the cluster.
                    let cross = singularities(&i, &fs).unwrap();
                    assert_eq!(
                        cross.faces, ss.faces,
                        "{} n {n_sym}: the signpost field on integers",
                        shape.name
                    );
                    assert!((cross.max_rounding - ss.max_rounding).abs() < 1e-8);
                    // And the 3D vectors at the input vertices.
                    let (vs, vi) = (
                        fs.input_vectors_3d(&s, &shape.positions).unwrap(),
                        fs.input_vectors_3d(&i, &shape.positions).unwrap(),
                    );
                    assert_eq!(vs, vi);
                    n_compared += 1;
                }
            }
        }
    }
    assert!(n_compared > 25, "{n_compared}");
}

#[test]
fn harmonic_fields_and_the_maximum_principle_agree() {
    let mut n_compared = 0;
    for shape in shapes() {
        for refine in [false, true] {
            let Some((s, i)) = pair(&shape, refine) else {
                continue;
            };
            if (0..s.n_vertices()).all(|v| !s.is_boundary_vertex(v)) {
                continue;
            }
            let boundary_input: Vec<u32> = (0..shape.mesh.n_vertices() as u32)
                .filter(|&v| shape.mesh.is_boundary_vertex(v as usize))
                .collect();
            for n_sym in [1u32, 2] {
                for laplacian in [
                    ConnectionLaplacian::FiniteDifference,
                    ConnectionLaplacian::FiniteElement,
                ] {
                    let given: Vec<(u32, C64)> = boundary_input
                        .iter()
                        .map(|&v| (v, C64::from_angle(0.37 * v as f64 + PI / 7.0)))
                        .collect();
                    for boundary in [
                        BoundaryValues::InwardNormal,
                        BoundaryValues::Tangent,
                        BoundaryValues::AtInputVertices(given),
                    ] {
                        let options = HarmonicFieldOptions {
                            n_sym,
                            laplacian,
                            boundary: boundary.clone(),
                        };
                        let (fs, _) = harmonic_field(&s, &options).unwrap();
                        let (fi, _) = harmonic_field(&i, &options).unwrap();
                        let scale = fs.values.iter().fold(0.0f64, |m, x| m.max(x.abs()));
                        let du = fs
                            .values
                            .iter()
                            .zip(&fi.values)
                            .map(|(a, b)| (*a - *b).abs())
                            .fold(0.0, f64::max)
                            / scale;
                        let ps = check_maximum_principle(&s, &fs, 1e-9);
                        let pi = check_maximum_principle_of(&i, &fi, 1e-9).unwrap();
                        eprintln!(
                            "{} refine={refine} n {n_sym} {laplacian:?} {}: field {du:.1e}; \
                             outside {} vs {}, margin {:.3e} vs {:.3e}",
                            shape.name,
                            match boundary {
                                BoundaryValues::InwardNormal => "normal",
                                BoundaryValues::Tangent => "tangent",
                                _ => "given",
                            },
                            ps.n_outside_cone,
                            pi.n_outside_cone,
                            ps.min_margin,
                            pi.min_margin
                        );
                        assert!(du < 1e-8, "{}: {du}", shape.name);
                        assert_eq!(ps.n_checked, pi.n_checked);
                        assert_eq!(ps.n_outside_cone, pi.n_outside_cone);
                        // The margin jumps where a neighbour lies at the angle `±π` (the
                        // branch cut of `arg`), which rounding decides: only its sign must agree.
                        assert!(
                            (ps.min_margin - pi.min_margin).abs() < 1e-8
                                || (ps.min_margin >= 0.0) == (pi.min_margin >= 0.0),
                            "{}: margins {} vs {}",
                            shape.name,
                            ps.min_margin,
                            pi.min_margin
                        );
                        n_compared += 1;
                    }
                }
            }
        }
    }
    assert!(n_compared > 20, "{n_compared}");
}

#[test]
fn bad_integers_are_frame_errors() {
    let shape = jittered_icosphere(2, 0.12, 5);
    let (_, i) = pair(&shape, true).unwrap();
    let mut parts = i.into_parts();
    let e = parts.normal.iter().position(|&n| n > 0).unwrap();
    parts.normal[e] += 1;
    let bad = IntegerTriangulation::from_parts(parts).unwrap();
    let options = DirectionFieldOptions::default();
    match smoothest_field(&bad, &options) {
        Err(VectorError::Frames { reason }) => assert!(!reason.is_empty()),
        other => panic!("{other:?}"),
    }
    assert!(matches!(
        transport_angles(&bad),
        Err(VectorError::Frames { .. })
    ));
}
