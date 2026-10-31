//! Runs that stall (SPEC §4.1 "Floating point", F6): three pairs of the `tree` mesh whose
//! Dijkstra paths pass two nearly coincident vertices (`1.4e-4`, `1.2e-4` and `5.7e-4` apart
//! on a mesh of size about 1). Straightening the joint next to them needs a last flip that
//! makes a needle triangle (shape `1.2e-7` to `6.9e-7`), which IT N5's absolute shape test
//! (`FLIP_SHAPE_LIMIT = 1e-6`) refuses; the outer arc is then longer than the joint and the
//! joint stalls, still bent by `6e-4` to `4e-3` rad. The reference's flip test is relative
//! (`A ≥ 1e-6 (A1 + A2)`, `gc: signpost_intrinsic_triangulation.cpp:335`) and makes that flip.
//!
//! * With `ShortenOptions::flip_shape = FlipShapeTest::N5` these runs must not report
//!   `Geodesic` (their final state does not certify), and the two representations agree on
//!   them bitwise.
//! * With the default, the reference's relative test (IT N52), all three certify, end with
//!   status `Geodesic`, and have the reference's length (its driver's traces,
//!   `runs_tmp/fo-speed/gc-trace-tree-*.txt`), on both representations, bitwise alike.
//!
//! Pairs: 12804 → 1144 (the reference comparison of SPEC §7 layer 5), and pairs 372 and 934
//! of `random_pairs(n, 1000, 2)` (SPEC §7 test 9, F3's open trials).

use std::f64::consts::PI;

mod common;

use common::{Outcome, load_data, run};
use ddg_flip_geodesics::{ShortenOptions, ShortenStatus, dijkstra_path};
use ddg_integer_coords::IntegerTriangulation;
use ddg_intrinsic::{FlipShapeTest, Triangulation};

/// `(a, b, status, length bits, vertices of the failing joints)` as measured with IT N5's
/// shape test.
const PAIRS: [(usize, usize, ShortenStatus, u64, &[usize]); 3] = [
    (
        12804,
        1144,
        ShortenStatus::Stalled,
        0x3ffb_b93b_5060_8ee1,
        &[12753],
    ),
    (
        6878,
        4216,
        ShortenStatus::SweepLimit,
        0x3ff4_9b1e_aa23_12b6,
        &[5010, 4146],
    ),
    (
        4260,
        11433,
        ShortenStatus::Stalled,
        0x3ff9_79a3_1991_53de,
        &[4104],
    ),
];

fn check<T: Triangulation + Clone>(out: &Outcome<T>, pair: (usize, usize), label: &str) {
    let (a, b, status, bits, vertices) = PAIRS
        .iter()
        .copied()
        .find(|p| (p.0, p.1) == pair)
        .expect("a listed pair");
    let r = &out.report;
    assert_eq!(r.status, status, "{label} {a} -> {b}: {r:?}");
    assert_ne!(r.status, ShortenStatus::Geodesic);
    assert!(!out.certificate.geodesic, "{label} {a} -> {b}");
    assert!(r.refused_flips >= 1 && r.stalled >= 1, "{label}: {r:?}");
    assert_eq!(out.network.length().to_bits(), bits, "{label} {a} -> {b}");
    let mut failing: Vec<usize> = out.certificate.failures.iter().map(|f| f.vertex).collect();
    failing.sort_unstable();
    let mut expected = vertices.to_vec();
    expected.sort_unstable();
    assert_eq!(failing, expected, "{label} {a} -> {b}");
    for f in &out.certificate.failures {
        assert!(f.flexible, "{label}: {f:?}");
        assert!(f.angle > PI - 1e-2 && f.angle < PI - 1e-4, "{label}: {f:?}");
    }
}

#[test]
fn tree_pairs_with_a_refused_needle_flip_stall_under_n5() {
    let Some(shape) = load_data("optimal-cones/optimal-cones-data/data/tree.obj") else {
        return;
    };
    let signposts = shape.triangulation();
    let integers = IntegerTriangulation::new(&shape.mesh, &shape.geometry).unwrap();
    for &(a, b, ..) in &PAIRS {
        let path = dijkstra_path(&signposts, a, b).unwrap();
        let options = ShortenOptions {
            flip_shape: FlipShapeTest::N5,
            ..ShortenOptions::default()
        };
        let sp = run(
            signposts.clone(),
            std::slice::from_ref(&path),
            &[],
            &options,
        );
        check(&sp, (a, b), "signposts");
        let ic = run(integers.clone(), std::slice::from_ref(&path), &[], &options);
        check(&ic, (a, b), "integers");
        assert_eq!(sp.report, ic.report, "{a} -> {b}");
        assert_eq!(sp.network.paths(), ic.network.paths(), "{a} -> {b}");
    }
}

/// `(a, b, the reference's length)`: the final lengths of geometry-central's `FlipEdgeNetwork`
/// on these pairs (`runs_tmp/fo-speed/gc-trace-tree-{12804-1144,6878-4216,4260-11433}.txt`;
/// the first is also in `runs_tmp/fo-reference-driver/results/tree.gc`).
const REFERENCE: [(usize, usize, f64); 3] = [
    (12804, 1144, 1.732_722_581_670_872_5),
    (6878, 4216, 1.287_865_321_670_522_7),
    (4260, 11433, 1.592_196_559_669_969_4),
];

#[test]
fn tree_pairs_certify_with_the_reference_flip_test() {
    let Some(shape) = load_data("optimal-cones/optimal-cones-data/data/tree.obj") else {
        return;
    };
    let signposts = shape.triangulation();
    let integers = IntegerTriangulation::new(&shape.mesh, &shape.geometry).unwrap();
    let options = ShortenOptions::default();
    assert_eq!(options.flip_shape, FlipShapeTest::REFERENCE);
    for &(a, b, reference) in &REFERENCE {
        let path = dijkstra_path(&signposts, a, b).unwrap();
        let sp = run(
            signposts.clone(),
            std::slice::from_ref(&path),
            &[],
            &options,
        );
        let ic = run(integers.clone(), std::slice::from_ref(&path), &[], &options);
        for (out, label) in [(&sp.report, "signposts"), (&ic.report, "integers")] {
            assert_eq!(
                out.status,
                ShortenStatus::Geodesic,
                "{label} {a} -> {b}: {out:?}"
            );
            assert_eq!(out.refused_flips, 0, "{label} {a} -> {b}: {out:?}");
            assert_eq!(out.stalled, 0, "{label} {a} -> {b}: {out:?}");
        }
        for (cert, label) in [
            (&sp.certificate, "signposts"),
            (&ic.certificate, "integers"),
        ] {
            assert!(cert.geodesic, "{label} {a} -> {b}: {cert:?}");
            assert!(cert.failures.is_empty(), "{label} {a} -> {b}: {cert:?}");
        }
        let length = sp.network.length();
        assert!(
            (length - reference).abs() <= 1e-12 * reference,
            "{a} -> {b}: length {length:?}, the reference's {reference:?}"
        );
        assert_eq!(sp.report, ic.report, "{a} -> {b}");
        assert_eq!(sp.network.paths(), ic.network.paths(), "{a} -> {b}");
        assert_eq!(
            length.to_bits(),
            ic.network.length().to_bits(),
            "{a} -> {b}"
        );
    }
}
