//! The intrinsic constrained Delaunay triangulation of a network (SPEC §6.3, F20):
//! `GeodesicNetwork::make_delaunay` on shortened paths, networks and Bézier curves, on both
//! representations. Every path is still the same edge path, every edge that does not carry a
//! segment is Delaunay, both certificates pass (on `T` and on `T⁰`), the integer certificate
//! passes, the two representations agree bit for bit, and the run is deterministic.

mod common;

use common::*;
use ddg_flip_geodesics::{
    BezierOptions, CdtReport, CertifyOptions, EdgePath, ExtractOptions, FlipError, GeodesicNetwork,
    Network, ShortenOptions, ShortenStatus, dijkstra_path, geodesic_bezier,
};
use ddg_integer_coords::IntegerTriangulation;
use ddg_intrinsic::{DelaunayCriterion, DelaunayOptions, DelaunayStatus, Triangulation};

/// Checks the iCDT of `net` (already shortened); returns it, its report, and a fingerprint of
/// the report and the triangulation (faces, lengths as bits).
fn check_cdt<T: Triangulation + Sync>(
    mut net: GeodesicNetwork<T>,
    label: &str,
) -> (GeodesicNetwork<T>, CdtReport, String) {
    let paths = net.paths();
    let path_lengths: Vec<u64> = (0..net.network().n_paths())
        .map(|p| net.path_length(p).to_bits())
        .collect();
    let options = DelaunayOptions::default();
    let report = net
        .make_delaunay(&options)
        .unwrap_or_else(|e| panic!("{label}: {e}"));
    // The network is unchanged: same edge paths, lengths bitwise.
    assert_eq!(net.paths(), paths, "{label}: paths changed");
    let after: Vec<u64> = (0..net.network().n_paths())
        .map(|p| net.path_length(p).to_bits())
        .collect();
    assert_eq!(after, path_lengths, "{label}: lengths changed");
    net.check_invariants()
        .unwrap_or_else(|e| panic!("{label}: invariants: {e}"));
    // Every edge carrying a segment is constrained; every other interior edge is Delaunay.
    let tri = net.triangulation();
    let complex = tri.complex();
    let carrying = (0..complex.n_edges())
        .filter(|&e| net.network().carries_segment(e))
        .count();
    assert_eq!(report.n_network_edges, carrying, "{label}");
    for e in 0..complex.n_edges() {
        if net.network().carries_segment(e) {
            assert!(tri.is_constrained(e), "{label}: edge {e} not constrained");
        } else if let Some(excess) = tri.edge_angle_excess(e) {
            assert!(
                excess <= options.criterion_tolerance(),
                "{label}: edge {e} is not Delaunay (excess {excess:e})"
            );
        }
    }
    assert_eq!(
        report.delaunay.status,
        DelaunayStatus::Delaunay,
        "{label}: {report:?}"
    );
    assert_eq!(report.delaunay.n_not_delaunay, 0);
    assert_eq!(report.delaunay.n_constrained, complex.n_constrained());
    assert!(tri.is_delaunay(&options.criterion), "{label}");
    // Both certificates.
    let cert = net.certify(&CertifyOptions::default()).unwrap();
    assert!(cert.geodesic, "{label}: {:?}", cert.failures);
    let x = net
        .extract_with_report(&ExtractOptions::default())
        .unwrap_or_else(|e| panic!("{label}: extraction: {e}"));
    let c = x.certificate.as_ref().unwrap();
    assert!(c.passed, "{label}: certificate on T0: {c:#?}");
    let faces = format!("{:?}", tri.complex().faces());
    let lengths: Vec<u64> = tri.edge_lengths().iter().map(|x| x.to_bits()).collect();
    let fingerprint = format!("{report:?} {faces} {lengths:?}");
    (net, report, fingerprint)
}

trait CriterionTolerance {
    fn criterion_tolerance(&self) -> f64;
}

impl CriterionTolerance for DelaunayOptions {
    fn criterion_tolerance(&self) -> f64 {
        match self.criterion {
            DelaunayCriterion::AngleSum { tolerance } => tolerance,
            // A cotan weight `≥ −tol` allows `α + β − π` up to about `2 tol sin α sin β`; the
            // angle test below is the weaker statement `≤ 2 tol`.
            DelaunayCriterion::CotanWeight { tolerance } => 2.0 * tolerance,
        }
    }
}

fn shortened<T: Triangulation + Clone>(
    tri: T,
    paths: &[EdgePath],
    pinned: &[usize],
    label: &str,
) -> GeodesicNetwork<T> {
    let mut net =
        GeodesicNetwork::new(tri, paths, pinned, None).unwrap_or_else(|e| panic!("{label}: {e}"));
    let r = net.shorten(&ShortenOptions::default()).unwrap();
    assert_eq!(r.status, ShortenStatus::Geodesic, "{label}");
    net
}

fn shapes() -> Vec<Shape> {
    let mut shapes = generated_shapes();
    shapes.push(waist_tube(10, 3));
    shapes.extend(data_shapes());
    shapes
}

#[test]
fn icdt_of_shortened_networks_on_both_representations() {
    let mut cases = 0;
    let mut flips = 0;
    let mut constrained = 0;
    for shape in shapes() {
        let signposts = shape.triangulation();
        let integers = IntegerTriangulation::new(&shape.mesh, &shape.geometry).unwrap();
        let n = signposts.complex().n_vertices();
        let mut inputs: Vec<Vec<EdgePath>> = Vec::new();
        // One path, and four paths at once (a network; skipped where they cross).
        let pairs = random_pairs(n, 8, 2026);
        for &(a, b) in &pairs[..3] {
            inputs.push(vec![dijkstra_path(&signposts, a, b).unwrap()]);
        }
        inputs.push(
            pairs[3..7]
                .iter()
                .map(|&(a, b)| dijkstra_path(&signposts, a, b).unwrap())
                .collect(),
        );
        for (k, paths) in inputs.iter().enumerate() {
            if Network::new(signposts.complex(), paths, &[], None).is_err() {
                continue;
            }
            let label = format!("{} input {k}", shape.name);
            let a = shortened(signposts.clone(), paths, &[], &label);
            let b = shortened(integers.clone(), paths, &[], &label);
            let (_, report, fa) = check_cdt(a, &label);
            let (b, _, fb) = check_cdt(b, &format!("{label} (integers)"));
            assert_eq!(fa, fb, "{label}: the representations differ");
            b.triangulation()
                .certify()
                .unwrap_or_else(|e| panic!("{label}: integer certificate: {e}"));
            // Deterministic.
            let again = check_cdt(shortened(signposts.clone(), paths, &[], &label), &label).2;
            assert_eq!(fa, again, "{label}: not deterministic");
            cases += 1;
            flips += report.delaunay.n_flips;
            constrained += report.n_network_edges;
        }
    }
    eprintln!(
        "{cases} networks: iCDT identical on both, {flips} flips, {constrained} constrained edges"
    );
    assert!(cases >= 40, "only {cases}");
}

#[test]
fn icdt_of_bezier_curves() {
    let mut cases = 0;
    for shape in data_shapes().into_iter().chain([Shape::generated(
        "icosphere 3",
        ddg_core::generators::icosphere(3).unwrap(),
    )]) {
        let integers = IntegerTriangulation::new(&shape.mesh, &shape.geometry).unwrap();
        let n = shape.mesh.n_vertices();
        for (a, b) in random_pairs(n, 4, 5) {
            let tri = shape.triangulation();
            let vs = dijkstra_path(&tri, a, b).unwrap().vertices(tri.complex());
            if vs.len() < 6 {
                continue;
            }
            let controls = [
                vs[0],
                vs[vs.len() / 3],
                vs[2 * vs.len() / 3],
                vs[vs.len() - 1],
            ];
            let label = format!("{} {controls:?}", shape.name);
            let opts = BezierOptions::default();
            let (sa, ra) = match geodesic_bezier(tri, &controls, 3, &opts) {
                Ok(x) => x,
                Err(FlipError::Crossing { .. } | FlipError::Bezier { .. }) => continue,
                Err(e) => panic!("{label}: {e}"),
            };
            let (sb, rb) = geodesic_bezier(integers.clone(), &controls, 3, &opts).unwrap();
            assert_eq!(format!("{ra:?}"), format!("{rb:?}"), "{label}");
            let (_, _, fa) = check_cdt(sa, &label);
            let (net, _, fb) = check_cdt(sb, &format!("{label} (integers)"));
            assert_eq!(fa, fb, "{label}");
            net.triangulation().certify().unwrap();
            cases += 1;
        }
    }
    assert!(cases >= 6, "only {cases}");
}

#[test]
fn icdt_options_marks_and_later_shortening() {
    let shape = Shape::generated("icosphere 2", ddg_core::generators::icosphere(2).unwrap());
    let tri = shape.triangulation();
    let paths = [dijkstra_path(&tri, 0, 100).unwrap()];
    let mut net = shortened(tri, &paths, &[], "icosphere");
    // A bad tolerance changes nothing.
    let before = (
        net.triangulation().complex().faces().to_vec(),
        net.triangulation().complex().n_constrained(),
    );
    let bad = DelaunayOptions {
        criterion: DelaunayCriterion::AngleSum { tolerance: -1.0 },
        ..DelaunayOptions::default()
    };
    assert!(matches!(
        net.make_delaunay(&bad),
        Err(FlipError::InvalidOptions { .. })
    ));
    let after = (
        net.triangulation().complex().faces().to_vec(),
        net.triangulation().complex().n_constrained(),
    );
    assert_eq!(before, after);
    // Twice: the second finds everything constrained and Delaunay already.
    let first = net.make_delaunay(&DelaunayOptions::default()).unwrap();
    assert_eq!(first.n_newly_constrained, first.n_network_edges);
    let second = net.make_delaunay(&DelaunayOptions::default()).unwrap();
    assert_eq!(second.n_newly_constrained, 0);
    assert_eq!(second.delaunay.n_flips, 0);
    // After clearing the marks, a new network on the same triangulation shortens as usual.
    assert_eq!(net.clear_constraints(), first.n_network_edges);
    assert_eq!(net.triangulation().complex().n_constrained(), 0);
    let (tri, _) = net.into_parts();
    let more = [dijkstra_path(&tri, 5, 150).unwrap()];
    let net = shortened(tri, &more, &[], "after clearing");
    assert!(net.certify(&CertifyOptions::default()).unwrap().geodesic);
}
