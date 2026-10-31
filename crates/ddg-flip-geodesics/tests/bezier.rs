//! Geodesic Bézier curves (SPEC §6.1, §7 test 7, F19) and edge splits that carry the network
//! (`GeodesicNetwork::split_edge`): on a plane the pinned points after `r` rounds are the
//! planar de Casteljau subdivision points; on curved surfaces the result does not depend on
//! the direction of the control polygon; every result certifies on `T` and on `T⁰`, on both
//! representations; splits of edges with several segments keep the layering. The timings
//! (ignored, release) are at the end.

mod common;

use std::time::Instant;

use common::*;
use ddg_core::generators::{self, SplitMix64};
use ddg_flip_geodesics::{
    BezierOptions, BezierReport, CertifyOptions, EdgePath, ExtractOptions, FlipError,
    GeodesicNetwork, Layering, MidpointKind, Network, SegmentRef, ShortenOptions, ShortenStatus,
    geodesic_bezier, polyline_3d,
};
use ddg_integer_coords::IntegerTriangulation;
use ddg_intrinsic::{DelaunayOptions, PathOnInput, Triangulation};

// ----- helpers -----

/// The positions in space of vertices of `tri` (through their places on the input mesh).
fn positions_of<T: Triangulation>(tri: &T, shape: &Shape, vertices: &[usize]) -> Vec<[f64; 3]> {
    let path = PathOnInput {
        points: vertices.iter().map(|&v| tri.vertex_location(v)).collect(),
        faces: Vec::new(),
        mismatch: 0.0,
    };
    polyline_3d(&shape.mesh, &shape.positions, &path).unwrap()
}

/// A flat unit square, `n × n` cells, interior vertices jittered by up to `±0.3` of a cell
/// (so no straight line between two vertices passes through a third), diagonals alternating.
fn jittered_square(n: usize, seed: u64) -> Shape {
    let mut rng = SplitMix64::new(seed);
    let h = 1.0 / n as f64;
    let mut positions = Vec::new();
    for j in 0..=n {
        for i in 0..=n {
            let (mut x, mut y) = (i as f64 * h, j as f64 * h);
            if i > 0 && i < n && j > 0 && j < n {
                x += (rng.next_f64() - 0.5) * 0.6 * h;
                y += (rng.next_f64() - 0.5) * 0.6 * h;
            }
            positions.push([x, y, 0.0]);
        }
    }
    let id = |i: usize, j: usize| (j * (n + 1) + i) as u32;
    let mut faces = Vec::new();
    for j in 0..n {
        for i in 0..n {
            let (a, b, c, d) = (id(i, j), id(i + 1, j), id(i, j + 1), id(i + 1, j + 1));
            if (i + j) % 2 == 0 {
                faces.extend([[a, b, d], [a, d, c]]);
            } else {
                faces.extend([[a, b, c], [b, d, c]]);
            }
        }
    }
    Shape::new(
        &format!("jittered square {n} seed {seed}"),
        positions,
        faces,
    )
}

/// The vertex nearest to `p`.
fn nearest(shape: &Shape, p: [f64; 3]) -> usize {
    (0..shape.positions.len())
        .min_by(|&a, &b| dist(shape.positions[a], p).total_cmp(&dist(shape.positions[b], p)))
        .unwrap()
}

fn mid(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [0, 1, 2].map(|k| 0.5 * (a[k] + b[k]))
}

/// The two control polygons of de Casteljau's subdivision at `½`.
fn split_half(p: &[[f64; 3]]) -> (Vec<[f64; 3]>, Vec<[f64; 3]>) {
    let mut level = p.to_vec();
    let mut left = vec![level[0]];
    let mut right = vec![level[level.len() - 1]];
    while level.len() > 1 {
        level = level.windows(2).map(|w| mid(w[0], w[1])).collect();
        left.push(level[0]);
        right.push(level[level.len() - 1]);
    }
    right.reverse();
    (left, right)
}

/// The control points of the `2^r` pieces of planar de Casteljau subdivision, in curve order,
/// shared end points once (`fo-de-casteljau-midpoint`).
fn de_casteljau(p: &[[f64; 3]], rounds: usize) -> Vec<[f64; 3]> {
    if rounds == 0 {
        return p.to_vec();
    }
    let (left, right) = split_half(p);
    let mut out = de_casteljau(&left, rounds - 1);
    out.pop();
    out.extend(de_casteljau(&right, rounds - 1));
    out
}

fn options_checked() -> BezierOptions {
    BezierOptions {
        shorten: ShortenOptions {
            check_invariants: true,
            ..ShortenOptions::default()
        },
        ..BezierOptions::default()
    }
}

/// Checks a Bézier result: the report's counts, the pinned vertices are the control points,
/// the certificate on `T` (pinned joints exempt), extraction with the certificate on `T⁰`.
fn check_result<T: Triangulation + Sync>(
    net: &GeodesicNetwork<T>,
    report: &BezierReport,
    d: usize,
    rounds: usize,
    label: &str,
) {
    assert_eq!(
        report.control_points.len(),
        d * (1 << rounds) + 1,
        "{label}: number of control points"
    );
    assert_eq!(
        report.midpoints.len(),
        ((1 << rounds) - 1) * d * (d + 1) / 2,
        "{label}: midpoints made (d (d + 1) / 2 per subdivision, d of them kept)"
    );
    assert_eq!(report.n_splits + report.n_snapped, report.midpoints.len());
    assert_eq!(
        report.n_not_geodesic, 0,
        "{label}: {:?}",
        report.other_statuses
    );
    net.check_invariants().unwrap();
    let complex = net.triangulation().complex();
    let pinned: Vec<usize> = net
        .network()
        .path(0)
        .vertices(complex)
        .into_iter()
        .filter(|&v| net.network().is_pinned(v))
        .collect();
    assert_eq!(
        pinned, report.control_points,
        "{label}: pinned along the path"
    );
    let mut all = report.control_points.clone();
    all.sort_unstable();
    all.dedup();
    assert_eq!(all, net.network().pinned_vertices(), "{label}: pinned set");
    assert!(
        (net.length() - net.recompute_length()).abs() <= 1e-12 * net.length(),
        "{label}: running length"
    );
    let cert = net.certify(&CertifyOptions::default()).unwrap();
    assert!(cert.geodesic, "{label}: {:?}", cert.failures);
    let x = net
        .extract_with_report(&ExtractOptions::default())
        .unwrap_or_else(|e| panic!("{label}: extraction: {e}"));
    let c = x.certificate.as_ref().unwrap();
    assert!(c.passed, "{label}: certificate on T0: {c:#?}");
}

// ----- §7 test 7: the plane -----

/// The lengths of the parts of path 0 between consecutive pinned points, in order.
fn part_lengths<T: Triangulation>(net: &GeodesicNetwork<T>) -> Vec<f64> {
    let tri = net.triangulation();
    let path = net.network().path(0);
    let vs = path.vertices(tri.complex());
    let mut out = Vec::new();
    let mut len = 0.0;
    for (s, v) in path.segments.iter().zip(&vs[1..]) {
        len += tri.edge_length(s.edge as usize);
        if net.network().is_pinned(*v) {
            out.push(len);
            len = 0.0;
        }
    }
    out
}

/// The pinned points after `r` rounds are planar de Casteljau subdivision (SPEC §7 test 7).
///
/// Both representations make the same splits with bitwise the same lengths (the report and the
/// paths are compared), so the intrinsic result is one: its part lengths are the distances
/// between consecutive de Casteljau points to `1e-12`, and on integer coordinates the places
/// of the points on the input mesh are those points to `1e-12`. On signposts the place of a
/// new vertex is a trace over the input mesh from the split edge's signpost (IT N7, N8), whose
/// direction carries the rounding of every flip before; those places are within `1e-10`
/// (measured: `1.2e-11` on the 24 x 24 square at 4 rounds, `1.5e-13` on 12 x 12).
#[test]
fn planar_bezier_is_de_casteljau_subdivision() {
    let polygons: [&[[f64; 3]]; 4] = [
        &[[0.15, 0.2, 0.0], [0.85, 0.75, 0.0]],
        &[[0.1, 0.15, 0.0], [0.5, 0.9, 0.0], [0.9, 0.2, 0.0]],
        &[
            [0.1, 0.1, 0.0],
            [0.3, 0.9, 0.0],
            [0.8, 0.85, 0.0],
            [0.9, 0.15, 0.0],
        ],
        &[
            [0.12, 0.5, 0.0],
            [0.3, 0.88, 0.0],
            [0.62, 0.1, 0.0],
            [0.88, 0.6, 0.0],
            [0.7, 0.9, 0.0],
        ],
    ];
    // Worst: [integer places, part lengths, signpost places].
    let mut worst = [0.0f64; 3];
    let mut cases = 0;
    let mut skipped = 0;
    for (n, seed) in [(12, 3), (24, 5)] {
        let shape = jittered_square(n, seed);
        for polygon in polygons {
            let controls: Vec<usize> = polygon.iter().map(|&p| nearest(&shape, p)).collect();
            let exact: Vec<[f64; 3]> = controls.iter().map(|&v| shape.positions[v]).collect();
            let d = controls.len() - 1;
            for rounds in 0..=4 {
                let want = de_casteljau(&exact, rounds);
                for delaunay in [false, true] {
                    let label = format!(
                        "{}: d {d}, rounds {rounds}, delaunay {delaunay}",
                        shape.name
                    );
                    let mut integers =
                        IntegerTriangulation::new(&shape.mesh, &shape.geometry).unwrap();
                    let signposts = if delaunay {
                        integers
                            .flip_to_delaunay(&DelaunayOptions::default())
                            .unwrap();
                        shape.delaunay()
                    } else {
                        shape.triangulation()
                    };
                    let (inet, ireport) =
                        match geodesic_bezier(integers, &controls, rounds, &options_checked()) {
                            Ok(x) => x,
                            Err(e @ FlipError::Crossing { .. }) => {
                                // The waypoint legs touch: refused on both, never repaired.
                                let other = geodesic_bezier(
                                    signposts,
                                    &controls,
                                    rounds,
                                    &options_checked(),
                                );
                                assert!(matches!(other, Err(FlipError::Crossing { .. })));
                                eprintln!("SKIPPED {label}: {e}");
                                skipped += 1;
                                continue;
                            }
                            Err(e) => panic!("{label} (integers): {e}"),
                        };
                    check_result(&inet, &ireport, d, rounds, &format!("{label} (integers)"));
                    inet.triangulation().certify().unwrap();
                    assert!(
                        ireport
                            .midpoints
                            .iter()
                            .all(|m| m.kind == MidpointKind::Split),
                        "{label}: no midpoint lands near a vertex here"
                    );
                    let got = positions_of(inet.triangulation(), &shape, &ireport.control_points);
                    assert_eq!(got.len(), want.len(), "{label}");
                    for (i, (p, q)) in got.iter().zip(&want).enumerate() {
                        let err = dist(*p, *q);
                        worst[0] = worst[0].max(err);
                        assert!(err <= 1e-12, "{label}: point {i}: {p:?} vs {q:?} ({err:e})");
                    }
                    for (i, (l, w)) in part_lengths(&inet).iter().zip(want.windows(2)).enumerate() {
                        let err = (l - dist(w[0], w[1])).abs();
                        worst[1] = worst[1].max(err);
                        assert!(err <= 1e-12, "{label}: part {i}: {l} ({err:e})");
                    }
                    // Signposts: the same subdivision; the places carry the trace's rounding.
                    let (snet, sreport) =
                        geodesic_bezier(signposts, &controls, rounds, &options_checked())
                            .unwrap_or_else(|e| panic!("{label}: {e}"));
                    check_result(&snet, &sreport, d, rounds, &label);
                    assert_eq!(format!("{sreport:?}"), format!("{ireport:?}"), "{label}");
                    assert_eq!(snet.paths(), inet.paths(), "{label}");
                    let got = positions_of(snet.triangulation(), &shape, &sreport.control_points);
                    for (i, (p, q)) in got.iter().zip(&want).enumerate() {
                        let err = dist(*p, *q);
                        worst[2] = worst[2].max(err);
                        assert!(err <= 1e-10, "{label} (signposts): point {i}: {err:e}");
                    }
                    cases += 1;
                }
            }
        }
    }
    assert!(skipped * 10 <= cases, "{skipped} skipped, {cases} run");
    eprintln!(
        "{cases} planar cases ({skipped} refused inputs skipped); worst distance to de Casteljau: integers {:e}, part lengths \
         {:e}, signposts {:e}",
        worst[0], worst[1], worst[2]
    );
}

/// On a regular grid the midpoints of a line between lattice points are lattice points or
/// edge midpoints: some land on vertices and are snapped (F19).
#[test]
fn midpoints_on_vertices_are_snapped() {
    let shape = grid(8, 8);
    let tri = shape.triangulation();
    // (0, 0) → (8, 8) along the diagonal edges of the grid's cells? Use a straight row.
    let controls = [gv(8, 0, 4), gv(8, 8, 4)];
    let (net, report) = geodesic_bezier(tri, &controls, 3, &options_checked()).unwrap();
    check_result(&net, &report, 1, 3, "grid row");
    assert_eq!(report.n_splits, 0, "{:?}", report.midpoints);
    assert_eq!(report.n_snapped, 7);
    let want: Vec<usize> = (0..=8).map(|i| gv(8, i, 4)).collect();
    assert_eq!(report.control_points, want);
    // Snapped midpoints are never split, and the triangulation is unchanged.
    assert_eq!(
        net.triangulation().complex().n_vertices(),
        shape.mesh.n_vertices()
    );
}

// ----- §7 test 7: curved surfaces, reversal -----

fn curved_shapes() -> Vec<Shape> {
    let mut bumpy = generators::icosphere(3).unwrap();
    generators::perturb(&mut bumpy.positions, [0.03; 3], 11);
    let mut shapes = vec![
        Shape::generated("icosphere 3", generators::icosphere(3).unwrap()),
        Shape::generated("jittered icosphere 3", bumpy),
        Shape::generated("torus 24x12", generators::torus(24, 12, 1.0, 0.4).unwrap()),
        Shape::generated(
            "bumpy grid 16x16",
            generators::bumpy_grid_disk(16, 16, 0.2).unwrap(),
        ),
    ];
    shapes.extend(data_shapes());
    shapes
}

/// `count` control polygons of `d + 1` distinct vertices: the ends `0.15`–`0.35` bounding-box
/// diagonals apart, the inner points random vertices within `0.06` diagonals of the points
/// at `i/d` of the shortest edge path between them (gentle polygons: the curve does not cross
/// itself).
fn local_polygons(shape: &Shape, d: usize, count: usize, seed: u64) -> Vec<Vec<usize>> {
    let tri = shape.triangulation();
    let scale = bbox_diagonal(shape);
    let mut rng = SplitMix64::new(seed);
    let n = shape.positions.len();
    let mut out = Vec::new();
    let mut tries = 0;
    while out.len() < count {
        tries += 1;
        assert!(tries < 100_000, "{}: no polygons", shape.name);
        let u = (rng.next_u64() % n as u64) as usize;
        let w = (rng.next_u64() % n as u64) as usize;
        let r = dist(shape.positions[u], shape.positions[w]) / scale;
        if !(0.15..0.35).contains(&r) {
            continue;
        }
        let Ok(path) = ddg_flip_geodesics::dijkstra_path(&tri, u, w) else {
            continue;
        };
        let vs = path.vertices(tri.complex());
        let mut poly = vec![u];
        for i in 1..d {
            let c = shape.positions[vs[i * (vs.len() - 1) / d]];
            let near: Vec<usize> = (0..n)
                .filter(|&v| dist(shape.positions[v], c) < 0.06 * scale)
                .filter(|v| !poly.contains(v) && *v != w)
                .collect();
            if near.is_empty() {
                break;
            }
            poly.push(near[(rng.next_u64() % near.len() as u64) as usize]);
        }
        poly.push(w);
        if poly.len() == d + 1 {
            out.push(poly);
        }
    }
    out
}

fn bbox_diagonal(shape: &Shape) -> f64 {
    let mut lo = [f64::INFINITY; 3];
    let mut hi = [f64::NEG_INFINITY; 3];
    for p in &shape.positions {
        for k in 0..3 {
            lo[k] = lo[k].min(p[k]);
            hi[k] = hi[k].max(p[k]);
        }
    }
    dist(lo, hi)
}

/// One curve and its reversal on `tri`: the largest distance between corresponding pinned
/// points over the bounding box diagonal, and whether a midpoint was nudged; `None` for an
/// input refused in either direction (a waypoint path that crosses itself or meets a control
/// point twice; Dijkstra's legs may differ in the two directions).
fn reversal_error<T: Triangulation + Sync + Clone>(
    tri: &T,
    shape: &Shape,
    controls: &[usize],
    options: &BezierOptions,
    label: &str,
) -> Option<(f64, bool)> {
    let rounds = 3;
    let reversed: Vec<usize> = controls.iter().rev().copied().collect();
    let run = |c: &[usize]| geodesic_bezier(tri.clone(), c, rounds, options);
    let refused =
        |e: &FlipError| matches!(e, FlipError::Bezier { .. } | FlipError::Crossing { .. });
    let (fwd, rf, bwd, rb) = match (run(controls), run(&reversed)) {
        (Ok((fwd, rf)), Ok((bwd, rb))) => (fwd, rf, bwd, rb),
        (Err(e), _) | (_, Err(e)) if refused(&e) => {
            eprintln!("SKIPPED {label}: {e}");
            return None;
        }
        (Err(e), _) | (_, Err(e)) => panic!("{label}: {e}"),
    };
    check_result(&fwd, &rf, controls.len() - 1, rounds, label);
    check_result(
        &bwd,
        &rb,
        controls.len() - 1,
        rounds,
        &format!("{label} reversed"),
    );
    let a = positions_of(fwd.triangulation(), shape, &rf.control_points);
    let mut b = positions_of(bwd.triangulation(), shape, &rb.control_points);
    b.reverse();
    assert_eq!(a.len(), b.len());
    let scale = bbox_diagonal(shape);
    let err = a
        .iter()
        .zip(&b)
        .map(|(p, q)| dist(*p, *q) / scale)
        .fold(0.0, f64::max);
    let nudged = rf.n_nudged + rb.n_nudged > 0;
    let length_err = (fwd.length() - bwd.length()).abs() / fwd.length();
    assert!(
        length_err <= if nudged { 1e-8 } else { 1e-12 },
        "{label}: lengths {} vs {}",
        fwd.length(),
        bwd.length()
    );
    Some((err, nudged))
}

/// SPEC §7 test 7 on curved surfaces: reversing the control points gives the same points.
///
/// Signposts run with `unlocated_nudge = 0` (they never need it): `1e-12`. The integers run
/// with the default: a midpoint that lies on an input edge is moved by `2e-9 ℓ` towards the
/// middle of the edge it falls on, and the two directions may split different edges there, so
/// such curves agree to `1e-8` only (counted); the others to `1e-12`.
#[test]
fn bezier_is_symmetric_under_reversal() {
    let mut worst = [0.0f64; 3];
    let mut cases = [0usize; 3];
    for shape in curved_shapes() {
        let integers = IntegerTriangulation::new(&shape.mesh, &shape.geometry).unwrap();
        let signposts = shape.triangulation();
        for (k, controls) in local_polygons(&shape, 3, 4, 41).into_iter().enumerate() {
            let label = format!("{} polygon {k} {controls:?}", shape.name);
            let no_nudge = BezierOptions {
                unlocated_nudge: 0.0,
                ..options_checked()
            };
            if let Some((err, _)) = reversal_error(&signposts, &shape, &controls, &no_nudge, &label)
            {
                assert!(err <= 1e-12, "{label}: {err:e}");
                worst[0] = worst[0].max(err);
                cases[0] += 1;
            }
            let label = format!("{label} (integers)");
            if let Some((err, nudged)) =
                reversal_error(&integers, &shape, &controls, &options_checked(), &label)
            {
                let (i, tol) = if nudged { (2, 1e-8) } else { (1, 1e-12) };
                assert!(err <= tol, "{label}: {err:e} (nudged {nudged})");
                worst[i] = worst[i].max(err);
                cases[i] += 1;
            }
        }
    }
    eprintln!(
        "reversal: signposts {} cases, worst {:e}; integers {} cases, worst {:e}, and {} with a \
         nudged midpoint, worst {:e} (of the bounding box diagonal)",
        cases[0], worst[0], cases[1], worst[1], cases[2], worst[2]
    );
    assert!(cases[0] >= 12 && cases[1] >= 10, "only {cases:?} cases");
}

// ----- extraction, both representations, determinism -----

#[test]
fn bezier_on_both_representations_certifies_and_is_deterministic() {
    let mut cases = 0;
    for shape in curved_shapes() {
        let integers = IntegerTriangulation::new(&shape.mesh, &shape.geometry).unwrap();
        for controls in local_polygons(&shape, 2, 2, 7) {
            let label = format!("{} {controls:?}", shape.name);
            let Ok((a, ra)) =
                geodesic_bezier(shape.triangulation(), &controls, 4, &options_checked())
            else {
                continue;
            };
            let (b, rb) = geodesic_bezier(integers.clone(), &controls, 4, &options_checked())
                .unwrap_or_else(|e| panic!("{label} (integers): {e}"));
            check_result(&a, &ra, 2, 4, &label);
            check_result(&b, &rb, 2, 4, &format!("{label} (integers)"));
            b.triangulation().certify().unwrap();
            assert_eq!(format!("{ra:?}"), format!("{rb:?}"), "{label}");
            assert_eq!(a.paths(), b.paths(), "{label}");
            // Determinism: the same run again, and the extraction at several thread counts.
            let (a2, ra2) =
                geodesic_bezier(shape.triangulation(), &controls, 4, &options_checked()).unwrap();
            assert_eq!(format!("{ra:?}"), format!("{ra2:?}"), "{label}");
            assert_eq!(a.paths(), a2.paths());
            let run = |threads: usize| {
                rayon::ThreadPoolBuilder::new()
                    .num_threads(threads)
                    .build()
                    .unwrap()
                    .install(|| b.extract_with_report(&ExtractOptions::default()).unwrap())
            };
            assert_eq!(run(1), run(8), "{label}");
            cases += 1;
        }
    }
    assert!(cases >= 8, "only {cases} cases");
}

// ----- refused inputs -----

#[test]
fn bad_inputs_are_refused() {
    let shape = grid(6, 6);
    let tri = shape.triangulation();
    let opts = BezierOptions::default();
    let e = geodesic_bezier(tri.clone(), &[3], 2, &opts).unwrap_err();
    assert!(matches!(e, FlipError::TooFewVertices { .. }), "{e}");
    let e = geodesic_bezier(tri.clone(), &[3, 3, 10], 2, &opts).unwrap_err();
    assert!(matches!(e, FlipError::WaypointLeg { .. }), "{e}");
    let e = geodesic_bezier(tri.clone(), &[3, 30, 10], 25, &opts).unwrap_err();
    assert!(matches!(e, FlipError::InvalidOptions { .. }), "{e}");
    let bad = BezierOptions {
        snap_tolerance: 0.5,
        ..opts
    };
    let e = geodesic_bezier(tri.clone(), &[3, 30, 10], 2, &bad).unwrap_err();
    assert!(matches!(e, FlipError::InvalidOptions { .. }), "{e}");
    // A control point repeated (not first = last).
    let e = geodesic_bezier(
        tri.clone(),
        &[gv(6, 1, 1), gv(6, 5, 1), gv(6, 1, 1), gv(6, 5, 5)],
        2,
        &opts,
    )
    .unwrap_err();
    assert!(
        matches!(e, FlipError::Bezier { .. } | FlipError::Crossing { .. }),
        "{e}"
    );
    // Not one open path: nothing changes.
    let loop_path = vpath(
        &tri,
        &[gv(6, 2, 2), gv(6, 2, 3), gv(6, 3, 3), gv(6, 3, 2)],
        true,
    );
    let mut net = GeodesicNetwork::new(tri.clone(), &[loop_path], &[], None).unwrap();
    let before = (net.paths(), net.triangulation().complex().faces().to_vec());
    let e = net.bezier_subdivide(2, &opts).unwrap_err();
    assert!(matches!(e, FlipError::Bezier { .. }), "{e}");
    assert_eq!(
        before,
        (net.paths(), net.triangulation().complex().faces().to_vec())
    );
    let two = [
        vpath(&tri, &[gv(6, 1, 1), gv(6, 1, 2)], false),
        vpath(&tri, &[gv(6, 4, 4), gv(6, 4, 5)], false),
    ];
    let mut net = GeodesicNetwork::new(tri.clone(), &two, &[], None).unwrap();
    assert!(matches!(
        net.bezier_subdivide(1, &opts),
        Err(FlipError::Bezier { .. })
    ));
    // A closed-looking curve (P₀ = P_d) is fine.
    let ring = [
        gv(6, 2, 2),
        gv(6, 2, 4),
        gv(6, 4, 4),
        gv(6, 4, 2),
        gv(6, 2, 2),
    ];
    let (net, report) = geodesic_bezier(tri, &ring, 2, &options_checked()).unwrap();
    check_result(&net, &report, 4, 2, "closed-looking");
    assert_eq!(report.control_points.first(), report.control_points.last());
}

// ----- splits of edges with several segments -----

/// Three copies of a path through the middle of a grid side by side (two in one direction,
/// one in the other): splitting their edges keeps every segment's place across the edge.
#[test]
fn split_edge_keeps_the_layering_of_every_segment() {
    let shape = grid(6, 6);
    for delaunay in [false, true] {
        let tri = if delaunay {
            shape.delaunay()
        } else {
            shape.triangulation()
        };
        let route = [
            gv(6, 1, 3),
            gv(6, 2, 3),
            gv(6, 3, 3),
            gv(6, 4, 3),
            gv(6, 5, 3),
        ];
        let Ok(a) = ddg_flip_geodesics::path_from_vertices(&tri, &route, false) else {
            continue;
        };
        let back = a.reversed();
        let paths = vec![a.clone(), a.clone(), back];
        // Left to right of path 0: path 0, path 1, path 2 (which runs the other way).
        let complex = tri.complex();
        let mut layering: Layering = vec![Vec::new(); complex.n_edges()];
        let n = a.segments.len();
        for (i, s) in a.segments.iter().enumerate() {
            let order = [
                SegmentRef { path: 0, index: i },
                SegmentRef { path: 1, index: i },
                SegmentRef {
                    path: 2,
                    index: n - 1 - i,
                },
            ];
            let list = &mut layering[s.edge as usize];
            if s.forward {
                list.extend(order);
            } else {
                list.extend(order.iter().rev());
            }
        }
        let mut net = GeodesicNetwork::new(tri.clone(), &paths, &[], Some(&layering)).unwrap();
        for (k, s) in a.segments.iter().enumerate() {
            let e = s.edge as usize;
            let before: Vec<_> = net
                .network()
                .layer(e)
                .map(|id| net.network().segment(id).unwrap().path)
                .collect();
            let t = 0.3 + 0.1 * k as f64;
            let split = net.split_edge(e, t).unwrap();
            let [h0, h1] = split.insertion.halves.unwrap();
            let m = split.insertion.vertex as usize;
            net.check_invariants()
                .unwrap_or_else(|err| panic!("delaunay {delaunay}, split {k}: {err}"));
            for h in [h0, h1] {
                let after: Vec<_> = net
                    .network()
                    .layer(h as usize)
                    .map(|id| net.network().segment(id).unwrap().path)
                    .collect();
                assert_eq!(after, before, "split {k}: layer of half {h}");
            }
            assert_eq!(split.segments.len(), 3);
            for &(first, second) in &split.segments {
                let f = net.network().segment(first).unwrap();
                let g = net.network().segment(second).unwrap();
                assert_eq!(f.next, Some(second));
                assert_eq!(g.prev, Some(first));
                let complex = net.triangulation().complex();
                assert_eq!(f.segment.head(complex), m);
                assert_eq!(g.segment.tail(complex), m);
            }
            // Every path now passes through the new vertex, with the joint straight.
            for p in 0..3 {
                assert!(
                    net.network()
                        .path(p)
                        .vertices(net.triangulation().complex())
                        .contains(&m)
                );
            }
            let cert = net.certify(&CertifyOptions::default()).unwrap();
            assert!(cert.geodesic, "{:?}", cert.failures);
        }
        for p in 0..3 {
            assert_eq!(net.network().n_path_segments(p), 2 * n);
        }
        assert!(
            (net.length() - 3.0 * a.length(&tri)).abs() <= 1e-12 * net.length(),
            "the halves add up"
        );
        // The network as rebuilt from its paths and layering is the same.
        let rebuilt = Network::new(
            net.triangulation().complex(),
            &net.paths(),
            &net.network().pinned_vertices(),
            Some(&net.network().layering()),
        )
        .unwrap();
        assert_eq!(rebuilt.paths(), net.paths());
        assert_eq!(rebuilt.layering(), net.network().layering());
        // Shortening the three copies afterwards keeps them side by side.
        let report = net
            .shorten(&ShortenOptions {
                check_invariants: true,
                ..ShortenOptions::default()
            })
            .unwrap();
        assert_eq!(report.status, ShortenStatus::Geodesic);
    }
    // A loop of one segment along a self-edge, split: two segments, still a loop.
    let shape = needle_tetrahedron(8.0, 4.0);
    let tri = shape.delaunay();
    let self_edges = tri.irregularities().self_edges.clone();
    for s in self_edges {
        let path = EdgePath::closed(vec![ddg_flip_geodesics::Segment::new(s as usize, true)]);
        let mut net = GeodesicNetwork::new(tri.clone(), &[path], &[], None).unwrap();
        let split = net.split_edge(s as usize, 0.5).unwrap();
        net.check_invariants().unwrap();
        assert_eq!(net.network().n_path_segments(0), 2);
        assert!(net.network().is_closed(0));
        assert_eq!(split.segments.len(), 1);
    }
    // Refused splits change nothing.
    let tri = grid(3, 3).triangulation();
    let p = vpath(&tri, &[0, 1, 2], false);
    let e = p.segments[0].edge as usize;
    let mut net = GeodesicNetwork::new(tri, &[p], &[], None).unwrap();
    let before = (net.paths(), net.network().layering());
    for t in [0.0, 1.0, -0.5, f64::NAN] {
        assert!(net.split_edge(e, t).is_err());
    }
    assert!(matches!(
        net.split_edge(10_000, 0.5),
        Err(FlipError::IndexOutOfRange { .. })
    ));
    assert_eq!(before, (net.paths(), net.network().layering()));
}

// ----- timings (SPEC §9.3) -----

/// Bézier curves with 4 rounds and the iCDT after them on spot and camelbad, both
/// representations, single thread. Run alone:
/// `taskset -c 5 cargo test --release -p ddg-flip-geodesics --test bezier -- --ignored
/// --nocapture --test-threads=1`.
#[test]
#[ignore]
fn timings() {
    rayon::ThreadPoolBuilder::new()
        .num_threads(1)
        .build_global()
        .ok();
    let load = std::fs::read_to_string("/proc/loadavg").unwrap_or_default();
    println!("load before: {}", load.trim());
    for shape in data_shapes() {
        let polygons = local_polygons(&shape, 3, 10, 2026);
        for rep in ["signposts", "integers"] {
            let mut bezier_ms = Vec::new();
            let mut cdt_ms = Vec::new();
            let mut stats = (0usize, 0usize, 0usize, 0usize, 0usize);
            for controls in &polygons {
                let opts = BezierOptions::default();
                macro_rules! go {
                    ($tri:expr) => {{
                        let tri = $tri;
                        let t0 = Instant::now();
                        let r = geodesic_bezier(tri, controls, 4, &opts);
                        let t1 = t0.elapsed().as_secs_f64() * 1e3;
                        match r {
                            Ok((mut net, report)) => {
                                bezier_ms.push(t1);
                                stats.0 += report.n_splits;
                                stats.1 += report.flips;
                                stats.2 += report.n_shortenings;
                                let t0 = Instant::now();
                                let cdt = net.make_delaunay(&DelaunayOptions::default()).unwrap();
                                cdt_ms.push(t0.elapsed().as_secs_f64() * 1e3);
                                stats.3 += cdt.delaunay.n_flips;
                                stats.4 += cdt.n_network_edges;
                            }
                            Err(e) => println!("{}: {controls:?}: {e}", shape.name),
                        }
                    }};
                }
                if rep == "signposts" {
                    go!(shape.triangulation());
                } else {
                    go!(IntegerTriangulation::new(&shape.mesh, &shape.geometry).unwrap());
                }
            }
            let med = |v: &mut Vec<f64>| {
                v.sort_by(f64::total_cmp);
                (v[v.len() / 2], v[v.len() - 1])
            };
            // For comparison: plain Delaunay flipping of the input, nothing constrained.
            let t0 = Instant::now();
            let plain = if rep == "signposts" {
                shape
                    .triangulation()
                    .flip_to_delaunay(&DelaunayOptions::default())
                    .unwrap()
            } else {
                IntegerTriangulation::new(&shape.mesh, &shape.geometry)
                    .unwrap()
                    .flip_to_delaunay(&DelaunayOptions::default())
                    .unwrap()
            };
            let plain_ms = t0.elapsed().as_secs_f64() * 1e3;
            println!(
                "{} ({rep}): plain Delaunay flipping of the input (with building it): {plain_ms:.2} \
                 ms, {} flips",
                shape.name, plain.n_flips
            );
            let n = bezier_ms.len();
            let (bm, bx) = med(&mut bezier_ms);
            let (cm, cx) = med(&mut cdt_ms);
            println!(
                "{} ({rep}): {n} curves (d = 3, 4 rounds): Bézier median {bm:.2} ms, max \
                 {bx:.2}; iCDT median {cm:.2} ms, max {cx:.2}; per curve {} splits, {} flips, \
                 {} shortenings; iCDT {} flips, {} constrained",
                shape.name,
                stats.0 / n,
                stats.1 / n,
                stats.2 / n,
                stats.3 / n,
                stats.4 / n
            );
        }
    }
    let load = std::fs::read_to_string("/proc/loadavg").unwrap_or_default();
    println!("load after: {}", load.trim());
}
