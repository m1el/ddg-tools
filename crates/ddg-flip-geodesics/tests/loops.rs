//! Loops (SPEC §4.3, App. B) and the stop rules (F13, F14): a loop slides to the waist of a
//! tube; a loop around a cone vertex contracts to a point and is reported (F12); a
//! single-segment loop around a vertex of degree one is moved across its face (App. B) and
//! then contracts; loops around a torus stay and become geodesic; early stopping.

mod common;

use std::f64::consts::PI;

use common::{long_torus, needle_tetrahedron, path_vertices, pyramid, run, vpath, waist_tube};
use ddg_core::generators::{icosphere, torus};
use ddg_flip_geodesics::{
    EdgePath, FlipError, GeodesicNetwork, Segment, ShortenOptions, ShortenStatus,
};
use ddg_intrinsic::{IntrinsicTriangulation, Triangulation};

/// The link of `v` (its neighbours counter-clockwise) as a closed vertex sequence.
fn link<T: Triangulation>(tri: &T, v: usize) -> Vec<usize> {
    let c = tri.complex();
    c.outgoing_halfedges(v).map(|h| c.to_vertex(h)).collect()
}

/// On a tube whose waist is ring `w`, a loop on another ring, plain or zig-zag, slides to the
/// waist: length `2 n sin(π/n)` (radius 1 there), every vertex on ring `w`.
#[test]
fn a_loop_slides_to_the_waist() {
    let (n, w) = (10, 3);
    let shape = waist_tube(n, w);
    let v = |j: usize, k: usize| j * n + k % n;
    let exact = 2.0 * n as f64 * (PI / n as f64).sin();
    let plain: Vec<usize> = (0..n).map(|k| v(w + 2, k)).collect();
    let low: Vec<usize> = (0..n).map(|k| v(w - 3, k)).collect();
    let mut zigzag = Vec::new();
    for k in (0..n).step_by(2) {
        zigzag.extend([v(w + 1, k), v(w + 2, k + 1), v(w + 1, k + 1)]);
    }
    for tri in [shape.triangulation(), shape.delaunay()] {
        for vertices in [&plain, &low, &zigzag] {
            let path = vpath(&tri, vertices, true);
            let out = run(tri.clone(), &[path], &[], &ShortenOptions::default());
            assert_eq!(out.report.status, ShortenStatus::Geodesic);
            assert!(out.certificate.geodesic, "{:?}", out.certificate.failures);
            let e = (out.network.length() - exact).abs() / exact;
            assert!(e <= 1e-12, "{} vs {exact} ({e:e})", out.network.length());
            let rings: Vec<usize> = path_vertices(&out.network, 0)
                .iter()
                .map(|&x| x / n)
                .collect();
            assert!(rings.iter().all(|&r| r == w), "not on the waist: {rings:?}");
        }
    }
}

/// A loop around a vertex of a closed convex surface contracts to a point (Fig. 20): the path
/// has no segment left, it is reported in `contracted`, and the network stays valid.
#[test]
fn a_loop_around_a_cone_vertex_contracts() {
    let sphere = common::Shape::generated("icosphere 2", icosphere(2).unwrap());
    let pyr = pyramid(0.5);
    let mut moves = 0;
    let mut contractions = 0;
    for (shape, vertices) in [(&sphere, vec![0, 7, 40]), (&pyr, vec![4, 0])] {
        for tri in [shape.triangulation(), shape.delaunay()] {
            for &v in &vertices {
                let path = vpath(&tri, &link(&tri, v), true);
                let out = run(tri.clone(), &[path], &[], &ShortenOptions::default());
                assert_eq!(
                    out.report.status,
                    ShortenStatus::Geodesic,
                    "{} vertex {v}",
                    shape.name
                );
                assert_eq!(out.report.contracted, vec![0], "{} vertex {v}", shape.name);
                assert_eq!(out.network.network().n_path_segments(0), 0);
                assert_eq!(out.network.length(), 0.0);
                assert!(out.certificate.geodesic);
                moves += out.report.loop_perturbations;
                contractions += out.report.contractions;
                assert_eq!(out.loop_moves, out.report.loop_perturbations);
            }
        }
    }
    assert!(contractions >= 10);
    eprintln!("contracting loops: {moves} App. B moves, {contractions} contractions");
}

/// The Delaunay triangulation of a tilted needle has a vertex `d` of degree one inside a
/// self-edge at `x`. A loop along the self-edge is a single-segment loop: on the side of `d`
/// its wedge is below `π`, it is moved across the face `(x, x, d)` (App. B), which makes it
/// out and back to `d`, and that contracts.
#[test]
fn a_single_segment_loop_around_a_degree_one_vertex_contracts() {
    let shape = needle_tetrahedron(8.0, 4.0);
    let tri = shape.delaunay();
    let irregular = tri.irregularities();
    let s = irregular.self_edges[0] as usize;
    for forward in [true, false] {
        let lp = EdgePath::closed(vec![Segment::new(s, forward)]);
        let out = run(tri.clone(), &[lp], &[], &ShortenOptions::default());
        assert_eq!(out.report.contracted, vec![0]);
        assert!(out.report.loop_perturbations >= 1, "{:?}", out.report);
        assert!(out.report.contractions >= 1, "{:?}", out.report);
        eprintln!("needle loop (forward {forward}): {:?}", out.report);
    }
}

/// Loops around a torus (the short way and the long way, plain and through the Delaunay
/// triangulation with self-edges) cannot contract: they become geodesic and stay.
#[test]
fn loops_around_a_torus_stay() {
    let shapes = [
        common::Shape::generated("torus 16x8", torus(16, 8, 1.0, 0.35).unwrap()),
        long_torus(8, 3, 4.0, 0.2, 0.37),
        long_torus(24, 6, 2.0, 0.5, 0.0),
    ];
    for shape in &shapes {
        for tri in [shape.triangulation(), shape.delaunay()] {
            let c = tri.complex();
            // Two loops by Dijkstra-free construction: walk along "ring" edges from vertex 0,
            // the first and the last outgoing halfedge of every vertex, until back at 0.
            for pick_last in [false, true] {
                let mut vertices = vec![0usize];
                let mut seen = vec![false; c.n_vertices()];
                seen[0] = true;
                loop {
                    let v = *vertices.last().unwrap();
                    let mut hs = c.outgoing_halfedges(v);
                    let h = if pick_last {
                        hs.last().unwrap()
                    } else {
                        hs.next().unwrap()
                    };
                    let u = c.to_vertex(h);
                    if u == 0 || seen[u] {
                        break;
                    }
                    seen[u] = true;
                    vertices.push(u);
                }
                if vertices.len() < 3 {
                    continue;
                }
                let Ok(path) = ddg_flip_geodesics::path_from_vertices(&tri, &vertices, true) else {
                    continue;
                };
                let Ok(_) =
                    GeodesicNetwork::new(tri.clone(), std::slice::from_ref(&path), &[], None)
                else {
                    continue;
                };
                let initial = path.length(&tri);
                let out = run(tri.clone(), &[path], &[], &ShortenOptions::default());
                assert!(
                    out.certificate.geodesic,
                    "{}: {:?}",
                    shape.name, out.certificate.failures
                );
                assert!(out.network.length() <= initial);
                eprintln!(
                    "{}: loop of {} vertices, {initial:.4} -> {:.4}, contracted {:?}, {:?}",
                    shape.name,
                    vertices.len(),
                    out.network.length(),
                    out.report.contracted,
                    out.report.status
                );
            }
        }
    }
}

/// Stop rules (F13, F14) and refused options.
#[test]
fn stop_rules_and_refused_options() {
    let (n, w) = (12, 4);
    let shape = waist_tube(n, w);
    let tri: IntrinsicTriangulation = shape.triangulation();
    let ring: Vec<usize> = (0..n).map(|k| (2 * w) * n + k).collect();
    let path = vpath(&tri, &ring, true);
    let initial = path.length(&tri);
    let net = GeodesicNetwork::new(tri, &[path], &[], None).unwrap();

    let mut a = net.clone();
    let r = a
        .shorten(&ShortenOptions {
            max_iterations: Some(3),
            ..Default::default()
        })
        .unwrap();
    assert_eq!((r.status, r.iterations), (ShortenStatus::IterationLimit, 3));

    let mut b = net.clone();
    let r = b
        .shorten(&ShortenOptions {
            stop_length_ratio: Some(0.95),
            ..Default::default()
        })
        .unwrap();
    assert_eq!(r.status, ShortenStatus::LengthRatio);
    assert!(b.length() < 0.95 * initial && r.final_length == b.length());
    assert!(b.recompute_length() < 0.95 * initial);

    // A zig-zag loop: its tips have small wedges, the rest are near pi.
    {
        let tri = net.triangulation();
        let mut zigzag = Vec::new();
        for k in (0..n).step_by(2) {
            zigzag.extend([(w + 1) * n + k, (w + 2) * n + k + 1, (w + 1) * n + k + 1]);
        }
        let path = vpath(tri, &zigzag, true);
        let mut c = GeodesicNetwork::new(tri.clone(), &[path], &[], None).unwrap();
        let r = c
            .shorten(&ShortenOptions {
                stop_min_angle: Some(2.0),
                ..Default::default()
            })
            .unwrap();
        assert_eq!(r.status, ShortenStatus::MinAngle);
        assert!(r.iterations > 0);
        let full = c.shorten(&ShortenOptions::default()).unwrap();
        assert_eq!(full.status, ShortenStatus::Geodesic);
        assert!(full.iterations > 0);
    }

    let mut d = net.clone();
    let r = d.shorten(&ShortenOptions::default()).unwrap();
    assert_eq!(r.status, ShortenStatus::Geodesic);
    assert!(r.iterations > 3);

    // Refused options change nothing.
    for bad in [
        ShortenOptions {
            straight_tolerance: -1.0,
            ..Default::default()
        },
        ShortenOptions {
            stop_length_ratio: Some(f64::NAN),
            ..Default::default()
        },
        ShortenOptions {
            stop_min_angle: Some(f64::INFINITY),
            ..Default::default()
        },
    ] {
        let mut e = net.clone();
        assert!(matches!(
            e.shorten(&bad),
            Err(FlipError::InvalidOptions { .. })
        ));
        assert_eq!(e.paths(), net.paths());
        assert_eq!(
            e.triangulation().edge_lengths(),
            net.triangulation().edge_lengths()
        );
        assert_eq!(e.length().to_bits(), net.length().to_bits());
    }
}
