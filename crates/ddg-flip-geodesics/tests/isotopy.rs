//! Isotopy (SPEC §7 test 6) and the Δ-complex cases of App. A (§7 layer 3): every run here
//! checks the network (non-crossing included) after every iteration (`check_invariants`),
//! every replacement strictly shortens and the result is certified (the `run` helper).
//! Multiply covered input as in Fig. 8 top (an open path spiralling twice around a tube; two
//! parallel loops given as two paths), a path pulled tight around its own end (Fig. 8
//! bottom), networks with shared junctions, and random paths on triangulations with
//! self-edges and vertices of degree one.

mod common;

use std::f64::consts::PI;

use common::{
    Shape, generated_shapes, long_torus, needle_tetrahedron, path_vertices, random_pairs, run,
    vpath, waist_tube,
};
use ddg_core::generators::icosphere;
use ddg_flip_geodesics::{
    EdgePath, GeodesicNetwork, Network, Segment, ShortenOptions, dijkstra_path,
};
use ddg_intrinsic::{IntrinsicTriangulation, Triangulation};

const N: usize = 10;
const W: usize = 3;

fn v(j: usize, k: usize) -> usize {
    j * N + k % N
}

/// The largest number of segments on one edge.
fn max_cover<T: Triangulation>(net: &GeodesicNetwork<T>) -> usize {
    let n = net.triangulation().complex().n_edges();
    (0..n)
        .map(|e| net.network().layer(e).len())
        .max()
        .unwrap_or(0)
}

/// Fig. 8 top as an open path: once around ring `W − 1`, across the waist, once around ring
/// `W + 1`. Both windings slide towards the waist and end up on the same edges, kept apart by
/// the layering.
#[test]
fn an_open_spiral_winding_twice() {
    let shape = waist_tube(N, W);
    let mut vertices: Vec<usize> = (0..N).map(|k| v(W - 1, k)).collect();
    vertices.push(v(W, 0));
    vertices.extend((1..=N).map(|k| v(W + 1, k)));
    for tri in [shape.triangulation(), shape.delaunay()] {
        let path = vpath(&tri, &vertices, false);
        let initial = path.length(&tri);
        let out = run(tri, &[path], &[], &ShortenOptions::default());
        assert!(out.certificate.geodesic, "{:?}", out.certificate.failures);
        assert!(out.network.length() < initial);
        let cover = max_cover(&out.network);
        assert!(cover >= 2, "the two windings do not share an edge");
        eprintln!(
            "spiral: {initial:.4} -> {:.4}, {} iterations, {} flips, max cover {cover}",
            out.network.length(),
            out.report.iterations,
            out.report.flips
        );
    }
}

/// Fig. 8 top as two loops given as two paths, on rings `W − 2` and `W + 2`: both slide to
/// the waist and cover its edges twice; each has the waist's length.
#[test]
fn two_parallel_loops() {
    let shape = waist_tube(N, W);
    let exact = 2.0 * N as f64 * (PI / N as f64).sin();
    for tri in [shape.triangulation(), shape.delaunay()] {
        let a = vpath(&tri, &(0..N).map(|k| v(W - 2, k)).collect::<Vec<_>>(), true);
        let b = vpath(&tri, &(0..N).map(|k| v(W + 2, k)).collect::<Vec<_>>(), true);
        let out = run(tri, &[a, b], &[], &ShortenOptions::default());
        assert!(out.certificate.geodesic, "{:?}", out.certificate.failures);
        for p in 0..2 {
            let e = (out.network.path_length(p) - exact).abs() / exact;
            assert!(
                e <= 1e-12,
                "loop {p}: {} vs {exact}",
                out.network.path_length(p)
            );
            assert!(path_vertices(&out.network, p).iter().all(|&x| x / N == W));
        }
        assert_eq!(max_cover(&out.network), 2);
    }
}

/// Fig. 8 bottom, on a flat annulus (hole of radius 0.4): from `s` (radius 0.7, angle 0)
/// once counter-clockwise around the hole, then on to `e` (radius 0.9, angle 30°). Pulled
/// taut, the last stretch (from the hole to `e`) would cross the first (from `s` to the hole)
/// transversally, which the isotopy forbids: the path is pulled tight around its own end `s`,
/// passing through it (a joint at a pinned vertex) or past it with a wedge blocked by the end.
#[test]
fn a_path_pulled_tight_around_its_own_end() {
    let (n_around, n_radial) = (24, 6);
    let shape = Shape::generated(
        "annulus 24x6",
        ddg_core::generators::annulus(n_around, n_radial, 0.4, 1.0).unwrap(),
    );
    let at = |i: usize, j: usize| (i % n_around) * (n_radial + 1) + j;
    let (s, e) = (at(0, 3), at(2, 1));
    let mut accepted = 0;
    for tri in [shape.triangulation(), shape.delaunay()] {
        for via in [at(1, 2), at(0, 2), at(1, 3)] {
            let around: Vec<usize> = (0..n_around).map(|i| at(i, 3)).collect();
            let mut path = vpath(&tri, &around, false);
            for (a, b) in [(at(n_around - 1, 3), via), (via, e)] {
                if a != b {
                    path.segments
                        .extend(dijkstra_path(&tri, a, b).unwrap().segments);
                }
            }
            if Network::new(tri.complex(), std::slice::from_ref(&path), &[], None).is_err() {
                continue;
            }
            accepted += 1;
            let initial = path.length(&tri);
            let out = run(tri.clone(), &[path], &[], &ShortenOptions::default());
            assert!(out.certificate.geodesic, "{:?}", out.certificate.failures);
            let caught = out.certificate.n_blocked + out.certificate.n_pinned;
            assert!(caught >= 1, "not caught on its end: {:?}", out.certificate);
            assert!(out.network.length() < initial);
            let vs = path_vertices(&out.network, 0);
            assert_eq!((vs[0], *vs.last().unwrap()), (s, e));
            eprintln!(
                "around its end: {initial:.4} -> {:.4}, {} joints at the end vertex, {} blocked \
                 by an end, {} iterations",
                out.network.length(),
                out.certificate.n_pinned,
                out.certificate.n_blocked,
                out.report.iterations
            );
        }
    }
    assert!(
        accepted >= 2,
        "only {accepted} initial paths are non-crossing"
    );
}

/// Networks: Dijkstra paths from one junction to several targets (the junction and the
/// targets are path ends, so pinned); shortened together in one queue with one layering.
#[test]
fn networks_with_a_junction() {
    let mut runs = 0;
    for shape in [
        Shape::generated("icosphere 3", icosphere(3).unwrap()),
        waist_tube(12, 4),
        common::grid(12, 12),
    ] {
        for tri in [shape.triangulation(), shape.delaunay()] {
            let n = tri.complex().n_vertices();
            for (k, (hub, _)) in random_pairs(n, 10, 3).into_iter().enumerate() {
                let targets: Vec<usize> = random_pairs(n, 5, 100 + k as u64)
                    .into_iter()
                    .map(|(t, _)| t)
                    .filter(|&t| t != hub)
                    .collect();
                let paths: Vec<EdgePath> = targets
                    .iter()
                    .map(|&t| dijkstra_path(&tri, hub, t).unwrap())
                    .collect();
                if Network::new(tri.complex(), &paths, &[], None).is_err() {
                    continue;
                }
                let initial: f64 = paths.iter().map(|p| p.length(&tri)).sum();
                let out = run(tri.clone(), &paths, &[], &ShortenOptions::default());
                assert!(
                    out.certificate.geodesic,
                    "{}: {:?}",
                    shape.name, out.certificate.failures
                );
                assert!(out.network.length() <= initial);
                for (p, &t) in targets.iter().enumerate() {
                    let vs = path_vertices(&out.network, p);
                    assert_eq!((vs[0], *vs.last().unwrap()), (hub, t), "ends are kept");
                }
                runs += 1;
            }
        }
    }
    assert!(runs >= 8, "only {runs} networks were non-crossing");
}

/// An icosphere with vertex 0 pushed out by `factor`: a cone point whose angle sum drops below
/// `π` for large enough factors.
fn spike(factor: f64) -> Shape {
    let mut g = icosphere(2).unwrap();
    for c in g.positions[0].iter_mut() {
        *c *= factor;
    }
    Shape::generated(&format!("spike {factor}"), g)
}

/// Flips edges at `v` until it has degree one; `false` if no edge at `v` can be flipped
/// before that.
fn flip_to_degree_one(tri: &mut IntrinsicTriangulation, v: usize) -> bool {
    loop {
        let c = tri.complex();
        if c.vertex_degree(v) <= 1 {
            return true;
        }
        let edges: Vec<usize> = c.outgoing_halfedges(v).map(|h| c.edge(h)).collect();
        if !edges.into_iter().any(|e| tri.flip_edge(e).is_ok()) {
            return false;
        }
    }
}

/// Random paths and networks on triangulations with self-edges and vertices of degree one:
/// Delaunay triangulations of a needle and of thin tori, after random flips, and spikes
/// whose tip (angle sum between `0.4π` and `π`) is flipped down to degree one, where a self-edge
/// at its neighbour `x` encloses it; there, joints `u → x → w` through every pair of
/// neighbours of `x` are shortened too, and those whose wedge holds the self-edge meet App. A's
/// cases. Every run keeps the network valid after every iteration and certifies; the cases are
/// counted and must occur.
#[test]
fn delta_complex_cases_are_exercised() {
    let shapes = [
        needle_tetrahedron(8.0, 4.0),
        needle_tetrahedron(20.0, 9.0),
        long_torus(8, 3, 4.0, 0.2, 0.37),
        long_torus(12, 3, 6.0, 0.15, 0.21),
        long_torus(6, 3, 5.0, 0.1, 0.5),
    ];
    let mut cases: Vec<(String, IntrinsicTriangulation, Vec<Vec<EdgePath>>)> = Vec::new();
    for shape in &shapes {
        let mut tris: Vec<IntrinsicTriangulation> = vec![shape.delaunay()];
        for seed in 0..4 {
            let mut t = shape.delaunay();
            let n_flips = 4 * t.complex().n_edges();
            common::random_flips(&mut t, n_flips, seed);
            tris.push(t);
        }
        for tri in tris {
            let n = tri.complex().n_vertices();
            let mut inputs: Vec<Vec<EdgePath>> = Vec::new();
            for a in 0..n {
                for b in 0..n {
                    if a != b {
                        inputs.push(vec![dijkstra_path(&tri, a, b).unwrap()]);
                    }
                }
            }
            for &s in &tri.irregularities().self_edges {
                for forward in [true, false] {
                    inputs.push(vec![EdgePath::closed(vec![Segment::new(
                        s as usize, forward,
                    )])]);
                }
            }
            cases.push((shape.name.clone(), tri, inputs));
        }
    }
    let mut n_spikes = 0;
    for factor in [1.6, 1.8, 2.0, 2.3, 2.7, 3.2] {
        let shape = spike(factor);
        let mut tri = shape.triangulation();
        let theta = tri.vertex_angle_sum(0);
        if !(theta > 0.4 * PI && theta < PI) || !flip_to_degree_one(&mut tri, 0) {
            continue;
        }
        n_spikes += 1;
        let c = tri.complex();
        let x = c.to_vertex(c.outgoing_halfedge(0));
        let mut ring: Vec<usize> = c.outgoing_halfedges(x).map(|h| c.to_vertex(h)).collect();
        ring.sort_unstable();
        ring.dedup();
        ring.retain(|&u| u != 0 && u != x);
        let mut inputs = Vec::new();
        for &u in &ring {
            for &w in &ring {
                if u != w {
                    inputs.push(vec![vpath(&tri, &[u, x, w], false)]);
                }
            }
        }
        for (a, b) in random_pairs(c.n_vertices(), 60, factor.to_bits()) {
            inputs.push(vec![dijkstra_path(&tri, a, b).unwrap()]);
        }
        cases.push((shape.name.clone(), tri, inputs));
    }
    assert!(n_spikes >= 2, "only {n_spikes} spikes reached degree one");

    let (mut twice, mut degree_one, mut contractions, mut moves, mut runs) = (0, 0, 0, 0, 0);
    for (name, tri, inputs) in &cases {
        for paths in inputs {
            if Network::new(tri.complex(), paths, &[], None).is_err() {
                continue;
            }
            let out = run(tri.clone(), paths, &[], &ShortenOptions::default());
            assert!(
                out.certificate.geodesic,
                "{name}: {:?}",
                out.certificate.failures
            );
            twice += out.report.twice_incident_wedges;
            degree_one += out.report.degree_one_wedges;
            contractions += out.report.contractions;
            moves += out.report.loop_perturbations;
            runs += 1;
        }
    }
    eprintln!(
        "Δ-complex cases over {runs} runs ({n_spikes} spikes): {twice} wedges with a \
         twice-incident edge, {degree_one} outer arcs through a degree-one vertex, \
         {contractions} out-and-back contractions, {moves} App. B moves"
    );
    assert!(twice > 0 && degree_one > 0 && contractions > 0 && moves > 0);
}

/// Every generated shape: random Dijkstra paths with the invariants checked after every
/// iteration.
#[test]
fn random_paths_keep_the_invariants() {
    let mut runs = 0;
    for shape in generated_shapes() {
        for tri in [shape.triangulation(), shape.delaunay()] {
            for (a, b) in random_pairs(tri.complex().n_vertices(), 8, 21) {
                let path = dijkstra_path(&tri, a, b).unwrap();
                let out = run(tri.clone(), &[path], &[], &ShortenOptions::default());
                assert!(
                    out.certificate.geodesic,
                    "{}: {a} -> {b}: {:?}",
                    shape.name, out.certificate.failures
                );
                runs += 1;
            }
        }
    }
    eprintln!("{runs} random paths on the generated shapes");
}
