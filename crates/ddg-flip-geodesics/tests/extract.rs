//! Initial paths (SPEC §2.2) and extraction on the input mesh (SPEC §5, §7 layer 1b, layer 2,
//! tests 9 and 10): waypoint paths and edge sets; the polylines of FlipOut's results on signposts
//! and on integer coordinates, the certificate on `T⁰` for both, the two compared; determinism
//! at one and many threads. The evaluation (ignored, release) prints the timings.

mod common;

use std::collections::{BTreeSet, VecDeque};
use std::time::Instant;

use common::*;
use ddg_flip_geodesics::{
    CertifyOptions, EdgePath, ExtractOptions, Extraction, FlipError, GeodesicNetwork, Network,
    Segment, ShortenOptions, ShortenStatus, crossing_sines, dijkstra_path, network_from_edges,
    polyline_3d, waypoint_path,
};
use ddg_integer_coords::IntegerTriangulation;
use ddg_intrinsic::{Complex, IntrinsicTriangulation, SurfacePoint, Triangulation};

// ===== initial paths =====

#[test]
fn waypoint_paths_are_the_legs_concatenated() {
    let mut shapes = generated_shapes();
    shapes.extend(data_shapes());
    for shape in &shapes {
        let tri = shape.triangulation();
        let c = tri.complex();
        let n = c.n_vertices();
        for (k, w) in random_pairs(n, 8, 5).chunks(2).enumerate() {
            let pts = [w[0].0, w[0].1, w[1].0, w[1].1];
            for closed in [false, true] {
                // Consecutive equal waypoints (or a repeated first) are an error; skip them.
                let cyc = |i: usize| (pts[i], pts[(i + 1) % 4]);
                let legs = if closed { 4 } else { 3 };
                if (0..legs).any(|i| cyc(i).0 == cyc(i).1) {
                    continue;
                }
                let Ok(legs_paths) = (0..legs)
                    .map(|i| dijkstra_path(&tri, cyc(i).0, cyc(i).1))
                    .collect::<Result<Vec<_>, _>>()
                else {
                    continue;
                };
                for pin in [false, true] {
                    let got = waypoint_path(&tri, &pts, closed, pin).unwrap();
                    assert_eq!(got.paths.len(), 1);
                    let want: Vec<Segment> =
                        legs_paths.iter().flat_map(|p| p.segments.clone()).collect();
                    assert_eq!(got.paths[0].segments, want, "{} #{k}", shape.name);
                    assert_eq!(got.paths[0].closed, closed);
                    let mut expect: Vec<usize> = if pin { pts.to_vec() } else { vec![] };
                    expect.sort_unstable();
                    expect.dedup();
                    assert_eq!(got.pinned, expect);
                    let vs = got.paths[0].vertices(c);
                    assert_eq!(vs[0], pts[0]);
                    if !closed {
                        assert_eq!(*vs.last().unwrap(), pts[3]);
                    }
                }
            }
        }
    }
}

#[test]
fn waypoint_errors_name_the_leg() {
    let shape = grid(4, 4);
    let tri = shape.triangulation();
    let e = waypoint_path(&tri, &[3], false, true).unwrap_err();
    assert!(matches!(e, FlipError::TooFewVertices { n: 1, .. }), "{e}");
    let e = waypoint_path(&tri, &[0, 6, 6, 24], false, true).unwrap_err();
    assert!(
        matches!(e, FlipError::WaypointLeg { leg: 1, from: 6, to: 6, .. }),
        "{e:?}"
    );
    let msg = e.to_string();
    assert!(msg.contains("leg 1") && msg.contains("vertex 6"), "{msg}");
    let e = waypoint_path(&tri, &[0, 99], false, false).unwrap_err();
    assert!(matches!(e, FlipError::WaypointLeg { leg: 0, .. }));
    assert!(e.to_string().contains("99"), "{e}");
    // A closed path closes with the leg last -> first.
    let e = waypoint_path(&tri, &[0, 6, 0], true, false).unwrap_err();
    assert!(matches!(e, FlipError::WaypointLeg { leg: 2, .. }), "{e:?}");
    // Two components.
    let two = Shape::new(
        "two triangles",
        vec![
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [5.0, 0.0, 0.0],
            [6.0, 0.0, 0.0],
            [5.0, 1.0, 0.0],
        ],
        vec![[0, 1, 2], [3, 4, 5]],
    );
    let tri = two.triangulation();
    let e = waypoint_path(&tri, &[0, 1, 4], false, false).unwrap_err();
    assert!(e.to_string().contains("not connected"), "{e}");
    assert!(e.to_string().contains("leg 1"), "{e}");
}

/// A tree-cotree cut graph of a closed surface: the edges outside a spanning tree of the dual
/// graph that avoids a spanning tree of the vertex graph, with leaves pruned.
fn cut_graph(c: &Complex) -> Vec<bool> {
    let (nv, ne, nf) = (c.n_vertices(), c.n_edges(), c.n_faces());
    let mut tree = vec![false; ne];
    let mut seen = vec![false; nv];
    let mut queue = VecDeque::from([0usize]);
    seen[0] = true;
    while let Some(v) = queue.pop_front() {
        for h in c.outgoing_halfedges(v) {
            let w = c.to_vertex(h);
            if !seen[w] {
                seen[w] = true;
                tree[c.edge(h)] = true;
                queue.push_back(w);
            }
        }
    }
    let mut dual = vec![false; ne];
    let mut fseen = vec![false; nf];
    let mut queue = VecDeque::from([0usize]);
    fseen[0] = true;
    while let Some(f) = queue.pop_front() {
        for h in 3 * f..3 * f + 3 {
            if let Some(t) = c.twin(h)
                && !tree[c.edge(h)]
                && !fseen[c.face(t)]
            {
                fseen[c.face(t)] = true;
                dual[c.edge(h)] = true;
                queue.push_back(c.face(t));
            }
        }
    }
    let mut mask: Vec<bool> = (0..ne).map(|e| !dual[e]).collect();
    prune_leaves(c, &mut mask);
    mask
}

fn degrees(c: &Complex, mask: &[bool]) -> Vec<usize> {
    let mut d = vec![0; c.n_vertices()];
    for e in (0..c.n_edges()).filter(|&e| mask[e]) {
        let [a, b] = c.edge_vertices(e);
        d[a] += 1;
        d[b] += 1;
    }
    d
}

fn prune_leaves(c: &Complex, mask: &mut [bool]) {
    loop {
        let d = degrees(c, mask);
        let mut changed = false;
        for e in 0..c.n_edges() {
            let [a, b] = c.edge_vertices(e);
            if mask[e] && (d[a] == 1 || d[b] == 1) {
                mask[e] = false;
                changed = true;
            }
        }
        if !changed {
            return;
        }
    }
}

/// The breadth-first spanning tree from vertex 0 (a seam network with leaves and junctions).
fn bfs_tree(c: &Complex) -> Vec<bool> {
    let mut tree = vec![false; c.n_edges()];
    let mut seen = vec![false; c.n_vertices()];
    let mut queue = VecDeque::from([0usize]);
    seen[0] = true;
    while let Some(v) = queue.pop_front() {
        for h in c.outgoing_halfedges(v) {
            let w = c.to_vertex(h);
            if !seen[w] {
                seen[w] = true;
                tree[c.edge(h)] = true;
                queue.push_back(w);
            }
        }
    }
    tree
}

/// What `network_from_edges` must satisfy, checked against the mask alone.
fn check_edge_set(tri: &IntrinsicTriangulation, mask: &[bool], extra: &[usize], what: &str) {
    let c = tri.complex();
    let out = network_from_edges(tri, mask, extra).unwrap_or_else(|e| panic!("{what}: {e}"));
    let deg = degrees(c, mask);
    let is_end = |v: usize| extra.contains(&v) || (deg[v] != 0 && deg[v] != 2);
    // Every masked edge exactly once.
    let mut count = vec![0; c.n_edges()];
    for p in &out.paths {
        for s in &p.segments {
            count[s.edge as usize] += 1;
        }
    }
    for e in 0..c.n_edges() {
        assert_eq!(count[e], usize::from(mask[e]), "{what}: edge {e}");
    }
    // Paths: open ones run between ends through vertices of degree two that are not ends;
    // closed ones consist of such vertices.
    let mut open = 0;
    for (i, p) in out.paths.iter().enumerate() {
        let vs = p.vertices(c);
        if p.closed {
            assert!(vs.iter().all(|&v| !is_end(v)), "{what}: loop {i}");
            assert_eq!(p.segments[0].tail(c), p.segments.last().unwrap().head(c));
        } else {
            open += 1;
            assert!(is_end(vs[0]) && is_end(*vs.last().unwrap()), "{what}: path {i}");
            assert!(vs[1..vs.len() - 1].iter().all(|&v| !is_end(v)), "{what}: path {i}");
        }
        for w in p.segments.windows(2) {
            assert_eq!(w[0].head(c), w[1].tail(c));
        }
    }
    let end_incidences: usize = (0..c.n_vertices()).filter(|&v| is_end(v)).map(|v| deg[v]).sum();
    assert_eq!(2 * open, end_incidences, "{what}: number of open paths");
    let pinned: Vec<usize> = (0..c.n_vertices()).filter(|&v| is_end(v) && deg[v] != 0).collect();
    assert_eq!(out.pinned, pinned, "{what}");
    // A valid network: ends pair nothing, other joints are at degree-two vertices.
    Network::new(c, &out.paths, &out.pinned, None).unwrap_or_else(|e| panic!("{what}: {e}"));
}

#[test]
fn edge_sets_become_paths_and_loops() {
    let mut shapes = generated_shapes();
    shapes.push(waist_tube(10, 3));
    shapes.extend(data_shapes());
    let mut rng = ddg_core::generators::SplitMix64::new(3);
    for shape in &shapes {
        let tri = shape.triangulation();
        let c = tri.complex();
        let closed_surface = (0..c.n_edges()).all(|e| !c.is_boundary_edge(e));
        // Trees and cut graphs.
        check_edge_set(&tri, &bfs_tree(c), &[], &format!("{} tree", shape.name));
        if closed_surface && c.euler_characteristic() == 0 {
            let cut = cut_graph(c);
            assert!(cut.iter().any(|&b| b), "{}: empty cut graph", shape.name);
            check_edge_set(&tri, &cut, &[], &format!("{} cut graph", shape.name));
            // Extra pinned vertices on it split its loops.
            let on: Vec<usize> = (0..c.n_edges()).filter(|&e| cut[e]).map(|e| c.edge_vertices(e)[0]).collect();
            check_edge_set(&tri, &cut, &on[..on.len().min(3)], &format!("{} cut graph pinned", shape.name));
        }
        // Random subsets, sparse and dense (every degree occurs).
        for (k, p) in [0.05, 0.3, 0.7].into_iter().enumerate() {
            let mask: Vec<bool> = (0..c.n_edges())
                .map(|_| (rng.next_u64() % 1000) as f64 / 1000.0 < p)
                .collect();
            check_edge_set(&tri, &mask, &[], &format!("{} random {k}", shape.name));
        }
        // The boundary, if any: closed loops along boundary edges walked against their halfedge.
        let boundary: Vec<bool> = (0..c.n_edges()).map(|e| c.is_boundary_edge(e)).collect();
        if boundary.iter().any(|&b| b) {
            check_edge_set(&tri, &boundary, &[], &format!("{} boundary", shape.name));
            let one = (0..c.n_edges()).find(|&e| boundary[e]).unwrap();
            let v = c.edge_vertices(one)[0];
            check_edge_set(&tri, &boundary, &[v], &format!("{} boundary pinned", shape.name));
        }
    }
}

#[test]
fn edge_set_of_a_path_or_loop_is_that_path_or_loop() {
    for shape in generated_shapes().iter().chain(data_shapes().iter()) {
        let tri = shape.triangulation();
        let c = tri.complex();
        for (a, b) in random_pairs(c.n_vertices(), 6, 8) {
            let path = dijkstra_path(&tri, a, b).unwrap();
            let mut mask = vec![false; c.n_edges()];
            for s in &path.segments {
                mask[s.edge as usize] = true;
            }
            let got = network_from_edges(&tri, &mask, &[]).unwrap();
            // A shortest path is simple: one open path, the same edges, either direction.
            assert_eq!(got.paths.len(), 1, "{}", shape.name);
            let g = &got.paths[0];
            assert!(!g.closed);
            assert!(*g == path || *g == path.reversed(), "{}: {g:?} vs {path:?}", shape.name);
            assert_eq!(got.pinned.len(), 2);
        }
    }
    // A ring of a tube is a closed loop; the same edges, some rotation, either direction.
    let (n, j) = (10, 4);
    let tube = waist_tube(n, 3);
    let tri = tube.triangulation();
    let ring = vpath(&tri, &(0..n).map(|k| j * n + k).collect::<Vec<_>>(), true);
    let mut mask = vec![false; tri.complex().n_edges()];
    for s in &ring.segments {
        mask[s.edge as usize] = true;
    }
    let got = network_from_edges(&tri, &mask, &[]).unwrap();
    assert_eq!(got.paths.len(), 1);
    assert!(got.paths[0].closed && got.pinned.is_empty());
    let a: BTreeSet<u32> = got.paths[0].segments.iter().map(|s| s.edge).collect();
    let b: BTreeSet<u32> = ring.segments.iter().map(|s| s.edge).collect();
    assert_eq!(a, b);
}

#[test]
fn edge_set_errors_name_what_and_where() {
    let shape = grid(3, 3);
    let tri = shape.triangulation();
    let ne = tri.complex().n_edges();
    let e = network_from_edges(&tri, &vec![false; ne - 1], &[]).unwrap_err();
    assert!(e.to_string().contains(&format!("{} entries", ne - 1)), "{e}");
    assert!(e.to_string().contains(&format!("{ne} edges")), "{e}");
    let e = network_from_edges(&tri, &vec![true; ne], &[1000]).unwrap_err();
    assert!(e.to_string().contains("1000"), "{e}");
}

// ===== extraction =====

fn shorten<T: Triangulation + Clone + Sync>(
    tri: T,
    paths: &[EdgePath],
    pinned: &[usize],
    label: &str,
) -> GeodesicNetwork<T> {
    let mut net = GeodesicNetwork::new(tri, paths, pinned, None)
        .unwrap_or_else(|e| panic!("{label}: network: {e}"));
    let report = net
        .shorten(&ShortenOptions::default())
        .unwrap_or_else(|e| panic!("{label}: shorten: {e}"));
    assert_eq!(report.status, ShortenStatus::Geodesic, "{label}");
    let cert = net.certify(&CertifyOptions::default()).unwrap();
    assert!(cert.geodesic, "{label}: {:?}", cert.failures);
    net
}

/// Worst values over a test.
#[derive(Default, Debug)]
struct Stats {
    paths: usize,
    pieces: usize,
    crossings: usize,
    compared_points: usize,
    /// Trials where the signpost extraction failed: `label`.
    signpost_failures: Vec<String>,
    worst_cross_dev: [f64; 2],
    worst_joint: [f64; 2],
    worst_wedge: [f64; 2],
    worst_length: [f64; 2],
    worst_raw: f64,
    worst_conditioned: f64,
    worst_mismatch: [f64; 2],
    max_clamp: f64,
    joints_checked: usize,
    joints_exempt: usize,
}

impl Stats {
    fn note(&mut self, i: usize, x: &Extraction) {
        let c = x.certificate.as_ref().unwrap();
        self.worst_cross_dev[i] = self.worst_cross_dev[i].max(c.max_crossing_deviation);
        if c.min_joint_angle.is_finite() {
            self.worst_joint[i] =
                self.worst_joint[i].max(std::f64::consts::PI - c.min_joint_angle);
        }
        self.joints_checked += c.n_joints;
        self.joints_exempt += c.n_exempt;
        self.worst_wedge[i] = self.worst_wedge[i].max(c.max_wedge_mismatch);
        self.worst_length[i] = self.worst_length[i].max(c.max_length_error);
        self.worst_mismatch[i] = self.worst_mismatch[i].max(x.report.max_relative_mismatch);
        if i == 1 {
            self.max_clamp = self.max_clamp.max(x.report.max_clamp);
        }
    }
}

fn require(x: &Extraction, label: &str) {
    let c = x.certificate.as_ref().unwrap();
    assert!(c.passed, "{label}: certificate on T0 failed: {c:#?}");
}

/// Runs `paths` on both representations of `shape` and checks everything (module doc).
fn check_both(
    shape: &Shape,
    integers: &IntegerTriangulation,
    label: &str,
    paths: &[EdgePath],
    pinned: &[usize],
    stats: &mut Stats,
) {
    let signposts = shape.triangulation();
    if Network::new(signposts.complex(), paths, pinned, None).is_err() {
        return;
    }
    let a = shorten(signposts, paths, pinned, label);
    let b = shorten(integers.clone(), paths, pinned, label);
    assert_eq!(a.paths(), b.paths(), "{label}: FlipOut differs");
    let options = ExtractOptions::default();
    // Integers: exact; must always succeed and certify.
    let xb = b
        .extract_with_report(&options)
        .unwrap_or_else(|e| panic!("{label}: integer extraction: {e}"));
    require(&xb, &format!("{label} (integers)"));
    stats.note(1, &xb);
    stats.paths += paths.len();
    stats.pieces += xb.report.n_pieces;
    stats.crossings += xb.report.n_crossings;
    // Signposts: may fail (IT N8); where they do not, the same checks and the comparison.
    let xa = match a.extract_with_report(&options) {
        Ok(x) => x,
        Err(e) => {
            assert!(matches!(e, FlipError::Extraction { .. }), "{label}: {e}");
            eprintln!("SIGNPOST TRACE FAILURE: {label}: {e}");
            stats.signpost_failures.push(label.to_string());
            return;
        }
    };
    require(&xa, &format!("{label} (signposts)"));
    stats.note(0, &xa);
    let lengths = integers.input().geometry().edge_lengths();
    for (p, (pa, pb)) in xa.paths.iter().zip(&xb.paths).enumerate() {
        assert_eq!(pa.points.len(), pb.points.len(), "{label}: path {p} points");
        assert_eq!(pa.faces.len(), pb.faces.len(), "{label}: path {p} faces");
        let sines = crossing_sines(integers, pb);
        for (i, (qa, qb)) in pa.points.iter().zip(&pb.points).enumerate() {
            match (qa, qb) {
                (SurfacePoint::Vertex(u), SurfacePoint::Vertex(v)) => {
                    assert_eq!(u, v, "{label}: path {p} point {i}")
                }
                (
                    SurfacePoint::Edge { edge: ea, t: ta },
                    SurfacePoint::Edge { edge: eb, t: tb },
                ) => {
                    assert_eq!(ea, eb, "{label}: path {p} point {i}: crossed input edge");
                    let raw = (ta - tb).abs() * lengths[*ea as usize];
                    let cond = raw * sines[i].unwrap_or(1.0);
                    stats.worst_raw = stats.worst_raw.max(raw);
                    stats.worst_conditioned = stats.worst_conditioned.max(cond);
                    assert!(cond <= 1e-9, "{label}: path {p} point {i}: |Δ| sin θ = {cond:e}");
                }
                _ => panic!("{label}: path {p} point {i}: {qa:?} vs {qb:?}"),
            }
            stats.compared_points += 1;
        }
    }
}

fn model_shapes() -> Vec<Shape> {
    let mut shapes = generated_shapes();
    shapes.push(waist_tube(10, 3));
    shapes.extend(data_shapes());
    shapes.extend(load_data("models/bob/bob/bob_tri.obj"));
    shapes.extend(load_data("models/blub/blub/blub_triangulated.obj"));
    shapes.extend(load_data("models/crumpleddevelopable/CrumpledDevelopable.obj"));
    shapes
}

#[test]
fn extraction_certifies_and_representations_agree() {
    let mut stats = Stats::default();
    for shape in &model_shapes() {
        let integers = IntegerTriangulation::new(&shape.mesh, &shape.geometry)
            .unwrap_or_else(|e| panic!("{}: {e}", shape.name));
        let tri = shape.triangulation();
        let n = tri.complex().n_vertices();
        let seed = 41;
        let pairs = random_pairs(n, if n > 3000 { 4 } else { 10 }, seed);
        for (k, &(a, b)) in pairs.iter().enumerate() {
            let path = dijkstra_path(&tri, a, b).unwrap();
            let label = format!("{} seed {seed} pair {k} ({a}->{b})", shape.name);
            check_both(shape, &integers, &label, &[path], &[], &mut stats);
        }
        // Several paths at once, and one through waypoints.
        let paths: Vec<EdgePath> = pairs
            .iter()
            .take(3)
            .map(|&(a, b)| dijkstra_path(&tri, a, b).unwrap())
            .collect();
        check_both(shape, &integers, &format!("{} three paths", shape.name), &paths, &[], &mut stats);
        let w = [pairs[0].0, pairs[1].0, pairs[2].0];
        if let Ok(init) = waypoint_path(&tri, &w, true, true) {
            check_both(
                shape,
                &integers,
                &format!("{} waypoint loop", shape.name),
                &init.paths,
                &init.pinned,
                &mut stats,
            );
        }
        // Loops along self-edges.
        for &s in &tri.irregularities().self_edges {
            let path = EdgePath::closed(vec![Segment::new(s as usize, true)]);
            check_both(shape, &integers, &format!("{} self-edge {s}", shape.name), &[path], &[], &mut stats);
        }
    }
    // Loops around a tube with a waist, and two parallel loops.
    let (n, w) = (10, 3);
    let tube = waist_tube(n, w);
    let integers = IntegerTriangulation::new(&tube.mesh, &tube.geometry).unwrap();
    let tri = tube.triangulation();
    let ring = |j: usize| vpath(&tri, &(0..n).map(|k| j * n + k).collect::<Vec<_>>(), true);
    check_both(&tube, &integers, "tube loop", &[ring(w + 2)], &[], &mut stats);
    check_both(&tube, &integers, "tube two loops", &[ring(w - 2), ring(w + 2)], &[], &mut stats);
    eprintln!("{stats:#?}");
    assert!(stats.paths > 100 && stats.crossings > 500 && stats.joints_checked > 200, "{stats:?}");
}

#[test]
fn cut_graph_and_seam_networks_extract() {
    let mut stats = Stats::default();
    let mut shapes = vec![
        Shape::generated("torus 12x6", ddg_core::generators::torus(12, 6, 1.0, 0.3).unwrap()),
        Shape::generated("icosphere 2", ddg_core::generators::icosphere(2).unwrap()),
        long_torus(8, 3, 4.0, 0.2, 0.37),
    ];
    shapes.extend(data_shapes().into_iter().take(1));
    for shape in &shapes {
        let integers = IntegerTriangulation::new(&shape.mesh, &shape.geometry).unwrap();
        let tri = shape.triangulation();
        let c = tri.complex();
        let mask = if c.euler_characteristic() == 0 { cut_graph(c) } else { bfs_tree(c) };
        let init = network_from_edges(&tri, &mask, &[]).unwrap();
        assert!(!init.paths.is_empty());
        check_both(shape, &integers, &format!("{} seams", shape.name), &init.paths, &init.pinned, &mut stats);
    }
    eprintln!("{stats:#?}");
    assert!(stats.paths >= 4);
}

#[test]
fn extraction_is_deterministic_across_thread_counts() {
    for shape in model_shapes().iter().filter(|s| s.mesh.n_vertices() < 6000) {
        let integers = IntegerTriangulation::new(&shape.mesh, &shape.geometry).unwrap();
        let tri = shape.triangulation();
        let paths: Vec<EdgePath> = random_pairs(tri.complex().n_vertices(), 6, 77)
            .into_iter()
            .map(|(a, b)| dijkstra_path(&tri, a, b).unwrap())
            .collect();
        let Ok(_) = Network::new(tri.complex(), &paths, &[], None) else { continue };
        let net = shorten(integers, &paths, &[], &shape.name);
        let options = ExtractOptions::default();
        let run = |threads: usize| {
            rayon::ThreadPoolBuilder::new()
                .num_threads(threads)
                .build()
                .unwrap()
                .install(|| net.extract_with_report(&options).unwrap())
        };
        let one = run(1);
        for threads in [3, 16] {
            let many = run(threads);
            assert_eq!(one, many, "{}: {threads} threads", shape.name);
        }
        // Bitwise, not just equal: the positions' bits.
        let bits = |x: &Extraction| -> Vec<u64> {
            x.paths
                .iter()
                .flat_map(|p| &p.points)
                .filter_map(|q| match q {
                    SurfacePoint::Edge { t, .. } => Some(t.to_bits()),
                    _ => None,
                })
                .collect()
        };
        assert_eq!(bits(&one), bits(&run(16)));
    }
}

#[test]
fn polylines_3d_interpolate_the_input() {
    let shape = grid(6, 4);
    let integers = IntegerTriangulation::new(&shape.mesh, &shape.geometry).unwrap();
    let tri = shape.triangulation();
    let path = dijkstra_path(&tri, gv(4, 0, 0), gv(4, 6, 4)).unwrap();
    let net = shorten(integers, &[path], &[], "grid");
    let lines = net
        .polylines_3d(&shape.positions, &ExtractOptions::default())
        .unwrap();
    // On a flat grid the geodesic is the straight segment between the corners: collinear, and its
    // length is the distance.
    let line = &lines[0];
    let (p, q) = (line[0], *line.last().unwrap());
    assert_eq!(p, shape.positions[gv(4, 0, 0)]);
    assert_eq!(q, shape.positions[gv(4, 6, 4)]);
    let total: f64 = line.windows(2).map(|w| dist(w[0], w[1])).sum();
    assert!((total - dist(p, q)).abs() <= 1e-12 * total, "{total} vs {}", dist(p, q));
    assert!(line.len() > 3);
    for r in line {
        let (u, v) = ([q[0] - p[0], q[1] - p[1]], [r[0] - p[0], r[1] - p[1]]);
        assert!((u[0] * v[1] - u[1] * v[0]).abs() < 1e-12);
    }
    let xs = net.extract(&ExtractOptions::default()).unwrap();
    let e = polyline_3d(net.triangulation().input().mesh(), &shape.positions[1..], &xs[0]);
    assert!(matches!(e, Err(FlipError::IndexOutOfRange { .. })));
}

#[test]
fn bad_options_are_refused() {
    let shape = grid(3, 3);
    let tri = shape.triangulation();
    let path = dijkstra_path(&tri, 0, 15).unwrap();
    let net = GeodesicNetwork::new(tri, &[path], &[], None).unwrap();
    let options = ExtractOptions { crossing_tolerance: f64::NAN, ..ExtractOptions::default() };
    let e = net.extract(&options).unwrap_err();
    assert!(e.to_string().contains("crossing_tolerance"), "{e}");
}

// ===== the evaluation (release, alone) =====

/// Per mesh: many random trials on both representations; timings of Dijkstra, FlipOut and
/// extraction on one thread; signpost trace failures (mesh, seed). Run alone:
/// `cargo test -p ddg-flip-geodesics --release --test extract -- --ignored --nocapture evaluation`.
#[test]
#[ignore]
fn evaluation() {
    let trials: usize = std::env::var("FO_TRIALS").ok().and_then(|s| s.parse().ok()).unwrap_or(200);
    let pool = rayon::ThreadPoolBuilder::new().num_threads(1).build().unwrap();
    println!("| mesh | V | trials | Dijkstra ms | FlipOut ms | extract ms signposts | extract ms integers | pieces median | signpost failures | worst crossing dev | worst joint | worst |Δ|sinθ | worst raw |Δ| | worst length err | worst int mismatch |");
    for shape in &model_shapes() {
        if shape.mesh.n_vertices() < 250 {
            continue;
        }
        let integers = IntegerTriangulation::new(&shape.mesh, &shape.geometry).unwrap();
        let tri = shape.triangulation();
        let n = tri.complex().n_vertices();
        let (mut t_dij, mut t_fo, mut t_xs, mut t_xi) = (vec![], vec![], vec![], vec![]);
        let mut pieces = vec![];
        let mut stats = Stats::default();
        let options = ExtractOptions { certify_input: false, ..ExtractOptions::default() };
        let cert_options = ExtractOptions::default();
        for (k, &(a, b)) in random_pairs(n, trials, 2024).iter().enumerate() {
            let t0 = Instant::now();
            let path = dijkstra_path(&tri, a, b).unwrap();
            t_dij.push(t0.elapsed().as_secs_f64() * 1e3);
            let label = format!("{} seed 2024 pair {k} ({a}->{b})", shape.name);
            pool.install(|| {
                let mut net = GeodesicNetwork::new(integers.clone(), &[path.clone()], &[], None).unwrap();
                let t0 = Instant::now();
                net.shorten(&ShortenOptions::default()).unwrap();
                t_fo.push(t0.elapsed().as_secs_f64() * 1e3);
                let mut sp = GeodesicNetwork::new(tri.clone(), &[path.clone()], &[], None).unwrap();
                sp.shorten(&ShortenOptions::default()).unwrap();
                let t0 = Instant::now();
                let xi = net.extract(&options).unwrap();
                t_xi.push(t0.elapsed().as_secs_f64() * 1e3);
                pieces.push(xi[0].faces.len());
                let t0 = Instant::now();
                let xs = sp.extract(&options);
                t_xs.push(t0.elapsed().as_secs_f64() * 1e3);
                // Certificates and comparison (untimed).
                let cb = net.extract_with_report(&cert_options).unwrap();
                stats.note(1, &cb);
                if !cb.certificate.as_ref().unwrap().passed {
                    println!("INTEGER CERTIFICATE FAILED: {label}: {:?}", cb.certificate);
                }
                match sp.extract_with_report(&cert_options) {
                    Err(e) => {
                        println!("SIGNPOST TRACE FAILURE: {label}: {e}");
                        stats.signpost_failures.push(label.clone());
                    }
                    Ok(ca) => {
                        stats.note(0, &ca);
                        let _ = xs;
                        let lengths = integers.input().geometry().edge_lengths();
                        let sines = crossing_sines(&integers, &cb.paths[0]);
                        if ca.paths[0].points.len() != cb.paths[0].points.len() {
                            println!("DIFFERENT POINT COUNT: {label}");
                        } else {
                            for (i, (qa, qb)) in ca.paths[0].points.iter().zip(&cb.paths[0].points).enumerate() {
                                if let (
                                    SurfacePoint::Edge { edge: ea, t: ta },
                                    SurfacePoint::Edge { edge: eb, t: tb },
                                ) = (qa, qb)
                                {
                                    if ea != eb {
                                        println!("DIFFERENT CROSSED EDGE: {label} point {i}");
                                        continue;
                                    }
                                    let raw = (ta - tb).abs() * lengths[*ea as usize];
                                    stats.worst_raw = stats.worst_raw.max(raw);
                                    stats.worst_conditioned =
                                        stats.worst_conditioned.max(raw * sines[i].unwrap_or(1.0));
                                }
                            }
                        }
                    }
                }
            });
        }
        let med = |v: &mut Vec<f64>| {
            v.sort_by(f64::total_cmp);
            v[v.len() / 2]
        };
        pieces.sort_unstable();
        println!(
            "| {} | {} | {} | {:.3} | {:.3} | {:.3} | {:.3} | {} | {} | {:.1e} / {:.1e} | {:.1e} / {:.1e} | {:.1e} | {:.1e} | {:.1e} / {:.1e} | {:.1e} |",
            shape.name, n, trials, med(&mut t_dij), med(&mut t_fo), med(&mut t_xs), med(&mut t_xi),
            pieces[pieces.len() / 2], stats.signpost_failures.len(),
            stats.worst_cross_dev[0], stats.worst_cross_dev[1],
            stats.worst_joint[0], stats.worst_joint[1],
            stats.worst_conditioned, stats.worst_raw,
            stats.worst_length[0], stats.worst_length[1], stats.worst_mismatch[1],
        );
    }
}
