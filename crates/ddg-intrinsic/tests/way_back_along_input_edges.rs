//! The way back of a cut segment that runs exactly along an input edge (SPEC N13, N28).
//!
//! A straight segment in a face of the triangulation whose two ends lie on the same input edge
//! runs along that edge. This happens wherever a cut follows an input edge: on bob's control
//! mesh, the start plane of the cut flow is the mirror plane, which contains input edges, and
//! the level set of a symmetric function vanishes exactly on them. The ends are found on the
//! input mesh by tracing (SPEC N11), and such a trace ends a rounding error inside one of the
//! two input faces of the edge, chosen by rounding. Before the fix the segment's trace stayed
//! in one face and its end had been placed in the other ("the place of its end is in none of
//! the input faces the path went through"), at every tolerance.
//!
//! The meshes here are strips along a straight line of input edges (`y = 0`), with a row of
//! vertices above and, except for the boundary case, below; the strip is flat or folded along
//! the line like a mirror-symmetric surface (`z = fold · |y|`), and also mirrored (`y → −y`)
//! so that the roles of the two faces of every line edge are swapped. The triangulation has
//! vertices inserted on the line edges (splits) and inside the faces, then Delaunay flips.
//! The segments are found from the edge traces:
//!
//! * between the two points where two sides of a face cross the same line edge;
//! * from a corner on the line (an input vertex or an inserted vertex on a line edge) to the
//!   point where the opposite side crosses the line edge, and back (ending at an input vertex);
//! * along a side of the triangulation that lies along a line edge (interior or boundary):
//!   between two points of it, from its corners and to them.

mod common;

use common::{Shape, dist, inner_bary, position};
use ddg_core::generators::SplitMix64;
use ddg_intrinsic::formulas::triangle_layout;
use ddg_intrinsic::{
    DelaunayOptions, EdgeTraceOptions, IntrinsicTriangulation, PathsOnInput, Segment, SurfacePoint,
};

/// A strip of `n` cells along the line `y = 0` from `x = 0` to `x = n`, with rows of vertices
/// at `y ≈ 1` and, if `below`, at `y ≈ −1`. The line vertices stay exactly on the line (their
/// spacing is jittered); the others are jittered in both directions. `fold` lifts the rows off
/// the plane (`z = fold · |y|`, a crease along the line); `mirror` flips `y`.
fn strip(n: usize, below: bool, fold: f64, mirror: bool, seed: u64) -> Shape {
    let mut rng = SplitMix64::new(seed);
    let mut positions = Vec::new();
    let row = |y: f64, positions: &mut Vec<[f64; 3]>, rng: &mut SplitMix64| {
        for i in 0..=n {
            let inner = i > 0 && i < n;
            let x = i as f64
                + if inner {
                    0.3 * rng.uniform(-1.0, 1.0)
                } else {
                    0.0
                };
            let y = if y == 0.0 {
                0.0
            } else {
                y * (1.0 + 0.25 * rng.uniform(-1.0, 1.0))
            };
            let y = if mirror { -y } else { y };
            positions.push([x, y, fold * y.abs()]);
        }
    };
    row(0.0, &mut positions, &mut rng);
    row(1.0, &mut positions, &mut rng);
    if below {
        row(-1.0, &mut positions, &mut rng);
    }
    let line = |i: usize| i as u32;
    let top = |i: usize| (n + 1 + i) as u32;
    let bottom = |i: usize| (2 * (n + 1) + i) as u32;
    let mut faces = Vec::new();
    for i in 0..n {
        // Alternate diagonals, so that the line vertices have different degrees.
        if i % 2 == 0 {
            faces.push([line(i), line(i + 1), top(i + 1)]);
            faces.push([line(i), top(i + 1), top(i)]);
        } else {
            faces.push([line(i), line(i + 1), top(i)]);
            faces.push([line(i + 1), top(i + 1), top(i)]);
        }
        if below {
            if i % 2 == 0 {
                faces.push([bottom(i), bottom(i + 1), line(i + 1)]);
                faces.push([bottom(i), line(i + 1), line(i)]);
            } else {
                faces.push([bottom(i), bottom(i + 1), line(i)]);
                faces.push([bottom(i + 1), line(i + 1), line(i)]);
            }
        }
    }
    if mirror {
        // y → −y reverses the orientation: reverse the faces back.
        for f in &mut faces {
            f.swap(1, 2);
        }
    }
    Shape::new(
        &format!("strip {n} below {below} fold {fold} mirror {mirror} seed {seed}"),
        positions,
        faces,
        fold == 0.0,
    )
}

/// The input edges on the line `y = 0`.
fn line_edges(shape: &Shape, n: usize) -> Vec<usize> {
    (0..n)
        .map(|i| shape.mesh.find_edge(i, i + 1).expect("a line edge"))
        .collect()
}

/// Splits every line edge at one or two random places, inserts vertices inside random faces,
/// and flips to Delaunay.
fn triangulation(shape: &Shape, line: &[usize], rng: &mut SplitMix64) -> IntrinsicTriangulation {
    let mut tri = shape.triangulation();
    for &e in line {
        tri.split_edge(e, rng.uniform(0.3, 0.7)).unwrap();
        if rng.next_f64() < 0.5 {
            // The half from the tail kept the id `e`.
            tri.split_edge(e, rng.uniform(0.3, 0.7)).unwrap();
        }
    }
    for _ in 0..3 * shape.mesh.n_faces() {
        let f = (rng.next_u64() % tri.n_faces() as u64) as usize;
        // Thin triangles are refused; that is fine.
        let _ = tri.insert_vertex(f, inner_bary(rng));
    }
    tri.flip_to_delaunay(&DelaunayOptions::default()).unwrap();
    tri
}

/// The input edge of the line that a vertex of the triangulation lies on (as an inserted
/// vertex), or that it ends (as an input vertex); several for an input vertex.
fn line_edges_at(tri: &IntrinsicTriangulation, line: &[usize], v: usize) -> Vec<usize> {
    let mesh = tri.input().mesh();
    match tri.vertex_location(v) {
        SurfacePoint::Vertex(u) => line
            .iter()
            .copied()
            .filter(|&e| mesh.edge_vertices(e).contains(&(u as usize)))
            .collect(),
        SurfacePoint::Edge { edge, .. } if line.contains(&(edge as usize)) => vec![edge as usize],
        _ => Vec::new(),
    }
}

#[derive(Default, Debug)]
struct Counts {
    crossing_to_crossing: usize,
    corner_to_crossing: usize,
    from_input_vertex: usize,
    along_a_side: usize,
}

/// The segments of `tri` that run along a line edge, as nodes and segments for
/// `segments_on_input`, with the line edge of every segment.
fn segments_along_the_line(
    tri: &IntrinsicTriangulation,
    line: &[usize],
    counts: &mut Counts,
) -> (Vec<SurfacePoint>, Vec<Segment>, Vec<usize>) {
    let traces = tri.trace_edges(&EdgeTraceOptions::default()).unwrap();
    let mut nodes = Vec::new();
    let mut segments = Vec::new();
    let mut along = Vec::new();
    let node = |p: SurfacePoint, nodes: &mut Vec<SurfacePoint>| {
        nodes.push(p);
        (nodes.len() - 1) as u32
    };
    let mut both_ways = |face: usize, a: u32, b: u32, e: usize| {
        for pair in [[a, b], [b, a]] {
            segments.push(Segment {
                face: face as u32,
                nodes: pair,
            });
            along.push(e);
        }
    };
    for f in 0..tri.n_faces() {
        // Where the sides of the face cross line edges: (slot, line edge, node).
        let mut crossings = Vec::new();
        for k in 0..3 {
            let e = tri.edge(3 * f + k);
            for c in traces.crossings(e) {
                if line.contains(&(c.edge as usize)) {
                    let t = c.distance / tri.edge_length(e);
                    assert!(t > 0.0 && t < 1.0);
                    crossings.push((k, c.edge as usize, SurfacePoint::Edge { edge: e as u32, t }));
                }
            }
        }
        for (i, &(k1, e1, p1)) in crossings.iter().enumerate() {
            for &(k2, e2, p2) in &crossings[i + 1..] {
                if k1 != k2 && e1 == e2 {
                    let (a, b) = (node(p1, &mut nodes), node(p2, &mut nodes));
                    both_ways(f, a, b, e1);
                    counts.crossing_to_crossing += 1;
                }
            }
        }
        for k in 0..3 {
            let v = tri.faces()[f][k] as usize;
            let at = line_edges_at(tri, line, v);
            // The side opposite corner k is halfedge k + 1.
            for &(slot, e, p) in &crossings {
                if slot == (k + 1) % 3 && at.contains(&e) {
                    let (a, b) = (
                        node(SurfacePoint::Vertex(v as u32), &mut nodes),
                        node(p, &mut nodes),
                    );
                    both_ways(f, a, b, e);
                    counts.corner_to_crossing += 1;
                    if matches!(tri.vertex_location(v), SurfacePoint::Vertex(_)) {
                        counts.from_input_vertex += 1;
                    }
                }
            }
            // A side along a line edge: two points of it, and its corners.
            let e = tri.edge(3 * f + k);
            if let Some(input) = tri.edge_on_input(e).filter(|i| line.contains(i)) {
                let [a, b] = tri.edge_vertices(e);
                let p = node(
                    SurfacePoint::Edge {
                        edge: e as u32,
                        t: 0.3,
                    },
                    &mut nodes,
                );
                let q = node(
                    SurfacePoint::Edge {
                        edge: e as u32,
                        t: 0.8,
                    },
                    &mut nodes,
                );
                let va = node(SurfacePoint::Vertex(a as u32), &mut nodes);
                let vb = node(SurfacePoint::Vertex(b as u32), &mut nodes);
                both_ways(f, p, q, input);
                both_ways(f, va, q, input);
                both_ways(f, p, vb, input);
                counts.along_a_side += 1;
            }
        }
    }
    (nodes, segments, along)
}

/// Whether a point of the input mesh lies in the closed input face `f`.
fn in_input_face(shape: &Shape, point: SurfacePoint, f: usize) -> bool {
    match point {
        SurfacePoint::Vertex(v) => shape.mesh.faces()[f].contains(&v),
        SurfacePoint::Edge { edge, .. } => {
            (0..3).any(|k| shape.mesh.edge(3 * f + k) == edge as usize)
        }
        SurfacePoint::Face { face, .. } => face as usize == f,
    }
}

/// Whether a point of the input mesh lies on the closed input edge `e`.
fn on_input_edge(shape: &Shape, point: SurfacePoint, e: usize) -> bool {
    match point {
        SurfacePoint::Vertex(v) => shape.mesh.edge_vertices(e).contains(&(v as usize)),
        SurfacePoint::Edge { edge, .. } => edge as usize == e,
        SurfacePoint::Face { .. } => false,
    }
}

/// The checks on the way back of segments along line edges: every end is placed on its line
/// edge (a point of both its faces), every piece lies in the input face it names, and the path
/// has the length of the segment and lies on the line.
fn check(
    shape: &Shape,
    tri: &IntrinsicTriangulation,
    nodes: &[SurfacePoint],
    segments: &[Segment],
    along: &[usize],
    back: &PathsOnInput,
) {
    let name = &shape.name;
    // Relative to the segment the mismatch can reach 5e-10: some segments are 1e-3 of the
    // strip's cells, in thin faces; relative to the cells (of size 1) it is a rounding error.
    assert!(
        back.max_relative_mismatch < 1e-9,
        "{name}: mismatch {}",
        back.max_relative_mismatch
    );
    for path in &back.paths {
        assert!(path.mismatch < 1e-12, "{name}: {path:?}");
    }
    for (s, (segment, path)) in segments.iter().zip(&back.paths).enumerate() {
        let context = format!(
            "{name}: segment {s} {segment:?} along input edge {}",
            along[s]
        );
        let [from, to] = segment.nodes.map(|n| n as usize);
        assert_eq!(path.points[0], back.nodes[from], "{context}");
        assert_eq!(*path.points.last().unwrap(), back.nodes[to], "{context}");
        // A segment along an input edge crosses nothing.
        assert_eq!(path.points.len(), 2, "{context}: {path:?}");
        for &p in &path.points {
            assert!(on_input_edge(shape, p, along[s]), "{context}: {p} {path:?}");
            assert!(
                in_input_face(shape, p, path.faces[0] as usize),
                "{context}: {path:?}"
            );
            assert!(position(shape, p)[1].abs() < 1e-12, "{context}: {p}");
        }
        // The length of the segment in its face.
        let f = segment.face as usize;
        let v = triangle_layout(tri.face_lengths(f)).unwrap();
        let at = |p: SurfacePoint| -> [f64; 2] {
            match p {
                SurfacePoint::Vertex(u) => v[tri.faces()[f].iter().position(|&w| w == u).unwrap()],
                SurfacePoint::Edge { edge, t } => {
                    let k = (0..3)
                        .find(|&k| tri.edge(3 * f + k) == edge as usize)
                        .unwrap();
                    let th = if tri.edge_halfedge(edge as usize) == 3 * f + k {
                        t
                    } else {
                        1.0 - t
                    };
                    let (a, b) = (v[k], v[(k + 1) % 3]);
                    [a[0] + th * (b[0] - a[0]), a[1] + th * (b[1] - a[1])]
                }
                SurfacePoint::Face { .. } => unreachable!(),
            }
        };
        let (p, q) = (at(nodes[from]), at(nodes[to]));
        let expected = (q[0] - p[0]).hypot(q[1] - p[1]);
        let length = dist(
            position(shape, path.points[0]),
            position(shape, path.points[1]),
        );
        assert!(
            (length - expected).abs() < 1e-11 * expected.max(1.0),
            "{context}: the path has length {length}, the segment {expected}"
        );
    }
}

fn run_case(shape: &Shape, n: usize, seed: u64, counts: &mut Counts) {
    let line = line_edges(shape, n);
    let mut rng = SplitMix64::new(seed);
    let tri = triangulation(shape, &line, &mut rng);
    let (nodes, segments, along) = segments_along_the_line(&tri, &line, counts);
    // Every segment alone, so that a failure names it.
    let options = EdgeTraceOptions::default();
    let mut failed = Vec::new();
    for (s, segment) in segments.iter().enumerate() {
        if let Err(e) = tri.segments_on_input(&nodes, &[*segment], &options) {
            failed.push(format!("segment {s} {segment:?}: {e}"));
        }
    }
    assert!(
        failed.is_empty(),
        "{}: {} of {} segments along input edges fail:\n{}",
        shape.name,
        failed.len(),
        segments.len(),
        failed.join("\n")
    );
    let back = tri.segments_on_input(&nodes, &segments, &options).unwrap();
    check(shape, &tri, &nodes, &segments, &along, &back);
    // The places of the nodes are those of the point queries.
    for (node, place) in nodes.iter().zip(&back.nodes) {
        assert_eq!(tri.point_on_input(*node).unwrap(), *place, "{}", shape.name);
    }
}

/// Segments between two points where sides of a face cross the same input edge, and from a
/// corner on an input edge (inserted or input vertex) to such a point and back, with the input
/// edge between two faces: flat and folded, both orientations (so both faces of every edge
/// are the one the trace starts in, in some case).
#[test]
fn segments_along_an_interior_input_edge_come_back() {
    let n = 6;
    let mut counts = Counts::default();
    for seed in 0..6u64 {
        for fold in [0.0, 0.7] {
            for mirror in [false, true] {
                let shape = strip(n, true, fold, mirror, 100 + seed);
                run_case(&shape, n, 200 + seed, &mut counts);
            }
        }
    }
    assert!(counts.crossing_to_crossing >= 20, "{counts:?}");
    assert!(counts.corner_to_crossing >= 20, "{counts:?}");
    assert!(counts.from_input_vertex >= 1, "{counts:?}");
    assert!(counts.along_a_side >= 20, "{counts:?}");
}

/// The same on a boundary edge: the line is the boundary of the strip, so the only segments
/// along it run along sides of the triangulation that lie on it.
#[test]
fn segments_along_a_boundary_input_edge_come_back() {
    let n = 6;
    let mut counts = Counts::default();
    for seed in 0..4u64 {
        for mirror in [false, true] {
            let shape = strip(n, false, 0.0, mirror, 300 + seed);
            run_case(&shape, n, 400 + seed, &mut counts);
        }
    }
    assert_eq!(counts.crossing_to_crossing, 0);
    assert!(counts.along_a_side >= 20, "{counts:?}");
}

/// A point of the triangulation on an input edge is placed on that edge by the point query,
/// not a rounding error inside one of its faces (SPEC N11, N28).
#[test]
fn points_on_an_input_edge_are_placed_on_it() {
    let n = 6;
    let mut placed = 0;
    for seed in 0..4u64 {
        let shape = strip(n, true, 0.7, seed % 2 == 1, 500 + seed);
        let line = line_edges(&shape, n);
        let mut rng = SplitMix64::new(600 + seed);
        let tri = triangulation(&shape, &line, &mut rng);
        let traces = tri.trace_edges(&EdgeTraceOptions::default()).unwrap();
        for e in 0..tri.n_edges() {
            for c in traces.crossings(e) {
                let t = c.distance / tri.edge_length(e);
                let point = SurfacePoint::Edge { edge: e as u32, t };
                let place = tri.point_on_input(point).unwrap();
                assert!(
                    matches!(place, SurfacePoint::Edge { edge, .. } if edge == c.edge),
                    "{}: {point} crosses input edge {} and is placed at {place}",
                    shape.name,
                    c.edge
                );
                let gap = dist(
                    position(&shape, place),
                    position(
                        &shape,
                        SurfacePoint::Edge {
                            edge: c.edge,
                            t: c.t,
                        },
                    ),
                );
                assert!(gap < 1e-12, "{}: {gap}", shape.name);
                placed += 1;
            }
        }
    }
    assert!(placed >= 100, "{placed}");
}
