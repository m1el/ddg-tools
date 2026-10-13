//! SPEC §5 item 6: the overlay has the input's area, is a valid mesh, and each of its faces
//! lies in one input face and one face of the triangulation. And the operation the cut flow
//! needs: a straight segment of the triangulation as a path over the input mesh.
//!
//! The reference program builds the overlay in its viewer only, so there is nothing to
//! compare with; the checks are the invariants, on triangulations made by flips, insertions,
//! splits and refinement, and round trips through point queries.

mod common;

use common::{Shape, dist, position, standard_shapes};
use ddg_core::generators::SplitMix64;
use ddg_core::{Geometry, TriMesh};
use ddg_intrinsic::{
    DelaunayOptions, EdgeTraceOptions, IntrinsicError, IntrinsicTriangulation, Overlay,
    RefineOptions, Segment, SurfacePoint,
};

fn sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
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

/// Whether a point of the triangulation lies in its closed face `f`.
fn in_face(tri: &IntrinsicTriangulation, point: SurfacePoint, f: usize) -> bool {
    match point {
        SurfacePoint::Vertex(v) => tri.faces()[f].contains(&v),
        SurfacePoint::Edge { edge, .. } => (0..3).any(|k| tri.edge(3 * f + k) == edge as usize),
        SurfacePoint::Face { face, .. } => face as usize == f,
    }
}

/// Item 6 for the overlay of `tri`. Returns the overlay.
fn assert_overlay(shape: &Shape, tri: &IntrinsicTriangulation, what: &str) -> Overlay {
    let context = format!("{} ({what})", shape.name);
    let traces = tri
        .trace_edges(&EdgeTraceOptions::default())
        .unwrap_or_else(|e| panic!("{context}: {e}"));
    let overlay = tri
        .overlay(&traces)
        .unwrap_or_else(|e| panic!("{context}: {e}"));
    let mesh = &shape.mesh;

    // Vertices: those of the triangulation, then the crossings.
    assert_eq!(
        overlay.vertices.len(),
        tri.n_vertices() + traces.report.n_crossings,
        "{context}"
    );
    for (v, vertex) in overlay.vertices.iter().enumerate().take(tri.n_vertices()) {
        assert_eq!(
            vertex.on_triangulation,
            SurfacePoint::Vertex(v as u32),
            "{context}"
        );
        assert_eq!(vertex.on_input, tri.vertex_location(v), "{context}");
    }
    for vertex in &overlay.vertices[tri.n_vertices()..] {
        assert!(
            matches!(vertex.on_input, SurfacePoint::Edge { .. }),
            "{context}"
        );
        assert!(
            matches!(vertex.on_triangulation, SurfacePoint::Edge { .. }),
            "{context}"
        );
        // The same point of the surface, by a point query.
        let found = tri.point_on_input(vertex.on_triangulation).unwrap();
        let gap = dist(position(shape, found), position(shape, vertex.on_input));
        assert!(
            gap < 1e-9 * shape.mean_edge_length(),
            "{context}: {vertex:?}: {gap}"
        );
    }
    let positions: Vec<[f64; 3]> = overlay
        .vertices
        .iter()
        .map(|v| position(shape, v.on_input))
        .collect();

    // Faces: in order of the faces of the triangulation, each in one face of either mesh.
    let mut area_in_input = vec![0.0; mesh.n_faces()];
    let mut area_in_tri = vec![0.0; tri.n_faces()];
    let mut pieces_of = vec![0usize; tri.n_faces()];
    for f in 0..overlay.n_faces() {
        let corners = overlay.face(f);
        let (g, t) = (
            overlay.face_input[f] as usize,
            overlay.face_triangulation[f] as usize,
        );
        assert!(
            (3..=6).contains(&corners.len()),
            "{context}: face {f}: {corners:?}"
        );
        if f > 0 {
            assert!(overlay.face_triangulation[f - 1] <= overlay.face_triangulation[f]);
        }
        pieces_of[t] += 1;
        for &c in corners {
            let vertex = overlay.vertices[c as usize];
            assert!(
                in_input_face(shape, vertex.on_input, g),
                "{context}: face {f}: {} is not in input face {g}",
                vertex.on_input
            );
            assert!(
                in_face(tri, vertex.on_triangulation, t),
                "{context}: face {f}: {} is not in face {t} of the triangulation",
                vertex.on_triangulation
            );
        }
        // A convex polygon in the plane of its input face, counter-clockwise.
        let normal = {
            let [a, b, c] = mesh.faces()[g].map(|v| shape.positions[v as usize]);
            cross(sub(b, a), sub(c, a))
        };
        let scale = dot(normal, normal).sqrt();
        let mut area = 0.0;
        let n = corners.len();
        for k in 0..n {
            let (a, b, c) = (
                positions[corners[k] as usize],
                positions[corners[(k + 1) % n] as usize],
                positions[corners[(k + 2) % n] as usize],
            );
            let turn = dot(cross(sub(b, a), sub(c, b)), normal) / scale;
            assert!(
                turn > 0.0,
                "{context}: face {f} is not convex at corner {k}: {turn}"
            );
            if k + 2 < n {
                let (o, p, q) = (positions[corners[0] as usize], b, c);
                area += 0.5 * dot(cross(sub(p, o), sub(q, o)), normal) / scale;
            }
        }
        assert!(area > 0.0, "{context}: face {f} has area {area}");
        area_in_input[g] += area;
        area_in_tri[t] += area;
    }
    assert!(pieces_of.iter().all(|&n| n > 0), "{context}");
    // The pieces in an input face make up that face; those in a face of the triangulation
    // make up that face; all together the surface.
    let total = shape.geometry.total_area();
    for (g, &area) in area_in_input.iter().enumerate() {
        let expected = shape.geometry.face_areas()[g];
        assert!(
            (area - expected).abs() < 1e-9 * total,
            "{context}: the pieces in input face {g} have area {area}, the face {expected}"
        );
    }
    for (t, &area) in area_in_tri.iter().enumerate() {
        let expected = tri.face_area(t);
        assert!(
            (area - expected).abs() < 1e-9 * total,
            "{context}: the pieces in face {t} have area {area}, the face {expected}"
        );
    }
    let sum: f64 = area_in_input.iter().sum();
    assert!(
        (sum - total).abs() < 1e-11 * total,
        "{context}: area {sum}, input {total}"
    );

    // A valid mesh: its faces cut into triangles are a manifold mesh with the topology of
    // the input, and every triangle has positive area.
    let (triangles, parents) = overlay.triangles();
    assert_eq!(triangles.len(), parents.len());
    let as_mesh = TriMesh::new(overlay.vertices.len(), triangles)
        .unwrap_or_else(|e| panic!("{context}: the overlay is not a valid mesh: {e}"));
    assert_eq!(
        as_mesh.euler_characteristic(),
        mesh.euler_characteristic(),
        "{context}"
    );
    assert_eq!(
        as_mesh.boundary_loops().len(),
        mesh.boundary_loops().len(),
        "{context}"
    );
    assert_eq!(as_mesh.n_components(), mesh.n_components(), "{context}");
    Geometry::from_positions(&as_mesh, &positions)
        .unwrap_or_else(|e| panic!("{context}: the overlay has a degenerate triangle: {e}"));
    overlay
}

#[test]
fn the_overlay_of_a_fresh_triangulation_is_the_input_mesh() {
    for shape in standard_shapes() {
        let tri = shape.triangulation();
        let overlay = assert_overlay(&shape, &tri, "fresh");
        assert_eq!(overlay.vertices.len(), shape.mesh.n_vertices());
        assert_eq!(overlay.corners, shape.mesh.faces().concat());
        let identity: Vec<u32> = (0..shape.mesh.n_faces() as u32).collect();
        assert_eq!(overlay.face_input, identity);
        assert_eq!(overlay.face_triangulation, identity);
    }
}

#[test]
fn the_overlay_is_the_common_subdivision() {
    for (s, shape) in standard_shapes().iter().enumerate() {
        let mut tri = shape.triangulation();
        tri.flip_to_delaunay(&DelaunayOptions::default()).unwrap();
        assert_overlay(shape, &tri, "Delaunay");

        // Edges flipped at random: long edges that cross many input faces.
        let mut rng = SplitMix64::new(910 + s as u64);
        for _ in 0..3 * tri.n_edges() {
            let e = (rng.next_u64() % tri.n_edges() as u64) as usize;
            let _ = tri.flip_edge(e);
        }
        let overlay = assert_overlay(shape, &tri, "random flips");
        assert!(overlay.n_faces() >= tri.n_faces());

        // Vertices inside input faces and on input edges, the boundary included.
        for _ in 0..tri.n_faces() / 2 {
            let f = (rng.next_u64() % tri.n_faces() as u64) as usize;
            let (a, b) = (rng.uniform(0.1, 0.45), rng.uniform(0.1, 0.45));
            let _ = tri.insert_vertex(f, [a, b, 1.0 - a - b]);
            let e = (rng.next_u64() % tri.n_edges() as u64) as usize;
            let _ = tri.split_edge(e, rng.uniform(0.2, 0.8));
        }
        assert_overlay(shape, &tri, "insertions");
        tri.flip_to_delaunay(&DelaunayOptions::default()).unwrap();
        assert_overlay(shape, &tri, "insertions, then Delaunay");
    }
}

#[test]
fn the_overlay_of_a_refined_triangulation() {
    for shape in standard_shapes() {
        let mut tri = shape.triangulation();
        let report = tri.refine(&RefineOptions::default()).unwrap();
        let overlay = assert_overlay(&shape, &tri, "refined");
        if report.n_inserted > 0 {
            assert!(overlay.vertices.len() > shape.mesh.n_vertices());
        }
    }
}

#[test]
fn traces_of_another_triangulation_are_refused() {
    let shape = common::jittered_grid(4, 0.3, 8);
    let mut tri = shape.triangulation();
    let traces = tri.trace_edges(&EdgeTraceOptions::default()).unwrap();
    tri.insert_vertex(0, [0.3, 0.3, 0.4]).unwrap();
    let error = tri.overlay(&traces).unwrap_err();
    assert!(
        matches!(error, IntrinsicError::InvalidOptions { .. }),
        "{error}"
    );
    assert!(
        error
            .to_string()
            .contains("the traces are of a triangulation with"),
        "{error}"
    );
}

fn random_bary(rng: &mut SplitMix64) -> [f64; 3] {
    let (mut a, mut b) = (rng.next_f64(), rng.next_f64());
    if a + b > 1.0 {
        (a, b) = (1.0 - a, 1.0 - b);
    }
    let shrink = |x: f64| 0.03 + 0.91 * x;
    [shrink(a), shrink(b), shrink(1.0 - a - b)]
}

/// The length in space of a path over the input mesh, with the check that every piece lies in
/// the input face it names.
fn checked_length(shape: &Shape, points: &[SurfacePoint], faces: &[u32], context: &str) -> f64 {
    assert_eq!(faces.len() + 1, points.len(), "{context}");
    let mut length = 0.0;
    for (i, pair) in points.windows(2).enumerate() {
        for &p in pair {
            assert!(
                in_input_face(shape, p, faces[i] as usize),
                "{context}: {p} is not in input face {}",
                faces[i]
            );
        }
        length += dist(position(shape, pair[0]), position(shape, pair[1]));
    }
    for p in &points[1..points.len() - 1] {
        assert!(matches!(p, SurfacePoint::Edge { .. }), "{context}: {p}");
    }
    length
}

/// A straight segment inside a face of the triangulation is a path of straight pieces over the
/// input mesh with the same length, between the places of its ends.
#[test]
fn segments_of_the_triangulation_are_paths_on_the_input_mesh() {
    let options = EdgeTraceOptions::default();
    for (s, shape) in standard_shapes().iter().enumerate() {
        let mut tri = shape.triangulation();
        tri.refine(&RefineOptions::default()).unwrap();
        let mut rng = SplitMix64::new(6000 + s as u64);
        let scale = shape.mean_edge_length();
        let mut crossings = 0;
        for query in 0..300 {
            let f = (rng.next_u64() % tri.n_faces() as u64) as usize;
            // Ends inside the face, on its sides, at a corner.
            let mut end = |kind: usize| -> [f64; 3] {
                let b = random_bary(&mut rng);
                match kind % 4 {
                    0 | 1 => b,
                    2 => {
                        let t = b[0] / (b[0] + b[1]);
                        [t, 1.0 - t, 0.0]
                    }
                    _ => [0.0, 0.0, 1.0],
                }
            };
            let (from, to) = (end(query), end(query / 4 + 1));
            if from == to {
                continue;
            }
            let context = format!(
                "{}: query {query}: face {f}, {from:?} to {to:?}",
                shape.name
            );
            let path = tri
                .segment_on_input(f, from, to, &options)
                .unwrap_or_else(|e| panic!("{context}: {e}"));
            // Its length, from the layout of the face.
            let l = tri.face_lengths(f);
            let v = ddg_intrinsic::formulas::triangle_layout(l).unwrap();
            let at = |b: [f64; 3]| [0, 1].map(|i| b[1] * v[1][i] + b[2] * v[2][i]);
            let (p, q) = (at(from), at(to));
            let expected = (q[0] - p[0]).hypot(q[1] - p[1]);
            let length = checked_length(shape, &path.points, &path.faces, &context);
            assert!(
                (length - expected).abs() < 1e-9 * expected.max(scale),
                "{context}: the path has length {length}, the segment {expected}"
            );
            assert!(path.mismatch <= options.tolerance * expected, "{context}");
            crossings += path.points.len() - 2;
            // Its ends are where point queries put them.
            let face = f as u32;
            for (bary, found) in [(from, path.points[0]), (to, *path.points.last().unwrap())] {
                let place = tri
                    .point_on_input(SurfacePoint::Face { face, bary })
                    .unwrap();
                let gap = dist(position(shape, place), position(shape, found));
                assert!(
                    gap < 1e-9 * scale,
                    "{context}: an end is {gap} from its place"
                );
            }
            // On a flat shape the path is the straight segment between them.
            if shape.flat {
                let (a, b) = (
                    position(shape, path.points[0]),
                    position(shape, *path.points.last().unwrap()),
                );
                let mut travelled = 0.0;
                for pair in path.points.windows(2) {
                    travelled += dist(position(shape, pair[0]), position(shape, pair[1]));
                    let x = position(shape, pair[1]);
                    let s = travelled / length;
                    for k in 0..3 {
                        assert!(
                            (x[k] - (a[k] + s * (b[k] - a[k]))).abs() < 1e-10,
                            "{context}"
                        );
                    }
                }
            }
        }
        assert!(
            crossings > 0 || tri.n_vertices() == shape.mesh.n_vertices(),
            "{}",
            shape.name
        );
    }
}

/// Segments that share a node meet in one point of the input mesh: a polyline through the
/// faces of the triangulation maps to a connected path.
#[test]
fn segments_that_share_a_node_meet_on_the_input_mesh() {
    let shape = common::jittered_icosphere(2, 0.35, 5);
    let mut tri = shape.triangulation();
    tri.refine(&RefineOptions::default()).unwrap();
    // The fan of faces around a vertex: a closed polyline through the midpoints of the
    // edges at the vertex, and the three segments from them to a point inside one face.
    let v = 40;
    let spokes: Vec<usize> = tri.outgoing_halfedges(v).collect();
    let mut nodes: Vec<SurfacePoint> = spokes
        .iter()
        .map(|&h| SurfacePoint::Edge {
            edge: tri.edge(h) as u32,
            t: 0.5,
        })
        .collect();
    let mut segments: Vec<Segment> = (0..spokes.len())
        .map(|i| Segment {
            face: tri.face(spokes[i]) as u32,
            nodes: [i as u32, ((i + 1) % spokes.len()) as u32],
        })
        .collect();
    // A triple point in the first face, joined to its three sides.
    let f = tri.face(spokes[0]);
    let centre = nodes.len() as u32;
    nodes.push(SurfacePoint::Face {
        face: f as u32,
        bary: [0.3, 0.3, 0.4],
    });
    for k in 0..3 {
        let side = nodes.len() as u32;
        nodes.push(SurfacePoint::Edge {
            edge: tri.edge(3 * f + k) as u32,
            t: 0.3,
        });
        segments.push(Segment {
            face: f as u32,
            nodes: [side, centre],
        });
    }
    let result = tri
        .segments_on_input(&nodes, &segments, &EdgeTraceOptions::default())
        .unwrap();
    assert_eq!(
        (result.nodes.len(), result.paths.len()),
        (nodes.len(), segments.len())
    );
    assert!(
        result.max_relative_mismatch < 1e-10,
        "{}",
        result.max_relative_mismatch
    );
    for (segment, path) in segments.iter().zip(&result.paths) {
        // Each path starts and ends exactly at the places of its nodes.
        assert_eq!(path.points[0], result.nodes[segment.nodes[0] as usize]);
        assert_eq!(
            *path.points.last().unwrap(),
            result.nodes[segment.nodes[1] as usize]
        );
        checked_length(&shape, &path.points, &path.faces, "fan");
    }
    for (node, place) in nodes.iter().zip(&result.nodes) {
        assert_eq!(tri.point_on_input(*node).unwrap(), *place);
    }
    // Errors: a node that is not in the face of its segment, an end of an edge, a bad index.
    let wrong = Segment {
        face: tri.face(spokes[2]) as u32,
        nodes: [centre, 0],
    };
    assert!(matches!(
        tri.segments_on_input(&nodes, &[wrong], &EdgeTraceOptions::default())
            .unwrap_err(),
        IntrinsicError::InvalidPoint { .. }
    ));
    let mut ends = nodes.clone();
    ends[0] = SurfacePoint::Edge {
        edge: tri.edge(spokes[0]) as u32,
        t: 1.0,
    };
    let error = tri
        .segments_on_input(&ends, &segments, &EdgeTraceOptions::default())
        .unwrap_err();
    assert!(
        error.to_string().contains("must be given as the vertex"),
        "{error}"
    );
    let out_of_range = Segment {
        face: 0,
        nodes: [0, 1000],
    };
    assert!(matches!(
        tri.segments_on_input(&nodes, &[out_of_range], &EdgeTraceOptions::default())
            .unwrap_err(),
        IntrinsicError::IndexOutOfRange { what: "node", .. }
    ));
}
