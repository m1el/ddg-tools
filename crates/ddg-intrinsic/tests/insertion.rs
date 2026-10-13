//! Vertex insertion, edge split and point queries (SPEC §4), with the invariants of SPEC §5
//! items 1 to 3 for triangulations that have vertices on input edges and inside input faces.

mod common;

use std::f64::consts::{PI, TAU};

use common::{Shape, dist, path_length, position, standard_shapes};
use ddg_core::generators::{self, SplitMix64};
use ddg_intrinsic::formulas::FLIP_SHAPE_LIMIT;
use ddg_intrinsic::{
    DelaunayOptions, EdgeTraceOptions, InsertRefusal, IntrinsicError, IntrinsicTriangulation,
    SurfacePoint,
};

fn random_bary(rng: &mut SplitMix64) -> [f64; 3] {
    let (mut a, mut b) = (rng.next_f64(), rng.next_f64());
    if a + b > 1.0 {
        (a, b) = (1.0 - a, 1.0 - b);
    }
    let shrink = |x: f64| 0.03 + 0.91 * x;
    [shrink(a), shrink(b), shrink(1.0 - a - b)]
}

/// Inserts vertices in random faces and splits random edges; returns the numbers of
/// insertions, interior splits and boundary splits made.
fn random_insertions(
    shape: &Shape,
    tri: &mut IntrinsicTriangulation,
    per_face: f64,
    rng: &mut SplitMix64,
) -> [usize; 3] {
    let mut done = [0usize; 3];
    for _ in 0..(per_face * tri.n_faces() as f64) as usize + 2 {
        let f = (rng.next_u64() % tri.n_faces() as u64) as usize;
        let (nv, ne, nf) = (tri.n_vertices(), tri.n_edges(), tri.n_faces());
        match tri.insert_vertex(f, random_bary(rng)) {
            Ok(insertion) => {
                done[0] += 1;
                assert_eq!(insertion.vertex as usize, nv);
                assert_eq!(insertion.new_edges, [ne, ne + 1, ne + 2].map(|e| e as u32));
                assert_eq!(insertion.new_faces, [nf, nf + 1].map(|g| g as u32));
                assert_eq!(tri.vertex_degree(nv), 3);
                assert_eq!(tri.vertex_angle_sum(nv), TAU);
                assert_eq!(tri.vertex_location(nv), insertion.location);
            }
            Err(IntrinsicError::InsertionRefused { .. }) => {}
            Err(other) => panic!("{}: {other}", shape.name),
        }
        let e = (rng.next_u64() % tri.n_edges() as u64) as usize;
        let boundary = tri.is_boundary_edge(e);
        let (nv, ne, nf) = (tri.n_vertices(), tri.n_edges(), tri.n_faces());
        let length = tri.edge_length(e);
        let [a, b] = tri.edge_vertices(e);
        let t = rng.uniform(0.15, 0.85);
        match tri.split_edge(e, t) {
            Ok(insertion) => {
                done[1 + usize::from(boundary)] += 1;
                assert_eq!(insertion.vertex as usize, nv);
                let (n_edges, n_faces) = if boundary { (2, 1) } else { (3, 2) };
                assert_eq!((tri.n_edges(), tri.n_faces()), (ne + n_edges, nf + n_faces));
                assert_eq!(tri.vertex_degree(nv), if boundary { 3 } else { 4 });
                assert_eq!(tri.vertex_angle_sum(nv), if boundary { PI } else { TAU });
                assert_eq!(tri.is_boundary_vertex(nv), boundary);
                // The half from the tail keeps the id; the two halves add up.
                assert_eq!(tri.edge_vertices(e), [a, nv]);
                let other = insertion.new_edges[0] as usize;
                assert_eq!(tri.edge_vertices(other), [nv, b]);
                assert!((tri.edge_length(e) - t * length).abs() < 1e-15 * length);
                assert!(
                    (tri.edge_length(e) + tri.edge_length(other) - length).abs() < 1e-15 * length
                );
                assert_eq!(tri.edge_on_input(e), tri.edge_on_input(other));
            }
            Err(IntrinsicError::InsertionRefused { .. }) => {}
            Err(other) => panic!("{}: {other}", shape.name),
        }
    }
    done
}

/// Items 1 to 3 on a triangulation with inserted vertices. Returns the largest mismatch of
/// an edge trace over the edge's length.
fn assert_consistent(shape: &Shape, tri: &IntrinsicTriangulation, what: &str) -> f64 {
    let report = tri
        .check()
        .unwrap_or_else(|e| panic!("{} ({what}): {e}", shape.name));
    let context = format!("{} ({what}): {report:?}", shape.name);
    assert!(report.min_triangle_slack > 0.0, "{context}");
    assert!(report.relative_area_error.abs() < 1e-11, "{context}");
    assert!(report.max_angle_sum_error < 1e-9, "{context}");
    assert!(report.max_signpost_error < 1e-9, "{context}");
    assert!(report.gauss_bonnet_residual.abs() < 1e-9, "{context}");
    assert_eq!(
        tri.euler_characteristic(),
        shape.mesh.euler_characteristic(),
        "{context}"
    );
    // Inserted vertices are flat, and are where their places say.
    for v in shape.mesh.n_vertices()..tri.n_vertices() {
        let flat = if tri.is_boundary_vertex(v) { PI } else { TAU };
        assert_eq!(tri.vertex_angle_sum(v), flat, "{context}: vertex {v}");
        match tri.vertex_location(v) {
            SurfacePoint::Vertex(_) => panic!("{context}: inserted vertex {v} at an input vertex"),
            SurfacePoint::Edge { edge, t } => {
                assert!(t > 0.0 && t < 1.0, "{context}");
                assert_eq!(
                    tri.is_boundary_vertex(v),
                    shape.mesh.is_boundary_edge(edge as usize),
                    "{context}: vertex {v}"
                );
            }
            SurfacePoint::Face { bary, .. } => {
                assert!(
                    bary.iter().all(|&b| b > 0.0),
                    "{context}: vertex {v}: {bary:?}"
                );
                assert!(!tri.is_boundary_vertex(v), "{context}");
            }
        }
    }
    // Every edge traces to the place of its far vertex; its path has its length in space.
    let traces = tri
        .trace_edges(&EdgeTraceOptions::default())
        .unwrap_or_else(|e| panic!("{} ({what}): {e}", shape.name));
    assert!(
        traces.report.max_relative_mismatch < 1e-8,
        "{context}: {:?}",
        traces.report
    );
    for e in 0..tri.n_edges() {
        let [a, b] = tri.edge_vertices(e);
        let length = tri.edge_length(e);
        let in_space = path_length(
            shape,
            tri.vertex_location(a),
            traces.crossings(e),
            tri.vertex_location(b),
        );
        assert!(
            (in_space - length).abs() <= 1e-8 * length,
            "{context}: edge {e} = ({a}, {b}) has length {length}, its path {in_space}"
        );
        if shape.flat {
            let straight = dist(
                position(shape, tri.vertex_location(a)),
                position(shape, tri.vertex_location(b)),
            );
            assert!((length - straight).abs() < 1e-10, "{context}: edge {e}");
        }
    }
    traces.report.max_relative_mismatch
}

#[test]
fn insertions_and_splits_keep_the_surface() {
    for (s, shape) in standard_shapes().iter().enumerate() {
        let mut tri = shape.triangulation();
        let mut rng = SplitMix64::new(700 + s as u64);
        let mut total = [0usize; 3];
        for round in 0..3 {
            let done = random_insertions(shape, &mut tri, 0.3, &mut rng);
            for k in 0..3 {
                total[k] += done[k];
            }
            assert_consistent(shape, &tri, &format!("insertions, round {round}"));
            let report = tri.flip_to_delaunay(&DelaunayOptions::default()).unwrap();
            assert_eq!(report.n_not_delaunay, 0, "{}", shape.name);
            assert_consistent(shape, &tri, &format!("insertions and flips, round {round}"));
        }
        assert!(total[0] > 0 && total[1] > 0, "{}: {total:?}", shape.name);
        if shape.mesh.n_boundary_edges() > 0 {
            assert!(total[2] > 0, "{}: no boundary edge was split", shape.name);
        }
        // Boundary edges stay boundary edges along input boundary edges.
        for e in 0..tri.n_edges() {
            if tri.is_boundary_edge(e) {
                let input_edge = tri
                    .edge_on_input(e)
                    .expect("a boundary edge lies on the input");
                assert!(shape.mesh.is_boundary_edge(input_edge), "{}", shape.name);
            }
        }
    }
}

/// On a flat grid the places and lengths are known exactly.
#[test]
fn a_vertex_is_placed_where_its_coordinates_say() {
    let shape = Shape::generated("grid", generators::grid_disk(4, 4).unwrap(), true);
    let mut tri = shape.triangulation();
    // In face 5, by barycentric coordinates.
    let bary = [0.2, 0.5, 0.3];
    let corners = shape.mesh.faces()[5];
    let insertion = tri.insert_vertex(5, bary).unwrap();
    let expected = position(&shape, SurfacePoint::Face { face: 5, bary });
    let placed = position(&shape, insertion.location);
    assert!(dist(expected, placed) < 1e-15, "{:?}", insertion.location);
    assert!(matches!(
        insertion.location,
        SurfacePoint::Face { face: 5, .. }
    ));
    assert_eq!(insertion.snapped, None);
    for (k, &e) in insertion.new_edges.iter().enumerate() {
        let [from, to] = tri.edge_vertices(e as usize);
        assert_eq!((from, to), (corners[k] as usize, insertion.vertex as usize));
        let length = dist(shape.positions[from], expected);
        assert!((tri.edge_length(e as usize) - length).abs() < 1e-15);
    }
    // The signposts of the new vertex are angles in the frame of input face 5.
    let origin = shape.positions[corners[0] as usize];
    let axis = shape.positions[corners[1] as usize];
    let zero = (axis[1] - origin[1]).atan2(axis[0] - origin[0]);
    for h in tri.outgoing_halfedges(insertion.vertex as usize) {
        let to = shape.positions[tri.to_vertex(h)];
        let world = (to[1] - expected[1]).atan2(to[0] - expected[0]);
        let d = (tri.signpost(h) - (world - zero)).rem_euclid(TAU);
        assert!(d.min(TAU - d) < 1e-14, "halfedge {h}: {}", tri.signpost(h));
    }
    // An edge that lies along an input edge: the new vertex is on that input edge, exactly.
    let e = shape.mesh.find_edge(6, 7).unwrap();
    let split = tri.split_edge(e, 0.25).unwrap();
    let [tail, _] = tri.edge_vertices(e);
    let t = if tail == 6 { 0.25 } else { 0.75 };
    assert_eq!(split.location, SurfacePoint::Edge { edge: e as u32, t });
    // Split again: still on the input edge, by interpolation.
    let again = tri.split_edge(e, 0.5).unwrap();
    let t = if tail == 6 { 0.125 } else { 0.875 };
    assert_eq!(again.location, SurfacePoint::Edge { edge: e as u32, t });
    assert_eq!(tri.edge_on_input(again.new_edges[0] as usize), Some(e));
    // A boundary edge: the vertex is a boundary vertex with directions from 0 to π.
    let b = (0..shape.mesh.n_edges())
        .find(|&b| shape.mesh.is_boundary_edge(b))
        .unwrap();
    let on_boundary = tri.split_edge(b, 0.4).unwrap();
    let v = on_boundary.vertex as usize;
    assert!(tri.is_boundary_vertex(v));
    let signposts: Vec<f64> = tri.outgoing_halfedges(v).map(|h| tri.signpost(h)).collect();
    assert_eq!(signposts.len(), 2);
    assert_eq!(signposts[0], 0.0);
    assert!(signposts[1] > 0.0 && signposts[1] < PI);
    tri.check().unwrap();
    tri.trace_edges(&EdgeTraceOptions::default()).unwrap();
}

/// A traced place that is on an input edge up to rounding is stored as a point of that edge.
#[test]
fn a_place_on_an_input_edge_is_snapped_to_it() {
    // The centre of a square cell lies on the cell's diagonal, an input edge. After the
    // diagonal is flipped, the centre is the midpoint of the new edge, which crosses it.
    let shape = Shape::generated("grid", generators::grid_disk(3, 3).unwrap(), true);
    let mut tri = shape.triangulation();
    let diagonal = shape.mesh.find_edge(5, 10).unwrap();
    assert!(!shape.mesh.is_boundary_edge(diagonal));
    tri.flip_edge(diagonal).unwrap();
    assert_eq!(tri.edge_on_input(diagonal), None);
    let insertion = tri.split_edge(diagonal, 0.5).unwrap();
    match insertion.location {
        SurfacePoint::Edge { edge, t } => {
            assert_eq!(edge as usize, diagonal);
            assert!((t - 0.5).abs() < 1e-15);
        }
        other => panic!("{other}"),
    }
    let moved = insertion
        .snapped
        .expect("the place was moved onto the edge");
    assert!(moved.abs() <= 1e-12);
    tri.check().unwrap();
    let traces = tri.trace_edges(&EdgeTraceOptions::default()).unwrap();
    assert!(traces.report.max_relative_mismatch < 1e-14);
    // The two halves of the input diagonal are edges of the triangulation again, and are
    // recognized by their traces.
    let v = insertion.vertex as usize;
    let along: Vec<u32> = tri
        .outgoing_halfedges(v)
        .filter_map(|h| traces.along_input_edge(tri.edge(h)))
        .collect();
    assert_eq!(along, vec![diagonal as u32, diagonal as u32]);
}

#[test]
fn insertions_that_would_make_thin_triangles_are_refused() {
    let shape = common::jittered_grid(4, 0.2, 5);
    let mut tri = shape.triangulation();
    let before = tri.clone();
    // Next to side 1 of face 3 (the weight of corner 0 is tiny).
    let near = [1e-7, 0.4, 0.6 - 1e-7];
    match tri.insert_vertex(3, near).unwrap_err() {
        IntrinsicError::InsertionRefused {
            at,
            reason:
                InsertRefusal::TooThin {
                    shape: number,
                    nearest_edge,
                },
        } => {
            assert_eq!(
                at,
                SurfacePoint::Face {
                    face: 3,
                    bary: near
                }
            );
            assert!(number < FLIP_SHAPE_LIMIT && number > 0.0);
            assert_eq!(nearest_edge, Some(tri.edge(3 * 3 + 1) as u32));
        }
        other => panic!("{other}"),
    }
    let message = tri.insert_vertex(3, near).unwrap_err().to_string();
    let at = SurfacePoint::Face {
        face: 3,
        bary: near,
    };
    assert!(
        message.starts_with(&format!(
            "no vertex can be inserted at {at}: it would create a triangle too thin"
        )),
        "{message}"
    );
    assert!(message.contains("the point is next to edge "), "{message}");
    // Next to an end of an edge.
    let e = (0..tri.n_edges())
        .find(|&e| !tri.is_boundary_edge(e))
        .unwrap();
    assert!(matches!(
        tri.split_edge(e, 1e-9).unwrap_err(),
        IntrinsicError::InsertionRefused {
            reason: InsertRefusal::TooThin {
                nearest_edge: None,
                ..
            },
            ..
        }
    ));
    // Not points at all.
    for bary in [[0.5, 0.6, -0.1], [0.3, 0.3, 0.3], [f64::NAN, 0.5, 0.5]] {
        assert!(matches!(
            tri.insert_vertex(3, bary).unwrap_err(),
            IntrinsicError::InvalidPoint { .. }
        ));
    }
    assert!(matches!(
        tri.insert_vertex(tri.n_faces(), [0.3, 0.3, 0.4])
            .unwrap_err(),
        IntrinsicError::InvalidPoint { .. }
    ));
    for t in [0.0, 1.0, -0.5, 1.5, f64::NAN] {
        assert!(matches!(
            tri.split_edge(e, t).unwrap_err(),
            IntrinsicError::InvalidPoint { .. }
        ));
    }
    assert!(matches!(
        tri.split_edge(tri.n_edges(), 0.5).unwrap_err(),
        IntrinsicError::IndexOutOfRange { what: "edge", .. }
    ));
    assert_eq!(tri, before, "a refused insertion changed the triangulation");
    // An edge to a vertex of degree one has the same face on both sides.
    let needle = common::needle_tetrahedron(8.0, 4.0);
    let mut tri = needle.triangulation();
    tri.flip_to_delaunay(&DelaunayOptions::default()).unwrap();
    let apex = tri.irregularities().degree_one_vertices[0] as usize;
    let e = tri.edge(tri.outgoing_halfedge(apex));
    let error = tri.split_edge(e, 0.5).unwrap_err();
    assert!(matches!(
        error,
        IntrinsicError::InsertionRefused {
            reason: InsertRefusal::SameFaceOnBothSides,
            ..
        }
    ));
    assert!(
        error
            .to_string()
            .contains("both sides of the edge are the same triangle")
    );
}

/// Point queries: a point of the input mesh, found on the triangulation and found back.
#[test]
fn point_queries_are_inverse_to_each_other() {
    for (s, shape) in standard_shapes().iter().enumerate() {
        let mut tri = shape.triangulation();
        let mut rng = SplitMix64::new(4100 + s as u64);
        tri.flip_to_delaunay(&DelaunayOptions::default()).unwrap();
        random_insertions(shape, &mut tri, 0.2, &mut rng);
        tri.flip_to_delaunay(&DelaunayOptions::default()).unwrap();
        let scale = shape.mean_edge_length();
        let mesh = &shape.mesh;
        for query in 0..200 {
            // A point of the input mesh: in a face, on an edge, a vertex.
            let on_input = match query % 3 {
                0 => SurfacePoint::Face {
                    face: (rng.next_u64() % mesh.n_faces() as u64) as u32,
                    bary: random_bary(&mut rng),
                },
                1 => SurfacePoint::Edge {
                    edge: (rng.next_u64() % mesh.n_edges() as u64) as u32,
                    t: rng.uniform(0.05, 0.95),
                },
                _ => SurfacePoint::Vertex((rng.next_u64() % mesh.n_vertices() as u64) as u32),
            };
            let on_triangulation = tri
                .point_on_triangulation(on_input)
                .unwrap_or_else(|e| panic!("{}: {on_input}: {e}", shape.name));
            let back = tri
                .point_on_input(on_triangulation)
                .unwrap_or_else(|e| panic!("{}: {on_triangulation}: {e}", shape.name));
            let gap = dist(position(shape, on_input), position(shape, back));
            assert!(
                gap < 1e-9 * scale,
                "{}: {on_input} -> {on_triangulation} -> {back}: {gap}",
                shape.name
            );
            if let SurfacePoint::Vertex(v) = on_input {
                assert_eq!(on_triangulation, SurfacePoint::Vertex(v));
                assert_eq!(back, on_input);
            }
        }
        // A point of the triangulation, found on the input mesh and found back: the same
        // place, measured on the input mesh.
        for query in 0..200 {
            let on_triangulation = match query % 3 {
                0 => SurfacePoint::Face {
                    face: (rng.next_u64() % tri.n_faces() as u64) as u32,
                    bary: random_bary(&mut rng),
                },
                1 => SurfacePoint::Edge {
                    edge: (rng.next_u64() % tri.n_edges() as u64) as u32,
                    t: rng.uniform(0.05, 0.95),
                },
                _ => SurfacePoint::Vertex((rng.next_u64() % tri.n_vertices() as u64) as u32),
            };
            let on_input = tri
                .point_on_input(on_triangulation)
                .unwrap_or_else(|e| panic!("{}: {on_triangulation}: {e}", shape.name));
            if let SurfacePoint::Vertex(v) = on_triangulation {
                assert_eq!(on_input, tri.vertex_location(v as usize));
                continue;
            }
            let there = tri
                .point_on_triangulation(on_input)
                .unwrap_or_else(|e| panic!("{}: {on_input}: {e}", shape.name));
            let again = tri.point_on_input(there).unwrap();
            let gap = dist(position(shape, on_input), position(shape, again));
            assert!(
                gap < 1e-9 * scale,
                "{}: {on_triangulation} -> {on_input} -> {there} -> {again}: {gap}",
                shape.name
            );
        }
    }
}

#[test]
fn invalid_point_queries_are_errors() {
    let shape = common::jittered_grid(3, 0.2, 2);
    let tri = shape.triangulation();
    for point in [
        SurfacePoint::Vertex(1000),
        SurfacePoint::Edge { edge: 1000, t: 0.5 },
        SurfacePoint::Edge { edge: 1, t: 2.0 },
        SurfacePoint::Face {
            face: 1000,
            bary: [0.3, 0.3, 0.4],
        },
        SurfacePoint::Face {
            face: 1,
            bary: [0.3, 0.3, 0.3],
        },
    ] {
        for result in [tri.point_on_input(point), tri.point_on_triangulation(point)] {
            assert!(matches!(
                result.unwrap_err(),
                IntrinsicError::InvalidPoint { .. }
            ));
        }
    }
    // Ends of edges and corners of faces are the vertices.
    let [a, b] = tri.edge_vertices(4);
    assert_eq!(
        tri.point_on_input(SurfacePoint::Edge { edge: 4, t: 0.0 })
            .unwrap(),
        SurfacePoint::Vertex(a as u32)
    );
    assert_eq!(
        tri.point_on_input(SurfacePoint::Edge { edge: 4, t: 1.0 })
            .unwrap(),
        SurfacePoint::Vertex(b as u32)
    );
    let corner = tri.faces()[2][1];
    assert_eq!(
        tri.point_on_triangulation(SurfacePoint::Face {
            face: 2,
            bary: [0.0, 1.0, 0.0]
        })
        .unwrap(),
        SurfacePoint::Vertex(corner)
    );
}
