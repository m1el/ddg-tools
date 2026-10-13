//! The edge flip (SPEC §4): what it changes, that flipping twice restores lengths and angles,
//! and every refusal.

mod common;

use std::f64::consts::TAU;

use common::{Shape, dist, needle_tetrahedron, standard_shapes};
use ddg_intrinsic::formulas::unfold_diamond;
use ddg_intrinsic::{DelaunayOptions, FlipRefusal, IntrinsicError, IntrinsicTriangulation};

fn angle_distance(a: f64, b: f64) -> f64 {
    let d = (a - b).abs() % TAU;
    d.min(TAU - d)
}

/// The signposts of the two directions of edge `e`, keyed by the vertex they leave.
fn edge_signposts(tri: &IntrinsicTriangulation, e: usize) -> Vec<(usize, f64)> {
    let h = tri.edge_halfedge(e);
    let mut out = vec![(tri.from_vertex(h), tri.signpost(h))];
    if let Some(t) = tri.twin(h) {
        out.push((tri.from_vertex(t), tri.signpost(t)));
    }
    out
}

#[test]
fn a_flip_replaces_the_diagonal() {
    // Two triangles on the diagonal (0, 2) of a planar quadrilateral.
    let positions = vec![
        [0.0, 0.0, 0.0],
        [2.0, -1.0, 0.0],
        [3.0, 0.5, 0.0],
        [1.0, 1.5, 0.0],
    ];
    let shape = Shape::new("quad", positions.clone(), vec![[0, 1, 2], [0, 2, 3]], true);
    let mut tri = shape.triangulation();
    let e = shape.mesh.find_edge(0, 2).unwrap();
    assert_eq!(tri.can_flip(e).unwrap(), Ok(()));
    let flip = tri.flip_edge(e).unwrap();
    assert_eq!(flip.edge as usize, e);
    assert_eq!(flip.old_length, dist(positions[0], positions[2]));
    assert!((flip.new_length - dist(positions[1], positions[3])).abs() < 1e-15);
    // The halfedge 0 → 2 had vertex 3 on its left and 1 on its right; turned counter-clockwise
    // it runs from 1 to 3.
    let old = [
        flip.old_vertices[0].min(flip.old_vertices[1]),
        flip.old_vertices[0].max(flip.old_vertices[1]),
    ];
    assert_eq!(old, [0, 2]);
    let h = tri.edge_halfedge(e);
    assert_eq!(
        [tri.from_vertex(h), tri.to_vertex(h)],
        flip.new_vertices.map(|v| v as usize)
    );
    let mut ends = flip.new_vertices;
    ends.sort_unstable();
    assert_eq!(ends, [1, 3]);
    assert_eq!(tri.edge_length(e), flip.new_length);
    assert_eq!(tri.edge_on_input(e), None);
    assert_eq!((tri.n_vertices(), tri.n_edges(), tri.n_faces()), (4, 5, 2));
    let report = tri.check().unwrap();
    assert!(report.max_signpost_error < 1e-15 && report.max_angle_sum_error < 1e-15);
    // The faces are now (0, 1, 3) and (1, 2, 3), up to rotation.
    let mut faces: Vec<[u32; 3]> = tri
        .faces()
        .iter()
        .map(|f| {
            let k = (0..3).min_by_key(|&k| f[k]).unwrap();
            [f[k], f[(k + 1) % 3], f[(k + 2) % 3]]
        })
        .collect();
    faces.sort_unstable();
    assert_eq!(faces, vec![[0, 1, 3], [1, 2, 3]]);
    // The four boundary edges are untouched.
    for other in (0..5).filter(|&o| o != e) {
        assert_eq!(tri.edge_on_input(other), Some(other));
        assert_eq!(tri.edge_length(other), shape.geometry.edge_lengths()[other]);
    }
}

/// The five lengths of the two triangles on interior edge `e`, the longest of them, and the
/// smaller of the two heights over the edge.
fn diamond(tri: &IntrinsicTriangulation, e: usize) -> (f64, f64) {
    let h = tri.edge_halfedge(e);
    let t = tri.twin(h).unwrap();
    let l = |g: usize| tri.edge_length(tri.edge(g));
    let lengths = [
        l(h),
        l(tri.next(h)),
        l(tri.prev(h)),
        l(tri.next(t)),
        l(tri.prev(t)),
    ];
    let d = unfold_diamond(lengths[0], lengths[1], lengths[2], lengths[3], lengths[4]).unwrap();
    (
        lengths.iter().fold(0.0, |a: f64, &b| a.max(b)),
        d.h_k.min(d.h_l),
    )
}

/// SPEC §4: flipping twice restores lengths and angles to rounding. How much rounding there
/// is depends on the two triangles: an error `ε ℓ` in the length of the flipped edge moves the
/// height of a triangle over it by about `ε ℓ² / h`, which is large for a cap (`h ≪ ℓ`). With
/// `c = ℓ_max / h_min` over the edge before and after the first flip, the length must come
/// back within `1e-14 · c · ℓ_max` and the signposts within `1e-13 · c²` radians; on a
/// quadrilateral of good shape (`c < 10`) that is `1e-13 ℓ` and `1e-11`.
#[test]
fn flipping_twice_restores_lengths_and_signposts() {
    let mut well_shaped = 0;
    for shape in standard_shapes() {
        let mut tri = shape.triangulation();
        tri.flip_to_delaunay(&DelaunayOptions::default()).unwrap();
        let mut flipped = 0;
        for e in 0..tri.n_edges() {
            if tri.can_flip(e).unwrap().is_err() {
                continue;
            }
            let before = tri.clone();
            let (scale_before, height_before) = diamond(&tri, e);
            let first = tri.flip_edge(e).unwrap();
            let (scale_after, height_after) = diamond(&tri, e);
            // The reverse flip can be refused only if the first one made a degenerate
            // quadrilateral, which a strictly convex one does not.
            let second = tri.flip_edge(e).unwrap();
            flipped += 1;
            let scale = scale_before.max(scale_after);
            let condition = scale / height_before.min(height_after);
            well_shaped += usize::from(condition < 10.0);
            let context = format!(
                "{}: edge {e}: {first:?}, {second:?}, condition {condition:.1e}",
                shape.name
            );
            assert_eq!(
                second.new_vertices,
                [first.old_vertices[1], first.old_vertices[0]],
                "{context}"
            );
            assert_eq!(second.old_vertices, first.new_vertices, "{context}");
            assert!(
                (second.new_length - first.old_length).abs() <= 1e-14 * condition * scale,
                "{context}"
            );
            let angle_tolerance = 1e-13 * condition * condition;
            let mut was = edge_signposts(&before, e);
            let mut is = edge_signposts(&tri, e);
            was.sort_by(|a, b| a.partial_cmp(b).unwrap());
            is.sort_by(|a, b| a.partial_cmp(b).unwrap());
            if was[0].0 == was[1].0 {
                // A self-edge: both signposts leave the same vertex; compare as a set.
                let direct =
                    angle_distance(was[0].1, is[0].1).max(angle_distance(was[1].1, is[1].1));
                let swapped =
                    angle_distance(was[0].1, is[1].1).max(angle_distance(was[1].1, is[0].1));
                assert!(
                    direct.min(swapped) < angle_tolerance,
                    "{context}: {was:?} {is:?}"
                );
            } else {
                for (a, b) in was.iter().zip(&is) {
                    assert_eq!(a.0, b.0, "{context}");
                    assert!(
                        angle_distance(a.1, b.1) < angle_tolerance,
                        "{context}: {was:?} {is:?}"
                    );
                }
            }
            // Nothing else changed.
            for other in (0..tri.n_edges()).filter(|&o| o != e) {
                assert_eq!(
                    tri.edge_length(other),
                    before.edge_length(other),
                    "{context}"
                );
                assert_eq!(
                    edge_signposts(&tri, other),
                    edge_signposts(&before, other),
                    "{context}"
                );
            }
            // Restore exactly, so that rounding does not accumulate over the edges.
            tri = before;
        }
        assert!(flipped > 0, "{}", shape.name);
    }
    assert!(
        well_shaped > 1000,
        "only {well_shaped} well-shaped quadrilaterals were tested"
    );
}

#[test]
fn four_flips_restore_the_direction_too() {
    let shape = common::jittered_grid(6, 0.3, 1);
    let mut tri = shape.triangulation();
    for e in 0..tri.n_edges() {
        if tri.can_flip(e).unwrap().is_err() {
            continue;
        }
        let before = tri.clone();
        for _ in 0..4 {
            tri.flip_edge(e).unwrap();
        }
        assert_eq!(tri.edge_vertices(e), before.edge_vertices(e));
        assert!((tri.edge_length(e) - before.edge_length(e)).abs() < 1e-14);
        tri.check().unwrap();
        tri = before;
    }
}

#[test]
fn boundary_edges_are_refused() {
    for shape in standard_shapes() {
        let mut tri = shape.triangulation();
        for e in 0..tri.n_edges() {
            if !tri.is_boundary_edge(e) {
                continue;
            }
            let before = tri.clone();
            assert_eq!(tri.can_flip(e).unwrap(), Err(FlipRefusal::BoundaryEdge));
            let error = tri.flip_edge(e).unwrap_err();
            let [a, b] = tri.edge_vertices(e);
            assert_eq!(
                error,
                IntrinsicError::FlipRefused {
                    edge: e as u32,
                    vertices: [a as u32, b as u32],
                    reason: FlipRefusal::BoundaryEdge
                }
            );
            assert_eq!(
                error.to_string(),
                format!("edge {e} = ({a}, {b}) cannot be flipped: it is a boundary edge")
            );
            assert_eq!(tri, before, "a refused flip changed the triangulation");
        }
    }
}

#[test]
fn a_quadrilateral_that_is_not_strictly_convex_is_refused() {
    // Vertex 3 is moved until the quadrilateral 0 1 2 3 on the diagonal (0, 2) has a
    // nearly straight, a straight, then a reflex corner at 0.
    for (x3, flippable) in [
        (1.0, true),
        (1e-4, true),
        (1e-7, false),
        (0.0, false),
        (-1.0, false),
    ] {
        let positions = vec![
            [0.0, 0.0, 0.0],
            [0.0, -1.0, 0.0],
            [2.0, 0.0, 0.0],
            [x3, 1.0, 0.0],
        ];
        let shape = Shape::new("quad", positions, vec![[0, 1, 2], [0, 2, 3]], true);
        let mut tri = shape.triangulation();
        let e = shape.mesh.find_edge(0, 2).unwrap();
        let before = tri.clone();
        let result = tri.flip_edge(e);
        if flippable {
            result.unwrap();
            continue;
        }
        match result.unwrap_err() {
            IntrinsicError::FlipRefused {
                edge,
                vertices,
                reason:
                    FlipRefusal::NotConvex {
                        area_at_tail,
                        area_at_head,
                    },
            } => {
                assert_eq!(edge as usize, e);
                // The corner that is not convex is at vertex 0, where twice the signed area
                // of the triangle a flip would create is x3.
                let (at_zero, at_two) = if vertices[0] == 0 {
                    (area_at_tail, area_at_head)
                } else {
                    (area_at_head, area_at_tail)
                };
                assert!(
                    (at_zero - x3).abs() < 1e-14 && at_two > 0.5,
                    "{at_zero} {at_two}"
                );
            }
            // With the corner exactly straight the signed area is zero up to rounding; if it
            // comes out positive the flip is refused because the new triangle is too thin.
            IntrinsicError::FlipRefused {
                reason: FlipRefusal::TooThin { new_length, shape },
                ..
            } => {
                // Vertex 0 would lie x3 / 2 from the new edge of length 2: the triangle
                // (0, 1, 3) would have the shape number x3 / 4.
                assert!(x3 == 0.0 || x3 == 1e-7);
                assert!((new_length - 2.0).abs() < 1e-13, "{new_length}");
                assert!((shape - x3 / 4.0).abs() < 1e-15, "{shape}");
                assert!(shape < ddg_intrinsic::formulas::FLIP_SHAPE_LIMIT);
            }
            other => panic!("{other}"),
        }
        assert_eq!(tri, before);
        let message = tri.flip_edge(e).unwrap_err().to_string();
        assert!(message.contains(&format!("edge {e} = ")), "{message}");
        if x3 < 0.0 {
            assert!(message.contains("strictly convex"), "{message}");
        } else if x3 > 0.0 {
            assert!(message.contains("too thin"), "{message}");
        }
    }
}

/// On the tilted needle tetrahedron Delaunay flipping leaves the apex with a single edge: a
/// vertex of degree one inside a triangle that is glued to itself, a self-edge around it, and
/// two edges between the same pair of vertices.
#[test]
fn an_edge_with_an_end_of_degree_one_is_refused() {
    let shape = needle_tetrahedron(8.0, 4.0);
    let mut tri = shape.triangulation();
    tri.flip_to_delaunay(&DelaunayOptions::default()).unwrap();
    let irregular = tri.irregularities();
    assert_eq!(irregular.degree_one_vertices, vec![3], "{irregular:?}");
    assert_eq!(irregular.self_edges.len(), 1, "{irregular:?}");
    assert!(!irregular.faces_with_repeated_vertex.is_empty());
    let apex = 3;
    assert_eq!(tri.vertex_degree(apex), 1);
    let h = tri.outgoing_halfedge(apex);
    let e = tri.edge(h);
    // Both sides of the edge are the same face.
    assert_eq!(tri.face(h), tri.face(tri.twin(h).unwrap()));
    let before = tri.clone();
    assert_eq!(
        tri.can_flip(e).unwrap(),
        Err(FlipRefusal::DegreeOneEndpoint {
            vertex: apex as u32
        })
    );
    let error = tri.flip_edge(e).unwrap_err();
    assert!(matches!(
        error,
        IntrinsicError::FlipRefused {
            reason: FlipRefusal::DegreeOneEndpoint { vertex: 3 },
            ..
        }
    ));
    assert!(
        error.to_string().contains("end point 3 has degree one"),
        "{error}"
    );
    assert_eq!(tri, before);
    // The self-edge goes around the apex; the opposite angles of the edge to the apex lie in
    // one triangle and sum to less than π.
    let loop_edge = irregular.self_edges[0] as usize;
    let [a, b] = tri.edge_vertices(loop_edge);
    assert_eq!(a, b);
    assert!(tri.edge_angle_excess(e).unwrap() < 0.0);
    tri.check().unwrap();
    // The loop itself can be flipped back (its two triangles are distinct).
    assert_ne!(
        tri.face(tri.edge_halfedge(loop_edge)),
        tri.face(tri.twin(tri.edge_halfedge(loop_edge)).unwrap())
    );
}

#[test]
fn an_edge_that_does_not_exist_is_an_error() {
    let shape = common::jittered_grid(3, 0.2, 9);
    let mut tri = shape.triangulation();
    let n = tri.n_edges();
    let expected = IntrinsicError::IndexOutOfRange {
        what: "edge",
        index: n,
        bound: n,
    };
    assert_eq!(tri.flip_edge(n).unwrap_err(), expected);
    assert_eq!(tri.can_flip(n).unwrap_err(), expected);
    assert_eq!(
        expected.to_string(),
        format!("edge {n} is out of range (must be < {n})")
    );
}
