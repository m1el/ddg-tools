//! Vertex repositioning (paper Sec. 3.3.3, Alg. 10; SPEC N19, N20) with the invariants of
//! SPEC §5 items 1 to 3, 6 and 9: the surface does not change, the moved vertex is where its
//! lengths say, point queries and the overlay stay consistent, input vertices never move, and
//! what cannot be done is refused with nothing changed.

mod common;

use std::f64::consts::{PI, TAU};

use common::{Shape, dist, path_length, position, standard_shapes};
use ddg_core::generators::{self, SplitMix64};
use ddg_intrinsic::{
    DelaunayOptions, EdgeTraceOptions, IntrinsicError, IntrinsicTriangulation, MoveRefusal,
    SurfacePoint,
};

fn random_bary(rng: &mut SplitMix64) -> [f64; 3] {
    let (mut a, mut b) = (rng.next_f64(), rng.next_f64());
    if a + b > 1.0 {
        (a, b) = (1.0 - a, 1.0 - b);
    }
    let shrink = |x: f64| 0.05 + 0.85 * x;
    [shrink(a), shrink(b), shrink(1.0 - a - b)]
}

/// Inserts vertices in random faces and splits random edges (boundary edges among them).
fn insert_some(shape: &Shape, tri: &mut IntrinsicTriangulation, n: usize, rng: &mut SplitMix64) {
    for _ in 0..n {
        let f = (rng.next_u64() % tri.n_faces() as u64) as usize;
        match tri.insert_vertex(f, random_bary(rng)) {
            Ok(_) | Err(IntrinsicError::InsertionRefused { .. }) => {}
            Err(other) => panic!("{}: {other}", shape.name),
        }
        let e = (rng.next_u64() % tri.n_edges() as u64) as usize;
        match tri.split_edge(e, rng.uniform(0.2, 0.8)) {
            Ok(_) | Err(IntrinsicError::InsertionRefused { .. }) => {}
            Err(other) => panic!("{}: {other}", shape.name),
        }
    }
}

/// SPEC §5 items 1 to 3: the triangulation describes the input surface, every edge traces to
/// the place of its far vertex and its path has its length in space, and on a flat surface
/// the length is the distance of the places. Returns the largest relative mismatch of a trace.
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
    // The angle sums at the input vertices are those of the input, computed from the lengths.
    let sums = tri.corner_angle_sums();
    for (v, sum) in sums.iter().enumerate().take(shape.mesh.n_vertices()) {
        let input = tri.input().vertex_angle_sum(v);
        assert!((sum - input).abs() < 1e-9, "{context}: vertex {v}: {sum} != {input}");
        assert_eq!(tri.vertex_location(v), SurfacePoint::Vertex(v as u32), "{context}");
    }
    for v in shape.mesh.n_vertices()..tri.n_vertices() {
        let flat = if tri.is_boundary_vertex(v) { PI } else { TAU };
        assert_eq!(tri.vertex_angle_sum(v), flat, "{context}: vertex {v}");
        match tri.vertex_location(v) {
            SurfacePoint::Vertex(u) => panic!("{context}: vertex {v} is at input vertex {u}"),
            SurfacePoint::Edge { edge, t } => {
                assert!(t > 0.0 && t < 1.0, "{context}");
                assert_eq!(
                    tri.is_boundary_vertex(v),
                    shape.mesh.is_boundary_edge(edge as usize),
                    "{context}: vertex {v}"
                );
            }
            SurfacePoint::Face { bary, .. } => {
                assert!(bary.iter().all(|&b| b >= 0.0), "{context}: {bary:?}");
                assert!(!tri.is_boundary_vertex(v), "{context}");
            }
        }
    }
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
    // The overlay is a subdivision of the input with its area.
    let overlay = tri
        .overlay(&traces)
        .unwrap_or_else(|e| panic!("{context}: {e}"));
    let mut area = 0.0;
    for f in 0..overlay.n_faces() {
        let corners: Vec<[f64; 3]> = overlay
            .face(f)
            .iter()
            .map(|&k| position(shape, overlay.vertices[k as usize].on_input))
            .collect();
        for k in 1..corners.len() - 1 {
            let (u, w) = (
                [0, 1, 2].map(|i| corners[k][i] - corners[0][i]),
                [0, 1, 2].map(|i| corners[k + 1][i] - corners[0][i]),
            );
            let c = [
                u[1] * w[2] - u[2] * w[1],
                u[2] * w[0] - u[0] * w[2],
                u[0] * w[1] - u[1] * w[0],
            ];
            area += 0.5 * (c[0] * c[0] + c[1] * c[1] + c[2] * c[2]).sqrt();
        }
    }
    let total = shape.geometry.total_area();
    assert!((area - total).abs() < 1e-9 * total, "{context}: overlay area {area}");
    traces.report.max_relative_mismatch
}

/// Point queries are inverse to each other after moves (SPEC N11).
fn assert_point_queries(shape: &Shape, tri: &IntrinsicTriangulation, rng: &mut SplitMix64) {
    let scale = shape.geometry.edge_lengths().iter().sum::<f64>()
        / shape.geometry.edge_lengths().len() as f64;
    for _ in 0..40 {
        let f = (rng.next_u64() % tri.n_faces() as u64) as usize;
        let point = SurfacePoint::Face {
            face: f as u32,
            bary: random_bary(rng),
        };
        let on_input = tri
            .point_on_input(point)
            .unwrap_or_else(|e| panic!("{}: {e}", shape.name));
        let back = tri
            .point_on_triangulation(on_input)
            .unwrap_or_else(|e| panic!("{}: {e}", shape.name));
        let again = tri.point_on_input(back).unwrap();
        let d = dist(position(shape, on_input), position(shape, again));
        assert!(d < 1e-9 * scale, "{}: {point} -> {on_input} -> {back}: {d}", shape.name);
    }
}

/// A random displacement of the vertex inside its star (about `share` of the way to the
/// border in a random direction), along the boundary for a boundary vertex.
fn random_displacement(
    tri: &IntrinsicTriangulation,
    v: usize,
    share: f64,
    rng: &mut SplitMix64,
) -> [f64; 2] {
    let shortest = tri
        .outgoing_halfedges(v)
        .map(|h| tri.edge_length(tri.edge(h)))
        .fold(f64::INFINITY, f64::min);
    let r = share * shortest;
    if tri.is_boundary_vertex(v) {
        let s = if rng.next_f64() < 0.5 { -r } else { r };
        return [s, 0.0];
    }
    let angle = rng.uniform(0.0, TAU);
    [r * angle.cos(), r * angle.sin()]
}

#[test]
fn moves_keep_the_surface() {
    for (s, shape) in standard_shapes().iter().enumerate() {
        let mut tri = shape.triangulation();
        let mut rng = SplitMix64::new(4100 + s as u64);
        let n = tri.n_faces() / 4 + 2;
        insert_some(shape, &mut tri, n, &mut rng);
        tri.flip_to_delaunay(&DelaunayOptions::default()).unwrap();
        let mut counts = [0usize; 4]; // moved, beyond the paper's height bound, refused, boundary
        for round in 0..3 {
            for v in shape.mesh.n_vertices()..tri.n_vertices() {
                let displacement = random_displacement(&tri, v, 0.8, &mut rng);
                let before = tri.clone();
                match tri.move_vertex(v, displacement) {
                    Ok(moved) => {
                        counts[0] += 1;
                        counts[1] += usize::from(!moved.within_height_bound);
                        counts[3] += usize::from(tri.is_boundary_vertex(v));
                        assert!(moved.min_shape >= 1e-6, "{}", shape.name);
                        assert_eq!(tri.vertex_location(v), moved.to);
                    }
                    Err(IntrinsicError::MoveRefused {
                        reason: MoveRefusal::NotInStar { fraction, .. },
                        ..
                    }) => {
                        counts[2] += 1;
                        assert!(fraction <= 1.0 && fraction > 0.0, "{fraction}");
                        assert_eq!(tri, before, "{}: a refusal changed something", shape.name);
                    }
                    Err(IntrinsicError::MoveRefused {
                        reason: MoveRefusal::TooThin { .. } | MoveRefusal::SelfEdge { .. },
                        ..
                    }) => {
                        counts[2] += 1;
                        assert_eq!(tri, before, "{}", shape.name);
                    }
                    Err(other) => panic!("{}: vertex {v}: {other}", shape.name),
                }
            }
            assert_consistent(shape, &tri, &format!("moves, round {round}"));
            assert_point_queries(shape, &tri, &mut rng);
            let report = tri.flip_to_delaunay(&DelaunayOptions::default()).unwrap();
            assert_eq!(report.n_not_delaunay, 0, "{}", shape.name);
        }
        assert_consistent(shape, &tri, "moves and flips");
        eprintln!("{}: moved, beyond the height bound, refused, on the boundary: {counts:?}", shape.name);
        assert!(counts[0] > 0, "{}: {counts:?}", shape.name);
        let inserted_on_boundary = (shape.mesh.n_vertices()..tri.n_vertices())
            .filter(|&v| tri.is_boundary_vertex(v))
            .count();
        if inserted_on_boundary > 1 {
            assert!(counts[3] > 0, "{}: no boundary vertex was moved", shape.name);
        }
    }
}

#[test]
fn a_move_is_undone_by_moving_back() {
    for (s, shape) in standard_shapes().iter().enumerate().step_by(3) {
        let mut tri = shape.triangulation();
        let mut rng = SplitMix64::new(77 + s as u64);
        insert_some(shape, &mut tri, 12, &mut rng);
        tri.flip_to_delaunay(&DelaunayOptions::default()).unwrap();
        for v in shape.mesh.n_vertices()..tri.n_vertices() {
            let before = tri.clone();
            let displacement = random_displacement(&tri, v, 0.3, &mut rng);
            let Ok(moved) = tri.move_vertex(v, displacement) else {
                continue;
            };
            let back = [
                moved.distance * moved.back_direction.cos(),
                moved.distance * moved.back_direction.sin(),
            ];
            let back = if tri.is_boundary_vertex(v) { [back[0], 0.0] } else { back };
            tri.move_vertex(v, back)
                .unwrap_or_else(|e| panic!("{}: {e}", shape.name));
            let scale = before.edge_lengths().iter().copied().fold(0.0, f64::max);
            for e in 0..tri.n_edges() {
                assert!(
                    (tri.edge_length(e) - before.edge_length(e)).abs() < 1e-12 * scale,
                    "{}: edge {e}",
                    shape.name
                );
            }
            for h in 0..tri.n_halfedges() {
                let d = (tri.signpost(h) - before.signpost(h)).abs();
                assert!(d.min((d - TAU).abs()) < 1e-10, "{}: halfedge {h}: {d}", shape.name);
            }
            let (p, q) = (
                position(shape, tri.vertex_location(v)),
                position(shape, before.vertex_location(v)),
            );
            assert!(dist(p, q) < 1e-12 * scale, "{}: vertex {v}", shape.name);
        }
    }
}

/// On a flat grid every place, length and direction is known in coordinates.
#[test]
fn on_a_flat_grid_the_moved_vertex_is_where_its_coordinates_say() {
    let shape = Shape::generated("grid", generators::grid_disk(5, 5).unwrap(), true);
    let mut tri = shape.triangulation();
    let face = 17;
    let bary = [0.3, 0.3, 0.4];
    let v = tri.insert_vertex(face, bary).unwrap().vertex as usize;
    tri.flip_to_delaunay(&DelaunayOptions::default()).unwrap();
    // The frame of a place: the world angle of its zero direction.
    let zero_of = |point: SurfacePoint| match point {
        SurfacePoint::Face { face, .. } => {
            let f = shape.mesh.faces()[face as usize];
            let (o, x) = (shape.positions[f[0] as usize], shape.positions[f[1] as usize]);
            (x[1] - o[1]).atan2(x[0] - o[0])
        }
        SurfacePoint::Edge { edge, .. } => {
            let [a, b] = shape.mesh.edge_vertices(edge as usize);
            let (o, x) = (shape.positions[a], shape.positions[b]);
            (x[1] - o[1]).atan2(x[0] - o[0])
        }
        SurfacePoint::Vertex(_) => unreachable!(),
    };
    let mut rng = SplitMix64::new(5);
    for step in 0..40 {
        let old = position(&shape, tri.vertex_location(v));
        let zero = zero_of(tri.vertex_location(v));
        let displacement = random_displacement(&tri, v, 0.6, &mut rng);
        let moved = match tri.move_vertex(v, displacement) {
            Ok(m) => m,
            Err(IntrinsicError::MoveRefused { .. }) => continue,
            Err(e) => panic!("{e}"),
        };
        let (c, s) = (zero.cos(), zero.sin());
        let expected = [
            old[0] + c * displacement[0] - s * displacement[1],
            old[1] + s * displacement[0] + c * displacement[1],
            0.0,
        ];
        let new = position(&shape, moved.to);
        assert!(dist(new, expected) < 1e-14, "step {step}: {new:?} != {expected:?}");
        // Lengths and signposts against coordinates, in the frame of the new place.
        let zero = zero_of(moved.to);
        for h in tri.outgoing_halfedges(v) {
            let to = position(&shape, tri.vertex_location(tri.to_vertex(h)));
            assert!((dist(to, new) - tri.edge_length(tri.edge(h))).abs() < 1e-14);
            let world = (to[1] - new[1]).atan2(to[0] - new[0]);
            let d = (tri.signpost(h) - (world - zero)).rem_euclid(TAU);
            assert!(d.min(TAU - d) < 1e-11, "step {step}, halfedge {h}: {:?} {d} {moved:?}", tri.vertex_location(v));
        }
        tri.flip_to_delaunay(&DelaunayOptions::default()).unwrap();
    }
    assert_consistent(&shape, &tri, "40 moves on a grid");
}

/// The update of the optimal Delaunay triangulation, computed in the frame of the vertex,
/// against coordinates, and against a finite difference of the energy of Chen and Holst (2011,
/// Lemma 4.2): `Σ_T |T| (x − c_T)` is `3/2` the gradient of `Σ_T (|x|² + |b|² + |c|²)|T| / 3`.
#[test]
fn the_odt_target_is_the_weighted_mean_of_the_circumcentres() {
    let shape = Shape::generated("grid", generators::grid_disk(5, 5).unwrap(), true);
    let mut tri = shape.triangulation();
    let v = tri.insert_vertex(12, [0.2, 0.5, 0.3]).unwrap().vertex as usize;
    tri.flip_to_delaunay(&DelaunayOptions::default()).unwrap();
    let x = position(&shape, tri.vertex_location(v));
    let SurfacePoint::Face { face, .. } = tri.vertex_location(v) else {
        panic!()
    };
    let f = shape.mesh.faces()[face as usize];
    let (o, a) = (shape.positions[f[0] as usize], shape.positions[f[1] as usize]);
    let zero = (a[1] - o[1]).atan2(a[0] - o[0]);
    let target = tri.odt_target(v, false).unwrap();
    let world = [
        zero.cos() * target[0] - zero.sin() * target[1],
        zero.sin() * target[0] + zero.cos() * target[1],
    ];
    // In coordinates: the star's link vertices, counter-clockwise.
    let link: Vec<[f64; 2]> = tri
        .outgoing_halfedges(v)
        .map(|h| {
            let p = position(&shape, tri.vertex_location(tri.to_vertex(h)));
            [p[0], p[1]]
        })
        .collect();
    let n = link.len();
    let star = |x: [f64; 2]| -> Vec<(f64, [f64; 2], [f64; 2])> {
        (0..n)
            .map(|k| {
                let (b, c) = (link[k], link[(k + 1) % n]);
                let (u, w) = ([b[0] - x[0], b[1] - x[1]], [c[0] - x[0], c[1] - x[1]]);
                (0.5 * (u[0] * w[1] - u[1] * w[0]), b, c)
            })
            .collect()
    };
    let x2 = [x[0], x[1]];
    let (mut sum, mut weighted) = (0.0, [0.0, 0.0]);
    for (area, b, c) in star(x2) {
        let (u, w) = ([b[0] - x2[0], b[1] - x2[1]], [c[0] - x2[0], c[1] - x2[1]]);
        let d = 2.0 * (u[0] * w[1] - u[1] * w[0]);
        let (uu, ww) = (u[0] * u[0] + u[1] * u[1], w[0] * w[0] + w[1] * w[1]);
        let centre = [(w[1] * uu - u[1] * ww) / d, (u[0] * ww - w[0] * uu) / d];
        sum += area;
        weighted[0] += area * centre[0];
        weighted[1] += area * centre[1];
    }
    let expected = [weighted[0] / sum, weighted[1] / sum];
    assert!(
        (world[0] - expected[0]).hypot(world[1] - expected[1]) < 1e-14,
        "{world:?} != {expected:?}"
    );
    // The energy and its gradient by central differences.
    let energy = |x: [f64; 2]| -> f64 {
        star(x)
            .iter()
            .map(|(area, b, c)| {
                let q = |p: [f64; 2]| p[0] * p[0] + p[1] * p[1];
                (q(x) + q(*b) + q(*c)) * area / 3.0
            })
            .sum()
    };
    let h = 1e-6;
    let gradient = [0, 1].map(|k| {
        let (mut p, mut m) = (x2, x2);
        p[k] += h;
        m[k] -= h;
        (energy(p) - energy(m)) / (2.0 * h)
    });
    // x_new − x = −(3 / 2|ω|) ∇E.
    let predicted = [-1.5 * gradient[0] / sum, -1.5 * gradient[1] / sum];
    assert!(
        (predicted[0] - expected[0]).hypot(predicted[1] - expected[1]) < 1e-7,
        "{predicted:?} != {expected:?}"
    );
}

#[test]
fn what_cannot_be_moved_is_refused_and_nothing_changes() {
    let shape = Shape::generated("grid", generators::grid_disk(4, 4).unwrap(), true);
    let mut tri = shape.triangulation();
    let v = tri.insert_vertex(5, [0.3, 0.3, 0.4]).unwrap().vertex as usize;
    let b = (0..shape.mesh.n_edges())
        .find(|&e| shape.mesh.is_boundary_edge(e))
        .unwrap();
    let w = tri.split_edge(b, 0.5).unwrap().vertex as usize;
    let before = tri.clone();
    let reason = |r: Result<_, IntrinsicError>| match r {
        Err(IntrinsicError::MoveRefused { reason, .. }) => reason,
        other => panic!("{other:?}"),
    };
    assert_eq!(reason(tri.move_vertex(0, [0.01, 0.0])), MoveRefusal::InputVertex);
    assert_eq!(
        reason(tri.move_vertex(v, [f64::NAN, 0.0])),
        MoveRefusal::NotFinite
    );
    match reason(tri.move_vertex(v, [10.0, 3.0])) {
        MoveRefusal::NotInStar { fraction, .. } => assert!(fraction > 0.0 && fraction < 0.1),
        other => panic!("{other:?}"),
    }
    assert_eq!(
        reason(tri.move_vertex(w, [0.01, 0.01])),
        MoveRefusal::LeavesBoundary { normal: 0.01 }
    );
    assert!(matches!(
        tri.move_vertex(tri.n_vertices(), [0.0, 0.0]),
        Err(IntrinsicError::IndexOutOfRange { .. })
    ));
    assert!(matches!(
        tri.odt_target(3, false),
        Err(IntrinsicError::MoveRefused {
            reason: MoveRefusal::InputVertex,
            ..
        })
    ));
    assert_eq!(tri, before);
    // The message names the vertex.
    let message = tri.move_vertex(v, [10.0, 3.0]).unwrap_err().to_string();
    assert!(message.contains(&format!("vertex {v}")), "{message}");
    // Just inside the star: accepted, beyond the paper's height bound.
    let MoveRefusal::NotInStar { fraction, .. } = reason(tri.move_vertex(v, [0.2, 0.05])) else {
        // Inside already: scale up until it is not.
        panic!("expected the move to leave the star");
    };
    let inside = [0.2 * 0.9 * fraction, 0.05 * 0.9 * fraction];
    let moved = tri.move_vertex(v, inside).unwrap();
    assert!(!moved.within_height_bound);
    assert_consistent(&shape, &tri, "a move beyond the height bound");
    // A zero move changes nothing.
    let now = tri.clone();
    tri.move_vertex(v, [0.0, 0.0]).unwrap();
    assert_eq!(tri, now);
    // Along the boundary: interpolated on the input edge.
    let SurfacePoint::Edge { edge, t } = tri.vertex_location(w) else {
        panic!()
    };
    let first = tri.outgoing_halfedge(w);
    let length = tri.edge_length(tri.edge(first));
    let moved = tri.move_vertex(w, [0.25 * length, 0.0]).unwrap();
    let SurfacePoint::Edge { edge: e2, t: t2 } = moved.to else {
        panic!()
    };
    assert_eq!(e2, edge);
    assert!((t2 - t).abs() > 0.0);
    assert_consistent(&shape, &tri, "a move along the boundary");
}
