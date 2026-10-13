//! SPEC §5 items 1, 2 and 3: after any sequence of flips the triangulation describes the same
//! surface (strict triangle inequalities, area, angle sums, Gauss–Bonnet), and its signposts
//! lead to the neighbouring vertices, on the intrinsic triangulation and on the input mesh.

mod common;

use std::f64::consts::TAU;

use common::{Shape, path_length, standard_shapes};
use ddg_core::generators::SplitMix64;
use ddg_core::operators::{gauss_bonnet_residual, vertex_angle_sums};
use ddg_intrinsic::{
    ConsistencyReport, DelaunayOptions, EdgeTraceOptions, IntrinsicError, IntrinsicTriangulation,
    SurfacePoint, TraceOptions, TraceStop,
};

/// How closely the identities of the surface must hold.
#[derive(Clone, Copy)]
struct Tolerances {
    /// Relative error of the total area.
    area: f64,
    /// Angle sums and signposts, in radians.
    angle: f64,
}

/// After Delaunay flipping, whose flips make triangles of good shape (shape numbers above
/// `4e-4` on every mesh here, so each flip moves an angle by `1e-12` at most).
const TIGHT: Tolerances = Tolerances {
    area: 1e-12,
    angle: 1e-11,
};

/// After flips of randomly chosen edges, which make triangles as thin as the flip allows
/// (`FLIP_SHAPE_LIMIT = 1e-6`: each such flip may move an angle by `1e-10`).
const AFTER_RANDOM_FLIPS: Tolerances = Tolerances {
    area: 1e-9,
    angle: 1e-9,
};

/// Items 1 and 2: triangle inequalities, area, angle sums, Gauss–Bonnet (to `1e-10`).
fn assert_same_surface(
    shape: &Shape,
    tri: &IntrinsicTriangulation,
    what: &str,
    tolerances: Tolerances,
) -> ConsistencyReport {
    let report = tri
        .check()
        .unwrap_or_else(|e| panic!("{} ({what}): {e}", shape.name));
    let context = format!("{} ({what}): {report:?}", shape.name);
    assert!(report.min_triangle_slack > 0.0, "{context}");
    assert!(
        report.relative_area_error.abs() < tolerances.area,
        "{context}"
    );
    assert!(report.max_angle_sum_error < tolerances.angle, "{context}");
    assert!(report.max_signpost_error < tolerances.angle, "{context}");
    assert!(report.gauss_bonnet_residual.abs() < 1e-10, "{context}");
    // The same against quantities computed by ddg_core on the input mesh.
    let input_sums = vertex_angle_sums(&shape.mesh, &shape.geometry).unwrap();
    for (v, (ours, theirs)) in tri.corner_angle_sums().iter().zip(&input_sums).enumerate() {
        assert!(
            (ours - theirs).abs() < tolerances.angle,
            "{context}: angle sum of vertex {v}"
        );
        // Stored when the triangulation was made; summed in another order than ddg_core does.
        assert!(
            (tri.vertex_angle_sum(v) - theirs).abs() < 1e-13,
            "{context}: stored angle sum of {v}"
        );
    }
    let total: f64 = tri.angle_defects().iter().sum::<f64>()
        + tri.boundary_geodesic_curvatures().iter().sum::<f64>();
    let input_residual = gauss_bonnet_residual(&shape.mesh, &shape.geometry).unwrap();
    let chi = shape.mesh.euler_characteristic();
    assert_eq!(tri.euler_characteristic(), chi, "{context}");
    assert!(
        (total - TAU * chi as f64 - input_residual).abs() < 1e-10,
        "{context}"
    );
    for &l in tri.edge_lengths() {
        assert!(l.is_finite() && l > 0.0, "{context}");
    }
    for &phi in tri.signposts() {
        assert!(phi.is_finite(), "{context}");
    }
    report
}

/// Flips `per_edge * n_edges` randomly chosen edges that can be flipped; returns the number
/// of flips made.
fn random_flips(tri: &mut IntrinsicTriangulation, per_edge: f64, rng: &mut SplitMix64) -> usize {
    let n = (per_edge * tri.n_edges() as f64) as usize;
    let mut done = 0;
    for _ in 0..20 * n {
        if done == n {
            break;
        }
        let e = (rng.next_u64() % tri.n_edges() as u64) as usize;
        match tri.flip_edge(e) {
            Ok(_) => done += 1,
            Err(IntrinsicError::FlipRefused { .. }) => {}
            Err(other) => panic!("{other}"),
        }
    }
    done
}

#[test]
fn a_fresh_triangulation_is_the_input_mesh() {
    for shape in standard_shapes() {
        let tri = shape.triangulation();
        let mesh = &shape.mesh;
        assert_eq!(tri.n_vertices(), mesh.n_vertices());
        assert_eq!(tri.n_edges(), mesh.n_edges());
        assert_eq!(tri.faces(), mesh.faces());
        for h in 0..mesh.n_halfedges() {
            assert_eq!(tri.twin(h), mesh.twin(h));
            assert_eq!(tri.edge(h), mesh.edge(h));
            assert_eq!(tri.from_vertex(h), mesh.from_vertex(h));
            assert_eq!(tri.signpost(h), tri.input().signpost(h));
        }
        for e in 0..mesh.n_edges() {
            assert_eq!(tri.edge_length(e), shape.geometry.edge_lengths()[e]);
            assert_eq!(tri.edge_on_input(e), Some(e));
            assert_eq!(tri.is_boundary_edge(e), mesh.is_boundary_edge(e));
        }
        for v in 0..mesh.n_vertices() {
            assert_eq!(tri.vertex_location(v), SurfacePoint::Vertex(v as u32));
            assert_eq!(tri.is_boundary_vertex(v), mesh.is_boundary_vertex(v));
            assert_eq!(tri.outgoing_halfedge(v), mesh.outgoing_halfedge(v));
            let ours: Vec<usize> = tri.outgoing_halfedges(v).collect();
            let theirs: Vec<usize> = mesh.outgoing_halfedges(v).collect();
            assert_eq!(ours, theirs);
            let neighbours = mesh.vertex_neighbors(v).count();
            assert_eq!(
                tri.vertex_degree(v),
                neighbours,
                "{}: vertex {v}",
                shape.name
            );
        }
        assert!(tri.irregularities().is_empty());
        let report = assert_same_surface(&shape, &tri, "fresh", TIGHT);
        assert!(report.max_angle_sum_error < 1e-13 && report.max_signpost_error < 1e-13);
        // The exact way back.
        let (mesh2, geometry2) = tri.to_tri_mesh().unwrap();
        assert_eq!(&mesh2, mesh);
        assert_eq!(geometry2, shape.geometry);
    }
}

#[test]
fn random_flips_keep_the_surface() {
    for (s, shape) in standard_shapes().iter().enumerate() {
        let mut tri = shape.triangulation();
        let mut rng = SplitMix64::new(100 + s as u64);
        let mut total = 0;
        for round in 0..6 {
            total += random_flips(&mut tri, 0.5, &mut rng);
            assert_same_surface(
                shape,
                &tri,
                &format!("random flips, round {round}"),
                AFTER_RANDOM_FLIPS,
            );
        }
        assert!(total > 0, "{}: no flip was possible", shape.name);
    }
}

#[test]
fn delaunay_flipping_keeps_the_surface() {
    for (s, shape) in standard_shapes().iter().enumerate() {
        // From the input, and from a triangulation made worse by random flips.
        let mut tri = shape.triangulation();
        tri.flip_to_delaunay(&DelaunayOptions::default()).unwrap();
        assert_same_surface(shape, &tri, "Delaunay", TIGHT);

        let mut tri = shape.triangulation();
        let mut rng = SplitMix64::new(7 + s as u64);
        random_flips(&mut tri, 2.0, &mut rng);
        tri.flip_to_delaunay(&DelaunayOptions::default()).unwrap();
        assert_same_surface(
            shape,
            &tri,
            "random flips, then Delaunay",
            AFTER_RANDOM_FLIPS,
        );
    }
}

/// Signposts on the intrinsic triangulation: tracing from the tail of a halfedge along its
/// signpost for the edge's length arrives at its head, exactly and without a crossing.
fn assert_signposts_lead_along_edges(shape: &Shape, tri: &IntrinsicTriangulation) {
    let options = TraceOptions {
        record_crossings: true,
        max_crossings: None,
    };
    for h in 0..tri.n_halfedges() {
        let (a, b) = (tri.from_vertex(h), tri.to_vertex(h));
        let e = tri.edge(h);
        let trace = tri
            .trace(
                SurfacePoint::Vertex(a as u32),
                tri.signpost(h),
                tri.edge_length(e),
                &options,
            )
            .unwrap_or_else(|err| panic!("{}: halfedge {h}: {err}", shape.name));
        let context = format!("{}: halfedge {h} = ({a}, {b}): {trace:?}", shape.name);
        assert_eq!(trace.end, SurfacePoint::Vertex(b as u32), "{context}");
        assert_eq!(trace.stop, TraceStop::Completed, "{context}");
        assert_eq!(trace.n_crossings, 0, "{context}");
        assert_eq!(trace.along_edge, Some(e as u32), "{context}");
        assert_eq!(trace.length, tri.edge_length(e), "{context}");
        // The way back is the signpost of the twin (π for the boundary).
        let expected = tri
            .twin(h)
            .map_or(std::f64::consts::PI, |t| tri.signpost(t));
        assert_eq!(trace.back_direction, expected, "{context}");
    }
}

/// Signposts on the input mesh: every edge traced over the input mesh arrives at its far
/// vertex within `1e-9` of its length, and the traced path has the length of the edge
/// (relative `1e-10`).
fn assert_edges_trace_to_their_ends(shape: &Shape, tri: &IntrinsicTriangulation, what: &str) {
    let traces = tri
        .trace_edges(&EdgeTraceOptions {
            tolerance: 1e-9,
            max_crossings: None,
        })
        .unwrap_or_else(|e| panic!("{} ({what}): {e}", shape.name));
    assert_eq!(traces.report.n_edges, tri.n_edges());
    assert_eq!(
        traces.report.n_on_input + traces.report.n_traced,
        tri.n_edges()
    );
    assert!(traces.report.max_relative_mismatch <= 1e-9);
    assert!(
        traces.report.max_clamp < 1e-9,
        "{} ({what}): {:?}",
        shape.name,
        traces.report
    );
    for e in 0..tri.n_edges() {
        let [a, b] = tri.edge_vertices(e);
        let crossings = traces.crossings(e);
        let context = format!("{} ({what}): edge {e} = ({a}, {b})", shape.name);
        if let Some(input_edge) = traces.along_input_edge(e) {
            assert!(crossings.is_empty(), "{context}");
            assert_eq!(traces.first_face(e), None, "{context}");
            let [p, q] = shape.mesh.edge_vertices(input_edge as usize);
            assert_eq!((a.min(b), a.max(b)), (p, q), "{context}");
        } else {
            assert!(traces.first_face(e).is_some(), "{context}");
        }
        let length = path_length(
            shape,
            tri.vertex_location(a),
            crossings,
            tri.vertex_location(b),
        );
        let expected = tri.edge_length(e);
        assert!(
            (length - expected).abs() <= 1e-10 * expected,
            "{context}: path length {length}, edge length {expected}"
        );
        let mut previous = 0.0;
        for c in crossings {
            assert!((0.0..=1.0).contains(&c.t), "{context}: {c:?}");
            assert!(
                c.distance >= previous - 1e-12 * expected,
                "{context}: {c:?}"
            );
            assert_eq!(
                shape.mesh.edge(c.halfedge as usize),
                c.edge as usize,
                "{context}"
            );
            previous = c.distance;
        }
        // Consecutive crossings lie on edges of one input face, the one the path is in.
        for pair in crossings.windows(2) {
            let entered = shape.mesh.twin(pair[0].halfedge as usize).unwrap() / 3;
            assert_eq!(entered, pair[1].halfedge as usize / 3, "{context}");
        }
        if let (Some(first), Some(face)) = (crossings.first(), traces.first_face(e)) {
            assert_eq!(first.halfedge / 3, face, "{context}");
        }
    }
}

#[test]
fn signposts_are_consistent() {
    for (s, shape) in standard_shapes().iter().enumerate() {
        let mut tri = shape.triangulation();
        assert_signposts_lead_along_edges(shape, &tri);
        assert_edges_trace_to_their_ends(shape, &tri, "fresh");

        let mut rng = SplitMix64::new(31 + s as u64);
        random_flips(&mut tri, 1.0, &mut rng);
        assert_signposts_lead_along_edges(shape, &tri);
        assert_edges_trace_to_their_ends(shape, &tri, "random flips");

        tri.flip_to_delaunay(&DelaunayOptions::default()).unwrap();
        assert_signposts_lead_along_edges(shape, &tri);
        assert_edges_trace_to_their_ends(shape, &tri, "Delaunay");
    }
}

/// On a flat surface every intrinsic edge is the straight segment between its end points:
/// its length is their Euclidean distance and its crossings lie on that segment.
#[test]
fn on_flat_shapes_intrinsic_edges_are_straight_segments() {
    for (s, shape) in standard_shapes().iter().enumerate() {
        if !shape.flat {
            continue;
        }
        let mut tri = shape.triangulation();
        let mut rng = SplitMix64::new(900 + s as u64);
        random_flips(&mut tri, 1.0, &mut rng);
        for delaunay in [false, true] {
            if delaunay {
                tri.flip_to_delaunay(&DelaunayOptions::default()).unwrap();
            }
            let traces = tri.trace_edges(&EdgeTraceOptions::default()).unwrap();
            for e in 0..tri.n_edges() {
                let [a, b] = tri.edge_vertices(e);
                let (pa, pb) = (shape.positions[a], shape.positions[b]);
                let straight = common::dist(pa, pb);
                assert!(
                    (tri.edge_length(e) - straight).abs() < 1e-12,
                    "{}: edge {e} = ({a}, {b}) has length {}, its end points are {straight} apart",
                    shape.name,
                    tri.edge_length(e)
                );
                for c in traces.crossings(e) {
                    let p = common::crossing_position(shape, c);
                    let along = c.distance / straight;
                    for k in 0..3 {
                        let on_segment = pa[k] + along * (pb[k] - pa[k]);
                        assert!(
                            (p[k] - on_segment).abs() < 1e-11,
                            "{}: edge {e}: {c:?}",
                            shape.name
                        );
                    }
                }
            }
        }
    }
}
