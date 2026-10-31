//! Tracing over a triangulation without stored signposts (`AngleFrames`, SPEC N50, IC-SPEC I9):
//! after random flips, insertions and splits, the signposts computed from the corner angles
//! agree with the stored ones (in the stored frame turned so that the outgoing halfedge of
//! each vertex points at angle 0), and traces over the triangulation from vertices and from
//! points of faces end at the same places.

mod common;

use std::f64::consts::TAU;

use common::{Shape, dist, inner_bary, position, standard_shapes};
use ddg_core::generators::SplitMix64;
use ddg_intrinsic::{AngleFrames, IntrinsicTriangulation, SurfacePoint, TraceOptions};

fn scramble(tri: &mut IntrinsicTriangulation, rng: &mut SplitMix64) {
    for _ in 0..2 * tri.n_edges() {
        let e = (rng.next_u64() % tri.n_edges() as u64) as usize;
        let _ = tri.flip_edge(e);
    }
    for _ in 0..10 {
        let f = (rng.next_u64() % tri.n_faces() as u64) as usize;
        let _ = tri.insert_vertex(f, inner_bary(rng));
        let e = (rng.next_u64() % tri.n_edges() as u64) as usize;
        let _ = tri.split_edge(e, rng.uniform(0.2, 0.8));
    }
    for _ in 0..tri.n_edges() {
        let e = (rng.next_u64() % tri.n_edges() as u64) as usize;
        let _ = tri.flip_edge(e);
    }
}

fn wrap(a: f64) -> f64 {
    a.rem_euclid(TAU)
}

/// The largest difference between computed and stored signposts (in the turned frame), and
/// between computed and stored angle sums relative to `2π`.
fn signpost_agreement(tri: &IntrinsicTriangulation) -> (f64, f64) {
    let frames = AngleFrames::of(tri);
    let (mut signposts, mut sums) = (0.0f64, 0.0f64);
    for v in 0..tri.n_vertices() {
        let base = tri.signpost(tri.outgoing_halfedge(v));
        let boundary = tri.is_boundary_vertex(v);
        for h in tri.outgoing_halfedges(v) {
            let stored = if boundary {
                tri.signpost(h)
            } else {
                wrap(tri.signpost(h) - base)
            };
            let d = (frames.signpost(h) - stored).abs();
            signposts = signposts.max(if boundary { d } else { d.min(TAU - d) });
        }
        sums = sums.max((frames.angle_sum(v) - tri.vertex_angle_sum(v)).abs() / TAU);
    }
    (signposts, sums)
}

fn end_position(shape: &Shape, tri: &IntrinsicTriangulation, p: SurfacePoint) -> [f64; 3] {
    position(shape, tri.point_on_input(p).unwrap())
}

#[test]
fn computed_signposts_and_traces_agree_with_the_stored_ones() {
    let options = TraceOptions::default();
    let (mut worst_signpost, mut worst_sum, mut worst_trace) = (0.0f64, 0.0f64, 0.0f64);
    let (mut n_traces, mut n_same_end) = (0usize, 0usize);
    for (i, shape) in standard_shapes().iter().enumerate() {
        let mut rng = SplitMix64::new(17 + i as u64);
        let mut tri = shape.triangulation();
        let fresh = signpost_agreement(&tri);
        scramble(&mut tri, &mut rng);
        tri.check().unwrap();
        let (signposts, sums) = signpost_agreement(&tri);
        println!(
            "{}: fresh {:.1e} / {:.1e}, after operations {signposts:.1e} / {sums:.1e}",
            shape.name, fresh.0, fresh.1
        );
        worst_signpost = worst_signpost.max(signposts).max(fresh.0);
        worst_sum = worst_sum.max(sums).max(fresh.1);

        let frames = AngleFrames::of(&tri);
        let scale = shape.mean_edge_length();
        for _ in 0..200 {
            let length = rng.uniform(0.1, 2.0) * scale;
            // From a vertex, in a direction given in the stored frame.
            let v = (rng.next_u64() % tri.n_vertices() as u64) as usize;
            let boundary = tri.is_boundary_vertex(v);
            let stored = if boundary {
                rng.uniform(0.0, std::f64::consts::PI)
            } else {
                rng.uniform(0.0, TAU)
            };
            let base = tri.signpost(tri.outgoing_halfedge(v));
            let computed = if boundary {
                stored
            } else {
                wrap(stored - base)
            };
            let start = SurfacePoint::Vertex(v as u32);
            let a = tri.trace(start, stored, length, &options).unwrap();
            let b = frames.trace(start, computed, length, &options).unwrap();
            // From a point of a face, in the frame of the face (no signpost at the start).
            let f = (rng.next_u64() % tri.n_faces() as u64) as usize;
            let start = SurfacePoint::Face {
                face: f as u32,
                bary: inner_bary(&mut rng),
            };
            let direction = rng.uniform(0.0, TAU);
            let c = tri.trace(start, direction, length, &options).unwrap();
            let d = frames.trace(start, direction, length, &options).unwrap();
            for (x, y) in [(&a, &b), (&c, &d)] {
                n_traces += 1;
                n_same_end += usize::from(x.end == y.end && x.stop == y.stop);
                assert_eq!(x.stop, y.stop, "{}", shape.name);
                let gap = dist(
                    end_position(shape, &tri, x.end),
                    end_position(shape, &tri, y.end),
                ) / scale;
                worst_trace = worst_trace.max(gap);
            }
        }
    }
    println!(
        "signposts {worst_signpost:.2e} rad, angle sums {worst_sum:.2e} of 2π, trace ends \
         {worst_trace:.2e} of the mean edge; {n_same_end} of {n_traces} ends bitwise equal"
    );
    // Measured (2026-10-31): signposts 6.2e-12 rad and angle sums 2.0e-12 of 2π at most, both
    // on the needle tetrahedron, whose stored angle sums are those of the input while the
    // corner angles of its flipped needles drift by rounding (`check().max_angle_sum_error`);
    // below 4e-12 on every other shape, 0 before any operation. Trace ends 6.2e-12 of the mean
    // edge, about half of them bitwise equal.
    assert!(
        worst_signpost < 1e-10,
        "signposts differ by {worst_signpost}"
    );
    assert!(worst_sum < 1e-10, "angle sums differ by {worst_sum}");
    assert!(worst_trace < 1e-9, "trace ends differ by {worst_trace}");
}

#[test]
fn bad_lengths_are_refused() {
    let shape = &standard_shapes()[0];
    let tri = shape.triangulation();
    let short = &tri.edge_lengths()[1..];
    assert!(AngleFrames::new(tri.complex(), short).is_err());
    let frames = AngleFrames::new(tri.complex(), tri.edge_lengths()).unwrap();
    let start = SurfacePoint::Vertex(0);
    assert!(
        frames
            .trace(start, 7.0, 1.0, &TraceOptions::default())
            .is_err()
    );
    assert!(
        frames
            .trace(start, 1.0, -1.0, &TraceOptions::default())
            .is_err()
    );
}
