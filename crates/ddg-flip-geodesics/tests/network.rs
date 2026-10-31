//! Networks: construction, layering, pinning, validity, refused inputs (SPEC §2.1, §2.2, §7
//! test 6).

mod common;

use common::*;
use ddg_flip_geodesics::{
    EdgePath, FlipError, GeodesicNetwork, Network, Segment, SegmentId, SegmentRef, dijkstra_path,
};

#[test]
fn construction_ids_layering_and_round_trip() {
    let shape = grid(4, 4);
    let tri = shape.triangulation();
    let g = |i, j| gv(4, i, j);
    let a = vpath(&tri, &[g(0, 2), g(1, 2), g(2, 2), g(3, 2)], false);
    let b = vpath(
        &tri,
        &[g(1, 1), g(2, 1), g(3, 1), g(3, 0), g(2, 0), g(1, 0)],
        true,
    );
    let net = Network::new(tri.complex(), &[a.clone(), b.clone()], &[g(4, 4)], None).unwrap();
    assert_eq!(net.n_paths(), 2);
    assert_eq!(net.paths(), vec![a.clone(), b.clone()]);
    assert_eq!(net.n_segments(), 3 + 6);
    // Ids in path order.
    assert_eq!(
        net.path_segment_ids(1),
        (3..9).map(SegmentId).collect::<Vec<_>>()
    );
    let info = net.segment(SegmentId(3)).unwrap();
    assert_eq!(info.prev, Some(SegmentId(8)), "a loop links last to first");
    assert!(!info.is_path_end());
    assert!(net.segment(SegmentId(0)).unwrap().is_path_end());
    assert!(net.segment(SegmentId(9)).is_none());
    // Joints: two in the open path, six in the loop.
    assert_eq!(net.path_joints(0).len(), 2);
    assert_eq!(net.path_joints(1).len(), 6);
    assert_eq!(net.joints().len(), 8);
    assert_eq!(net.path_joints(1)[0].incoming, SegmentId(8));
    // Pinned: the explicit vertex and the ends of the open path.
    assert_eq!(net.pinned_vertices(), {
        let mut v = vec![g(0, 2), g(3, 2), g(4, 4)];
        v.sort();
        v
    });
    assert_eq!(net.path_ends(tri.complex(), 0), Some([g(0, 2), g(3, 2)]));
    assert_eq!(net.path_ends(tri.complex(), 1), None);
    // Every segment is the only one on its edge; outermost on both sides.
    for p in 0..2 {
        for id in net.path_segment_ids(p) {
            let s = net.segment(id).unwrap().segment;
            assert_eq!(net.layer(s.edge as usize).collect::<Vec<_>>(), vec![id]);
            let h = tri.complex().edge_halfedge(s.edge as usize);
            assert_eq!(net.outermost(tri.complex(), h), Some(id));
            assert!(net.carries_segment(s.edge as usize));
        }
    }
    // Round trip through the exported layering.
    let again = Network::new(
        tri.complex(),
        &net.paths(),
        &net.pinned_vertices(),
        Some(&net.layering()),
    )
    .unwrap();
    assert_eq!(again, net);

    let gn = GeodesicNetwork::new(tri.clone(), &[a.clone(), b.clone()], &[], None).unwrap();
    let expected = a.length(&tri) + b.length(&tri);
    assert_eq!(gn.length(), expected);
    assert_eq!(gn.recompute_length(), expected);
    assert!((gn.length() - 2.25).abs() < 1e-12, "{}", gn.length());
    gn.check_invariants().unwrap();
    assert_eq!(gn.paths(), vec![a, b]);
}

#[test]
fn default_layering_and_out_and_back() {
    // Out and back along one edge, then on: the backward segment goes to the back.
    let shape = grid(4, 4);
    let tri = shape.triangulation();
    let g = |i, j| gv(4, i, j);
    let path = vpath(&tri, &[g(1, 1), g(2, 1), g(1, 1), g(1, 2)], false);
    let net = Network::new(tri.complex(), std::slice::from_ref(&path), &[], None).unwrap();
    let e = path.segments[0].edge as usize;
    assert_eq!(path.segments[1].edge as usize, e);
    let layer: Vec<SegmentId> = net.layer(e).collect();
    let (fwd, bwd) = if path.segments[0].forward {
        (SegmentId(0), SegmentId(1))
    } else {
        (SegmentId(1), SegmentId(0))
    };
    assert_eq!(layer, vec![fwd, bwd]);
    // The same edge traversed by two paths in the same direction: the later one outermost
    // on the side of the halfedge.
    let p1 = vpath(&tri, &[g(1, 2), g(2, 2), g(3, 2)], false);
    let p2 = vpath(&tri, &[g(2, 3), g(2, 2), g(3, 2)], false);
    let net = Network::new(tri.complex(), &[p1.clone(), p2.clone()], &[], None).unwrap();
    let e = p1.segments[1].edge as usize;
    let layer: Vec<SegmentId> = net.layer(e).collect();
    if p1.segments[1].forward {
        assert_eq!(layer, vec![SegmentId(3), SegmentId(1)]);
    } else {
        assert_eq!(layer, vec![SegmentId(1), SegmentId(3)]);
    }
}

#[test]
fn crossing_inputs_are_refused() {
    let shape = grid(4, 4);
    let tri = shape.triangulation();
    let c = tri.complex();
    let g = |i, j| gv(4, i, j);
    let centre = g(2, 2);

    // Two paths that cross at the centre.
    let h = vpath(&tri, &[g(1, 2), g(2, 2), g(3, 2)], false);
    let v = vpath(&tri, &[g(2, 1), g(2, 2), g(2, 3)], false);
    let err = Network::new(c, &[h.clone(), v.clone()], &[], None).unwrap_err();
    assert!(
        matches!(err, FlipError::Crossing { vertex, paths } if vertex == centre
            && (paths == [0, 1] || paths == [1, 0])),
        "{err}"
    );
    assert!(
        err.to_string().contains(&format!("vertex {centre}")),
        "{err}"
    );
    // Touching without crossing is fine: a path that turns back at the centre.
    let touch = vpath(&tri, &[g(2, 1), g(2, 2), g(1, 2)], false);
    let other = vpath(&tri, &[g(2, 3), g(2, 2), g(3, 2)], false);
    Network::new(c, &[touch, other], &[], None).unwrap();

    // A loop that crosses itself at the centre (a figure eight whose lobes cross).
    let eight = vpath(
        &tri,
        &[
            centre,
            g(3, 2),
            g(3, 3),
            g(2, 3),
            centre,
            g(2, 1),
            g(1, 1),
            g(1, 2),
        ],
        true,
    );
    let err = Network::new(c, std::slice::from_ref(&eight), &[], None).unwrap_err();
    assert!(
        matches!(err, FlipError::Crossing { vertex, paths: [0, 0] } if vertex == centre),
        "{err}"
    );
    assert!(err.to_string().contains("crosses itself"), "{err}");
    // The same figure eight with lobes that only touch.
    let touching = vpath(
        &tri,
        &[
            centre,
            g(3, 2),
            g(3, 3),
            g(2, 3),
            centre,
            g(1, 2),
            g(1, 1),
            g(2, 1),
        ],
        true,
    );
    Network::new(c, &[touching], &[], None).unwrap();
}

#[test]
fn shared_edges_need_a_consistent_layering() {
    let shape = grid(4, 4);
    let tri = shape.triangulation();
    let c = tri.complex();
    let g = |i, j| gv(4, i, j);
    // Both paths run along edge (2,2)-(3,2); one comes from the west, one from the north.
    let p1 = vpath(&tri, &[g(1, 2), g(2, 2), g(3, 2)], false);
    let p2 = vpath(&tri, &[g(2, 3), g(2, 2), g(3, 2)], false);
    let net = Network::new(c, &[p1.clone(), p2.clone()], &[], None).unwrap();
    let mut layering = net.layering();
    let e = p1.segments[1].edge as usize;
    layering[e].reverse();
    let err = Network::new(c, &[p1.clone(), p2.clone()], &[], Some(&layering)).unwrap_err();
    assert!(
        matches!(err, FlipError::Crossing { vertex, .. } if vertex == g(2, 2)),
        "{err}"
    );
    // Swap the roles: from the south instead of the north, the reversed order is the right one.
    let p3 = vpath(&tri, &[g(2, 1), g(2, 2), g(3, 2)], false);
    Network::new(c, &[p1.clone(), p3.clone()], &[], None).unwrap_err();
    let first = SegmentRef { path: 1, index: 0 };
    layering[p2.segments[0].edge as usize].retain(|r| *r != first);
    layering[p3.segments[0].edge as usize].push(first);
    Network::new(c, &[p1, p3], &[], Some(&layering)).unwrap();
}

#[test]
fn layering_errors_name_the_segment() {
    let shape = grid(3, 3);
    let tri = shape.triangulation();
    let c = tri.complex();
    let g = |i, j| gv(3, i, j);
    let p = vpath(&tri, &[g(0, 1), g(1, 1), g(2, 1)], false);
    let paths = [p.clone()];
    let good = Network::new(c, &paths, &[], None).unwrap().layering();
    let (e0, e1) = (p.segments[0].edge as usize, p.segments[1].edge as usize);

    let short: Vec<Vec<SegmentRef>> = good[..good.len() - 1].to_vec();
    assert!(matches!(
        Network::new(c, &paths, &[], Some(&short)),
        Err(FlipError::LayeringSize { .. })
    ));
    let mut wrong = good.clone();
    wrong[e0].clear();
    wrong[e1].push(SegmentRef { path: 0, index: 0 });
    let err = Network::new(c, &paths, &[], Some(&wrong)).unwrap_err();
    assert_eq!(
        err,
        FlipError::LayeringWrongEdge {
            edge: e1,
            path: 0,
            index: 0,
            segment_edge: e0
        }
    );
    assert!(err.to_string().contains(&format!("edge {e1}")), "{err}");
    let mut dup = good.clone();
    dup[e0].push(SegmentRef { path: 0, index: 0 });
    assert!(matches!(
        Network::new(c, &paths, &[], Some(&dup)),
        Err(FlipError::LayeringDuplicate {
            path: 0,
            index: 0,
            ..
        })
    ));
    let mut missing = good.clone();
    missing[e1].clear();
    assert_eq!(
        Network::new(c, &paths, &[], Some(&missing)).unwrap_err(),
        FlipError::LayeringMissing {
            edge: e1,
            path: 0,
            index: 1
        }
    );
    let mut out = good.clone();
    out[e0].push(SegmentRef { path: 3, index: 0 });
    assert!(matches!(
        Network::new(c, &paths, &[], Some(&out)),
        Err(FlipError::IndexOutOfRange { what: "path", .. })
    ));
}

#[test]
fn bad_paths_are_refused_with_messages() {
    let shape = grid(3, 3);
    let tri = shape.triangulation();
    let c = tri.complex();
    let g = |i, j| gv(3, i, j);
    let p = vpath(&tri, &[g(0, 1), g(1, 1), g(2, 1)], false);

    let err = Network::new(c, &[p.clone(), EdgePath::open(vec![])], &[], None).unwrap_err();
    assert_eq!(err, FlipError::EmptyPath { path: 1 });

    let bad_edge = EdgePath::open(vec![Segment::new(c.n_edges() + 4, true)]);
    let err = Network::new(c, &[bad_edge], &[], None).unwrap_err();
    assert!(matches!(
        err,
        FlipError::IndexOutOfRange { what: "edge", .. }
    ));

    let mut gap = p.clone();
    gap.segments.swap(0, 1);
    let err = Network::new(c, &[gap], &[], None).unwrap_err();
    assert!(
        matches!(
            err,
            FlipError::PathNotConnected {
                path: 0,
                index: 0,
                ..
            }
        ),
        "{err}"
    );
    assert!(
        err.to_string().contains("segment 0 ends at vertex"),
        "{err}"
    );

    let mut not_closed = p.clone();
    not_closed.closed = true;
    let err = Network::new(c, &[not_closed], &[], None).unwrap_err();
    assert_eq!(
        err,
        FlipError::LoopNotClosed {
            path: 0,
            first: g(0, 1),
            last: g(2, 1)
        }
    );

    let err = Network::new(c, std::slice::from_ref(&p), &[c.n_vertices()], None).unwrap_err();
    assert!(matches!(
        err,
        FlipError::IndexOutOfRange { what: "vertex", .. }
    ));
}

#[test]
fn pinning() {
    let shape = grid(3, 3);
    let tri = shape.triangulation();
    let c = tri.complex();
    let g = |i, j| gv(3, i, j);
    let p = vpath(&tri, &[g(0, 1), g(1, 1), g(2, 1)], false);
    let mut net = Network::new(c, &[p], &[g(1, 1)], None).unwrap();
    assert!(net.is_pinned(g(1, 1)) && net.is_pinned(g(0, 1)) && net.is_pinned(g(2, 1)));
    let before = net.clone();
    let err = net.unpin(c, g(2, 1)).unwrap_err();
    assert_eq!(
        err,
        FlipError::UnpinnedEnd {
            path: 0,
            vertex: g(2, 1)
        }
    );
    assert_eq!(net, before, "a refused unpin changes nothing");
    net.unpin(c, g(1, 1)).unwrap();
    assert!(!net.is_pinned(g(1, 1)));
    net.pin(g(3, 3));
    assert!(net.is_pinned(g(3, 3)));
    net.pin(c.n_vertices() + 2); // a vertex inserted later
    assert!(net.is_pinned(c.n_vertices() + 2));
    net.check(c).unwrap();
}

/// Fig. 8 top: a path winding twice around a vertex, as an open spiral, is non-crossing with
/// the default layering (the second round inside the first); with the rounds swapped it
/// crosses. A closed path that runs around twice crosses itself under every layering: its
/// rounds must change places somewhere.
#[test]
fn multiply_covered_paths() {
    let shape = grid(4, 4);
    let tri = shape.triangulation();
    let c = tri.complex();
    let g = |i, j| gv(4, i, j);
    let ring = [
        g(1, 1),
        g(2, 1),
        g(3, 1),
        g(3, 2),
        g(3, 3),
        g(2, 3),
        g(1, 3),
        g(1, 2),
    ];
    let mut vertices = vec![g(1, 0)];
    vertices.extend_from_slice(&ring);
    vertices.extend_from_slice(&ring);
    vertices.push(g(2, 2));
    let spiral = vpath(&tri, &vertices, false);
    let net = Network::new(c, std::slice::from_ref(&spiral), &[], None).unwrap();
    for i in 0..8 {
        assert_eq!(
            net.layer(spiral.segments[1 + i].edge as usize).len(),
            2 - usize::from(i == 7)
        );
    }
    let mut swapped = net.layering();
    for l in swapped.iter_mut() {
        l.reverse();
    }
    assert!(matches!(
        Network::new(c, std::slice::from_ref(&spiral), &[], Some(&swapped)),
        Err(FlipError::Crossing { .. })
    ));

    // The closed double loop: all 2^4 layerings of a square ring cross.
    let square = [g(1, 1), g(2, 1), g(2, 2), g(1, 2)];
    let twice: Vec<usize> = square.iter().chain(square.iter()).copied().collect();
    let double = vpath(&tri, &twice, true);
    let paths = [double.clone()];
    let base = Network::new(c, &[vpath(&tri, &square, true)], &[], None)
        .unwrap()
        .layering();
    for mask in 0..16u32 {
        let mut layering = vec![Vec::new(); base.len()];
        for k in 0..4 {
            let e = double.segments[k].edge as usize;
            let (x, y) = (
                SegmentRef { path: 0, index: k },
                SegmentRef {
                    path: 0,
                    index: k + 4,
                },
            );
            layering[e] = if mask & (1 << k) == 0 {
                vec![x, y]
            } else {
                vec![y, x]
            };
        }
        let err = Network::new(c, &paths, &[], Some(&layering)).unwrap_err();
        assert!(
            matches!(err, FlipError::Crossing { paths: [0, 0], .. }),
            "{err}"
        );
    }
}

#[test]
fn construction_is_deterministic_on_meshes() {
    let mut shapes = generated_shapes();
    shapes.extend(data_shapes());
    for shape in &shapes {
        let tri = shape.delaunay();
        let n = tri.complex().n_vertices();
        let paths: Vec<EdgePath> = random_pairs(n, 6, 11)
            .into_iter()
            .filter_map(|(a, b)| dijkstra_path(&tri, a, b).ok())
            .collect();
        // Dijkstra paths are simple, but two of them may cross each other: build each alone.
        for path in &paths {
            let one = GeodesicNetwork::new(tri.clone(), std::slice::from_ref(path), &[], None)
                .unwrap_or_else(|e| panic!("{}: {e}", shape.name));
            let two =
                GeodesicNetwork::new(tri.clone(), std::slice::from_ref(path), &[], None).unwrap();
            assert_eq!(one.network(), two.network());
            assert_eq!(one.length().to_bits(), two.length().to_bits());
            one.check_invariants().unwrap();
        }
    }
}
