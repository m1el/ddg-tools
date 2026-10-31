//! Wedge angles against hand computations, the joint classification and the Δ-complex cases
//! (SPEC §3).

mod common;

use std::f64::consts::PI;

use common::*;
use ddg_flip_geodesics::{
    EdgePath, GeodesicNetwork, Joint, JointGeometry, Network, Segment, Side, Turn, VertexStar,
    dijkstra_path,
};
use ddg_intrinsic::IntrinsicTriangulation;

const TOL: f64 = 1e-5;

/// The geometry of the `k`-th joint of path 0 of a network made of `paths`.
fn joint_of(tri: &IntrinsicTriangulation, paths: &[EdgePath], k: usize) -> JointGeometry {
    let net = Network::new(tri.complex(), paths, &[], None).unwrap();
    JointGeometry::new(tri, &net, net.path_joints(0)[k]).unwrap()
}

fn close(a: f64, b: f64, tol: f64) -> bool {
    (a - b).abs() <= tol * (1.0 + b.abs())
}

/// The centroid of the face of corner `h`.
fn centroid(shape: &Shape, tri: &IntrinsicTriangulation, h: u32) -> [f64; 3] {
    let f = tri.complex().faces()[h as usize / 3];
    let p = |k: usize| shape.positions[f[k] as usize];
    [0, 1, 2].map(|c| (p(0)[c] + p(1)[c] + p(2)[c]) / 3.0)
}

#[test]
fn flat_grid_by_hand() {
    let shape = grid(4, 4);
    let tri = shape.triangulation();
    let g = |i, j| gv(4, i, j);
    let y0 = 0.5;

    // Straight through an interior vertex: π on both sides; the left wedge lies north.
    let j = joint_of(&tri, &[vpath(&tri, &[g(1, 2), g(2, 2), g(3, 2)], false)], 0);
    assert!(close(j.left.angle, PI, 1e-14) && close(j.right.angle, PI, 1e-14));
    assert_eq!(j.turn(TOL), Turn::Straight);
    assert!(j.sides_to_shorten(TOL).is_empty());
    assert!(
        j.left
            .corners
            .iter()
            .all(|&h| centroid(&shape, &tri, h)[1] > y0)
    );
    assert!(
        j.right
            .corners
            .iter()
            .all(|&h| centroid(&shape, &tri, h)[1] < y0)
    );
    assert!(!j.left.boundary && !j.right.boundary);
    assert_eq!(j.left.degree() + j.right.degree(), 6);

    // East then north: a left turn of π/2.
    let j = joint_of(&tri, &[vpath(&tri, &[g(1, 2), g(2, 2), g(2, 3)], false)], 0);
    assert!(close(j.left.angle, PI / 2.0, 1e-14), "{}", j.left.angle);
    assert!(close(j.right.angle, 1.5 * PI, 1e-14));
    assert_eq!(j.min_side(), Side::Left);
    assert_eq!(j.turn(TOL), Turn::Left);
    assert_eq!(j.sides_to_shorten(TOL), vec![(Side::Left, j.left.angle)]);

    // A boundary corner: west then north around (0, 0); the outside is infinite.
    let j = joint_of(&tri, &[vpath(&tri, &[g(1, 0), g(0, 0), g(0, 1)], false)], 0);
    assert!(j.left.boundary && j.left.angle == f64::INFINITY);
    assert!(close(j.right.angle, PI / 2.0, 1e-14));
    assert_eq!(j.turn(TOL), Turn::Right);
    // Along the boundary: π inside.
    let j = joint_of(&tri, &[vpath(&tri, &[g(0, 0), g(1, 0), g(2, 0)], false)], 0);
    assert!(close(j.left.angle, PI, 1e-14) && j.right.angle == f64::INFINITY);
    assert_eq!(j.turn(TOL), Turn::Straight);
}

/// A loop counter-clockwise along the boundary of a flat polygon: its inner (left) wedge
/// angles add up to `(n − 2) π`, the outer ones are infinite.
#[test]
fn flat_polygon_angles_sum() {
    let (nx, ny) = (3, 2);
    let shape = grid(nx, ny);
    let tri = shape.triangulation();
    let g = |i, j| gv(ny, i, j);
    let mut ring = Vec::new();
    ring.extend((0..nx).map(|i| g(i, 0)));
    ring.extend((0..ny).map(|j| g(nx, j)));
    ring.extend((1..=nx).rev().map(|i| g(i, ny)));
    ring.extend((1..=ny).rev().map(|j| g(0, j)));
    let n = ring.len();
    let path = vpath(&tri, &ring, true);
    let net = Network::new(tri.complex(), &[path], &[], None).unwrap();
    let mut sum = 0.0;
    for joint in net.path_joints(0) {
        let j = JointGeometry::new(&tri, &net, joint).unwrap();
        assert_eq!(j.right.angle, f64::INFINITY);
        sum += j.left.angle;
    }
    assert!(close(sum, (n - 2) as f64 * PI, 1e-14), "{sum}");
}

#[test]
fn cone_apex_of_a_pyramid() {
    for h in [0.5, 5.0] {
        let shape = pyramid(h);
        let tri = shape.triangulation();
        let p = &shape.positions;
        let theta = angle_at(p[0], p[4], p[1]);
        let big_theta = tri.vertex_angle_sum(4);
        assert!(close(big_theta, 4.0 * theta, 1e-14));
        // Opposite base corners: 2θ on each side.
        let j = joint_of(&tri, &[vpath(&tri, &[0, 4, 2], false)], 0);
        assert!(close(j.left.angle, 2.0 * theta, 1e-14));
        assert!(close(j.right.angle, 2.0 * theta, 1e-14));
        // Neighbouring base corners: θ on the side of face (0, 1, 4), 3θ on the other.
        let j = joint_of(&tri, &[vpath(&tri, &[0, 4, 1], false)], 0);
        assert!(
            close(j.right.angle, theta, 1e-14),
            "{} {theta}",
            j.right.angle
        );
        assert!(close(j.left.angle, 3.0 * theta, 1e-14));
        assert!(close(j.left.angle + j.right.angle, big_theta, 1e-14));
        assert_eq!(j.turn(TOL), Turn::Right);
        // A cone point never carries a geodesic: α^min ≤ Θ/2 < π.
        assert!(j.min_angle() < PI);
        let sides = j.sides_to_shorten(TOL);
        if big_theta < PI {
            // Both wedges below π: both are queued, the smaller first.
            assert_eq!(sides.len(), 2, "h {h}");
            assert_eq!(sides[0].0, Side::Right);
        } else {
            assert_eq!(
                sides.len(),
                1 + usize::from(3.0 * theta < PI - TOL),
                "h {h}"
            );
        }
    }
}

#[test]
fn saddle_with_three_pi() {
    let shape = saddle();
    let tri = shape.triangulation();
    assert!(close(tri.vertex_angle_sum(0), 3.0 * PI, 1e-14));
    let j = joint_of(&tri, &[vpath(&tri, &[1, 0, 4], false)], 0);
    assert!(close(j.left.angle, 1.5 * PI, 1e-14) && close(j.right.angle, 1.5 * PI, 1e-14));
    assert_eq!(j.turn(TOL), Turn::Straight);
    let j = joint_of(&tri, &[vpath(&tri, &[1, 0, 3], false)], 0);
    assert!(close(j.min_angle(), PI, 1e-14) && close(j.max_angle(), 2.0 * PI, 1e-14));
    assert_eq!(j.turn(TOL), Turn::Straight);
    let j = joint_of(&tri, &[vpath(&tri, &[1, 0, 2], false)], 0);
    assert!(close(j.right.angle, PI / 2.0, 1e-14) && close(j.left.angle, 2.5 * PI, 1e-14));
    assert_eq!((j.right.degree(), j.left.degree()), (1, 5));
    assert_eq!(j.turn(TOL), Turn::Right);
    // Corners in the order of the outer arc a = n₀ = 1, n₁, …, n₅ = c = 2: corner i spans
    // n_i, n_{i+1}; for the left wedge n_{i+1} is the head of its halfedge.
    let arc: Vec<usize> = j
        .corners_from_incoming(Side::Left)
        .iter()
        .map(|&h| tri.complex().to_vertex(h as usize))
        .collect();
    assert_eq!(arc, vec![6, 5, 4, 3, 2]);
}

/// On every shape, after random flips: for every joint of Dijkstra paths, the corners of the
/// two wedges are exactly the corners at the vertex and their angles add up to `Θ_b`; at a
/// boundary vertex exactly one side is infinite.
#[test]
fn wedges_partition_the_star() {
    let mut shapes = generated_shapes();
    shapes.extend(data_shapes());
    for shape in &shapes {
        for flips in [0usize, 1] {
            let mut tri = shape.triangulation();
            if flips > 0 {
                let n_flips = 3 * tri.complex().n_edges();
                random_flips(&mut tri, n_flips, 5);
            }
            let n = tri.complex().n_vertices();
            let mut n_joints = 0;
            for (a, b) in random_pairs(n, 8, 3) {
                let Ok(path) = dijkstra_path(&tri, a, b) else {
                    continue;
                };
                let net = Network::new(tri.complex(), &[path], &[], None).unwrap();
                for joint in net.path_joints(0) {
                    n_joints += 1;
                    let j = JointGeometry::new(&tri, &net, joint).unwrap();
                    let v = j.vertex();
                    let mut corners: Vec<u32> = j
                        .left
                        .corners
                        .iter()
                        .chain(&j.right.corners)
                        .copied()
                        .collect();
                    corners.sort_unstable();
                    let mut all: Vec<u32> = tri
                        .complex()
                        .outgoing_halfedges(v)
                        .map(|h| h as u32)
                        .collect();
                    all.sort_unstable();
                    assert_eq!(corners, all, "{}: vertex {v}", shape.name);
                    let theta = tri.vertex_angle_sum(v);
                    if tri.complex().is_boundary_vertex(v) {
                        assert!(j.left.boundary != j.right.boundary, "{}", shape.name);
                        assert!(j.min_angle() <= theta * (1.0 + 1e-12));
                    } else {
                        let sum = j.left.angle + j.right.angle;
                        // The corners summed in another order: equal to rounding. The stored
                        // angle sum of the input: equal to the accuracy of angles from lengths
                        // (looser on the needles of camelbad).
                        let star_sum = j.star.angle_sum(&tri);
                        assert!(close(sum, star_sum, 1e-14), "{}", shape.name);
                        assert!(
                            close(sum, theta, 1e-10),
                            "{}: vertex {v}: {sum} vs {theta}",
                            shape.name
                        );
                    }
                    // Bitwise reproducible.
                    assert_eq!(j, JointGeometry::new(&tri, &net, joint).unwrap());
                }
            }
            assert!(n_joints > 0 || n < 10, "{}", shape.name);
        }
    }
}

/// The degree-one vertex and the self-edge that Delaunay flipping makes on a tilted needle
/// (as `ddg-intrinsic`'s tests): out and back to the degree-one vertex has a wedge of degree
/// zero (angle 0) on the side between the two segments in the layering and the whole vertex
/// on the other; a loop along the self-edge has its two wedges add up to `Θ`, and the
/// out-and-back path inside one of them blocks it with its ends.
#[test]
fn degree_one_vertex_and_self_edge() {
    let shape = needle_tetrahedron(8.0, 4.0);
    let tri = shape.delaunay();
    let irregular = tri.irregularities();
    assert_eq!(irregular.degree_one_vertices.len(), 1);
    assert_eq!(irregular.self_edges.len(), 1);
    let d = irregular.degree_one_vertices[0] as usize;
    let c = tri.complex();
    assert_eq!(c.vertex_degree(d), 1);
    let x = c.to_vertex(c.outgoing_halfedge(d));
    let theta_d = tri.vertex_angle_sum(d);

    let out_and_back = vpath(&tri, &[x, d, x], false);
    let net = Network::new(c, std::slice::from_ref(&out_and_back), &[], None).unwrap();
    let joint = net.path_joints(0)[0];
    let j = JointGeometry::new(&tri, &net, joint).unwrap();
    let (zero, full) = if j.left.degree() == 0 {
        (Side::Left, Side::Right)
    } else {
        (Side::Right, Side::Left)
    };
    assert_eq!(j.angle(zero), 0.0);
    assert_eq!(j.wedge(full).degree(), 1);
    assert!(close(j.angle(full), theta_d, 1e-12));
    assert_eq!(
        j.turn(TOL),
        if zero == Side::Left {
            Turn::Left
        } else {
            Turn::Right
        }
    );
    assert!(j.is_flexible(zero) && j.is_flexible(full));
    // The other layering puts the empty wedge on the other side.
    let mut layering = net.layering();
    layering[out_and_back.segments[0].edge as usize].reverse();
    let net2 = Network::new(c, std::slice::from_ref(&out_and_back), &[], Some(&layering)).unwrap();
    let j2 = JointGeometry::new(&tri, &net2, net2.path_joints(0)[0]).unwrap();
    assert_eq!(j2.angle(full), 0.0);
    assert!(close(j2.angle(zero), theta_d, 1e-12));

    // The loop along the self-edge, at its vertex.
    let s = irregular.self_edges[0] as usize;
    let y = c.edge_vertices(s)[0];
    let lp = EdgePath::closed(vec![Segment::new(s, true)]);
    let net = Network::new(c, std::slice::from_ref(&lp), &[], None).unwrap();
    let j = JointGeometry::new(&tri, &net, net.path_joints(0)[0]).unwrap();
    assert_eq!(j.vertex(), y);
    assert_eq!(
        j.left.degree() + j.right.degree(),
        c.outgoing_halfedges(y).count()
    );
    if !c.is_boundary_vertex(y) {
        assert!(close(
            j.left.angle + j.right.angle,
            tri.vertex_angle_sum(y),
            1e-12
        ));
    }
    // The self-edge goes around the degree-one vertex: y = x. With the out-and-back path
    // ending at x inside one wedge of the loop.
    assert_eq!(y, x);
    {
        let both = [lp.clone(), out_and_back.clone()];
        let net = Network::new(c, &both, &[], None).unwrap();
        let j = JointGeometry::new(&tri, &net, net.path_joints(0)[0]).unwrap();
        let blocked: Vec<Side> = [Side::Left, Side::Right]
            .into_iter()
            .filter(|&s| j.is_blocked_by_end(&net, s))
            .collect();
        assert_eq!(blocked.len(), 1, "the path's ends lie in exactly one wedge");
        assert!(!j.is_flexible(blocked[0]));
        assert!(j.is_flexible(blocked[0].opposite()));
    }
}

/// A self-edge at a vertex of the Delaunay triangulation of a thin torus lies inside a wedge
/// of a joint through that vertex with both its ends (twice incident, App. A). With a loop on
/// it, that wedge is not flexible, and not blocked by an end either.
#[test]
fn twice_incident_self_edge() {
    let shape = long_torus(8, 3, 4.0, 0.2, 0.37);
    let tri = shape.delaunay();
    let c = tri.complex();
    let s = tri.irregularities().self_edges[0] as usize;
    let v = c.edge_vertices(s)[0];
    let lp = EdgePath::closed(vec![Segment::new(s, true)]);
    let net = Network::new(c, std::slice::from_ref(&lp), &[], None).unwrap();
    let star = VertexStar::new(c, &net, v);
    let n = star.ends.len();
    let ends_of_s: Vec<usize> = (0..n)
        .filter(|&i| star.ends[i].edge as usize == s)
        .collect();
    assert_eq!(ends_of_s.len(), 2);
    // Edge ends just outside the two ends of s, skipping other self-edges.
    let is_self = |i: usize| {
        let [a, b] = c.edge_vertices(star.ends[i].edge as usize);
        a == b
    };
    let mut p = (ends_of_s[0] + n - 1) % n;
    while is_self(p) {
        p = (p + n - 1) % n;
    }
    let mut q = (ends_of_s[1] + 1) % n;
    while is_self(q) {
        q = (q + 1) % n;
    }
    let (ep, eq) = (star.ends[p], star.ends[q]);
    let incoming = Segment::new(ep.edge as usize, !ep.at_tail);
    let outgoing = Segment::new(eq.edge as usize, eq.at_tail);
    assert_eq!(incoming.head(c), v);
    assert_eq!(outgoing.tail(c), v);
    let through = EdgePath::open(vec![incoming, outgoing]);
    let paths = [through, lp];
    let net = Network::new(c, &paths, &[], None).unwrap();
    let joint: Joint = net.path_joints(0)[0];
    let j = JointGeometry::new(&tri, &net, joint).unwrap();
    // The right wedge (counter-clockwise from ab to bc) holds both ends of s.
    let inside: Vec<_> = j.star.inside(&j.right).map(|e| e.segment).collect();
    assert_eq!(inside.len(), 2, "both ends of the loop's segment");
    let both_ends_counted = j
        .right
        .corners
        .iter()
        .filter(|&&h| c.edge(h as usize) == s)
        .count();
    assert_eq!(
        both_ends_counted, 2,
        "s bounds two corners of the right wedge"
    );
    assert!(close(
        j.left.angle + j.right.angle,
        tri.vertex_angle_sum(v),
        1e-12
    ));
    assert!(!j.is_flexible(Side::Right));
    assert!(!j.is_blocked_by_end(&net, Side::Right));
    assert!(j.is_flexible(Side::Left));
    let gn = GeodesicNetwork::new(tri.clone(), &paths, &[], None).unwrap();
    assert_eq!(gn.joint_geometry(joint).unwrap(), j);
}

#[test]
fn bad_joints_are_refused() {
    let shape = grid(3, 3);
    let tri = shape.triangulation();
    let g = |i, j| gv(3, i, j);
    let p = vpath(&tri, &[g(0, 1), g(1, 1), g(2, 1)], false);
    let net = Network::new(tri.complex(), &[p], &[], None).unwrap();
    let mut joint = net.path_joints(0)[0];
    std::mem::swap(&mut joint.incoming, &mut joint.outgoing);
    assert!(JointGeometry::new(&tri, &net, joint).is_err());
}
