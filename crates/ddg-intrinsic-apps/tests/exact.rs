//! The exact polyhedral distance (the ground truth of Sec. 5.3) against cases with known
//! answers, against other triangulations of the same surface, and against bounds.

mod common;

use std::collections::BinaryHeap;

use common::{Shape, cube, dist, jittered_grid, load_data, sheared_parallelogram};
use ddg_core::generators::{annulus, bumpy_grid_disk, icosphere, torus};
use ddg_intrinsic::{DelaunayOptions, ExportOptions, RefineOptions};
use ddg_intrinsic_apps::exact::{ExactError, ExactOptions, exact_distance};

fn exact(shape: &Shape, sources: &[usize]) -> Vec<f64> {
    exact_distance(&shape.mesh, &shape.geometry, sources, &ExactOptions::default())
        .unwrap_or_else(|e| panic!("{}: {e}", shape.name))
        .distances
}

#[test]
fn flat_convex_domains_give_euclidean_distances() {
    for shape in [
        sheared_parallelogram(10, 3.0),
        jittered_grid(12, 0.35, 7),
        Shape::generated("grid", ddg_core::generators::grid_disk(9, 5).unwrap()),
    ] {
        for source in [0, shape.mesh.n_vertices() / 2 + 3] {
            let d = exact(&shape, &[source]);
            for (v, &x) in d.iter().enumerate() {
                let e = dist(shape.positions[source], shape.positions[v]);
                assert!((x - e).abs() < 1e-10, "{}: vertex {v}: {x} != {e}", shape.name);
            }
        }
    }
}

/// Shortest paths in the plane region between two convex polygons (the outer one convex, the
/// inner one a hole) through the vertices of the hole: Dijkstra on the visibility graph.
fn annulus_reference(points: &[[f64; 2]], hole: &[[f64; 2]], source: usize) -> Vec<f64> {
    // Does the segment pq pass through the open interior of the convex polygon `hole`
    // (counter-clockwise)? Cyrus-Beck clipping against the strict half-planes.
    let blocked = |p: [f64; 2], q: [f64; 2]| {
        let (mut t0, mut t1) = (0.0f64, 1.0f64);
        let d = [q[0] - p[0], q[1] - p[1]];
        for k in 0..hole.len() {
            let (a, b) = (hole[k], hole[(k + 1) % hole.len()]);
            // Inside: left of a -> b, i.e. cross(b − a, x − a) > 0.
            let e = [b[0] - a[0], b[1] - a[1]];
            let f0 = e[0] * (p[1] - a[1]) - e[1] * (p[0] - a[0]);
            let fd = e[0] * d[1] - e[1] * d[0];
            if fd.abs() < 1e-15 {
                if f0 <= 1e-12 {
                    return false;
                }
            } else {
                let t = -f0 / fd;
                if fd > 0.0 {
                    t0 = t0.max(t);
                } else {
                    t1 = t1.min(t);
                }
            }
        }
        t1 - t0 > 1e-9
    };
    let n = points.len();
    let mut d = vec![f64::INFINITY; n];
    let mut done = vec![false; n];
    d[source] = 0.0;
    let mut heap = BinaryHeap::new();
    heap.push((std::cmp::Reverse(0u64), source));
    while let Some((_, v)) = heap.pop() {
        if done[v] {
            continue;
        }
        done[v] = true;
        for w in 0..n {
            if done[w] || blocked(points[v], points[w]) {
                continue;
            }
            let x = d[v] + (points[v][0] - points[w][0]).hypot(points[v][1] - points[w][1]);
            if x < d[w] {
                d[w] = x;
                heap.push((std::cmp::Reverse(x.to_bits()), w));
            }
        }
    }
    d
}

#[test]
fn paths_bend_around_a_hole_at_boundary_vertices() {
    let (n_around, n_radial) = (24, 4);
    let shape = Shape::generated("annulus", annulus(n_around, n_radial, 0.4, 1.0).unwrap());
    let points: Vec<[f64; 2]> = shape.positions.iter().map(|p| [p[0], p[1]]).collect();
    // The inner ring: j = n_radial, counter-clockwise in i.
    let hole: Vec<[f64; 2]> = (0..n_around)
        .map(|i| points[i * (n_radial + 1) + n_radial])
        .collect();
    for source in [0, 2, n_radial] {
        let d = exact(&shape, &[source]);
        let reference = annulus_reference(&points, &hole, source);
        let mut bent = 0;
        for v in 0..points.len() {
            assert!(
                (d[v] - reference[v]).abs() < 1e-10,
                "source {source}, vertex {v}: {} != {}",
                d[v],
                reference[v]
            );
            if reference[v] > dist(shape.positions[source], shape.positions[v]) + 1e-6 {
                bent += 1;
            }
        }
        assert!(bent > 10, "source {source}: only {bent} paths bend around the hole");
    }
}

#[test]
fn cube_distances_are_those_of_the_unfolding() {
    let n = 4;
    let shape = cube(n);
    let corner = |p: [f64; 3]| {
        shape
            .positions
            .iter()
            .position(|&q| q == p)
            .expect("a corner of the cube")
    };
    let source = corner([0.0, 0.0, 0.0]);
    let d = exact(&shape, &[source]);
    // Corners: along an edge 1, across a face √2, the opposite corner √5 (unfold two faces).
    for (p, expected) in [
        ([1.0, 0.0, 0.0], 1.0),
        ([1.0, 1.0, 0.0], 2f64.sqrt()),
        ([1.0, 1.0, 1.0], 5f64.sqrt()),
        ([0.0, 1.0, 1.0], 2f64.sqrt()),
    ] {
        let v = corner(p);
        assert!((d[v] - expected).abs() < 1e-12, "{p:?}: {} != {expected}", d[v]);
    }
    // On the three faces at the source the straight segment is shortest.
    for (v, p) in shape.positions.iter().enumerate() {
        if p.contains(&0.0) {
            let e = dist(*p, shape.positions[source]);
            assert!((d[v] - e).abs() < 1e-12, "vertex {v} {p:?}: {} != {e}", d[v]);
        }
    }
    // The far faces: the shortest of the unfoldings over one or two faces.
    for (v, p) in shape.positions.iter().enumerate() {
        if (p[0] - 1.0).abs() < 1e-15 {
            // On x = 1: unfold across the edge x = 1 of the face z = 0 or y = 0, or across two
            // faces (z = 0 then y = 1 ... ) — the first two are enough for points with y, z ≤ 1:
            // over z = 0: (1 + z, y); over y = 0: (1 + y, z).
            let a = (1.0 + p[2]).hypot(p[1]);
            let b = (1.0 + p[1]).hypot(p[2]);
            let expected = a.min(b);
            assert!(d[v] <= expected + 1e-12, "vertex {v} {p:?}: {} > {expected}", d[v]);
            assert!((d[v] - expected).abs() < 1e-12, "vertex {v} {p:?}: {} != {expected}", d[v]);
        }
    }
}

#[test]
fn the_distance_is_that_of_the_surface_not_of_the_triangulation() {
    // The same surface triangulated differently (intrinsic Delaunay, refined, exported as
    // meshes) has the same distances between input vertices.
    let mut shapes = vec![
        Shape::generated("icosphere 2", icosphere(2).unwrap()),
        Shape::generated("torus 16x6", torus(16, 6, 1.0, 0.35).unwrap()),
        Shape::generated("bumpy grid", bumpy_grid_disk(10, 10, 0.3).unwrap()),
        sheared_parallelogram(8, 4.0),
    ];
    shapes.extend(load_data("optimal-cones/optimal-cones-data/data/camelbad.obj"));
    for shape in &shapes {
        let n = shape.mesh.n_vertices();
        let sources = [0, n / 3];
        let reference = exact(shape, &sources[..1]);
        let multi = exact(shape, &sources);
        let other = exact(shape, &sources[1..]);
        for v in 0..n {
            assert!((multi[v] - reference[v].min(other[v])).abs() < 1e-12 * (1.0 + multi[v]));
        }
        for refine in [false, true] {
            let mut tri = shape.triangulation();
            tri.flip_to_delaunay(&DelaunayOptions::default()).unwrap();
            if refine {
                tri.refine(&RefineOptions::default()).unwrap();
            }
            let export = tri.export_simplicial(&ExportOptions::default()).unwrap();
            let d = exact_distance(&export.mesh, &export.geometry, &[0], &ExactOptions::default())
                .unwrap()
                .distances;
            let scale = reference.iter().fold(0.0f64, |m, &x| m.max(x));
            let worst = (0..n)
                .map(|v| (d[v] - reference[v]).abs())
                .fold(0.0f64, f64::max);
            assert!(
                worst < 1e-9 * scale,
                "{} (refine {refine}): largest difference {worst:e} of {scale}",
                shape.name
            );
        }
    }
}

#[test]
fn bounded_by_the_straight_line_and_by_the_edge_graph() {
    for shape in [
        Shape::generated("torus", torus(20, 7, 1.0, 0.4).unwrap()),
        Shape::generated("bumpy grid", bumpy_grid_disk(12, 12, 0.25).unwrap()),
    ] {
        let d = exact(&shape, &[5]);
        // Dijkstra over edges.
        let n = shape.mesh.n_vertices();
        let mut g = vec![f64::INFINITY; n];
        g[5] = 0.0;
        let mut heap = BinaryHeap::new();
        heap.push((std::cmp::Reverse(0u64), 5usize));
        while let Some((_, v)) = heap.pop() {
            for h in shape.mesh.outgoing_halfedges(v).collect::<Vec<_>>() {
                for (a, b) in [
                    (shape.mesh.to_vertex(h), h),
                    (shape.mesh.from_vertex(shape.mesh.prev(h)), shape.mesh.prev(h)),
                ] {
                    let x = g[v] + shape.geometry.edge_lengths()[shape.mesh.edge(b)];
                    if x < g[a] {
                        g[a] = x;
                        heap.push((std::cmp::Reverse(x.to_bits()), a));
                    }
                }
            }
        }
        let mut strictly_shorter = 0;
        for v in 0..n {
            let straight = dist(shape.positions[5], shape.positions[v]);
            assert!(d[v] >= straight - 1e-12, "{}: vertex {v}", shape.name);
            assert!(d[v] <= g[v] + 1e-12, "{}: vertex {v}", shape.name);
            if d[v] < g[v] - 1e-6 {
                strictly_shorter += 1;
            }
        }
        assert!(strictly_shorter > n / 2, "{}: {strictly_shorter} of {n}", shape.name);
    }
}

#[test]
fn bad_input_is_an_error() {
    let shape = Shape::generated("icosphere 1", icosphere(1).unwrap());
    let options = ExactOptions::default();
    assert_eq!(
        exact_distance(&shape.mesh, &shape.geometry, &[], &options),
        Err(ExactError::NoSource)
    );
    assert!(matches!(
        exact_distance(&shape.mesh, &shape.geometry, &[1000], &options),
        Err(ExactError::SourceOutOfRange { source: 1000, .. })
    ));
    let limited = ExactOptions {
        max_windows: Some(3),
        ..options
    };
    assert_eq!(
        exact_distance(&shape.mesh, &shape.geometry, &[0], &limited),
        Err(ExactError::WindowLimit { limit: 3 })
    );
}
