//! Steiner tree approximation over intrinsic triangulations (paper Sec. 5.1).

mod common;

use std::f64::consts::{PI, TAU};

use common::{Shape, cube, dist, jittered_grid, load_data, position, sheared_torus};
use ddg_core::generators::{icosphere, torus};
use ddg_intrinsic::{IntrinsicTriangulation, SurfacePoint};
use ddg_intrinsic_apps::steiner::{
    SteinerError, SteinerOptions, SteinerTree, TriangleCandidates, cut_to_disk, fermat_point,
    minimum_spanning_tree, steiner_tree,
};

/// The tree joins every vertex (with the vertices of every boundary loop counted as joined)
/// and has no cycle.
fn assert_spanning_tree(tri: &IntrinsicTriangulation, edges: &[u32], what: &str) {
    let n = tri.n_vertices();
    let mut parent: Vec<usize> = (0..n).collect();
    fn find(p: &mut [usize], mut x: usize) -> usize {
        while p[x] != x {
            p[x] = p[p[x]];
            x = p[x];
        }
        x
    }
    let mut n_sets = n;
    for e in 0..tri.n_edges() {
        if tri.is_boundary_edge(e) {
            let [a, b] = tri.edge_vertices(e);
            let (ra, rb) = (find(&mut parent, a), find(&mut parent, b));
            if ra != rb {
                parent[ra] = rb;
                n_sets -= 1;
            }
        }
    }
    for &e in edges {
        let [a, b] = tri.edge_vertices(e as usize);
        let (ra, rb) = (find(&mut parent, a), find(&mut parent, b));
        assert_ne!(ra, rb, "{what}: edge {e} closes a cycle");
        parent[ra] = rb;
        n_sets -= 1;
    }
    assert_eq!(n_sets, 1, "{what}: the tree leaves {n_sets} pieces");
}

fn check(shape: &Shape, tree: &SteinerTree, what: &str) {
    let tri = &tree.triangulation;
    assert_spanning_tree(tri, &tree.edges, what);
    let length: f64 = tree.edges.iter().map(|&e| tri.edge_length(e as usize)).sum();
    assert!((length - tree.length).abs() < 1e-12 * length);
    // The cut leaves a disk that is flat inside and unfolds.
    let cut = cut_to_disk(tri, &tree.edges).unwrap();
    assert!(cut.is_disk, "{what}: {cut:?}");
    assert_eq!(cut.euler_characteristic, 1);
    assert!(cut.max_interior_defect < 1e-9, "{what}: {}", cut.max_interior_defect);
    assert!((cut.boundary_turning - TAU).abs() < 1e-8, "{what}: {}", cut.boundary_turning);
    assert!(cut.max_layout_mismatch < 1e-9, "{what}: {}", cut.max_layout_mismatch);
    let boundary: f64 = (0..tri.n_edges())
        .filter(|&e| tri.is_boundary_edge(e))
        .map(|e| tri.edge_length(e))
        .sum();
    assert!((cut.boundary_length - 2.0 * tree.length - boundary).abs() < 1e-10 * cut.boundary_length);
    // The paths on the input mesh have the lengths of the edges.
    assert_eq!(tree.traces.len(), tree.edges.len());
    assert!(tree.report.max_trace_mismatch < 1e-8, "{what}");
    for t in &tree.traces {
        let [a, b] = t.vertices;
        let mut points = vec![position(shape, tri.vertex_location(a as usize))];
        for c in &t.crossings {
            points.push(position(shape, SurfacePoint::Edge { edge: c.edge, t: c.t }));
        }
        points.push(position(shape, tri.vertex_location(b as usize)));
        let along: f64 = points.windows(2).map(|w| dist(w[0], w[1])).sum();
        let l = tri.edge_length(t.edge as usize);
        assert!((along - l).abs() < 1e-8 * l, "{what}: edge {}: {along} != {l}", t.edge);
    }
}

#[test]
fn the_fermat_point_joins_the_corners_at_120_degrees() {
    for l in [[1.0, 1.0, 1.0], [3.0, 4.0, 5.0], [1.0, 1.9, 1.5], [2.0, 0.6, 1.7]] {
        let (bary, d, s) = fermat_point(l).unwrap();
        assert!(bary.iter().all(|&b| b > 0.0));
        // Corners in the plane and the point.
        let a = [0.0, 0.0];
        let b = [l[0], 0.0];
        let x = (l[0] * l[0] + l[2] * l[2] - l[1] * l[1]) / (2.0 * l[0]);
        let c = [x, (l[2] * l[2] - x * x).sqrt()];
        let p = [0, 1].map(|k| bary[0] * a[k] + bary[1] * b[k] + bary[2] * c[k]);
        let dirs = [a, b, c].map(|q| [q[0] - p[0], q[1] - p[1]]);
        for k in 0..3 {
            let (u, v) = (dirs[k], dirs[(k + 1) % 3]);
            let angle = (u[0] * v[1] - u[1] * v[0]).atan2(u[0] * v[0] + u[1] * v[1]).abs();
            assert!((angle - 2.0 * PI / 3.0).abs() < 1e-12, "{l:?}: {angle}");
            assert!((d[k] - dirs[k][0].hypot(dirs[k][1])).abs() < 1e-12);
        }
        // Length: s² = (a² + b² + c²)/2 + 2√3 A.
        let area = 0.5 * l[0] * c[1];
        let expected = ((l[0] * l[0] + l[1] * l[1] + l[2] * l[2]) / 2.0
            + 2.0 * 3f64.sqrt() * area)
            .sqrt();
        assert!((s - expected).abs() < 1e-12 * s, "{l:?}: {s} {expected}");
        // The paper's coordinates csc(θ + π/3) are trilinear: they differ.
    }
    // An angle of 120° or more: no Fermat point inside.
    assert!(fermat_point([1.0, 1.0, 1.9]).is_none());
}

#[test]
fn trees_span_all_vertices_and_cut_closed_surfaces_into_flat_disks() {
    let mut shapes = vec![
        Shape::generated("icosphere 2", icosphere(2).unwrap()),
        cube(3),
        sheared_torus(24, 8, 1.0, 0.4, 4),
    ];
    shapes.extend(load_data("models/spot/spot/spot_control_mesh.obj"));
    for shape in &shapes {
        let tri = shape.triangulation();
        let genus0 = tri.euler_characteristic() == 2;
        for (name, options) in [
            ("intrinsic steiner", SteinerOptions::default()),
            (
                "intrinsic mst",
                SteinerOptions {
                    steiner_points: false,
                    ..SteinerOptions::default()
                },
            ),
            (
                "extrinsic mst",
                SteinerOptions {
                    delaunay: None,
                    steiner_points: false,
                    ..SteinerOptions::default()
                },
            ),
            (
                "extrinsic steiner",
                SteinerOptions {
                    delaunay: None,
                    ..SteinerOptions::default()
                },
            ),
            (
                "intrinsic steiner, all triangles",
                SteinerOptions {
                    candidates: TriangleCandidates::All,
                    ..SteinerOptions::default()
                },
            ),
        ] {
            let what = format!("{} {name}", shape.name);
            let tree = steiner_tree(&tri, &options).unwrap();
            assert!(tree.length <= tree.report.mst_length * (1.0 + 1e-12) || name.contains("all"));
            if genus0 {
                check(shape, &tree, &what);
            } else {
                // A torus cut along a tree is not a disk: one boundary loop and genus one.
                assert_spanning_tree(&tree.triangulation, &tree.edges, &what);
                let cut = cut_to_disk(&tree.triangulation, &tree.edges).unwrap();
                assert!(!cut.is_disk);
                assert_eq!(cut.euler_characteristic, -1, "{what}");
            }
            if options.steiner_points {
                // The sheared torus has an angle above 120° in every input triangle: no
                // Steiner point on the input mesh.
                assert!(
                    tree.report.n_steiner_points > 0 || options.delaunay.is_none(),
                    "{what}"
                );
                assert_eq!(tree.steiner_vertices.len(), tree.report.n_steiner_points);
                for &v in &tree.steiner_vertices {
                    assert_eq!(tree.triangulation.vertex_degree(v as usize), 3);
                }
            }
        }
    }
}

#[test]
fn intrinsic_trees_are_shorter() {
    let mut shapes = vec![
        jittered_grid(10, 0.4, 9),
        sheared_torus(30, 8, 1.0, 0.4, 5),
        Shape::generated("torus", torus(20, 8, 1.0, 0.3).unwrap()),
    ];
    shapes.extend(load_data("optimal-cones/optimal-cones-data/data/camelbad.obj"));
    for shape in &shapes {
        let tri = shape.triangulation();
        let (_, extrinsic) = minimum_spanning_tree(&tri);
        let mst = steiner_tree(
            &tri,
            &SteinerOptions {
                steiner_points: false,
                trace: None,
                ..SteinerOptions::default()
            },
        )
        .unwrap();
        let steiner = steiner_tree(
            &tri,
            &SteinerOptions {
                trace: None,
                ..SteinerOptions::default()
            },
        )
        .unwrap();
        println!(
            "{}: extrinsic MST {extrinsic}, intrinsic MST {}, intrinsic Steiner {}",
            shape.name, mst.length, steiner.length
        );
        assert!(mst.length <= extrinsic * (1.0 + 1e-12), "{}", shape.name);
        assert!(steiner.length < mst.length, "{}", shape.name);
    }
}

#[test]
fn surfaces_with_boundary_count_the_boundary_as_cut() {
    let shape = jittered_grid(8, 0.3, 1);
    let tri = shape.triangulation();
    let tree = steiner_tree(&tri, &SteinerOptions::default()).unwrap();
    assert_eq!(tree.report.n_boundary_loops, 1);
    check(&shape, &tree, "jittered grid");
    // An annulus: two loops, joined by the tree; the cut is a disk.
    let annulus =
        Shape::generated("annulus", ddg_core::generators::annulus(12, 3, 0.4, 1.0).unwrap());
    let tree = steiner_tree(&annulus.triangulation(), &SteinerOptions::default()).unwrap();
    assert_eq!(tree.report.n_boundary_loops, 2);
    check(&annulus, &tree, "annulus");
}

#[test]
fn bad_input_is_an_error() {
    let shape = Shape::generated("icosphere 1", icosphere(1).unwrap());
    assert_eq!(
        cut_to_disk(&shape.triangulation(), &[100_000]),
        Err(SteinerError::EdgeOutOfRange {
            edge: 100_000,
            n_edges: shape.mesh.n_edges()
        })
    );
}
