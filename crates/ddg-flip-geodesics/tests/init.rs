//! Initial paths: Dijkstra on the current triangulation (F15) and paths through vertices
//! (SPEC §2.2).

mod common;

use common::*;
use ddg_flip_geodesics::{FlipError, Network, dijkstra_path, path_from_vertices};
use ddg_intrinsic::IntrinsicTriangulation;

/// Single-source distances by Bellman–Ford over the edges (an independent algorithm).
fn distances(tri: &IntrinsicTriangulation, from: usize) -> Vec<f64> {
    let c = tri.complex();
    let mut d = vec![f64::INFINITY; c.n_vertices()];
    d[from] = 0.0;
    loop {
        let mut changed = false;
        for e in 0..c.n_edges() {
            let [a, b] = c.edge_vertices(e);
            let l = tri.edge_length(e);
            for (x, y) in [(a, b), (b, a)] {
                if d[x] + l < d[y] {
                    d[y] = d[x] + l;
                    changed = true;
                }
            }
        }
        if !changed {
            return d;
        }
    }
}

#[test]
fn dijkstra_paths_are_shortest_edge_paths() {
    let mut shapes = generated_shapes();
    shapes.extend(
        data_shapes()
            .into_iter()
            .filter(|s| s.mesh.n_vertices() < 4000),
    );
    for shape in &shapes {
        for flipped in [false, true] {
            let mut tri = shape.triangulation();
            if flipped {
                let n = 2 * tri.complex().n_edges();
                random_flips(&mut tri, n, 9);
            }
            let c = tri.complex();
            for (a, b) in random_pairs(c.n_vertices(), 4, 21) {
                let d = distances(&tri, a);
                match dijkstra_path(&tri, a, b) {
                    Ok(path) => {
                        let v = path.vertices(c);
                        assert_eq!((v[0], *v.last().unwrap()), (a, b));
                        let length = path.length(&tri);
                        assert!(
                            (length - d[b]).abs() <= 1e-12 * d[b],
                            "{}: {a} -> {b}: {length} vs {}",
                            shape.name,
                            d[b]
                        );
                        // A shortest path is simple, so it is a valid network.
                        let mut sorted = v.clone();
                        sorted.sort_unstable();
                        sorted.dedup();
                        assert_eq!(sorted.len(), v.len());
                        Network::new(c, std::slice::from_ref(&path), &[], None).unwrap();
                        // Deterministic.
                        assert_eq!(dijkstra_path(&tri, a, b).unwrap(), path);
                    }
                    Err(e) => {
                        assert_eq!(e, FlipError::NotConnected { from: a, to: b });
                        assert!(d[b].is_infinite());
                    }
                }
            }
        }
    }
}

/// Ties by `(distance, halfedge id)`, a boundary edge walked backwards after every halfedge.
#[test]
fn dijkstra_ties() {
    // The unit square split along a–c: b → d is 2 both ways round.
    let (a, b, c, d) = (0, 1, 2, 3);
    let shape = Shape::new(
        "square",
        vec![
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [1.0, 1.0, 0.0],
            [0.0, 1.0, 0.0],
        ],
        vec![
            [a as u32, b as u32, c as u32],
            [a as u32, c as u32, d as u32],
        ],
    );
    let tri = shape.triangulation();
    let cx = tri.complex();
    // From b: b → c runs along halfedge 1 and b → a against boundary halfedge 0, so c is
    // settled first, and d is reached from c along halfedge 4 before from a.
    assert_eq!(
        dijkstra_path(&tri, b, d).unwrap().vertices(cx),
        vec![b, c, d]
    );
    // From d: d → a is halfedge 5, d → c is against boundary halfedge 4; a first, then a → b
    // along halfedge 0 before c → b against boundary halfedge 1.
    assert_eq!(
        dijkstra_path(&tri, d, b).unwrap().vertices(cx),
        vec![d, a, b]
    );
}

#[test]
fn dijkstra_errors() {
    // Two separate triangles.
    let shape = Shape::new(
        "two triangles",
        vec![
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [5.0, 0.0, 0.0],
            [6.0, 0.0, 0.0],
            [5.0, 1.0, 0.0],
        ],
        vec![[0, 1, 2], [3, 4, 5]],
    );
    let tri = shape.triangulation();
    let err = dijkstra_path(&tri, 1, 1).unwrap_err();
    assert_eq!(err, FlipError::SameVertex { vertex: 1 });
    assert!(err.to_string().contains("vertex 1"), "{err}");
    let err = dijkstra_path(&tri, 0, 4).unwrap_err();
    assert_eq!(err, FlipError::NotConnected { from: 0, to: 4 });
    assert!(err.to_string().contains("0 and 4"), "{err}");
    assert!(matches!(
        dijkstra_path(&tri, 0, 6),
        Err(FlipError::IndexOutOfRange {
            what: "vertex",
            index: 6,
            bound: 6
        })
    ));
    assert_eq!(dijkstra_path(&tri, 3, 5).unwrap().segments.len(), 1);
}

#[test]
fn paths_through_vertices() {
    let shape = grid(3, 3);
    let tri = shape.triangulation();
    let g = |i, j| gv(3, i, j);
    let p = path_from_vertices(&tri, &[g(0, 0), g(1, 0), g(1, 1)], false).unwrap();
    assert_eq!(p.vertices(tri.complex()), vec![g(0, 0), g(1, 0), g(1, 1)]);
    assert!(!p.closed);
    let q = path_from_vertices(&tri, &[g(1, 1), g(2, 1), g(2, 2), g(1, 2)], true).unwrap();
    assert_eq!(q.segments.len(), 4);
    assert_eq!(
        q.vertices(tri.complex()),
        vec![g(1, 1), g(2, 1), g(2, 2), g(1, 2)]
    );
    assert_eq!(q.reversed().reversed(), q);
    assert_eq!(
        path_from_vertices(&tri, &[g(0, 0), g(2, 2)], false).unwrap_err(),
        FlipError::NotAdjacent {
            index: 0,
            from: g(0, 0),
            to: g(2, 2)
        }
    );
    assert_eq!(
        path_from_vertices(&tri, &[g(0, 0)], false).unwrap_err(),
        FlipError::TooFewVertices {
            n: 1,
            closed: false
        }
    );
    assert_eq!(
        path_from_vertices(&tri, &[], true).unwrap_err(),
        FlipError::TooFewVertices { n: 0, closed: true }
    );

    // A Δ-complex: a loop along a self-edge from one vertex, and the shorter of two edges.
    let needle = needle_tetrahedron(8.0, 4.0).delaunay();
    let irregular = needle.irregularities();
    let s = irregular.self_edges[0] as usize;
    let y = needle.complex().edge_vertices(s)[0];
    let lp = path_from_vertices(&needle, &[y], true).unwrap();
    assert_eq!(lp.segments.len(), 1);
    assert_eq!(lp.segments[0].edge as usize, s);
    for &e in &irregular.multi_edges {
        let [u, w] = needle.complex().edge_vertices(e as usize);
        let chosen = path_from_vertices(&needle, &[u, w], false)
            .unwrap()
            .segments[0];
        let shortest = irregular
            .multi_edges
            .iter()
            .filter(|&&f| {
                let [x, z] = needle.complex().edge_vertices(f as usize);
                (x, z) == (u, w) || (x, z) == (w, u)
            })
            .map(|&f| needle.edge_length(f as usize))
            .fold(f64::INFINITY, f64::min);
        assert_eq!(needle.edge_length(chosen.edge as usize), shortest);
    }
}
