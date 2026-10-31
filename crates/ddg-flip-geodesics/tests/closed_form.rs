//! Closed forms (SPEC §7 test 3): the lateral surface of a pyramid (a polyhedral cone whose
//! apex is never on a geodesic), a cube corner to the opposite corner (`√5 a`), a loop around
//! a prism tube (the cross-section's perimeter), a path through a saddle with `Θ = 3π`
//! (unchanged, no flip).

mod common;

use std::f64::consts::PI;

use common::{cube, dist, open_pyramid, path_vertices, run, saddle, tube, vertex_at, vpath};
use ddg_flip_geodesics::{ShortenOptions, ShortenStatus, dijkstra_path};

fn rel(a: f64, b: f64) -> f64 {
    (a - b).abs() / b.abs()
}

/// From one base corner to the opposite one over the lateral faces: unfolded, two faces with
/// apex angle `θ` each (`sin(θ/2) = 1/s`, `s` the slant edge) give the chord `2 s sin θ`
/// when `2θ < π`. Dijkstra goes over the apex (`2s < 4`), whose angle sum `4θ < 2π` makes
/// both wedges smaller than `π`.
#[test]
fn pyramid_lateral_surface() {
    for h in [0.1, 0.2, 0.5, 1.0, 2.0] {
        let shape = open_pyramid(h);
        let s = (2.0 + h * h).sqrt();
        let theta = 2.0 * (1.0 / s).asin();
        assert!(2.0 * theta < PI);
        let exact = 2.0 * s * theta.sin();
        for tri in [shape.triangulation(), shape.delaunay()] {
            for (a, b) in [(0, 2), (2, 0), (1, 3), (3, 1)] {
                let path = dijkstra_path(&tri, a, b).unwrap();
                let out = run(tri.clone(), &[path], &[], &ShortenOptions::default());
                assert_eq!(out.report.status, ShortenStatus::Geodesic);
                assert!(
                    out.certificate.geodesic,
                    "h {h}: {:?}",
                    out.certificate.failures
                );
                assert!(
                    !path_vertices(&out.network, 0).contains(&4),
                    "h {h}: through the apex"
                );
                let e = rel(out.network.length(), exact);
                assert!(
                    e <= 1e-12,
                    "h {h}, {a} -> {b}: {} vs {exact} ({e:e})",
                    out.network.length()
                );
            }
        }
    }
}

/// Corner to opposite corner of the cube of side `a`: `√5 a` (two faces unfolded), on grids of
/// several sizes; the corners themselves (`Θ = 3π/2`) are never passed.
#[test]
fn cube_corner_to_opposite_corner() {
    for (n, a) in [(1, 1.0), (2, 1.0), (3, 2.5), (4, 1.0), (6, 0.7)] {
        let shape = cube(n, a);
        let corners: Vec<usize> = [[0.0, 0.0, 0.0], [a, a, a], [a, 0.0, 0.0], [0.0, a, a]]
            .iter()
            .map(|&p| vertex_at(&shape, p))
            .collect();
        let exact = 5f64.sqrt() * a;
        for tri in [shape.triangulation(), shape.delaunay()] {
            for (from, to) in [(corners[0], corners[1]), (corners[2], corners[3])] {
                let path = dijkstra_path(&tri, from, to).unwrap();
                let out = run(tri.clone(), &[path], &[], &ShortenOptions::default());
                assert!(
                    out.certificate.geodesic,
                    "cube {n}: {:?}",
                    out.certificate.failures
                );
                let e = rel(out.network.length(), exact);
                assert!(
                    e <= 1e-12,
                    "cube {n}: {} vs {exact} ({e:e})",
                    out.network.length()
                );
                let interior = &path_vertices(&out.network, 0);
                let interior = &interior[1..interior.len() - 1];
                for v in interior {
                    let p = shape.positions[*v];
                    let at_corner = p.iter().all(|&c| c.abs() < 1e-12 || (c - a).abs() < 1e-12);
                    assert!(!at_corner, "cube {n}: the path passes the corner {v}");
                }
            }
        }
    }
}

/// A zig-zag loop around a prism tube (`n`-gon of circumradius 1) shortens to a ring: the
/// cross-section's perimeter `2 n sin(π/n)`.
#[test]
fn loop_around_a_prism_tube() {
    for n in [3, 4, 6, 8, 12] {
        let shape = tube(n, &[1.0; 5], 0.4);
        let v = |j: usize, k: usize| j * n + k % n;
        // Blocks (1, k) -> (2, k + 1) -> (1, k + 1), k even; for odd n the last block is
        // (1, n - 1) -> back to (1, 0).
        let mut vertices = Vec::new();
        let mut k = 0;
        while k + 1 < n {
            vertices.extend([v(1, k), v(2, k + 1), v(1, k + 1)]);
            k += 2;
        }
        if k < n {
            vertices.push(v(1, k));
        }
        let exact = 2.0 * n as f64 * (PI / n as f64).sin();
        for tri in [shape.triangulation(), shape.delaunay()] {
            let path = vpath(&tri, &vertices, true);
            let out = run(tri, &[path], &[], &ShortenOptions::default());
            assert_eq!(out.report.status, ShortenStatus::Geodesic, "n {n}");
            assert!(
                out.certificate.geodesic,
                "n {n}: {:?}",
                out.certificate.failures
            );
            assert!(out.report.contracted.is_empty());
            let e = rel(out.network.length(), exact);
            assert!(
                e <= 1e-12,
                "n {n}: {} vs {exact} ({e:e})",
                out.network.length()
            );
            let rings: Vec<usize> = path_vertices(&out.network, 0)
                .iter()
                .map(|&x| x / n)
                .collect();
            assert!(
                rings.iter().all(|&r| r == rings[0]),
                "n {n}: not on one ring: {rings:?}"
            );
        }
    }
}

/// Through the centre of a fan of six right isosceles triangles (`Θ = 3π`): opposite ring
/// vertices have both wedges `3π/2 ≥ π`, so the path stays as it is with no flip and no
/// iteration. Two corners apart (`π` and `2π`) it is straight too; one corner apart (`π/2`) it
/// is shortened to the ring edge.
#[test]
fn saddle_with_three_pi_keeps_the_path() {
    let shape = saddle();
    for tri in [shape.triangulation(), shape.delaunay()] {
        for (a, b) in [(1, 4), (2, 5), (3, 6), (1, 3), (2, 6)] {
            let path = vpath(&tri, &[a, 0, b], false);
            let before = path.clone();
            let out = run(tri.clone(), &[path], &[], &ShortenOptions::default());
            assert_eq!(out.report.status, ShortenStatus::Geodesic);
            assert_eq!(
                (out.report.flips, out.report.iterations),
                (0, 0),
                "{a} 0 {b}"
            );
            assert_eq!(out.network.paths(), vec![before]);
            assert!(out.certificate.geodesic);
        }
        let path = vpath(&tri, &[1, 0, 2], false);
        let out = run(tri.clone(), &[path], &[], &ShortenOptions::default());
        assert!(out.certificate.geodesic);
        assert_eq!(path_vertices(&out.network, 0), vec![1, 2]);
        assert!(
            rel(
                out.network.length(),
                dist(shape.positions[1], shape.positions[2])
            ) < 1e-15
        );
    }
}
