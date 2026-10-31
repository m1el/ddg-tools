//! Flat domains (SPEC §7 tests 1, 2 and the boundary corners of test 3): in a convex planar
//! polygon the geodesic between two vertices is the straight segment; in an L the shortest
//! path bends at the reflex corner, which it keeps; a path around a convex boundary corner is
//! pushed off it.

mod common;

use common::{
    Outcome, Shape, dist, grid, gv, l_distance, l_shape, path_vertices, random_pairs, run,
    sheared_parallelogram, vpath,
};
use ddg_core::generators::{self, SplitMix64};
use ddg_flip_geodesics::{ShortenOptions, ShortenStatus, dijkstra_path};
use ddg_intrinsic::IntrinsicTriangulation;

/// A flat grid with its interior vertices moved in the plane: convex, flat, irregular.
fn jittered_grid(n: usize, amplitude: f64, seed: u64) -> Shape {
    let g = generators::grid_disk(n, n).unwrap();
    let h = 1.0 / n as f64;
    let mut rng = SplitMix64::new(seed);
    let mut positions = g.positions.clone();
    for i in 1..n {
        for j in 1..n {
            let p = &mut positions[i * (n + 1) + j];
            for c in p.iter_mut().take(2) {
                let u = (rng.next_u64() >> 11) as f64 / (1u64 << 53) as f64;
                *c += amplitude * h * (2.0 * u - 1.0);
            }
        }
    }
    Shape::new(&format!("jittered grid {n}"), positions, g.faces)
}

fn shortest(
    shape: &Shape,
    tri: &IntrinsicTriangulation,
    a: usize,
    b: usize,
) -> Outcome<IntrinsicTriangulation> {
    let path = dijkstra_path(tri, a, b).unwrap();
    let out = run(tri.clone(), &[path], &[], &ShortenOptions::default());
    assert_eq!(
        out.report.status,
        ShortenStatus::Geodesic,
        "{}: {a} -> {b}",
        shape.name
    );
    out
}

#[test]
fn flat_convex_domains_give_straight_segments() {
    let mut worst: f64 = 0.0;
    let mut n_runs = 0;
    let mut n_flips = 0;
    for shape in [
        grid(10, 10),
        grid(7, 4),
        sheared_parallelogram(8, 3.0),
        sheared_parallelogram(10, 6.0),
        jittered_grid(12, 0.35, 7),
    ] {
        for tri in [shape.triangulation(), shape.delaunay()] {
            for (a, b) in random_pairs(shape.mesh.n_vertices(), 25, 11) {
                let out = shortest(&shape, &tri, a, b);
                assert!(
                    out.certificate.geodesic,
                    "{}: {a} -> {b}: {:?}",
                    shape.name, out.certificate.failures
                );
                let (pa, pb) = (shape.positions[a], shape.positions[b]);
                let exact = dist(pa, pb);
                let length = out.network.length();
                let rel = (length - exact).abs() / exact;
                assert!(
                    rel <= 1e-12,
                    "{}: {a} -> {b}: length {length} vs {exact} (rel {rel:e})",
                    shape.name
                );
                // Every vertex of the path lies on the segment.
                for v in path_vertices(&out.network, 0) {
                    let p = shape.positions[v];
                    let off = dist(pa, p) + dist(p, pb) - exact;
                    assert!(
                        off <= 1e-12 * exact,
                        "{}: vertex {v} is off the segment by {off:e}",
                        shape.name
                    );
                }
                worst = worst.max(rel);
                n_runs += 1;
                n_flips += out.report.flips;
            }
        }
    }
    eprintln!(
        "flat convex: {n_runs} paths, {n_flips} flips, worst relative length error {worst:e}"
    );
}

#[test]
fn flat_l_shape_bends_at_the_reflex_corner() {
    let shape = l_shape(6);
    let corner = shape
        .positions
        .iter()
        .position(|p| dist(*p, [1.0, 1.0, 0.0]) < 1e-12)
        .unwrap();
    let (mut bent, mut worst) = (0, 0.0f64);
    for tri in [shape.triangulation(), shape.delaunay()] {
        for (a, b) in random_pairs(shape.mesh.n_vertices(), 150, 5) {
            let out = shortest(&shape, &tri, a, b);
            assert!(
                out.certificate.geodesic,
                "{a} -> {b}: {:?}",
                out.certificate.failures
            );
            let exact = l_distance(shape.positions[a], shape.positions[b]);
            let rel = (out.network.length() - exact).abs() / exact;
            assert!(
                rel <= 1e-12,
                "{a} -> {b}: {} vs {exact} (rel {rel:e})",
                out.network.length()
            );
            worst = worst.max(rel);
            if exact > dist(shape.positions[a], shape.positions[b]) + 1e-9 {
                assert!(
                    path_vertices(&out.network, 0).contains(&corner),
                    "{a} -> {b} misses the corner"
                );
                bent += 1;
            }
        }
    }
    assert!(bent >= 10, "only {bent} paths bend at the reflex corner");
    eprintln!(
        "L-shape: {bent} paths bend at the reflex corner (kept); worst relative error {worst:e}"
    );
}

#[test]
fn a_convex_boundary_corner_pushes_the_path_off() {
    let shape = grid(4, 4);
    let tri = shape.triangulation();
    // (0, 1) -> (0, 0) -> (1, 0) around the corner (0, 0), interior angle pi/2.
    let path = vpath(&tri, &[gv(4, 0, 1), gv(4, 0, 0), gv(4, 1, 0)], false);
    let out = run(tri, &[path], &[], &ShortenOptions::default());
    assert_eq!(out.report.status, ShortenStatus::Geodesic);
    assert!(out.certificate.geodesic);
    assert!(!path_vertices(&out.network, 0).contains(&gv(4, 0, 0)));
    let exact = 2f64.sqrt() / 4.0;
    assert!(
        (out.network.length() - exact).abs() <= 1e-15,
        "{}",
        out.network.length()
    );
}
