//! SPEC §5 item 5: geodesic distances between input vertices computed on the input mesh and
//! on the intrinsic triangulation agree to the accuracy of the distance solver.
//!
//! The solver is `ddg_core::distance::upwind_distance` (first order, single pass), run on the
//! input mesh and on the intrinsic Delaunay triangulation handed out as a `TriMesh` (the
//! shapes here have simplicial Delaunay triangulations). Both triangulate the same surface, so
//! both approximate the same distance; what differs is the error of the solver, which is of
//! the order of an edge length on a Delaunay triangulation and can be many edge lengths on
//! obtuse triangles.
//!
//! Run with `-- --nocapture` for the tables.

mod common;

use common::{Shape, dist, jittered_icosphere, sheared_parallelogram};
use ddg_core::distance::upwind_distance;
use ddg_core::generators;
use ddg_intrinsic::{DelaunayOptions, DelaunayStatus};

/// Mean and maximum of `|a − b|`.
fn difference(a: &[f64], b: &[f64]) -> (f64, f64) {
    let mut sum = 0.0;
    let mut max = 0.0f64;
    for (x, y) in a.iter().zip(b) {
        sum += (x - y).abs();
        max = max.max((x - y).abs());
    }
    (sum / a.len() as f64, max)
}

struct Distances {
    input: Vec<f64>,
    intrinsic: Vec<f64>,
    flips: usize,
}

fn distances(shape: &Shape, seed: u32) -> Distances {
    let input = upwind_distance(&shape.mesh, &shape.geometry, &[(seed, 0.0)], None).unwrap();
    assert_eq!(input.n_unreached, 0);
    let mut tri = shape.triangulation();
    let report = tri.flip_to_delaunay(&DelaunayOptions::default()).unwrap();
    assert_eq!(report.status, DelaunayStatus::Delaunay);
    let (mesh, geometry) = tri
        .to_tri_mesh()
        .unwrap_or_else(|e| panic!("{}: {e}", shape.name));
    let intrinsic = upwind_distance(&mesh, &geometry, &[(seed, 0.0)], None).unwrap();
    assert_eq!(intrinsic.n_unreached, 0);
    Distances {
        input: input.values,
        intrinsic: intrinsic.values,
        flips: report.n_flips,
    }
}

/// On a flat convex shape the exact distance is the Euclidean one, so the error of each
/// triangulation is known.
///
/// Measured, in units of the lattice width `h` (errors of the single-pass solver against the
/// exact distance; lattices of 16 to 48 cells, shear 0.5 to 6, seeds at two corners and in the
/// middle):
///
/// * intrinsic Delaunay triangulation: mean error 0.03 to 0.30, largest 0.06 to 0.96, the same
///   at every resolution: first order, as on any mesh without obtuse angles.
/// * input mesh, seed at the corner the shear leans away from: mean 0.004 to 0.07 (its long
///   edges happen to point along the paths from there, which makes it better than the
///   Delaunay triangulation).
/// * input mesh, other seeds: mean 0.4 to 12.6 and largest 1.3 to 32, growing with the
///   resolution: no convergence in units of `h`.
#[test]
fn on_sheared_lattices_the_intrinsic_distances_are_first_order_accurate() {
    println!(
        "{:<44} {:>6} {:>11} {:>10} {:>11} {:>10}",
        "shape (errors in units of h)",
        "flips",
        "mean input",
        "max input",
        "mean intr.",
        "max intr."
    );
    for (n, sigma) in [
        (16, 0.5),
        (16, 1.0),
        (32, 1.0),
        (16, 2.0),
        (24, 3.0),
        (24, 6.0),
    ] {
        let shape = sheared_parallelogram(n, sigma);
        let h = 1.0 / n as f64;
        // The corner the shear leans away from, the opposite corner of the bottom row, and
        // the middle.
        for (which, seed) in [0, n, n * (n + 1) / 2 + n / 2].into_iter().enumerate() {
            let d = distances(&shape, seed as u32);
            let exact: Vec<f64> = shape
                .positions
                .iter()
                .map(|&p| dist(p, shape.positions[seed]))
                .collect();
            let (mean_in, max_in) = difference(&d.input, &exact);
            let (mean_it, max_it) = difference(&d.intrinsic, &exact);
            println!(
                "{:<44} {:>6} {:>11.3} {:>10.3} {:>11.3} {:>10.3}",
                format!("{} seed {seed}", shape.name),
                d.flips,
                mean_in / h,
                max_in / h,
                mean_it / h,
                max_it / h
            );
            assert!(d.flips > 0);
            // The solver never undershoots the exact distance by more than rounding.
            for (v, (&x, &e)) in d.intrinsic.iter().zip(&exact).enumerate() {
                assert!(x >= e - 1e-12, "{}: vertex {v}: {x} < {e}", shape.name);
            }
            // First order on the Delaunay triangulation.
            assert!(max_it < 1.2 * h, "{}: {max_it} vs h = {h}", shape.name);
            assert!(mean_it < 0.4 * h, "{}: {mean_it} vs h = {h}", shape.name);
            // The two agree to the accuracy of the solver on each.
            let (_, apart) = difference(&d.input, &d.intrinsic);
            assert!(apart <= max_in + max_it, "{}", shape.name);
            if which > 0 {
                // Not on the input mesh.
                assert!(
                    mean_it < 0.5 * mean_in,
                    "{}: {mean_it} vs {mean_in}",
                    shape.name
                );
                assert!(
                    max_it < 0.5 * max_in,
                    "{}: {max_it} vs {max_in}",
                    shape.name
                );
            }
        }
    }
}

/// On curved shapes with triangles of good shape the two solutions agree to a fraction of an
/// edge length (measured: 0.02 to 0.04 `h` in the mean, 0.12 to 0.21 `h` at most), and on the
/// sphere both are within 0.4 `h` of the great-circle distance.
#[test]
fn on_curved_shapes_the_distances_agree() {
    println!(
        "{:<36} {:>6} {:>12} {:>12} {:>12}",
        "shape", "flips", "mean diff/h", "max diff/h", "max err/h"
    );
    let mut torus = generators::torus(40, 16, 1.0, 0.4).unwrap();
    generators::perturb(&mut torus.positions, [0.01; 3], 3);
    let shapes = [
        (jittered_icosphere(3, 0.25, 4), true),
        (jittered_icosphere(4, 0.25, 9), true),
        (
            Shape::generated("perturbed torus 40x16", torus, false),
            false,
        ),
        (
            Shape::generated(
                "bumpy grid",
                generators::bumpy_grid_disk(24, 24, 0.15).unwrap(),
                false,
            ),
            false,
        ),
    ];
    for (shape, sphere) in &shapes {
        let h = shape.mean_edge_length();
        let d = distances(shape, 0);
        let (mean, max) = difference(&d.input, &d.intrinsic);
        // On the unit sphere the distance on the polyhedron is the great-circle distance up
        // to O(h²).
        let mut error = 0.0;
        if *sphere {
            let p0 = shape.positions[0];
            let exact: Vec<f64> = shape
                .positions
                .iter()
                .map(|p| {
                    (p[0] * p0[0] + p[1] * p0[1] + p[2] * p0[2])
                        .clamp(-1.0, 1.0)
                        .acos()
                })
                .collect();
            let (_, max_in) = difference(&d.input, &exact);
            let (_, max_it) = difference(&d.intrinsic, &exact);
            assert!(max_it < 0.6 * h, "{}: {max_it} vs h = {h}", shape.name);
            assert!(max_in < 0.6 * h, "{}: {max_in} vs h = {h}", shape.name);
            error = max_it;
        }
        println!(
            "{:<36} {:>6} {:>12.3} {:>12.3} {:>12.3}",
            shape.name,
            d.flips,
            mean / h,
            max / h,
            error / h
        );
        assert!(d.flips > 0, "{}", shape.name);
        assert!(
            max < 0.5 * h && mean < 0.1 * h,
            "{}: {mean} {max} vs h = {h}",
            shape.name
        );
    }
}
