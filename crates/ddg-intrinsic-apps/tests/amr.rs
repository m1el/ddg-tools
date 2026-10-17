//! Tests of the adaptive refinement (SPEC §5.4): the solves, the estimator, the marking, the
//! hook in `ddg-intrinsic`, and convergence to the uniformly refined solution.

mod common;

use common::*;
use ddg_core::generators::SplitMix64;
use ddg_intrinsic::{DelaunayOptions, InsertionLimit, IntrinsicTriangulation, RefineOptions};
use ddg_intrinsic_apps::amr::*;

fn uniform(shape: &Shape, vertices: usize) -> IntrinsicTriangulation {
    let mut tri = shape.triangulation();
    let r = (shape.geometry.total_area() / (2.6 * vertices as f64)).sqrt();
    tri.refine(&RefineOptions {
        max_circumradius: Some(r),
        max_insertions: InsertionLimit::Unlimited,
        ..Default::default()
    })
    .unwrap();
    tri
}

/// A linear solution is exact: on a flat disk bounded by a regular polygon on the unit circle
/// the harmonic disk map is a rotation, reproduced exactly, and the estimator vanishes. On a
/// curved disk it does not.
#[test]
fn the_estimator_vanishes_on_an_exact_solution() {
    for twist in [0.0, 0.25] {
        let shape = polar_disk(6, twist);
        let mut idt = shape.triangulation();
        idt.flip_to_delaunay(&DelaunayOptions::default()).unwrap();
        // (Refinement would put boundary vertices on the chords of the polygon, which the map
        // sends to the circle: then the map is no longer linear.)
        for (what, tri) in [("input", shape.triangulation()), ("iDT", idt)] {
            let sol = solve(&tri, &Problem::HarmonicDisk).unwrap();
            // The map is a rotation (the arc length starts at the boundary vertex of smallest id).
            let start = (0..shape.positions.len())
                .find(|&v| shape.mesh.is_boundary_vertex(v))
                .unwrap();
            let a = shape.positions[start][1].atan2(shape.positions[start][0]);
            for v in 0..tri.n_input_vertices() {
                let p = shape.positions[v];
                let (s, c) = (-a).sin_cos();
                let expected = [c * p[0] - s * p[1], s * p[0] + c * p[1]];
                assert!((sol.components[0][v] - expected[0]).abs() < 1e-12, "{what} twist {twist}");
                assert!((sol.components[1][v] - expected[1]).abs() < 1e-12);
            }
            let est = estimate(&tri, &Problem::HarmonicDisk, &sol).unwrap();
            assert!(est.total < 1e-11, "{what} twist {twist}: {}", est.total);
            let q = map_quality(&tri, &sol);
            assert_eq!(q.n_flipped, 0);
            assert!((q.mean_conformal - 1.0).abs() < 1e-9, "{what} twist {twist}: {q:?}");
        }
    }
    let bumpy = bumpy_sheared_grid(8, 0.3, 0.0);
    let tri = bumpy.triangulation();
    let sol = solve(&tri, &Problem::HarmonicDisk).unwrap();
    assert!(estimate(&tri, &Problem::HarmonicDisk, &sol).unwrap().total > 1e-2);
}

/// The estimator of a smooth solution decreases under uniform refinement, roughly like the
/// mesh size; the oscillation of the heat kernel faster.
#[test]
fn the_estimator_decreases_under_uniform_refinement() {
    let shape = bumpy_sheared_grid(8, 0.3, 0.8);
    let mut last = f64::INFINITY;
    for vertices in [500, 2000, 8000] {
        let tri = uniform(&shape, vertices);
        let sol = solve(&tri, &Problem::HarmonicDisk).unwrap();
        let est = estimate(&tri, &Problem::HarmonicDisk, &sol).unwrap();
        eprintln!("{} vertices: estimate {:.4e}", tri.n_vertices(), est.total);
        assert!(est.total < 0.7 * last);
        last = est.total;
    }
}

/// The solves: the Green's function has mean zero and solves its system; one backward Euler
/// step conserves the integral; bad input is refused with a reason.
#[test]
fn the_model_problems_are_solved() {
    let shape = jittered_icosphere(2, 0.1, 3);
    let mut tri = shape.triangulation();
    tri.flip_to_delaunay(&DelaunayOptions::default()).unwrap();
    let green = solve(&tri, &Problem::GreenFunction { source: 5 }).unwrap();
    assert!(green.residual < 1e-12);
    let mass = tri.lumped_mass();
    let mean: f64 = green.components[0].iter().zip(&mass).map(|(u, m)| u * m).sum();
    assert!(mean.abs() < 1e-12);
    // The source is the maximum.
    let u = &green.components[0];
    assert!(u.iter().all(|&x| x <= u[5]));
    let heat = solve(&tri, &Problem::HeatKernel { source: 5, time: 0.01 }).unwrap();
    let total: f64 = heat.components[0].iter().zip(&mass).map(|(u, m)| u * m).sum();
    assert!((total - 1.0).abs() < 1e-12, "{total}");
    assert!(matches!(
        solve(&tri, &Problem::HeatKernel { source: 5, time: 0.0 }),
        Err(AmrError::InvalidOptions { .. })
    ));
    assert!(matches!(
        solve(&tri, &Problem::GreenFunction { source: 10_000 }),
        Err(AmrError::InvalidInput { .. })
    ));
    assert!(matches!(
        solve(&tri, &Problem::HarmonicDisk),
        Err(AmrError::InvalidInput { .. })
    ));
}

/// The bulk criterion: the marked faces carry the fraction of the total, and are the largest.
#[test]
fn marking_meets_the_bulk_criterion() {
    let mut rng = SplitMix64::new(9);
    let n = 500;
    let indicator: Vec<f64> = (0..n).map(|_| rng.uniform(0.0, 1.0).powi(4)).collect();
    let oscillation: Vec<f64> = (0..n).map(|_| rng.uniform(0.0, 1.0).powi(8)).collect();
    let est = Estimate {
        total: indicator.iter().map(|x| x * x).sum::<f64>().sqrt(),
        total_oscillation: oscillation.iter().map(|x| x * x).sum::<f64>().sqrt(),
        indicator: indicator.clone(),
        oscillation: oscillation.clone(),
    };
    for (sum, theta) in [(BulkSum::Indicators, 0.5), (BulkSum::Squares, 0.5), (BulkSum::Indicators, 0.8)] {
        let marked = mark(&est, theta, sum);
        let value = |x: f64| if sum == BulkSum::Squares { x * x } else { x };
        let frac = if sum == BulkSum::Squares { theta * theta } else { theta };
        for data in [&indicator, &oscillation] {
            let total: f64 = data.iter().map(|&x| value(x)).sum();
            let have: f64 = marked.iter().map(|&f| value(data[f as usize])).sum();
            assert!(have >= frac * total);
        }
        // Without the marked face of smallest indicator, the indicators alone would not reach it
        // unless the oscillation needed that face.
        assert!(marked.len() < n);
        let unmarked_max = (0..n)
            .filter(|f| !marked.contains(&(*f as u32)))
            .map(|f| indicator[f])
            .fold(0.0f64, f64::max);
        let n_above = marked.iter().filter(|&&f| indicator[f as usize] >= unmarked_max).count();
        assert!(n_above * 2 >= marked.len(), "the marked faces are mostly the largest");
    }
}

/// The hook in `ddg-intrinsic`: with no face marked, refinement is the plain one bit for bit;
/// every marked face is refined.
#[test]
fn refining_marked_faces() {
    let shape = jittered_icosphere(2, 0.1, 21);
    let mut plain = shape.triangulation();
    plain.refine(&RefineOptions::default()).unwrap();
    let mut hooked = shape.triangulation();
    refine_marked(&mut hooked, &RefineOptions::default(), &[]).unwrap();
    assert_eq!(plain.fingerprint(), hooked.fingerprint());

    let mut tri = plain.clone();
    let marked: Vec<u32> = (0..tri.n_faces() as u32).step_by(7).collect();
    let before: Vec<[usize; 3]> = marked
        .iter()
        .map(|&f| [0, 1, 2].map(|k| tri.edge(3 * f as usize + k)))
        .collect();
    let report = refine_marked(&mut tri, &RefineOptions::default(), &marked).unwrap();
    assert!(report.n_inserted >= marked.len() / 3, "{}", report.n_inserted);
    for (i, &f) in marked.iter().enumerate() {
        let now = [0, 1, 2].map(|k| tri.edge(3 * f as usize + k));
        assert_ne!(now, before[i], "marked face {f} is unchanged");
    }
    assert!(tri.check().unwrap().max_angle_sum_error < 1e-9);
    assert!(refine_marked(&mut tri, &RefineOptions::default(), &[u32::MAX]).is_err());
}

/// The adaptive solution converges to the uniformly refined one (relative L² error against a
/// fine uniform reference), and with fewer vertices than uniform refinement for the same error
/// on a localized problem.
#[test]
fn adaptive_refinement_converges_to_the_uniform_solution() {
    let shape = jittered_icosphere(2, 0.05, 2);
    let problem = Problem::HeatKernel { source: 3, time: 0.002 };
    // The reference: adaptive (fine at the source, where the kernel has a logarithmic
    // singularity), then uniformly fine.
    let mut reference_tri = shape.triangulation();
    let options = AmrOptions {
        max_vertices: 40_000,
        ..Default::default()
    };
    adaptive_solve(&mut reference_tri, &problem, &options).unwrap();
    let r = (shape.geometry.total_area() / (2.6 * 80_000.0)).sqrt();
    reference_tri
        .refine(&RefineOptions {
            max_circumradius: Some(r),
            max_insertions: InsertionLimit::Unlimited,
            ..Default::default()
        })
        .unwrap();
    let reference = solve(&reference_tri, &problem).unwrap();
    let mut errors = Vec::new();
    for max_vertices in [1000, 4000, 12_000] {
        let mut tri = shape.triangulation();
        let options = AmrOptions {
            max_vertices,
            ..Default::default()
        };
        let (sol, report) = adaptive_solve(&mut tri, &problem, &options).unwrap();
        assert_eq!(report.status, AmrStatus::VertexLimit);
        assert!(tri.n_vertices() <= max_vertices);
        let error = l2_error(&reference_tri, &reference, &tri, &sol).unwrap();
        eprintln!(
            "adaptive {} vertices ({} steps): L2 error {:.3e}, estimate {:.3e}",
            tri.n_vertices(),
            report.iterations.len(),
            error,
            report.iterations.last().unwrap().estimate
        );
        errors.push((tri.n_vertices(), error));
    }
    assert!(errors[1].1 < 0.6 * errors[0].1 && errors[2].1 < 0.6 * errors[1].1, "{errors:?}");
    assert!(errors[2].1 < 2e-2, "{errors:?}");
    // Uniform refinement with the vertices of the finest adaptive run is less accurate.
    let tri = uniform(&shape, errors[2].0);
    let sol = solve(&tri, &problem).unwrap();
    let error = l2_error(&reference_tri, &reference, &tri, &sol).unwrap();
    eprintln!("uniform {} vertices: L2 error {error:.3e}", tri.n_vertices());
    assert!(error > errors[2].1);
}

/// Paper Fig. 20: the harmonic map of a disk is injective on the intrinsic Delaunay
/// triangulation (non-negative weights, a convex boundary polygon: Tutte), and after adaptive
/// refinement.
#[test]
fn the_harmonic_disk_map_is_injective_on_the_intrinsic_delaunay_triangulation() {
    let mut shapes = vec![sheared_grid(10, 2.5, 0.45, 5), bumpy_sheared_grid(10, 0.3, 2.0)];
    shapes.extend(dataset("bump"));
    for shape in &shapes {
        if shape.mesh.boundary_loops().len() != 1 || shape.mesh.n_vertices() as i64 - shape.mesh.n_edges() as i64 + shape.mesh.n_faces() as i64 != 1 {
            continue;
        }
        let input = shape.triangulation();
        let sol = solve(&input, &Problem::HarmonicDisk).unwrap();
        let q_in = map_quality(&input, &sol);
        let mut tri = shape.triangulation();
        tri.flip_to_delaunay(&DelaunayOptions::default()).unwrap();
        let sol = solve(&tri, &Problem::HarmonicDisk).unwrap();
        let q = map_quality(&tri, &sol);
        eprintln!(
            "{}: input {} flipped (min area ratio {:.2e}), iDT {} flipped (min {:.2e})",
            shape.name, q_in.n_flipped, q_in.min_area_ratio, q.n_flipped, q.min_area_ratio
        );
        assert_eq!(q.n_flipped, 0, "{}", shape.name);
        if shape.name == "sheared_grid_10" {
            // The point of the example: the map of the input mesh folds over.
            assert!(q_in.n_flipped > 0);
        }
        let options = AmrOptions {
            max_vertices: 4 * shape.mesh.n_vertices().max(500),
            ..Default::default()
        };
        let (sol, _) = adaptive_solve(&mut tri, &Problem::HarmonicDisk, &options).unwrap();
        assert_eq!(map_quality(&tri, &sol).n_flipped, 0);
    }
}
