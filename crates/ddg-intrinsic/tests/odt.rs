//! The intrinsic optimal Delaunay triangulation (paper Sec. 4.3; SPEC N21 to N23, §5 item 11):
//! the surface does not change, the result is Delaunay after every iteration, input vertices
//! stay, the run is deterministic, and the options and limits do what they say.

mod common;

use common::{Shape, assert_consistent, assert_point_queries, bad_meshes, models, standard_shapes};
use ddg_core::generators::{self, SplitMix64};
use ddg_intrinsic::{
    BoundaryVertices, DelaunayOptions, InsertionLimit, IntrinsicError, IntrinsicTriangulation,
    OdtCenters, OptimalDelaunayOptions, OptimalDelaunayStatus, RefineOptions, SurfacePoint,
    TargetLength,
};

fn options(target: TargetLength, iterations: usize) -> OptimalDelaunayOptions {
    OptimalDelaunayOptions {
        target_length: target,
        max_iterations: iterations,
        ..OptimalDelaunayOptions::default()
    }
}

fn assert_input_vertices_stay(shape: &Shape, tri: &IntrinsicTriangulation) {
    for v in 0..shape.mesh.n_vertices() {
        assert_eq!(tri.vertex_location(v), SurfacePoint::Vertex(v as u32), "{}", shape.name);
        assert_eq!(
            tri.vertex_angle_sum(v),
            tri.input().vertex_angle_sum(v),
            "{}",
            shape.name
        );
    }
}

#[test]
fn the_optimal_delaunay_triangulation_keeps_the_surface() {
    for (s, shape) in standard_shapes().iter().enumerate() {
        let mut rng = SplitMix64::new(900 + s as u64);
        for (iterations, boundary) in [
            (1, BoundaryVertices::Fixed),
            (3, BoundaryVertices::AlongBoundary),
            (12, BoundaryVertices::Fixed),
        ] {
            let mut tri = shape.triangulation();
            let report = tri
                .optimal_delaunay(&OptimalDelaunayOptions {
                    boundary,
                    ..options(TargetLength::MeanInputEdge(1.0), iterations)
                })
                .unwrap_or_else(|e| panic!("{}: {e}", shape.name));
            let what = format!("iODT, {iterations} iterations, {boundary:?}");
            // Delaunay after the last of these iterations: by the report and by a fresh test.
            assert_eq!(report.delaunay.n_not_delaunay, 0, "{} ({what})", shape.name);
            assert!(tri.is_delaunay(&DelaunayOptions::default().criterion), "{}", shape.name);
            assert_eq!(report.n_iterations, report.iterations.len());
            assert!(report.n_iterations <= iterations);
            assert_eq!(
                tri.n_vertices(),
                shape.mesh.n_vertices() + report.n_inserted,
                "{}",
                shape.name
            );
            let stats = tri.shape_statistics();
            assert_eq!(stats, report.statistics);
            assert_input_vertices_stay(shape, &tri);
            assert_consistent(shape, &tri, &what);
            assert_point_queries(shape, &tri, &mut rng);
        }
    }
}

#[test]
fn models_and_bad_meshes() {
    for shape in models().iter().chain(bad_meshes().iter()) {
        let mut tri = shape.triangulation();
        let report = tri
            .optimal_delaunay(&options(TargetLength::MeanInputEdge(1.0), 8))
            .unwrap_or_else(|e| panic!("{}: {e}", shape.name));
        eprintln!(
            "{}: {:?} after {} iterations, {} inserted, last relative motion {:.2e}, \
             angles {:.2}..{:.2} degrees, area cv {:.3}",
            shape.name,
            report.status,
            report.n_iterations,
            report.n_inserted,
            report.last_max_relative_motion,
            report.statistics.angles.min.to_degrees(),
            report.statistics.angles.max.to_degrees(),
            report.statistics.area_cv
        );
        assert_eq!(report.delaunay.n_not_delaunay, 0, "{}", shape.name);
        assert_input_vertices_stay(shape, &tri);
        let check = tri.check().unwrap();
        assert!(check.relative_area_error.abs() < 1e-11, "{}: {check:?}", shape.name);
        assert!(check.max_angle_sum_error < 1e-8, "{}: {check:?}", shape.name);
        assert!(check.gauss_bonnet_residual.abs() < 1e-8, "{}: {check:?}", shape.name);
        let traces = tri.trace_edges(&Default::default()).unwrap();
        assert!(traces.report.max_relative_mismatch < 1e-8, "{}", shape.name);
    }
}

/// The run is deterministic, and the iterations of a longer run are those of a shorter one.
#[test]
fn the_iteration_is_deterministic() {
    let shape = &standard_shapes()[1]; // the jittered icosphere
    let run = |iterations| {
        let mut tri = shape.triangulation();
        let report = tri
            .optimal_delaunay(&options(TargetLength::MeanInputEdge(0.8), iterations))
            .unwrap();
        (tri, report)
    };
    let (a, ra) = run(6);
    let (b, rb) = run(6);
    assert_eq!(a, b);
    assert_eq!(ra, rb);
    let (c, rc) = run(3);
    assert_eq!(rc.iterations[..], ra.iterations[..3]);
    assert_ne!(a, c);
}

/// On a flat square the iteration converges, and it evens out the areas of the triangles
/// (paper Sec. 4.3: "a better distribution of areas").
#[test]
fn on_a_jittered_grid_the_areas_even_out() {
    let mut g = generators::grid_disk(12, 12).unwrap();
    let mut rng = SplitMix64::new(3);
    let h = 1.0 / 12.0;
    for (i, p) in g.positions.iter_mut().enumerate() {
        let (x, y) = (i / 13, i % 13);
        if x > 0 && x < 12 && y > 0 && y < 12 {
            p[0] += 0.3 * h * rng.uniform(-1.0, 1.0);
            p[1] += 0.3 * h * rng.uniform(-1.0, 1.0);
        }
    }
    let shape = Shape::generated("jittered grid", g, true);
    let mut idt = shape.triangulation();
    idt.flip_to_delaunay(&DelaunayOptions::default()).unwrap();
    let mut refined = shape.triangulation();
    refined.refine(&RefineOptions::default()).unwrap();
    let mut tri = shape.triangulation();
    let report = tri
        .optimal_delaunay(&options(TargetLength::Absolute(0.6 * h), 100))
        .unwrap();
    assert_eq!(report.status, OptimalDelaunayStatus::Converged, "{report:?}");
    assert!(report.last_max_relative_motion <= 1e-3);
    let (a, b, c) = (
        idt.shape_statistics(),
        refined.shape_statistics(),
        report.statistics,
    );
    eprintln!(
        "area cv: iDT {:.3}, iDR {:.3}, iODT {:.3}; min angle {:.1} {:.1} {:.1}",
        a.area_cv,
        b.area_cv,
        c.area_cv,
        a.angles.min.to_degrees(),
        b.angles.min.to_degrees(),
        c.angles.min.to_degrees()
    );
    assert!(c.area_cv < a.area_cv && c.area_cv < b.area_cv);
    assert_consistent(&shape, &tri, "iODT on a jittered grid");
}

#[test]
fn boundary_vertices_stay_or_slide() {
    let shape = Shape::generated("grid", generators::grid_disk(6, 4).unwrap(), true);
    // Splits of the boundary edges make the boundary vertices.
    let mut base = shape.triangulation();
    base.optimal_delaunay(&options(TargetLength::MeanInputEdge(0.6), 1))
        .unwrap();
    let boundary: Vec<usize> = (shape.mesh.n_vertices()..base.n_vertices())
        .filter(|&v| base.is_boundary_vertex(v))
        .collect();
    assert!(!boundary.is_empty());
    for (mode, centers) in [
        (BoundaryVertices::Fixed, OdtCenters::BarycentersAtBoundary),
        (BoundaryVertices::AlongBoundary, OdtCenters::BarycentersAtBoundary),
        (BoundaryVertices::AlongBoundary, OdtCenters::Circumcenters),
    ] {
        let mut tri = base.clone();
        let report = tri
            .optimal_delaunay(&OptimalDelaunayOptions {
                boundary: mode,
                centers,
                ..options(TargetLength::NoSplits, 10)
            })
            .unwrap();
        assert_eq!(report.n_inserted, 0);
        let moved = boundary
            .iter()
            .filter(|&&v| tri.vertex_location(v) != base.vertex_location(v))
            .count();
        match mode {
            BoundaryVertices::Fixed => assert_eq!(moved, 0),
            BoundaryVertices::AlongBoundary => assert!(moved > 0, "{centers:?}"),
        }
        for &v in &boundary {
            assert!(tri.is_boundary_vertex(v));
        }
        assert_consistent(&shape, &tri, &format!("{mode:?}, {centers:?}"));
    }
}

#[test]
fn limits_and_options() {
    let shape = Shape::generated("ico", generators::icosphere(1).unwrap(), false);
    // Nothing inserted and no splits: nothing to move, converged at once.
    let mut tri = shape.triangulation();
    let report = tri
        .optimal_delaunay(&options(TargetLength::NoSplits, 10))
        .unwrap();
    assert_eq!(report.status, OptimalDelaunayStatus::Converged);
    assert_eq!((report.n_iterations, report.n_inserted), (1, 0));
    // The insertion limit.
    let mut tri = shape.triangulation();
    let report = tri
        .optimal_delaunay(&OptimalDelaunayOptions {
            max_insertions: InsertionLimit::Vertices(5),
            ..options(TargetLength::MeanInputEdge(0.3), 10)
        })
        .unwrap();
    assert_eq!(report.status, OptimalDelaunayStatus::InsertionLimit);
    assert_eq!(report.n_inserted, 5);
    // The iteration limit.
    let mut tri = shape.triangulation();
    let report = tri
        .optimal_delaunay(&options(TargetLength::MeanInputEdge(0.3), 2))
        .unwrap();
    assert_eq!(report.status, OptimalDelaunayStatus::IterationLimit);
    assert_eq!(report.n_iterations, 2);
    // Invalid options change nothing.
    let mut tri = shape.triangulation();
    let fresh = tri.clone();
    for bad in [
        options(TargetLength::Absolute(-1.0), 3),
        options(TargetLength::MeanInputEdge(f64::NAN), 3),
        OptimalDelaunayOptions {
            step: 0.0,
            ..OptimalDelaunayOptions::default()
        },
        OptimalDelaunayOptions {
            step: 1.5,
            ..OptimalDelaunayOptions::default()
        },
        OptimalDelaunayOptions {
            motion_tolerance: f64::NAN,
            ..OptimalDelaunayOptions::default()
        },
    ] {
        assert!(matches!(
            tri.optimal_delaunay(&bad),
            Err(IntrinsicError::InvalidOptions { .. })
        ));
        assert_eq!(tri, fresh);
    }
}
