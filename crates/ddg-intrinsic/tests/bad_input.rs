//! SPEC §5 item 9: bad input and failed traces are errors with a message that names the
//! element; no NaN anywhere.
//!
//! Invalid tracing queries are in `tracing.rs`, refused flips in `flip.rs`, invalid Delaunay
//! options in `delaunay.rs`; traces that fail because the data is wrong are unit tests of the
//! crate (`edge_trace.rs`), since no triangulation built through the public interface is wrong.

mod common;

use common::{Shape, long_torus, needle_tetrahedron, sheared_parallelogram, standard_shapes};
use ddg_core::generators;
use ddg_core::{Geometry, GeometryError, TriMesh};
use ddg_intrinsic::{
    DelaunayOptions, EdgeTraceFailureReason, EdgeTraceOptions, IntrinsicError,
    IntrinsicTriangulation, SurfacePoint,
};

#[test]
fn a_geometry_of_another_mesh_is_refused() {
    let small = Shape::generated("grid", generators::grid_disk(2, 2).unwrap(), true);
    let large = Shape::generated("grid", generators::grid_disk(3, 2).unwrap(), true);
    let error = IntrinsicTriangulation::new(&small.mesh, &large.geometry).unwrap_err();
    assert_eq!(
        error,
        IntrinsicError::Geometry(GeometryError::MeshMismatch {
            mesh: (small.mesh.n_edges(), small.mesh.n_faces()),
            geometry: (large.mesh.n_edges(), large.mesh.n_faces()),
        })
    );
    assert!(
        error
            .to_string()
            .starts_with("geometry does not belong to this mesh")
    );
    // Degenerate input never gets this far: ddg_core refuses it with the face named.
    let mesh = TriMesh::new(3, vec![[0, 1, 2]]).unwrap();
    let flat = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [2.0, 0.0, 0.0]];
    assert!(matches!(
        Geometry::from_positions(&mesh, &flat),
        Err(GeometryError::DegenerateFaces { .. })
    ));
}

/// With a tolerance of zero the rounding error of a trace is a failure: the error lists
/// every edge concerned, with its end points and how far the trace was off.
#[test]
fn a_trace_that_misses_its_vertex_is_an_error_that_names_the_edge() {
    let shape = sheared_parallelogram(6, 2.0);
    let mut tri = shape.triangulation();
    tri.flip_to_delaunay(&DelaunayOptions::default()).unwrap();
    let good = tri.trace_edges(&EdgeTraceOptions::default()).unwrap();
    assert!(good.report.n_traced > 50 && good.report.max_relative_mismatch > 0.0);
    let strict = EdgeTraceOptions {
        tolerance: 0.0,
        max_crossings: None,
    };
    let error = tri.trace_edges(&strict).unwrap_err();
    let IntrinsicError::EdgeTraces { n_edges, failures } = &error else {
        panic!("{error}");
    };
    assert_eq!(*n_edges, tri.n_edges());
    // Exactly the edges whose trace does not land on the vertex to the last bit.
    let inexact: Vec<u32> = (0..tri.n_edges())
        .filter(|&e| good.mismatch(e) > 0.0)
        .map(|e| e as u32)
        .collect();
    assert!(!inexact.is_empty());
    assert_eq!(failures.iter().map(|f| f.edge).collect::<Vec<_>>(), inexact);
    for failure in failures {
        let e = failure.edge as usize;
        assert_eq!(failure.vertices.map(|v| v as usize), tri.edge_vertices(e));
        assert_eq!(
            (failure.length, failure.tolerance),
            (tri.edge_length(e), 0.0)
        );
        match &failure.reason {
            EdgeTraceFailureReason::Missed {
                mismatch: Some(d), ..
            } => {
                assert_eq!(*d, good.mismatch(e));
                assert!(*d < 1e-13 * failure.length);
            }
            // An edge that ends at a boundary vertex: a path that passes the vertex by the
            // last bit leaves the surface there.
            EdgeTraceFailureReason::HitBoundary { at, travelled } => {
                assert!(tri.is_boundary_vertex(failure.vertices[1] as usize));
                assert!(matches!(at, SurfacePoint::Edge { .. }));
                assert!((travelled - failure.length).abs() < 1e-13 * failure.length);
            }
            other => panic!("{other:?}"),
        }
        // One edge alone gives the same failure.
        assert_eq!(
            tri.trace_edge(e, &strict).unwrap_err(),
            IntrinsicError::EdgeTrace(failure.clone())
        );
    }
    let message = error.to_string();
    let first = &failures[0];
    assert!(
        message.starts_with(&format!(
            "{} of {} intrinsic edges could not be traced over the input mesh: the trace of \
             intrinsic edge {} = ({}, {}) of length {:?} over the input mesh ",
            failures.len(),
            tri.n_edges(),
            first.edge,
            first.vertices[0],
            first.vertices[1],
            first.length
        )),
        "{message}"
    );
    assert!(message.contains("the tolerance is 0.0"), "{message}");
    assert!(
        message.ends_with(&format!("... ({} more)", failures.len() - 4)),
        "{message}"
    );
}

#[test]
fn invalid_trace_options_and_edges_are_errors() {
    let shape = sheared_parallelogram(4, 1.0);
    let mut tri = shape.triangulation();
    tri.flip_to_delaunay(&DelaunayOptions::default()).unwrap();
    for tolerance in [-1e-9, f64::NAN, f64::INFINITY] {
        let options = EdgeTraceOptions {
            tolerance,
            max_crossings: None,
        };
        for error in [
            tri.trace_edges(&options).map(|_| ()).unwrap_err(),
            tri.trace_edge(0, &options).map(|_| ()).unwrap_err(),
        ] {
            assert!(
                matches!(error, IntrinsicError::InvalidOptions { .. }),
                "{error}"
            );
            assert!(
                error.to_string().contains("tolerance of an edge trace"),
                "{error}"
            );
        }
    }
    let n = tri.n_edges();
    assert_eq!(
        tri.trace_edge(n, &EdgeTraceOptions::default()).unwrap_err(),
        IntrinsicError::IndexOutOfRange {
            what: "edge",
            index: n,
            bound: n
        }
    );
}

#[test]
fn a_triangulation_that_is_not_simplicial_says_why() {
    let torus = long_torus(8, 3, 4.0, 0.2, 0.37);
    let mut tri = torus.triangulation();
    assert!(tri.to_tri_mesh().is_ok());
    tri.flip_to_delaunay(&DelaunayOptions::default()).unwrap();
    let irregular = tri.irregularities();
    let error = tri.to_tri_mesh().unwrap_err();
    let e = irregular.self_edges[0];
    let [a, _] = tri.edge_vertices(e as usize);
    assert_eq!(
        error.to_string(),
        format!(
            "the intrinsic triangulation is not a simplicial complex: edge {e} joins vertex {a} \
             to itself (6 self-edges in total)"
        )
    );
    // Without self-edges the first double edge is named.
    let needle = needle_tetrahedron(8.0, 4.0);
    let mut tri = needle.triangulation();
    tri.flip_edge(needle.mesh.find_edge(1, 3).unwrap()).unwrap();
    let irregular = tri.irregularities();
    assert!(irregular.self_edges.is_empty() && irregular.multi_edges.len() == 2);
    let message = tri.to_tri_mesh().unwrap_err().to_string();
    let e = irregular.multi_edges[0];
    let [a, b] = tri.edge_vertices(e as usize);
    assert_eq!(
        message,
        format!(
            "the intrinsic triangulation is not a simplicial complex: edge {e} = ({a}, {b}) is \
             one of several edges between its end points (2 such edges in total)"
        )
    );
}

/// Every number an operation hands out is finite, on the badly shaped inputs too.
#[test]
fn no_nan_anywhere() {
    let finite = |x: f64, what: &str, name: &str| assert!(x.is_finite(), "{name}: {what} is {x}");
    for shape in standard_shapes() {
        let name = &shape.name;
        let mut tri = shape.triangulation();
        let report = tri.flip_to_delaunay(&DelaunayOptions::default()).unwrap();
        finite(report.max_angle_excess, "max_angle_excess", name);
        finite(report.min_cotan_weight, "min_cotan_weight", name);
        for &l in tri.edge_lengths() {
            finite(l, "an edge length", name);
        }
        for &phi in tri.signposts() {
            finite(phi, "a signpost", name);
        }
        for h in 0..tri.n_halfedges() {
            finite(tri.corner_angle(h), "a corner angle", name);
        }
        for e in 0..tri.n_edges() {
            finite(tri.edge_cotan_weight(e), "a cotan weight", name);
        }
        for &v in tri.cotan_laplacian().values() {
            finite(v, "a Laplacian entry", name);
        }
        for m in tri.lumped_mass() {
            assert!(m.is_finite() && m > 0.0, "{name}: a lumped area is {m}");
        }
        let q = tri.quality(0.0);
        for x in [
            q.min_angle,
            q.max_angle,
            q.min_cotan_weight,
            q.negative_weight_fraction(),
        ] {
            finite(x, "a quality number", name);
        }
        let c = tri.check().unwrap();
        for x in [
            c.max_angle_sum_error,
            c.max_signpost_error,
            c.relative_area_error,
            c.gauss_bonnet_residual,
            c.min_triangle_slack,
        ] {
            finite(x, "a consistency number", name);
        }
        let traces = tri.trace_edges(&EdgeTraceOptions::default()).unwrap();
        finite(
            traces.report.max_relative_mismatch,
            "the largest mismatch",
            name,
        );
        for e in 0..tri.n_edges() {
            finite(traces.mismatch(e), "a mismatch", name);
            for c in traces.crossings(e) {
                finite(c.t, "a crossing parameter", name);
                finite(c.distance, "a crossing distance", name);
            }
        }
    }
}
