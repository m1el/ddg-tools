//! The invariants of SPEC §5 (items 1 to 4 and 9) on real meshes: the models of `data/models/`
//! (spot, bob, blub) and the deliberately bad meshes of the optimal-cones data set
//! (`camelbad`, `camelbadnonuniform`, `camel_noise1`, `camel_noise2`, `camelshotnoise`), and
//! how often an edge trace fails its check on them.
//!
//! The data is not part of the repository; a file that is absent is skipped with a notice.
//! Run with `cargo test -p ddg-intrinsic --test models -- --nocapture` for the tables.

mod common;

use common::{Shape, bad_meshes, models, path_length};
use ddg_core::generators::SplitMix64;
use ddg_intrinsic::{
    DelaunayCriterion, DelaunayOptions, DelaunayStatus, EdgeTraceOptions, IntrinsicError,
    IntrinsicTriangulation,
};

/// Traces all edges with the default tolerance (`1e-6` of the edge length) and returns, for
/// the edges that were traced, `mismatch / length`; the number of failures comes first.
fn trace_all(
    shape: &Shape,
    tri: &IntrinsicTriangulation,
    what: &str,
    path_tolerance: f64,
) -> (usize, Vec<f64>) {
    match tri.trace_edges(&EdgeTraceOptions::default()) {
        Ok(traces) => {
            let mut relative = Vec::new();
            for e in 0..tri.n_edges() {
                if tri.edge_on_input(e).is_some() {
                    continue;
                }
                let [a, b] = tri.edge_vertices(e);
                let length = tri.edge_length(e);
                relative.push(traces.mismatch(e) / length);
                // The traced path, measured in space, has the length of the edge.
                let in_space = path_length(
                    shape,
                    tri.vertex_location(a),
                    traces.crossings(e),
                    tri.vertex_location(b),
                );
                assert!(
                    (in_space - length).abs() <= path_tolerance * length,
                    "{} ({what}): edge {e} = ({a}, {b}) has length {length}, its path {in_space}",
                    shape.name
                );
            }
            (0, relative)
        }
        Err(IntrinsicError::EdgeTraces { failures, .. }) => {
            for failure in failures.iter().take(5) {
                println!("  {} ({what}): {failure}", shape.name);
            }
            (failures.len(), Vec::new())
        }
        Err(other) => panic!("{} ({what}): {other}", shape.name),
    }
}

fn above(values: &[f64], threshold: f64) -> usize {
    values.iter().filter(|&&x| x > threshold).count()
}

fn run(shapes: &[Shape], angle_tolerance: f64) {
    println!(
        "{:<20} {:>7} {:>6} {:>5} {:>6} {:>9} {:>9} {:>7} {:>6} {:>9} {:>6} {:>6} {:>6}",
        "mesh",
        "edges",
        "flips",
        "self",
        "multi",
        "area err",
        "angle err",
        "traced",
        "failed",
        "max mism.",
        ">1e-12",
        ">1e-10",
        ">1e-8"
    );
    for shape in shapes {
        let mut tri = shape.triangulation();
        let report = tri.flip_to_delaunay(&DelaunayOptions::default()).unwrap();
        let context = format!("{}: {report:?}", shape.name);
        // Item 4.
        assert_eq!(report.status, DelaunayStatus::Delaunay, "{context}");
        assert_eq!(
            (report.n_refused, report.n_not_delaunay),
            (0, 0),
            "{context}"
        );
        assert!(report.max_angle_excess <= 1e-10, "{context}");
        assert!(tri.is_delaunay(&DelaunayCriterion::default()), "{context}");
        // Items 1 and 2.
        let consistency = tri
            .check()
            .unwrap_or_else(|e| panic!("{}: {e}", shape.name));
        let context = format!("{context}\n{consistency:?}");
        assert!(consistency.min_triangle_slack > 0.0, "{context}");
        assert!(consistency.relative_area_error.abs() < 1e-11, "{context}");
        assert!(
            consistency.max_angle_sum_error < angle_tolerance,
            "{context}"
        );
        assert!(
            consistency.max_signpost_error < angle_tolerance,
            "{context}"
        );
        assert!(consistency.gauss_bonnet_residual.abs() < 1e-9, "{context}");
        // Item 3, and item 9: no failed trace at the default tolerance.
        // Measured: the path of an edge has its length to 2e-12 relative at worst.
        let (failed, relative) = trace_all(shape, &tri, "Delaunay", 1e-10);
        let irregular = tri.irregularities();
        let largest = relative.iter().fold(0.0f64, |a, &b| a.max(b));
        println!(
            "{:<20} {:>7} {:>6} {:>5} {:>6} {:>9.1e} {:>9.1e} {:>7} {:>6} {:>9.1e} {:>6} {:>6} {:>6}",
            shape.name,
            tri.n_edges(),
            report.n_flips,
            irregular.self_edges.len(),
            irregular.multi_edges.len(),
            consistency.relative_area_error.abs(),
            consistency.max_angle_sum_error,
            relative.len() + failed,
            failed,
            largest,
            above(&relative, 1e-12),
            above(&relative, 1e-10),
            above(&relative, 1e-8),
        );
        assert_eq!(failed, 0, "{context}");
        assert!(largest < 1e-9, "{context}: largest mismatch {largest}");
        // No NaN.
        assert!(tri.edge_lengths().iter().all(|l| l.is_finite() && *l > 0.0));
        assert!(tri.signposts().iter().all(|p| p.is_finite()));
        assert!(tri.cotan_laplacian().values().iter().all(|v| v.is_finite()));
    }
}

#[test]
fn models_flip_to_delaunay_and_every_edge_traces() {
    // Angle sums and signposts to 1e-12 rad on meshes of good shape.
    run(&models(), 1e-12);
}

#[test]
fn bad_meshes_flip_to_delaunay_and_every_edge_traces() {
    // Their thinnest triangles have angles of a few 1e-5 rad; an angle of such a triangle
    // follows from its lengths to about 1e-11 rad only.
    run(&bad_meshes(), 1e-9);
}

/// The irregular elements the reference program finds on the same inputs
/// (`data/intrinsic-reference/SUMMARY.md`): self-edges and vertex pairs with several edges.
#[test]
fn bad_meshes_have_the_irregular_elements_of_the_reference() {
    for (name, self_edges, pairs) in [
        ("camelbad", 0, 0),
        ("camel_noise1", 1, 5),
        ("camelshotnoise", 22, 24),
    ] {
        let Some(shape) =
            common::load_data(&format!("optimal-cones/optimal-cones-data/data/{name}.obj"))
        else {
            continue;
        };
        let mut tri = shape.triangulation();
        tri.flip_to_delaunay(&DelaunayOptions::default()).unwrap();
        let irregular = tri.irregularities();
        let mut multi: Vec<[usize; 2]> = irregular
            .multi_edges
            .iter()
            .map(|&e| {
                let [a, b] = tri.edge_vertices(e as usize);
                [a.min(b), a.max(b)]
            })
            .collect();
        multi.sort_unstable();
        multi.dedup();
        assert_eq!(
            (irregular.self_edges.len(), multi.len()),
            (self_edges, pairs),
            "{name}"
        );
        assert_eq!(irregular.degree_one_vertices.len(), self_edges, "{name}");
    }
}

/// A harder case for the traces than Delaunay flipping: edges flipped at random, which makes
/// long thin triangles and edges that cross many input faces.
#[test]
fn after_random_flips_every_edge_still_traces() {
    println!(
        "{:<20} {:>7} {:>7} {:>7} {:>6} {:>9} {:>9} {:>6} {:>6}",
        "mesh", "flips", "traced", "cross.", "failed", "max mism.", "angle err", ">1e-10", ">1e-8"
    );
    let mut shapes = models();
    shapes.extend(bad_meshes());
    for (s, shape) in shapes.iter().enumerate() {
        let mut tri = shape.triangulation();
        let mut rng = SplitMix64::new(4000 + s as u64);
        let mut flips = 0;
        for _ in 0..tri.n_edges() {
            let e = (rng.next_u64() % tri.n_edges() as u64) as usize;
            flips += usize::from(tri.flip_edge(e).is_ok());
        }
        let consistency = tri
            .check()
            .unwrap_or_else(|e| panic!("{}: {e}", shape.name));
        let (failed, relative) = trace_all(shape, &tri, "random flips", 1e-8);
        let crossings = tri
            .trace_edges(&EdgeTraceOptions::default())
            .map_or(0, |t| t.report.n_crossings);
        let largest = relative.iter().fold(0.0f64, |a, &b| a.max(b));
        println!(
            "{:<20} {:>7} {:>7} {:>7} {:>6} {:>9.1e} {:>9.1e} {:>6} {:>6}",
            shape.name,
            flips,
            relative.len() + failed,
            crossings,
            failed,
            largest,
            consistency.max_angle_sum_error,
            above(&relative, 1e-10),
            above(&relative, 1e-8),
        );
        assert_eq!(failed, 0, "{}", shape.name);
        // Each of these flips may move an angle by 1e-10 (SPEC N5).
        assert!(
            consistency.max_angle_sum_error < 1e-7,
            "{}: {consistency:?}",
            shape.name
        );
        assert!(
            consistency.relative_area_error.abs() < 1e-8,
            "{}",
            shape.name
        );
    }
}
