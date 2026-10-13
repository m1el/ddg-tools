//! SPEC §5 item 7: after refinement the angle bound holds everywhere except at the corners the
//! report lists, each with its reason; and item 8 for refinement: counts and angle statistics
//! comparable to the authors' program on the same inputs.
//!
//! Run with `cargo test --release -p ddg-intrinsic --test refinement -- --nocapture` for the
//! tables.

mod common;

use std::collections::HashMap;
use std::f64::consts::{FRAC_PI_2, PI};
use std::path::PathBuf;
use std::time::Instant;

use common::{Shape, dist, needle_tetrahedron, sheared_parallelogram, standard_shapes};
use ddg_core::distance::upwind_distance;
use ddg_intrinsic::{
    BelowReason, DelaunayCriterion, DelaunayOptions, EdgeTraceOptions, InsertionLimit,
    IntrinsicError, IntrinsicTriangulation, RefineOptions, RefineReport, RefineStatus,
    SurfacePoint,
};

fn options(degrees: f64) -> RefineOptions {
    RefineOptions {
        min_angle: degrees.to_radians(),
        ..RefineOptions::default()
    }
}

/// Item 7 and what else a refined triangulation must be.
fn assert_refined(
    shape: &Shape,
    tri: &IntrinsicTriangulation,
    report: &RefineReport,
    options: &RefineOptions,
) {
    let bound = options.min_angle;
    let context = format!(
        "{} at {:.0} degrees: {:?}, {} inserted",
        shape.name,
        bound.to_degrees(),
        report.status,
        report.n_inserted
    );
    // The counts.
    assert_eq!(
        report.n_inserted,
        report.n_face_insertions + report.n_edge_splits + report.n_boundary_splits,
        "{context}"
    );
    assert_eq!(
        tri.n_vertices(),
        shape.mesh.n_vertices() + report.n_inserted,
        "{context}"
    );
    assert_eq!(
        tri.euler_characteristic(),
        shape.mesh.euler_characteristic(),
        "{context}"
    );
    // Every corner below the bound is listed, with a reason that holds.
    let listed: HashMap<(u32, u32), &BelowReason> = report
        .below_bound
        .iter()
        .map(|b| ((b.face, b.vertex), &b.reason))
        .collect();
    let (mut smallest, mut largest, mut reachable) = (f64::INFINITY, 0.0f64, f64::INFINITY);
    for f in 0..tri.n_faces() {
        let mut face_listed = false;
        let mut face_largest = 0.0f64;
        for k in 0..3 {
            let vertex = tri.faces()[f][k];
            let angle = tri.corner_angle(3 * f + k);
            smallest = smallest.min(angle);
            largest = largest.max(angle);
            face_largest = face_largest.max(angle);
            let unreachable = tri.vertex_angle_sum(vertex as usize) < bound;
            if !unreachable {
                reachable = reachable.min(angle);
            }
            if angle >= bound {
                continue;
            }
            face_listed = true;
            let reason = listed
                .get(&(f as u32, vertex))
                .unwrap_or_else(|| panic!("{context}: the corner of face {f} at vertex {vertex} has {angle} and is not reported"));
            match reason {
                BelowReason::AngleSum { angle_sum } => {
                    assert!(unreachable, "{context}: face {f}, vertex {vertex}");
                    assert_eq!(
                        *angle_sum,
                        tri.vertex_angle_sum(vertex as usize),
                        "{context}"
                    );
                    // An input vertex: inserted vertices are flat.
                    assert!((vertex as usize) < shape.mesh.n_vertices(), "{context}");
                }
                BelowReason::NeedleVertex {
                    vertex: needle,
                    angle_sum,
                } => {
                    assert!(*angle_sum < PI / 3.0 && *angle_sum >= bound, "{context}");
                    assert!(!tri.is_boundary_vertex(*needle as usize), "{context}");
                }
                BelowReason::NotInsertable { detail } => {
                    assert!(!detail.is_empty(), "{context}");
                }
                BelowReason::InsertionLimit => {
                    assert_eq!(report.status, RefineStatus::InsertionLimit, "{context}");
                }
            }
        }
        // A triangle all of whose angles meet the bound has no angle above π − 2 · bound.
        if !face_listed {
            assert!(
                face_largest <= PI - 2.0 * bound + 1e-12,
                "{context}: face {f}"
            );
        }
    }
    assert_eq!(
        report.below_bound.len(),
        listed.len(),
        "{context}: a corner is listed twice"
    );
    for b in &report.below_bound {
        assert!(b.angle < bound, "{context}: {b:?}");
        assert_eq!(
            b.angle,
            tri.corner_angle(
                3 * b.face as usize
                    + (0..3)
                        .find(|&k| tri.faces()[b.face as usize][k] == b.vertex)
                        .unwrap()
            )
            .min(b.angle),
            "{context}: {b:?}"
        );
    }
    assert_eq!(report.min_angle, smallest, "{context}");
    assert_eq!(report.max_angle, largest, "{context}");
    assert_eq!(report.min_angle_reachable, reachable, "{context}");
    if report
        .below_bound
        .iter()
        .all(|b| matches!(b.reason, BelowReason::AngleSum { .. }))
    {
        assert!(report.min_angle_reachable >= bound, "{context}");
    }
    // Delaunay, and no boundary edge with an obtuse angle opposite: every cotan weight is
    // non-negative, the boundary edges included.
    assert_eq!(report.n_not_delaunay, 0, "{context}");
    assert!(tri.is_delaunay(&DelaunayCriterion::default()), "{context}");
    if report.status == RefineStatus::Done {
        assert_eq!(report.n_encroached, 0, "{context}");
        for e in 0..tri.n_edges() {
            let weight = tri.edge_cotan_weight(e);
            assert!(
                weight >= -1e-8,
                "{context}: edge {e} has cotan weight {weight}"
            );
            if tri.is_boundary_edge(e) {
                let opposite = tri.corner_angle(tri.prev(tri.edge_halfedge(e)));
                assert!(opposite <= FRAC_PI_2 + 1e-9, "{context}: boundary edge {e}");
            }
        }
    }
    // The same surface, and every vertex where it says.
    let consistency = tri.check().unwrap_or_else(|e| panic!("{context}: {e}"));
    assert!(
        consistency.relative_area_error.abs() < 1e-11,
        "{context}: {consistency:?}"
    );
    assert!(
        consistency.max_angle_sum_error < 1e-9,
        "{context}: {consistency:?}"
    );
    assert!(
        consistency.max_signpost_error < 1e-9,
        "{context}: {consistency:?}"
    );
    for v in 0..shape.mesh.n_vertices() {
        assert_eq!(
            tri.vertex_location(v),
            SurfacePoint::Vertex(v as u32),
            "{context}"
        );
    }
    let traces = tri
        .trace_edges(&EdgeTraceOptions::default())
        .unwrap_or_else(|e| panic!("{context}: {e}"));
    assert!(
        traces.report.max_relative_mismatch < 1e-8,
        "{context}: {:?}",
        traces.report
    );
}

#[test]
fn the_angle_bound_holds_except_where_reported() {
    for shape in standard_shapes() {
        for degrees in [25.0, 30.0] {
            let options = options(degrees);
            let mut tri = shape.triangulation();
            let report = tri.refine(&options).unwrap();
            assert_eq!(
                report.status,
                RefineStatus::Done,
                "{}: {degrees}",
                shape.name
            );
            assert_refined(&shape, &tri, &report, &options);
            // Nothing but angle sums stands in the way on these shapes.
            for b in &report.below_bound {
                assert!(
                    matches!(b.reason, BelowReason::AngleSum { .. }),
                    "{}: {degrees}: {b:?}",
                    shape.name
                );
            }
            // Refinement is deterministic.
            let mut again = shape.triangulation();
            assert_eq!(again.refine(&options).unwrap(), report, "{}", shape.name);
            assert_eq!(again, tri, "{}", shape.name);
            assert_eq!(again.fingerprint(), tri.fingerprint());
            // A second run has nothing to do.
            let second = tri.refine(&options).unwrap();
            assert_eq!(
                (second.n_inserted, second.n_flips),
                (0, 0),
                "{}",
                shape.name
            );
        }
    }
}

/// SPEC N1 on a shape whose boundary has corners of 9.5°, 18.4° and obtuse triangles along it.
#[test]
fn boundary_edges_are_split_where_they_are_encroached() {
    for (n, sigma, corner) in [(10, 6.0, 9.462f64), (8, 3.0, 18.435)] {
        let shape = sheared_parallelogram(n, sigma);
        let mut tri = shape.triangulation();
        let options = options(25.0);
        // Before: boundary edges with an obtuse angle opposite, i.e. a negative cotan weight,
        // which flipping cannot remove.
        tri.flip_to_delaunay(&DelaunayOptions::default()).unwrap();
        let negative = (0..tri.n_edges())
            .filter(|&e| tri.is_boundary_edge(e) && tri.edge_cotan_weight(e) < -1e-9)
            .count();
        assert!(negative > 0, "{}", shape.name);
        let report = tri.refine(&options).unwrap();
        assert_eq!(report.status, RefineStatus::Done);
        assert_refined(&shape, &tri, &report, &options);
        assert!(report.n_boundary_splits > 0 && report.n_face_insertions > 0);
        // The two acute corners of the parallelogram are all that is left below the bound.
        assert_eq!(report.below_bound.len(), 2, "{:?}", report.below_bound);
        for b in &report.below_bound {
            assert!((b.angle.to_degrees() - corner).abs() < 1e-3, "{b:?}");
            assert!(matches!(b.reason, BelowReason::AngleSum { .. }));
            // Its triangle is the corner cut off by two boundary edges of equal length (the
            // boundary edges at a small angle are split at the same distances).
            let h = (0..3)
                .map(|k| 3 * b.face as usize + k)
                .find(|&h| tri.from_vertex(h) == b.vertex as usize)
                .unwrap();
            let (leaving, arriving) = (tri.edge(h), tri.edge(tri.prev(h)));
            assert!(tri.is_boundary_edge(leaving) && tri.is_boundary_edge(arriving));
            let (a, c) = (tri.edge_length(leaving), tri.edge_length(arriving));
            assert!((a - c).abs() < 1e-12 * a, "{a} and {c}");
        }
        // Vertices inserted on the boundary lie on input boundary edges.
        let mut on_boundary = 0;
        for v in shape.mesh.n_vertices()..tri.n_vertices() {
            if tri.is_boundary_vertex(v) {
                on_boundary += 1;
                match tri.vertex_location(v) {
                    SurfacePoint::Edge { edge, .. } => {
                        assert!(shape.mesh.is_boundary_edge(edge as usize))
                    }
                    other => panic!("{other}"),
                }
            }
        }
        assert_eq!(on_boundary, report.n_boundary_splits);
        // With an angle bound of zero only the boundary rule acts.
        let mut only_boundary = shape.triangulation();
        let zero = options_with(0.0, InsertionLimit::PerVertex(10));
        let report = only_boundary.refine(&zero).unwrap();
        assert_eq!(report.status, RefineStatus::Done);
        assert_eq!((report.n_face_insertions, report.n_edge_splits), (0, 0));
        assert!(report.n_boundary_splits > 0 && report.below_bound.is_empty());
        assert_refined(&shape, &only_boundary, &report, &zero);
    }
}

fn options_with(degrees: f64, max_insertions: InsertionLimit) -> RefineOptions {
    RefineOptions {
        max_insertions,
        ..options(degrees)
    }
}

/// SPEC N4: vertices with a small angle sum. A tetrahedron over a base of radius 1 whose apex
/// has the given angle sum, for sums below the bound, between the bound and 60°, and above.
#[test]
fn corners_at_vertices_with_a_small_angle_sum_are_reported() {
    println!(
        "{:>9} {:>4} {:>9} {:>9} {:>9} {:>9} {:>6} {:>10}",
        "apex sum", "deg", "inserted", "min", "min rch", "max", "below", "apex deg."
    );
    for apex_degrees in [
        2.0f64, 5.0, 10.0, 20.0, 24.0, 26.0, 29.0, 31.0, 40.0, 50.0, 59.0, 61.0, 90.0,
    ] {
        // Each of the three faces at the apex has a third of the sum there.
        let half = (apex_degrees / 6.0).to_radians();
        let side = 0.75f64.sqrt() / half.sin();
        let shape = needle_tetrahedron((side * side - 1.0).sqrt(), 0.0);
        let apex = 3;
        for degrees in [25.0, 30.0] {
            let options = options(degrees);
            let mut tri = shape.triangulation();
            let angle_sum = tri.vertex_angle_sum(apex);
            assert!((angle_sum.to_degrees() - apex_degrees).abs() < 1e-9);
            let report = tri.refine(&options).unwrap();
            println!(
                "{:>9.1} {:>4.0} {:>9} {:>9.3} {:>9.3} {:>9.3} {:>6} {:>10}",
                apex_degrees,
                degrees,
                report.n_inserted,
                report.min_angle.to_degrees(),
                report.min_angle_reachable.to_degrees(),
                report.max_angle.to_degrees(),
                report.below_bound.len(),
                tri.vertex_degree(apex),
            );
            assert_eq!(report.status, RefineStatus::Done, "{apex_degrees}");
            assert_refined(&shape, &tri, &report, &options);
            if apex_degrees < degrees {
                // What is left: the corners at the apex, all of them, and nothing else.
                assert!(!report.below_bound.is_empty());
                let mut at_apex = 0.0;
                for b in &report.below_bound {
                    assert_eq!(b.vertex as usize, apex, "{b:?}");
                    assert_eq!(b.reason, BelowReason::AngleSum { angle_sum });
                    at_apex += b.angle;
                }
                assert!((at_apex - angle_sum).abs() < 1e-12);
                assert!(report.min_angle < options.min_angle);
            } else {
                // The bound is reached: around the apex too.
                assert!(report.below_bound.is_empty(), "{:?}", report.below_bound);
                assert!(report.min_angle >= options.min_angle);
            }
            assert!(report.min_angle_reachable >= options.min_angle);
        }
    }
    // An angle sum equal to the bound: the one triangle around the apex has the corner of the
    // bound, to rounding above it or below. Below, it is reported: as unreachable if the sum
    // came out below the bound too, otherwise as the loop a circumcentre would make there.
    for degrees in [20.0f64, 25.0, 30.0] {
        let half = (degrees / 6.0).to_radians();
        let side = 0.75f64.sqrt() / half.sin();
        let shape = needle_tetrahedron((side * side - 1.0).sqrt(), 0.0);
        let options = options(degrees);
        let mut tri = shape.triangulation();
        let report = tri.refine(&options).unwrap();
        assert_eq!(report.status, RefineStatus::Done);
        assert_refined(&shape, &tri, &report, &options);
        assert!(report.min_angle > options.min_angle - 1e-12);
        for b in &report.below_bound {
            assert_eq!(b.vertex, 3, "{b:?}");
            match &b.reason {
                BelowReason::AngleSum { .. } => {}
                BelowReason::NeedleVertex { vertex, angle_sum } => {
                    assert_eq!(*vertex, 3);
                    assert_eq!(*angle_sum, tri.vertex_angle_sum(3));
                }
                other => panic!("{other:?}"),
            }
        }
        println!(
            "apex sum = bound = {degrees}: {} inserted, {:?}",
            report.n_inserted, report.below_bound
        );
    }
}

#[test]
fn limits_and_options() {
    let shape = sheared_parallelogram(8, 3.0);
    // The insertion limit is reported, and what it left is listed.
    let limited = options_with(30.0, InsertionLimit::Vertices(12));
    let mut tri = shape.triangulation();
    let report = tri.refine(&limited).unwrap();
    assert_eq!(report.status, RefineStatus::InsertionLimit);
    assert_eq!(report.n_inserted, 12);
    assert!(
        report
            .below_bound
            .iter()
            .any(|b| b.reason == BelowReason::InsertionLimit)
            || report.n_encroached > 0
    );
    assert_refined(&shape, &tri, &report, &limited);
    // Carrying on reaches the bound.
    let unlimited = options_with(30.0, InsertionLimit::Unlimited);
    let rest = tri.refine(&unlimited).unwrap();
    assert_eq!(rest.status, RefineStatus::Done);
    assert_refined(
        &shape,
        &tri,
        &RefineReport {
            n_inserted: 12 + rest.n_inserted,
            n_face_insertions: report.n_face_insertions + rest.n_face_insertions,
            n_edge_splits: report.n_edge_splits + rest.n_edge_splits,
            n_boundary_splits: report.n_boundary_splits + rest.n_boundary_splits,
            ..rest.clone()
        },
        &unlimited,
    );
    // A bound on the circumradius: more vertices, smaller triangles.
    let radius = 0.04;
    let sized = RefineOptions {
        max_circumradius: Some(radius),
        ..options(25.0)
    };
    let mut fine = shape.triangulation();
    let fine_report = fine.refine(&sized).unwrap();
    let mut coarse = shape.triangulation();
    let coarse_report = coarse.refine(&options(25.0)).unwrap();
    assert_eq!(fine_report.status, RefineStatus::Done);
    assert!(fine_report.too_large.is_empty());
    assert!(fine_report.n_inserted > 2 * coarse_report.n_inserted);
    for f in 0..fine.n_faces() {
        let [a, b, c] = fine.face_lengths(f);
        assert!(a * b * c / (4.0 * fine.face_area(f)) <= radius * (1.0 + 1e-12));
    }
    assert_refined(&shape, &fine, &fine_report, &sized);
    // Invalid options.
    for bad in [
        options(-1.0),
        options(36.0),
        options(f64::NAN),
        RefineOptions {
            max_circumradius: Some(0.0),
            ..options(25.0)
        },
        RefineOptions {
            criterion: DelaunayCriterion::AngleSum { tolerance: -1.0 },
            ..options(25.0)
        },
    ] {
        let mut tri = shape.triangulation();
        let error = tri.refine(&bad).unwrap_err();
        assert!(
            matches!(error, IntrinsicError::InvalidOptions { .. }),
            "{error}"
        );
        assert_eq!(tri, shape.triangulation());
    }
    assert_eq!(RefineOptions::default().min_angle, 25f64.to_radians());
}

/// SPEC §5 item 5, second half: the distances of the first-order solver after refinement.
/// Sheared lattices, exact distance known; errors at the input vertices in units of the
/// input's lattice width, on the input mesh, on the intrinsic Delaunay triangulation, after
/// refinement to 25° and to 30°, and after refinement to 25° with the circumradius bounded by
/// half a lattice width.
#[test]
fn distances_after_refinement() {
    println!(
        "{:<34} {:>5} | {:>6} {:>6} | {:>6} {:>6} | {:>5} {:>6} {:>6} | {:>5} {:>6} {:>6} | {:>5} {:>6} {:>6}",
        "shape",
        "seed",
        "input",
        "max",
        "Del.",
        "max",
        "+V",
        "25",
        "max",
        "+V",
        "30",
        "max",
        "+V",
        "sized",
        "max"
    );
    for (n, sigma) in [(16, 1.5), (16, 2.5), (16, 3.7), (24, 3.7), (16, 6.3)] {
        let shape = sheared_parallelogram(n, sigma);
        let h = 1.0 / n as f64;
        for seed in [n, n * (n + 1) / 2 + n / 2] {
            let exact = |p: [f64; 3]| dist(p, shape.positions[seed]);
            let errors = |mesh: &ddg_core::TriMesh, geometry: &ddg_core::Geometry| {
                let d = upwind_distance(mesh, geometry, &[(seed as u32, 0.0)], None).unwrap();
                let mut sum = 0.0;
                let mut max = 0.0f64;
                for v in 0..shape.mesh.n_vertices() {
                    let e = (d.values[v] - exact(shape.positions[v])).abs();
                    sum += e;
                    max = max.max(e);
                }
                (sum / shape.mesh.n_vertices() as f64 / h, max / h)
            };
            let on = |tri: &IntrinsicTriangulation| {
                let (mesh, geometry) = tri.to_tri_mesh().unwrap();
                errors(&mesh, &geometry)
            };
            let input = errors(&shape.mesh, &shape.geometry);
            let mut delaunay = shape.triangulation();
            delaunay
                .flip_to_delaunay(&DelaunayOptions::default())
                .unwrap();
            let flipped = on(&delaunay);
            let mut refined = Vec::new();
            for options in [
                options(25.0),
                options(30.0),
                RefineOptions {
                    max_circumradius: Some(0.5 * h),
                    max_insertions: InsertionLimit::Unlimited,
                    ..options(25.0)
                },
            ] {
                let mut tri = shape.triangulation();
                let report = tri.refine(&options).unwrap();
                assert_eq!(report.status, RefineStatus::Done);
                refined.push((report.n_inserted, on(&tri)));
            }
            println!(
                "{:<34} {:>5} | {:>6.3} {:>6.3} | {:>6.3} {:>6.3} | {:>5} {:>6.3} {:>6.3} | {:>5} {:>6.3} {:>6.3} | {:>5} {:>6.3} {:>6.3}",
                shape.name,
                seed,
                input.0,
                input.1,
                flipped.0,
                flipped.1,
                refined[0].0,
                refined[0].1.0,
                refined[0].1.1,
                refined[1].0,
                refined[1].1.0,
                refined[1].1.1,
                refined[2].0,
                refined[2].1.0,
                refined[2].1.1,
            );
            // What holds: far better than on the input mesh, and of the size of the error on
            // the Delaunay triangulation. Not better than that one: the solver is of first
            // order in the mesh width, which an angle bound does not reduce, and refinement
            // leaves obtuse angles (up to 180° minus twice the bound), which the Delaunay
            // triangulation of a lattice does not have (SPEC N16).
            for (_, (mean, max)) in &refined {
                assert!(
                    *mean < input.0 / 3.0 && *max < input.1 / 3.0,
                    "{}",
                    shape.name
                );
                assert!(
                    *mean < 2.5 * flipped.0 && *max < 2.5 * flipped.1,
                    "{}",
                    shape.name
                );
            }
        }
    }
}

// ----- against the reference -----

fn reference_dir() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../data/intrinsic-reference")
}

/// The number after `"key":` in a flat JSON text.
fn json_number(text: &str, key: &str) -> f64 {
    let at = text
        .find(&format!("\"{key}\":"))
        .unwrap_or_else(|| panic!("no {key}"))
        + key.len()
        + 3;
    text[at..]
        .trim_start()
        .split([',', '\n', '}'])
        .next()
        .unwrap()
        .trim()
        .parse()
        .unwrap_or_else(|e| panic!("{key}: {e}"))
}

struct Row {
    name: &'static str,
    degrees: f64,
    reference_inserted: usize,
    reference_min: f64,
    reference_max: f64,
    report: RefineReport,
    seconds: f64,
    n_input: usize,
}

/// Counts and angle statistics against the reference's refinement of the same inputs
/// (`data/intrinsic-reference/<input>/refine25|refine30`). The rules differ where there is a
/// boundary or a vertex with a small angle sum (SPEC N1, N4), so the counts are comparable,
/// not equal: within 3 % on the closed surfaces without such vertices.
#[test]
fn refinement_is_comparable_to_the_reference_s() {
    const INPUTS: [&str; 15] = [
        "spot",
        "bob",
        "blub",
        "bunny",
        "camel",
        "camelbad",
        "camelbadnonuniform",
        "camel_noise1",
        "camelshotnoise",
        "bump",
        "icosphere",
        "icosphere_jitter",
        "thin_cylinder",
        "long_torus",
        "sheared_grid",
    ];
    let root = reference_dir();
    if !root.is_dir() {
        eprintln!(
            "SKIPPED: reference directory {} not found; nothing was checked",
            root.display()
        );
        return;
    }
    // One job per input and bound, in parallel.
    let rows: Vec<Option<Row>> = std::thread::scope(|scope| {
        let handles: Vec<_> = INPUTS
            .iter()
            .flat_map(|&name| [25.0, 30.0].map(|degrees| (name, degrees)))
            .map(|(name, degrees)| {
                let root = root.clone();
                scope.spawn(move || {
                    let dir = root.join(name);
                    let mut shape =
                        common::load_data(&format!("intrinsic-reference/{name}/input.obj"))?;
                    shape.name = name.to_string();
                    let analysis = std::fs::read_to_string(
                        dir.join(format!("refine{degrees:.0}"))
                            .join("analysis.json"),
                    )
                    .unwrap();
                    let options = options(degrees);
                    let mut tri = shape.triangulation();
                    let start = Instant::now();
                    let report = tri.refine(&options).unwrap();
                    let seconds = start.elapsed().as_secs_f64();
                    assert_refined(&shape, &tri, &report, &options);
                    Some(Row {
                        name,
                        degrees,
                        reference_inserted: json_number(&analysis, "inserted_vertices_net")
                            as usize,
                        reference_min: json_number(&analysis, "min_angle_deg"),
                        reference_max: json_number(&analysis, "max_angle_deg"),
                        report,
                        seconds,
                        n_input: shape.mesh.n_vertices(),
                    })
                })
            })
            .collect();
        handles.into_iter().map(|h| h.join().unwrap()).collect()
    });
    println!(
        "{:<20} {:>3} {:>7} | {:>8} {:>8} | {:>7} {:>7} {:>7} | {:>7} {:>7} | {:>5} {:>6} | {:>8}",
        "input",
        "deg",
        "V in",
        "inserted",
        "(ref)",
        "min",
        "min rch",
        "(ref)",
        "max",
        "(ref)",
        "below",
        "bdry",
        "ms"
    );
    for row in rows.iter().flatten() {
        let r = &row.report;
        println!(
            "{:<20} {:>3.0} {:>7} | {:>8} {:>8} | {:>7.3} {:>7.3} {:>7.3} | {:>7.2} {:>7.2} | {:>5} {:>6} | {:>8.1}",
            row.name,
            row.degrees,
            row.n_input,
            r.n_inserted,
            row.reference_inserted,
            r.min_angle.to_degrees(),
            r.min_angle_reachable.to_degrees(),
            row.reference_min,
            r.max_angle.to_degrees(),
            row.reference_max,
            r.below_bound.len(),
            r.n_boundary_splits,
            1e3 * row.seconds,
        );
        assert_eq!(r.status, RefineStatus::Done, "{}", row.name);
        // The bound is reached wherever a triangulation can reach it.
        assert!(
            r.min_angle_reachable >= row.degrees.to_radians(),
            "{}",
            row.name
        );
        for b in &r.below_bound {
            assert!(
                matches!(b.reason, BelowReason::AngleSum { .. }),
                "{}: {b:?}",
                row.name
            );
        }
        // Never worse than the reference in the smallest angle.
        assert!(
            r.min_angle.to_degrees() >= row.reference_min - 1e-6,
            "{}: {} against {}",
            row.name,
            r.min_angle.to_degrees(),
            row.reference_min
        );
        // Comparable counts: equal up to 3 % on closed surfaces without small angle sums.
        let closed =
            !["bump", "thin_cylinder", "sheared_grid", "camelshotnoise"].contains(&row.name);
        if closed && row.reference_inserted > 50 {
            let ratio = r.n_inserted as f64 / row.reference_inserted as f64;
            assert!(
                (0.97..=1.03).contains(&ratio),
                "{}: ratio {ratio}",
                row.name
            );
        }
    }
}
