//! SPEC §6: the simplicial export, the way back for a cut, and what identifies a state. And
//! the end-to-end run the export exists for: a badly shaped mesh is refined, exported, the cut
//! flow of `ddg-variational-cuts` runs on the export, and its cut is carried back to the input
//! mesh.
//!
//! Run with `cargo test --release -p ddg-intrinsic --test export -- --nocapture` for the
//! tables.

mod common;

use std::collections::HashMap;

use common::{Shape, dist, needles_and_caps, position, sheared_square, standard_shapes};
use ddg_core::{Geometry, TriMesh};
use ddg_intrinsic::{
    DelaunayCriterion, DelaunayOptions, EdgeTraceOptions, ExportOptions, IntrinsicError,
    IntrinsicTriangulation, PathsOnInput, RefineOptions, RefineStatus, Segment, SimplicialExport,
    SurfacePoint,
};
use ddg_variational_cuts::{Cut, CutProblem, Energy, Flow, FlowOptions, StepSize};

fn refine_options(degrees: f64) -> RefineOptions {
    RefineOptions {
        min_angle: degrees.to_radians(),
        ..RefineOptions::default()
    }
}

/// What an export must be, whatever the triangulation was.
fn assert_export(shape: &Shape, tri: &IntrinsicTriangulation, export: &SimplicialExport) {
    let context = &shape.name;
    // A simplicial complex: `TriMesh` accepted it, and the triangulation has nothing irregular.
    assert!(tri.irregularities().is_empty(), "{context}");
    assert_eq!(export.mesh.n_vertices(), tri.n_vertices(), "{context}");
    assert_eq!(export.mesh.faces(), tri.faces(), "{context}");
    assert_eq!(export.mesh.n_edges(), tri.n_edges(), "{context}");
    assert_eq!(
        export.mesh.euler_characteristic(),
        shape.mesh.euler_characteristic(),
        "{context}"
    );
    assert_eq!(
        export.mesh.boundary_loops().len(),
        shape.mesh.boundary_loops().len(),
        "{context}"
    );
    export.geometry.check_mesh(&export.mesh).unwrap();
    // The same surface.
    let area = shape.geometry.total_area();
    assert!(
        (export.geometry.total_area() - area).abs() < 1e-11 * area,
        "{context}"
    );
    // The edges of the mesh are those of the triangulation, with their lengths.
    let mut seen = vec![false; tri.n_edges()];
    for e in 0..export.mesh.n_edges() {
        let te = export.triangulation_edge(e);
        assert!(
            !std::mem::replace(&mut seen[te], true),
            "{context}: edge {te} twice"
        );
        let [a, b] = export.mesh.edge_vertices(e);
        let [tail, head] = tri.edge_vertices(te);
        assert_eq!([a, b], [tail.min(head), tail.max(head)], "{context}");
        assert_eq!(
            export.geometry.edge_lengths()[e],
            tri.edge_length(te),
            "{context}"
        );
        // A point of the mesh edge, measured from its lower vertex, as a point of the
        // triangulation's edge, measured from its tail.
        let expected = if tail < head { 0.25 } else { 0.75 };
        assert_eq!(
            export.edge_point(e, 0.25),
            SurfacePoint::Edge {
                edge: te as u32,
                t: expected
            },
            "{context}"
        );
    }
    // The angles it reports are its own.
    let (mut smallest, mut largest) = (f64::INFINITY, 0.0f64);
    for h in 0..tri.n_halfedges() {
        smallest = smallest.min(tri.corner_angle(h));
        largest = largest.max(tri.corner_angle(h));
    }
    assert!((export.min_angle - smallest).abs() < 1e-12, "{context}");
    assert!((export.max_angle - largest).abs() < 1e-12, "{context}");
    // The places of the vertices, and the state.
    assert_eq!(export.locations, tri.vertex_locations(), "{context}");
    for v in 0..shape.mesh.n_vertices() {
        assert_eq!(
            export.locations[v],
            SurfacePoint::Vertex(v as u32),
            "{context}"
        );
    }
    assert_eq!(export.fingerprint, tri.fingerprint(), "{context}");
    export
        .check_state(export.fingerprint, export.mesh.n_vertices())
        .unwrap();
    // Still a consistent triangulation whose edges can be traced.
    let consistency = tri.check().unwrap_or_else(|e| panic!("{context}: {e}"));
    assert!(
        consistency.max_angle_sum_error < 1e-9,
        "{context}: {consistency:?}"
    );
    tri.trace_edges(&EdgeTraceOptions::default())
        .unwrap_or_else(|e| panic!("{context}: {e}"));
}

#[test]
fn an_export_is_a_simplicial_mesh_of_the_same_surface() {
    println!(
        "{:<36} {:<10} {:>6} {:>6} {:>6} {:>7} {:>7} {:>6}",
        "shape", "state", "self", "multi", "faces", "splits", "flips", "V"
    );
    let mut n_irregular = 0;
    for shape in standard_shapes() {
        for state in ["as given", "Delaunay", "refined"] {
            let mut tri = shape.triangulation();
            match state {
                "Delaunay" => {
                    tri.flip_to_delaunay(&DelaunayOptions::default()).unwrap();
                }
                "refined" => {
                    tri.refine(&RefineOptions::default()).unwrap();
                }
                _ => {}
            }
            let before = tri.irregularities();
            let unchanged = tri.clone();
            // A triangulation as it is, flips or none, is exported without flips; then nothing
            // changes unless something has to be split.
            let options = ExportOptions {
                delaunay: (state != "as given").then(DelaunayOptions::default),
                ..ExportOptions::default()
            };
            let export = tri.export_simplicial(&options).unwrap();
            println!(
                "{:<36} {:<10} {:>6} {:>6} {:>6} {:>7} {:>7} {:>6}",
                shape.name,
                state,
                before.self_edges.len(),
                before.multi_edges.len(),
                before.faces_with_repeated_vertex.len(),
                export.n_splits,
                export.n_flips,
                export.mesh.n_vertices(),
            );
            assert_export(&shape, &tri, &export);
            // Splits exactly where the triangulation was not a simplicial complex, one new
            // vertex each.
            assert_eq!(export.n_splits > 0, !before.is_empty(), "{}", shape.name);
            assert_eq!(
                tri.n_vertices(),
                unchanged.n_vertices() + export.n_splits,
                "{}",
                shape.name
            );
            n_irregular += usize::from(!before.is_empty());
            assert_eq!(
                export.is_delaunay,
                tri.is_delaunay(&DelaunayCriterion::default()),
                "{}",
                shape.name
            );
            if state == "as given" {
                // The input mesh is a simplicial complex: it is handed out as it is.
                assert!(before.is_empty(), "{}", shape.name);
                assert_eq!(tri, unchanged, "{}", shape.name);
                assert_eq!(export.mesh.faces(), shape.mesh.faces());
                assert_eq!(
                    export.geometry.edge_lengths(),
                    shape.geometry.edge_lengths()
                );
                // With the default options it is flipped to Delaunay first.
                let mut flipped = shape.triangulation();
                let delaunay = flipped
                    .export_simplicial(&ExportOptions::default())
                    .unwrap();
                assert!(delaunay.is_delaunay, "{}", shape.name);
                assert_export(&shape, &flipped, &delaunay);
            } else {
                assert!(export.is_delaunay, "{}", shape.name);
                if before.is_empty() {
                    assert_eq!((export.n_splits, export.n_flips), (0, 0));
                    assert_eq!(tri, unchanged, "{}", shape.name);
                }
            }
            // Exporting again changes nothing.
            let again = tri.export_simplicial(&options).unwrap();
            assert_eq!((again.n_splits, again.n_flips), (0, 0), "{}", shape.name);
            assert_eq!(again.fingerprint, export.fingerprint, "{}", shape.name);
            assert_eq!(again.mesh.faces(), export.mesh.faces());
        }
    }
    assert!(
        n_irregular >= 4,
        "the shapes must include irregular triangulations"
    );
}

/// The smallest angle at a vertex where the bound can be reached (its angle sum is not below
/// the bound), and the number of such corners below the bound.
fn reachable_angles(tri: &IntrinsicTriangulation, bound: f64) -> (f64, usize) {
    let (mut reachable, mut below) = (f64::INFINITY, 0usize);
    for h in 0..tri.n_halfedges() {
        if tri.vertex_angle_sum(tri.from_vertex(h)) >= bound {
            let angle = tri.corner_angle(h);
            reachable = reachable.min(angle);
            below += usize::from(angle < bound);
        }
    }
    (reachable, below)
}

/// What the export does to the angles of a refined triangulation (SPEC N15): nothing where
/// nothing is split; where edges are split (thin tubes and handles, spikes), the result is
/// Delaunay but the angle bound of the refinement is lost. With refinement between the rounds
/// of splits the bound comes back on tubes and handles; at spikes it cannot. One job per mesh
/// and bound, in parallel; `--nocapture` prints the table.
#[test]
fn the_export_keeps_the_angles_only_where_nothing_is_split() {
    let mut shapes = standard_shapes();
    shapes.extend(common::models());
    shapes.extend(common::bad_meshes());
    struct Row {
        line: String,
        n_splits: usize,
    }
    let rows: Vec<Row> = std::thread::scope(|scope| {
        let handles: Vec<_> = shapes
            .iter()
            .flat_map(|shape| [25.0f64, 30.0].map(|degrees| (shape, degrees)))
            .map(|(shape, degrees)| {
                scope.spawn(move || {
                    let bound = degrees.to_radians();
                    let options = refine_options(degrees);
                    let mut tri = shape.triangulation();
                    let report = tri.refine(&options).unwrap();
                    assert_eq!(report.status, RefineStatus::Done, "{}", shape.name);
                    let irregular = tri.irregularities();
                    let refined = tri.clone();
                    // Splits and flips only.
                    let export = tri.export_simplicial(&ExportOptions::default()).unwrap();
                    assert!(export.is_delaunay, "{}", shape.name);
                    assert_eq!((export.n_refine_rounds, export.n_refined), (0, 0));
                    assert_eq!(export.refinement, None);
                    let (reachable, below) = reachable_angles(&tri, bound);
                    if export.n_splits == 0 {
                        // Handed out as it was: the report of the refinement holds for it.
                        assert_eq!(tri, refined, "{}", shape.name);
                        assert_eq!(export.min_angle, report.min_angle, "{}", shape.name);
                        assert_eq!(export.max_angle, report.max_angle, "{}", shape.name);
                        assert_eq!(below, 0, "{}", shape.name);
                    }
                    // With refinement between the rounds of splits.
                    let mut again = refined.clone();
                    let start = std::time::Instant::now();
                    let with = again
                        .export_simplicial(&ExportOptions {
                            refine: Some(options),
                            ..ExportOptions::default()
                        })
                        .unwrap();
                    let seconds = start.elapsed().as_secs_f64();
                    assert_export(shape, &again, &with);
                    assert!(with.is_delaunay, "{}", shape.name);
                    let (reachable_with, below_with) = reachable_angles(&again, bound);
                    match &with.refinement {
                        // It ended with a refinement: the bound holds as that one reports.
                        Some(last) => {
                            assert_eq!(last.status, RefineStatus::Done, "{}", shape.name);
                            assert_eq!(below_with, 0, "{}", shape.name);
                            assert_eq!(with.min_angle, last.min_angle, "{}", shape.name);
                        }
                        None => assert!(with.n_refine_rounds <= 6, "{}", shape.name),
                    }
                    if export.n_splits == 0 {
                        // Nothing to split: one refinement that finds nothing to do.
                        assert_eq!((with.n_refine_rounds, with.n_refined), (1, 0));
                        assert_eq!(again, refined, "{}", shape.name);
                    }
                    // Thin tubes and handles: the bound comes back.
                    if ["torus", "long torus", "tube"]
                        .iter()
                        .any(|kind| shape.name.starts_with(kind))
                    {
                        assert!(with.refinement.is_some(), "{}", shape.name);
                    }
                    Row {
                        line: format!(
                            "{:<38} {:>3.0} {:>6} {:>6} | {:>4} {:>5} | {:>7.3} {:>7.2} | {:>6} {:>5} {:>7.3} {:>7.2} {:>5} | {:>2} {:>3} {:>6} {:>6} {:>7.3} {:>7.2} {:>5} {:>7.1}",
                            shape.name,
                            degrees,
                            shape.mesh.n_vertices(),
                            report.n_inserted,
                            irregular.self_edges.len(),
                            irregular.multi_edges.len(),
                            report.min_angle_reachable.to_degrees(),
                            report.max_angle.to_degrees(),
                            export.n_splits,
                            export.n_flips,
                            reachable.to_degrees(),
                            export.max_angle.to_degrees(),
                            below,
                            with.n_refine_rounds,
                            if with.refinement.is_some() { "yes" } else { "no" },
                            with.n_splits,
                            with.n_refined,
                            reachable_with.to_degrees(),
                            with.max_angle.to_degrees(),
                            below_with,
                            1e3 * seconds,
                        ),
                        n_splits: export.n_splits,
                    }
                })
            })
            .collect();
        handles.into_iter().map(|h| h.join().unwrap()).collect()
    });
    println!(
        "refined: self-edges, edges in groups between two vertices, smallest reachable and largest angle;          export by splits and flips: splits, flips, angles, corners below the bound;          export with refinement between the rounds: rounds, ended refined, splits, inserted, angles, below, ms"
    );
    println!(
        "{:<38} {:>3} {:>6} {:>6} | {:>4} {:>5} | {:>7} {:>7} | {:>6} {:>5} {:>7} {:>7} {:>5} | {:>2} {:>3} {:>6} {:>6} {:>7} {:>7} {:>5} {:>7}",
        "mesh",
        "deg",
        "V",
        "+V",
        "self",
        "multi",
        "min",
        "max",
        "splits",
        "flips",
        "min",
        "max",
        "below",
        "r",
        "ref",
        "splits",
        "+V",
        "min",
        "max",
        "below",
        "ms"
    );
    let n_unsplit = rows.iter().filter(|r| r.n_splits == 0).count();
    for row in rows.iter().filter(|r| r.n_splits > 0) {
        println!("{}", row.line);
    }
    println!(
        "{n_unsplit} of {} runs needed no split and kept their angles",
        rows.len()
    );
    assert!(n_unsplit > rows.len() / 2);
}

#[test]
fn export_options_and_refusals() {
    let shape = common::long_torus(6, 4, 3.0, 0.2, 0.0);
    let mut tri = shape.triangulation();
    tri.flip_to_delaunay(&DelaunayOptions::default()).unwrap();
    let irregular = tri.irregularities();
    assert!(
        !irregular.is_empty(),
        "this torus must have double edges after flipping"
    );
    // A limit on the splits: an error that says so, and nothing was changed.
    let before = tri.clone();
    let error = tri
        .export_simplicial(&ExportOptions {
            max_splits: 0,
            ..ExportOptions::default()
        })
        .unwrap_err();
    assert!(
        matches!(error, IntrinsicError::ExportFailed { n_splits: 0, .. }),
        "{error}"
    );
    assert!(error.to_string().contains("more than 0 splits"), "{error}");
    assert_eq!(tri, before);
    // Without flips: simplicial, with the splits only.
    let mut plain = tri.clone();
    let export = plain
        .export_simplicial(&ExportOptions {
            delaunay: None,
            ..ExportOptions::default()
        })
        .unwrap();
    assert!(export.n_splits > 0 && export.n_flips == 0);
    assert_export(&shape, &plain, &export);
    assert_eq!(
        export.is_delaunay,
        plain.is_delaunay(&DelaunayCriterion::default()),
        "without flips the flag is the test's verdict"
    );
    // With flips: Delaunay as well.
    let export = tri.export_simplicial(&ExportOptions::default()).unwrap();
    assert!(export.n_splits > 0 && export.is_delaunay);
    assert_export(&shape, &tri, &export);
}

#[test]
fn the_export_is_reproducible_and_a_stale_state_is_refused() {
    for shape in [
        common::long_torus(6, 4, 3.0, 0.2, 0.0),
        sheared_square(12, 2.0),
        needles_and_caps(6),
    ] {
        let run = || {
            let mut tri = shape.triangulation();
            tri.refine(&RefineOptions::default()).unwrap();
            let export = tri.export_simplicial(&ExportOptions::default()).unwrap();
            (tri, export)
        };
        let (tri, export) = run();
        // The same input and options: the same export, bit for bit.
        let (tri_again, again) = run();
        assert_eq!(again, export, "{}", shape.name);
        assert_eq!(tri_again, tri, "{}", shape.name);
        again
            .check_state(export.fingerprint, export.mesh.n_vertices())
            .unwrap();
        // Another treatment (flips only): another export, and a state made for the first is
        // refused by it.
        let other = {
            let mut tri = shape.triangulation();
            tri.export_simplicial(&ExportOptions::default()).unwrap()
        };
        assert_ne!(other.fingerprint, export.fingerprint, "{}", shape.name);
        let error = other
            .check_state(export.fingerprint, export.mesh.n_vertices())
            .unwrap_err();
        assert!(
            matches!(error, IntrinsicError::StaleState { .. }),
            "{error}"
        );
        let message = error.to_string();
        assert!(
            message.contains(&format!("{:016x}", export.fingerprint))
                && message.contains(&format!("{:016x}", other.fingerprint)),
            "{message}"
        );
        // The right export, the wrong number of values.
        let error = export
            .check_state(export.fingerprint, export.mesh.n_vertices() + 1)
            .unwrap_err();
        assert!(
            matches!(error, IntrinsicError::IndexOutOfRange { .. }),
            "{error}"
        );
        // Any change of the triangulation changes its fingerprint: a flip (same vertices),
        // a split (same connectivity elsewhere).
        let mut flipped = tri.clone();
        let e = (0..flipped.n_edges())
            .find(|&e| flipped.can_flip(e).unwrap().is_ok())
            .unwrap();
        flipped.flip_edge(e).unwrap();
        assert_ne!(flipped.fingerprint(), export.fingerprint, "{}", shape.name);
        let mut split = tri.clone();
        split.split_edge(0, 0.5).unwrap();
        assert_ne!(split.fingerprint(), export.fingerprint, "{}", shape.name);
        // And the input mesh itself is part of it: the same connectivity scaled.
        let scaled: Vec<[f64; 3]> = shape.positions.iter().map(|p| p.map(|x| 1.5 * x)).collect();
        let bigger = Shape::new("scaled", scaled, shape.mesh.faces().to_vec(), shape.flat);
        assert_ne!(
            bigger.triangulation().fingerprint(),
            shape.triangulation().fingerprint()
        );
    }
}

// ----- the way back -----

/// Whether a point of the input mesh lies in the closed input face `f`.
fn in_input_face(shape: &Shape, point: SurfacePoint, f: usize) -> bool {
    match point {
        SurfacePoint::Vertex(v) => shape.mesh.faces()[f].contains(&v),
        SurfacePoint::Edge { edge, .. } => {
            (0..3).any(|k| shape.mesh.edge(3 * f + k) == edge as usize)
        }
        SurfacePoint::Face { face, .. } => face as usize == f,
    }
}

/// What carrying a cut back gave.
#[derive(Debug, Default)]
struct Back {
    n_nodes: usize,
    n_segments: usize,
    n_pieces: usize,
    /// Nodes at input vertices, on input edges, inside input faces.
    node_places: [usize; 3],
    /// Largest distance between a node's place on the input mesh and the same point found
    /// from the positions of the exported vertices (flat shapes).
    max_node_error: f64,
    /// Largest difference between the length of a segment on the export and the length of its
    /// pieces on the input mesh, relative to the segment.
    max_length_error: f64,
    max_mismatch: f64,
    /// The positions of the nodes on the input mesh.
    positions: Vec<[f64; 3]>,
}

/// Carries `cut`, found on `export`, back to the input mesh and checks what SPEC §6 promises:
/// every node located on the input mesh, every segment a chain of straight pieces that each
/// lie in one input face, the chains joined where the segments are.
fn carry_back(
    shape: &Shape,
    tri: &IntrinsicTriangulation,
    export: &SimplicialExport,
    cut: &Cut,
) -> (PathsOnInput, Back) {
    assert!(
        shape.flat,
        "the checks below compare positions in the plane"
    );
    // The cut as the triangulation is asked for it: nodes, then segments between them.
    let mut nodes: Vec<SurfacePoint> = cut
        .points
        .iter()
        .map(|p| export.edge_point(p.edge as usize, p.t))
        .collect();
    nodes.extend(cut.triple_points.iter().map(|t| SurfacePoint::Face {
        face: t.face,
        bary: t.barycentric,
    }));
    let segments: Vec<Segment> = cut
        .segments
        .iter()
        .map(|s| Segment {
            face: s.face,
            nodes: s.points,
        })
        .collect();
    let paths = tri
        .segments_on_input(&nodes, &segments, &EdgeTraceOptions::default())
        .unwrap_or_else(|e| panic!("{}: {e}", shape.name));
    assert_eq!(paths.nodes.len(), nodes.len());
    assert_eq!(paths.paths.len(), segments.len());

    // Where the nodes are, from the positions of the exported vertices: the export is a
    // triangulation of the same flat region, so its edges and faces are straight there.
    let vertex: Vec<[f64; 3]> = export
        .locations
        .iter()
        .map(|&l| position(shape, l))
        .collect();
    let mut expected: Vec<[f64; 3]> = cut
        .points
        .iter()
        .map(|p| {
            let [a, b] = export.mesh.edge_vertices(p.edge as usize);
            [0, 1, 2].map(|k| (1.0 - p.t) * vertex[a][k] + p.t * vertex[b][k])
        })
        .collect();
    expected.extend(cut.triple_points.iter().map(|t| {
        let f = export.mesh.faces()[t.face as usize];
        [0, 1, 2].map(|k| {
            (0..3)
                .map(|c| t.barycentric[c] * vertex[f[c] as usize][k])
                .sum()
        })
    }));
    let mut back = Back {
        n_nodes: nodes.len(),
        n_segments: segments.len(),
        max_mismatch: paths.max_relative_mismatch,
        ..Back::default()
    };
    for (i, &node) in paths.nodes.iter().enumerate() {
        let p = position(shape, node);
        back.max_node_error = back.max_node_error.max(dist(p, expected[i]));
        back.node_places[match node {
            SurfacePoint::Vertex(_) => 0,
            SurfacePoint::Edge { .. } => 1,
            SurfacePoint::Face { .. } => 2,
        }] += 1;
        back.positions.push(p);
    }
    for (segment, path) in segments.iter().zip(&paths.paths) {
        // From its first node to its second.
        assert_eq!(path.points[0], paths.nodes[segment.nodes[0] as usize]);
        assert_eq!(
            *path.points.last().unwrap(),
            paths.nodes[segment.nodes[1] as usize]
        );
        // Straight pieces, each in one input face.
        assert_eq!(path.faces.len() + 1, path.points.len());
        assert!(!path.faces.is_empty());
        let mut length = 0.0;
        for (k, &f) in path.faces.iter().enumerate() {
            let (a, b) = (path.points[k], path.points[k + 1]);
            assert!(
                in_input_face(shape, a, f as usize) && in_input_face(shape, b, f as usize),
                "{}: a piece from {a} to {b} is not in input face {f}",
                shape.name
            );
            length += dist(position(shape, a), position(shape, b));
        }
        // Between the pieces the path crosses an input edge (or passes an input vertex).
        for point in &path.points[1..path.points.len() - 1] {
            assert!(!matches!(point, SurfacePoint::Face { .. }), "{point}");
        }
        // The pieces make the straight segment: their lengths add up to its length.
        let straight = dist(
            expected[segment.nodes[0] as usize],
            expected[segment.nodes[1] as usize],
        );
        back.max_length_error = back
            .max_length_error
            .max((length - straight).abs() / straight);
        back.n_pieces += path.faces.len();
    }
    (paths, back)
}

/// How often every node is the end of a segment.
fn node_degrees(cut: &Cut) -> Vec<usize> {
    let mut degree = vec![0; cut.points.len() + cut.triple_points.len()];
    for s in &cut.segments {
        degree[s.points[0] as usize] += 1;
        degree[s.points[1] as usize] += 1;
    }
    degree
}

/// A cut with a triple point: three regions around a point of the shape, by three linear
/// functions.
fn three_region_cut(shape: &Shape, export: &SimplicialExport, center: [f64; 2]) -> Cut {
    let phis: Vec<Vec<f64>> = (0..3)
        .map(|k| {
            let angle = 0.3 + 2.0 * std::f64::consts::PI * k as f64 / 3.0;
            export
                .locations
                .iter()
                .map(|&l| {
                    let p = position(shape, l);
                    -((p[0] - center[0]) * angle.cos() + (p[1] - center[1]) * angle.sin())
                })
                .collect()
        })
        .collect();
    Cut::from_regions(&export.mesh, &phis, 1e-3).unwrap()
}

#[test]
fn a_cut_on_an_export_comes_back_as_pieces_on_input_faces() {
    println!(
        "{:<28} {:<9} {:>6} {:>6} {:>7} | {:>5} {:>5} {:>5} | {:>9} {:>9} {:>9}",
        "shape",
        "cut",
        "nodes",
        "segm.",
        "pieces",
        "at V",
        "on E",
        "in F",
        "node err",
        "length",
        "mismatch"
    );
    for shape in [
        sheared_square(16, 2.0),
        needles_and_caps(8),
        common::sheared_parallelogram(10, 6.0),
        common::jittered_grid(12, 0.3, 7),
    ] {
        let mut tri = shape.triangulation();
        tri.refine(&RefineOptions::default()).unwrap();
        let export = tri.export_simplicial(&ExportOptions::default()).unwrap();
        let center = {
            let n = shape.positions.len() as f64;
            let sum = shape
                .positions
                .iter()
                .fold([0.0; 2], |s, p| [s[0] + p[0], s[1] + p[1]]);
            [sum[0] / n + 0.011, sum[1] / n - 0.007]
        };
        let extent = shape
            .positions
            .iter()
            .map(|p| (p[1] - center[1]).abs())
            .fold(0.0, f64::max);
        // A closed curve: the circle of a level-set function.
        let phi: Vec<f64> = export
            .locations
            .iter()
            .map(|&l| {
                let p = position(&shape, l);
                (p[0] - center[0]).hypot(p[1] - center[1]) - 0.55 * extent
            })
            .collect();
        let circle = Cut::from_level_set(&export.mesh, &phi, 1e-3).unwrap();
        let network = three_region_cut(&shape, &export, center);
        assert_eq!(network.triple_points.len(), 1, "{}", shape.name);
        for (name, cut) in [("circle", &circle), ("network", &network)] {
            let (paths, back) = carry_back(&shape, &tri, &export, cut);
            println!(
                "{:<28} {:<9} {:>6} {:>6} {:>7} | {:>5} {:>5} {:>5} | {:>9.1e} {:>9.1e} {:>9.1e}",
                shape.name,
                name,
                back.n_nodes,
                back.n_segments,
                back.n_pieces,
                back.node_places[0],
                back.node_places[1],
                back.node_places[2],
                back.max_node_error,
                back.max_length_error,
                back.max_mismatch,
            );
            assert!(back.max_node_error < 1e-10, "{}: {back:?}", shape.name);
            assert!(back.max_length_error < 1e-9, "{}", shape.name);
            assert!(back.max_mismatch < 1e-9, "{}", shape.name);
            assert!(back.n_pieces >= back.n_segments);
            // The pieces are joined where the segments are: the path of one segment ends at
            // the very point where the next one starts (the nodes are located once).
            let degree = node_degrees(cut);
            let mut ends: HashMap<u32, usize> = HashMap::new();
            for (s, path) in cut.segments.iter().zip(&paths.paths) {
                for (node, point) in [
                    (s.points[0], path.points[0]),
                    (s.points[1], *path.points.last().unwrap()),
                ] {
                    assert_eq!(point, paths.nodes[node as usize]);
                    *ends.entry(node).or_insert(0) += 1;
                }
            }
            for (node, &d) in degree.iter().enumerate() {
                assert_eq!(ends.get(&(node as u32)).copied().unwrap_or(0), d);
            }
            // A curve ends only on the boundary of the input mesh.
            for (node, &d) in degree.iter().enumerate() {
                if d == 1 {
                    match paths.nodes[node] {
                        SurfacePoint::Edge { edge, .. } => {
                            assert!(shape.mesh.is_boundary_edge(edge as usize), "{}", shape.name)
                        }
                        other => panic!("{}: a curve ends at {other}", shape.name),
                    }
                } else if node < cut.points.len() {
                    assert_eq!(d, 2, "{}", shape.name);
                } else {
                    // The triple point joins three, and it lies inside an input face or on an
                    // input edge.
                    assert_eq!(d, 3, "{}", shape.name);
                    assert!(!matches!(paths.nodes[node], SurfacePoint::Vertex(_)));
                }
            }
            // Refinement put vertices inside input faces and on input edges, so nodes of the
            // cut lie inside input faces: what `ddg_uv::cutting` has to take.
            assert!(back.node_places[2] > 0, "{}", shape.name);
        }
        // A node that is not on its segment's face is refused, with its number.
        let mut nodes: Vec<SurfacePoint> = circle
            .points
            .iter()
            .map(|p| export.edge_point(p.edge as usize, p.t))
            .collect();
        let segments: Vec<Segment> = circle
            .segments
            .iter()
            .map(|s| Segment {
                face: s.face,
                nodes: s.points,
            })
            .collect();
        let first = segments[0];
        let foreign = (0..tri.n_edges())
            .find(|&e| (0..3).all(|k| tri.edge(3 * first.face as usize + k) != e))
            .unwrap();
        nodes[first.nodes[0] as usize] = SurfacePoint::Edge {
            edge: foreign as u32,
            t: 0.5,
        };
        let error = tri
            .segments_on_input(&nodes, &segments, &EdgeTraceOptions::default())
            .unwrap_err();
        assert!(
            matches!(error, IntrinsicError::InvalidPoint { .. }),
            "{error}"
        );
    }
}

// ----- the cut flow on an export -----

/// A mesh the cut flow can run on, with the positions of its vertices in the plane.
struct FlowMesh {
    label: String,
    mesh: TriMesh,
    geometry: Geometry,
    positions: Vec<[f64; 3]>,
}

/// The redistancing test of `ddg-variational-cuts/examples/bad_mesh_study.rs` (stage 2): start
/// from the exact signed distance `exact`, redistance, and compare. Returns the mean and the
/// largest `|φ − d|` over the vertices that are not next to the cut and have `|d| ≤ 0.25`,
/// and the number of vertices no front reached.
fn redistancing_error(m: &FlowMesh, exact: &dyn Fn([f64; 3]) -> f64) -> (f64, f64, usize) {
    let d: Vec<f64> = m.positions.iter().map(|&p| exact(p)).collect();
    let mut flow = Flow::new(&m.mesh, &m.geometry, &d).unwrap();
    let unreached = flow.redistance().unwrap();
    let phi = flow.level_set();
    let adjacent = Cut::from_level_set(&m.mesh, &d, 1e-3)
        .unwrap()
        .cut_adjacent_vertices(&m.mesh);
    let (mut sum, mut max, mut n) = (0.0, 0.0f64, 0usize);
    for v in 0..phi.len() {
        if adjacent[v] || d[v].abs() > 0.25 {
            continue;
        }
        let e = (phi[v] - d[v]).abs();
        sum += e;
        max = max.max(e);
        n += 1;
    }
    (sum / n.max(1) as f64, max, unreached)
}

/// What a run of the circle law gave.
struct Law {
    /// Relative error of the mean radius of the cut points against the exact radius.
    error: f64,
    /// Standard deviation of the radius over its mean.
    roundness: f64,
    clamped: usize,
    unreached: usize,
    phi: Vec<f64>,
}

const CIRCLE: (f64, f64, f64) = (0.3, 1.0, 0.03);

fn radius_statistics(points: &[[f64; 3]], center: [f64; 2]) -> (f64, f64) {
    let radii: Vec<f64> = points
        .iter()
        .map(|p| (p[0] - center[0]).hypot(p[1] - center[1]))
        .collect();
    let mean = radii.iter().sum::<f64>() / radii.len() as f64;
    let deviation =
        (radii.iter().map(|r| (r - mean).powi(2)).sum::<f64>() / radii.len() as f64).sqrt();
    (mean, deviation / mean)
}

/// The circle law of the same study (stage 3): on a flat shape of area 1 a circle of radius
/// `ρ₀ = 0.3` shrinks under the length term alone, `ρ' = −α/ρ`, to `√(ρ₀² − 2 α t)`.
fn circle_law(m: &FlowMesh, center: [f64; 2], paper: bool, steps: usize) -> Result<Law, String> {
    let (rho0, alpha, t_end) = CIRCLE;
    let exact = (rho0 * rho0 - 2.0 * alpha * t_end).sqrt();
    let phi0: Vec<f64> = m
        .positions
        .iter()
        .map(|p| (p[0] - center[0]).hypot(p[1] - center[1]) - rho0)
        .collect();
    let mut flow = Flow::new(&m.mesh, &m.geometry, &phi0).map_err(|e| e.to_string())?;
    let problem = CutProblem::new(Energy::Dirichlet, alpha);
    let step = StepSize::Fixed(t_end / steps as f64);
    let options = if paper {
        FlowOptions::paper(step)
    } else {
        FlowOptions::new(step)
    };
    let (mut clamped, mut unreached) = (0, 0);
    for _ in 0..steps {
        let report = flow.step(&problem, &options).map_err(|e| e.to_string())?;
        clamped = clamped.max(report.n_clamped_crossings);
        unreached = unreached.max(report.n_unreached);
    }
    let phi = flow.level_set();
    let cut = Cut::from_level_set(&m.mesh, &phi, 1e-3).map_err(|e| e.to_string())?;
    if cut.points.is_empty() {
        return Err("the cut vanished".to_string());
    }
    let points = cut
        .point_positions(&m.mesh, &m.positions)
        .map_err(|e| e.to_string())?;
    let (mean, roundness) = radius_statistics(&points, center);
    Ok(Law {
        error: mean / exact - 1.0,
        roundness,
        clamped,
        unreached,
        phi,
    })
}

/// A triangulation of the shape as a mesh for the flow.
fn flow_mesh(shape: &Shape, label: &str, export: &SimplicialExport) -> FlowMesh {
    FlowMesh {
        label: label.to_string(),
        mesh: export.mesh.clone(),
        geometry: export.geometry.clone(),
        positions: export
            .locations
            .iter()
            .map(|&l| position(shape, l))
            .collect(),
    }
}

/// One line of the end-to-end table, and what the test asserts on.
struct EndToEnd {
    line: String,
    refined: bool,
    law_default: Result<(f64, f64), String>,
    law_paper: Result<(f64, f64), String>,
    redistancing: [(f64, f64, usize); 2],
    back: Option<(Back, f64, f64)>,
}

/// The whole chain on one shape and one treatment: `None` runs the flow on the input mesh
/// itself, `Some(None)` on the export of the Delaunay-flipped triangulation, and
/// `Some(Some(bound))` on the export of the refined one.
fn end_to_end(
    shape: &Shape,
    center: [f64; 2],
    treatment: Option<Option<f64>>,
    steps: usize,
) -> EndToEnd {
    let mut tri = shape.triangulation();
    let (label, inserted, export) = match treatment {
        None => ("input mesh".to_string(), 0, None),
        Some(bound) => {
            let (label, inserted) = match bound {
                None => {
                    tri.flip_to_delaunay(&DelaunayOptions::default()).unwrap();
                    ("flipped".to_string(), 0)
                }
                Some(degrees) => {
                    let report = tri.refine(&refine_options(degrees)).unwrap();
                    assert_eq!(report.status, RefineStatus::Done, "{}", shape.name);
                    (format!("refined {degrees:.0}"), report.n_inserted)
                }
            };
            let export = tri.export_simplicial(&ExportOptions::default()).unwrap();
            assert_export(shape, &tri, &export);
            (label, inserted + export.n_splits, Some(export))
        }
    };
    let m = match &export {
        Some(export) => flow_mesh(shape, &label, export),
        None => FlowMesh {
            label,
            mesh: shape.mesh.clone(),
            geometry: shape.geometry.clone(),
            positions: shape.positions.clone(),
        },
    };
    // Errors in units of the input's mesh width, so that the rows compare.
    let h = (shape.geometry.total_area() / shape.mesh.n_vertices() as f64).sqrt();
    let line_distance = |p: [f64; 3]| p[0] - center[0];
    let circle_distance = |p: [f64; 3]| (p[0] - center[0]).hypot(p[1] - center[1]) - CIRCLE.0;
    let redistancing = [
        redistancing_error(&m, &line_distance),
        redistancing_error(&m, &circle_distance),
    ]
    .map(|(mean, max, unreached)| (mean / h, max / h, unreached));
    let default = circle_law(&m, center, false, steps);
    let paper = circle_law(&m, center, true, steps);
    // The way back, for the cut the default stages ended with.
    let back = match (&export, &default) {
        (Some(export), Ok(law)) => {
            export
                .check_state(export.fingerprint, law.phi.len())
                .unwrap();
            let cut = Cut::from_level_set(&export.mesh, &law.phi, 1e-3).unwrap();
            let (_, back) = carry_back(shape, &tri, export, &cut);
            assert!(node_degrees(&cut).iter().all(|&d| d == 2), "{}", shape.name);
            let exact = (CIRCLE.0 * CIRCLE.0 - 2.0 * CIRCLE.1 * CIRCLE.2).sqrt();
            let (mean, roundness) = radius_statistics(&back.positions, center);
            Some((back, mean / exact - 1.0, roundness))
        }
        _ => None,
    };
    let law_text = |law: &Result<Law, String>| match law {
        Ok(l) => format!(
            "{:>+9.2e} {:>8.1e} {:>3} {:>3}",
            l.error, l.roundness, l.unreached, l.clamped
        ),
        Err(e) => format!("{:>25}", format!("ERR[{}]", &e[..e.len().min(18)])),
    };
    let back_text = match &back {
        Some((b, error, roundness)) => format!(
            "{:>5} {:>6} {:>5}/{:<5} {:>8.1e} {:>8.1e} {:>+9.2e} {:>8.1e}",
            b.n_nodes,
            b.n_pieces,
            b.node_places[1],
            b.node_places[2],
            b.max_node_error,
            b.max_length_error,
            error,
            roundness
        ),
        None => String::new(),
    };
    let line = format!(
        "{:<26} {:<11} {:>5} {:>5} | {:>5.3} {:>6.3} {:>2} | {:>5.3} {:>6.3} {:>2} | {} | {} | {}",
        shape.name,
        m.label,
        m.mesh.n_vertices(),
        inserted,
        redistancing[0].0,
        redistancing[0].1,
        redistancing[0].2,
        redistancing[1].0,
        redistancing[1].1,
        redistancing[1].2,
        law_text(&default),
        law_text(&paper),
        back_text,
    );
    EndToEnd {
        line,
        refined: matches!(treatment, Some(Some(_))),
        law_default: default.map(|l| (l.error, l.roundness)),
        law_paper: paper.map(|l| (l.error, l.roundness)),
        redistancing,
        back,
    }
}

/// The end-to-end run: the two badly shaped flat meshes of the cut flow's bad-mesh study (the
/// sheared grid with angles of 153°, the needles-and-caps grid), refined, exported, the flow
/// run on the export, the cut carried back to the input mesh. For comparison the same on the
/// input mesh itself and on the Delaunay-flipped triangulation.
#[test]
fn the_cut_flow_runs_on_a_refined_export_and_its_cut_comes_back() {
    // The study's size (`FLAT_N = 32`, 96 steps) in a release build. A debug build of the
    // flow is some eighty times slower; there the same chain runs at half the size
    // (`E2E_SIZE=study` forces the study's size, `E2E_SIZE=small` the small one).
    let small = match std::env::var("E2E_SIZE").as_deref() {
        Ok("study") => false,
        Ok("small") => true,
        _ => cfg!(debug_assertions),
    };
    let (n, steps) = if small { (16, 24) } else { (32, 96) };
    let shapes = [
        common::jittered_grid(n, 0.0, 5),
        sheared_square(n, 2.0),
        needles_and_caps(n),
    ];
    // The centre of the study's circle, `(0.511, 0.493)` of the unit square, scaled with it.
    let centers: Vec<[f64; 2]> = shapes
        .iter()
        .map(|s| {
            let width = s.positions.iter().map(|p| p[0]).fold(0.0, f64::max);
            let height = s.positions.iter().map(|p| p[1]).fold(0.0, f64::max);
            // The sheared square was scaled to area 1 from the unit square.
            let scale = if s.name.starts_with("sheared") {
                height
            } else {
                1.0
            };
            assert!(width > 0.9);
            [0.511 * scale, 0.493 * scale]
        })
        .collect();
    let treatments = [None, Some(None), Some(Some(25.0)), Some(Some(30.0))];
    // One job per shape and treatment, in parallel.
    let rows: Vec<EndToEnd> = std::thread::scope(|scope| {
        let handles: Vec<_> = shapes
            .iter()
            .zip(&centers)
            .flat_map(|(shape, &center)| treatments.map(|t| (shape, center, t)))
            .map(|(shape, center, treatment)| {
                scope.spawn(move || end_to_end(shape, center, treatment, steps))
            })
            .collect();
        handles.into_iter().map(|h| h.join().unwrap()).collect()
    });
    println!(
        "redistancing: mean and max |phi - d| in input mesh widths, vertices not reached; circle law \
         ({steps} steps): relative error of the mean radius, roundness, unreached, clamped; \
         back: nodes, pieces on input faces, nodes on input edges / inside input faces, \
         node error, length error, radius error and roundness of the cut on the input mesh"
    );
    println!(
        "{:<26} {:<11} {:>5} {:>5} | {:^15} | {:^15} | {:^25} | {:^25} | back",
        "shape", "mesh", "V", "+V", "line", "circle", "law, default stages", "law, paper's stages"
    );
    for row in &rows {
        println!("{}", row.line);
    }
    // What the same flow does on a regular grid of the same size: the yardstick. (The paper's
    // stages have an error of this size in the radius there too; that is theirs, not the
    // mesh's.)
    let regular = &rows[0];
    let (base_default, base_paper) = (
        regular.law_default.clone().unwrap(),
        regular.law_paper.clone().unwrap(),
    );
    let base_circle = regular.redistancing[1];
    // "As well as on the regular grid": within a fifth at the study's size, within a factor
    // of two at the small one, where the circle has 36 cut points.
    let slack = if small { 2.0 } else { 1.2 };
    for row in rows.iter().filter(|r| r.refined) {
        // The flow ran with both stage sets, and the circle obeys its law as well as on the
        // regular grid.
        let (error, roundness) = row.law_default.clone().unwrap_or_else(|e| panic!("{e}"));
        assert!(
            error.abs() < slack * base_default.0.abs() && roundness < slack * base_default.1,
            "{}",
            row.line
        );
        let (error, roundness) = row.law_paper.clone().unwrap_or_else(|e| panic!("{e}"));
        assert!(
            error.abs() < slack * base_paper.0.abs() && roundness < slack * base_paper.1,
            "{}",
            row.line
        );
        // Redistancing: no vertex unreached; a linear function is kept, the distance to a
        // circle is as accurate as on the regular grid.
        let [line, circle] = row.redistancing;
        assert_eq!((line.2, circle.2), (0, 0), "{}", row.line);
        assert!(line.1 < 0.05, "{}", row.line);
        assert!(
            circle.0 < slack * base_circle.0 && circle.1 < slack * base_circle.1,
            "{}",
            row.line
        );
        // The cut came back: located to rounding, and the same circle on the input mesh.
        let (back, back_error, back_roundness) = row.back.as_ref().unwrap();
        assert!(
            back.max_node_error < 1e-10 && back.max_length_error < 1e-9,
            "{}",
            row.line
        );
        let (error, roundness) = row.law_default.clone().unwrap();
        assert!((back_error - error).abs() < 1e-12, "{}", row.line);
        assert!((back_roundness - roundness).abs() < 1e-12, "{}", row.line);
    }
    if small {
        return;
    }
    // And the meshes are bad: on the input mesh itself the same flow is worse.
    let of = |shape: usize, treatment: usize| &rows[4 * shape + treatment];
    // The sheared grid: redistancing is off by mesh widths, the circle is not round, and the
    // paper's stages lose the cut.
    let (sheared, sheared_refined) = (of(1, 0), of(1, 2));
    assert!(sheared.redistancing[1].0 > 10.0 * sheared_refined.redistancing[1].0);
    assert!(
        sheared.law_default.clone().unwrap().1
            > 5.0 * sheared_refined.law_default.clone().unwrap().1
    );
    assert!(sheared.law_paper.is_err(), "{}", sheared.line);
    // Needles and caps: redistancing and roundness are several times worse, the paper's
    // stages have more than twice the error.
    let (needles, needles_refined) = (of(2, 0), of(2, 2));
    assert!(needles.redistancing[1].0 > 2.0 * needles_refined.redistancing[1].0);
    assert!(
        needles.law_default.clone().unwrap().1
            > 5.0 * needles_refined.law_default.clone().unwrap().1
    );
    assert!(
        needles.law_paper.clone().unwrap().0.abs()
            > 2.0 * needles_refined.law_paper.clone().unwrap().0.abs()
    );
    // On the refined needles the nodes of the cut lie inside input faces.
    assert!(needles_refined.back.as_ref().unwrap().0.node_places[2] > 0);
}
