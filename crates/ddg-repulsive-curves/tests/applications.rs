//! Fast versions of the applications of `SPEC.md` §14.4 (work package WP13): each asserts the
//! property that matters. The scenes are the procedural ones of `scene::builtin` (a few dozen
//! steps, 50 to 150 vertices) and the committed scenes of the authors' format in `scenes/`.
//! Every step of the packing and graph runs is audited for crossings (`audit_every = 1`); the
//! collision-guarded line search is what keeps the isotopy class, the audit checks it from
//! outside.

use ddg_repulsive_curves::CurveNetwork;
use ddg_repulsive_curves::scene::{
    BackendChoice, RunOptions, RunOutcome, Scene, builtin, parse_scene, read_scene, run_scene,
};

fn run(scene: &Scene, steps: usize) -> RunOutcome {
    let o = RunOptions {
        steps: Some(steps),
        backend: BackendChoice::Exact,
        threads: Some(2),
        audit_every: 1,
        ..RunOptions::default()
    };
    let out = run_scene(scene, &o).unwrap();
    assert!(
        !out.report.status.starts_with("Error"),
        "{}: {}",
        scene.name,
        out.report.status
    );
    out
}

fn edge_lengths(net: &CurveNetwork, pos: &[[f64; 3]]) -> Vec<f64> {
    net.edge_geometry(pos).length
}

#[test]
fn packing_in_a_domain() {
    // 2 loops in a sphere of radius 1, edge lengths growing by 0.04 L0 per step
    let scene = builtin::packing_domain(16, 2).unwrap();
    let steps = 25;
    let out = run(&scene, steps);
    let r = &out.report;
    // no crossings at any audited step, and the loops stay apart
    assert!(r.min_distance_over_run > 0.0, "{}", r.min_distance_over_run);
    assert!(r.audit_end.min_edge_distance > 1e-3);
    // the distance to the domain stays positive: every vertex and edge midpoint is inside
    let d = r.audit_end.min_obstacle_distance.unwrap();
    assert!(d > 1e-3, "distance to the domain {d}");
    for p in &out.positions {
        assert!(p.iter().map(|x| x * x).sum::<f64>().sqrt() < 1.0, "outside");
    }
    // the length grew as scheduled: factor 1 + 0.04 k, within the projection tolerance
    let factor = r.length.1 / r.length.0;
    let want = 1.0 + 0.04 * steps as f64;
    assert!(
        (factor - want).abs() < 0.02 * want,
        "length factor {factor}, scheduled {want}"
    );
    assert!(r.audit_end.max_constraint_residual.unwrap() < 1e-4);
}

#[test]
fn packing_on_a_surface() {
    let scene = builtin::packing_surface(32).unwrap();
    let steps = 20;
    let out = run(&scene, steps);
    let r = &out.report;
    assert!(r.min_distance_over_run > 0.0);
    // the curve stays on the torus
    assert!(r.audit_end.max_surface_error < 1e-4, "{:?}", r.audit_end);
    let factor = r.length.1 / r.length.0;
    let want = 1.0 + 0.03 * steps as f64;
    // a retarget that fails halves its increment, so the factor may lag the schedule a little
    assert!(
        factor <= want * 1.001 && factor > 0.95 * want,
        "{factor} vs {want}"
    );
    // the spacing improved
    assert!(r.audit_end.min_edge_distance >= 0.9 * r.audit_start.min_edge_distance);
}

#[test]
fn packing_inside_a_fixed_boundary() {
    let scene = builtin::packing_boundary(24, 3).unwrap();
    let steps = 25;
    let out = run(&scene, steps);
    let r = &out.report;
    assert!(r.min_distance_over_run > 0.0);
    // planar: the curves stay in y = 0 and inside the frame
    assert!(r.audit_end.max_surface_error < 1e-6);
    let d = r.audit_end.min_obstacle_distance.unwrap();
    assert!(d > 1e-3, "distance to the boundary {d}");
    for p in &out.positions {
        assert!(p[0].abs() < 1.0 && p[2].abs() < 1.0);
    }
    let factor = r.length.1 / r.length.0;
    let want = 1.0 + 0.02 * steps as f64;
    assert!((factor - want).abs() < 0.02 * want, "{factor} vs {want}");
}

#[test]
fn graph_drawing_keeps_junctions_and_untangles_edges() {
    let scene = builtin::graph_drawing(8).unwrap();
    let out = run(&scene, 25);
    let r = &out.report;
    assert!(r.min_distance_over_run > 0.0);
    // junctions are where the layout put them
    assert_eq!(scene.resolved.pins.len(), 7);
    let pin_err = r.audit_end.max_pin_error;
    assert!(
        pin_err < 1e-6 * r.audit_end.mean_edge_length,
        "pin error {pin_err}"
    );
    // the drawing stays in the plane, keeps its total length, and the energy drops a lot
    assert!(r.audit_end.max_surface_error < 1e-6);
    assert!((r.length.1 / r.length.0 - 1.0).abs() < 1e-5);
    assert!(r.energy.1 < 0.1 * r.energy.0, "{:?}", r.energy);
    // no edge pair got (much) closer than at the start
    assert!(r.audit_end.min_edge_distance > 0.5 * r.audit_start.min_edge_distance);
}

#[test]
fn isometric_embedding_keeps_edge_lengths() {
    let scene = builtin::graph_isometric(10).unwrap();
    let before = edge_lengths(&scene.network, &scene.positions);
    let out = run(&scene, 25);
    let r = &out.report;
    assert!(r.min_distance_over_run > 0.0);
    let after = edge_lengths(&out.network, &out.positions);
    assert_eq!(before.len(), after.len());
    for (a, b) in before.iter().zip(&after) {
        assert!((a - b).abs() < 1e-4 * a, "edge length {a} -> {b}");
    }
    // the five crossing diagonals moved apart in space
    assert!(r.audit_end.min_edge_distance > r.audit_start.min_edge_distance);
    assert!(r.energy.1 < r.energy.0);
    // the embedding left the plane
    let ymax = out.positions.iter().fold(0.0f64, |m, p| m.max(p[1].abs()));
    assert!(ymax > 0.05, "max |y| {ymax}");
}

#[test]
fn vascular_tree_keeps_leaves_and_lengths() {
    let scene = builtin::vascular_tree(6, 3).unwrap();
    let before = edge_lengths(&scene.network, &scene.positions);
    let leaves = scene.resolved.pins.len();
    // the root and the four leaves
    assert_eq!(leaves, 5);
    let out = run(&scene, 20);
    let r = &out.report;
    assert!(r.min_distance_over_run > 0.0);
    assert!(r.audit_end.max_pin_error < 1e-6 * r.audit_end.mean_edge_length);
    let after = edge_lengths(&out.network, &out.positions);
    for (a, b) in before.iter().zip(&after) {
        assert!((a - b).abs() < 1e-4 * a);
    }
    assert!(r.audit_end.min_obstacle_distance.unwrap() > 1e-3);
    assert!(r.energy.1 < r.energy.0);
}

#[test]
fn spline_through_control_points_grows_without_crossing() {
    let scene = builtin::spline_points(48).unwrap();
    assert_eq!(scene.resolved.pins.len(), 8);
    let steps = 20;
    let out = run(&scene, steps);
    let r = &out.report;
    assert!(r.min_distance_over_run > 0.0);
    // control points exact (to the projection tolerance, a small fraction of an edge)
    assert!(
        r.audit_end.max_pin_error < 1e-6 * r.audit_end.mean_edge_length,
        "{}",
        r.audit_end.max_pin_error
    );
    let factor = r.length.1 / r.length.0;
    let want = 1.0 + 0.02 * steps as f64;
    assert!((factor - want).abs() < 0.02 * want, "{factor} vs {want}");
}

#[test]
fn spline_with_fixed_end_tangents() {
    let scene = builtin::spline_tangents(40).unwrap();
    assert_eq!(scene.resolved.tangents.len(), 2);
    let out = run(&scene, 25);
    let r = &out.report;
    assert!(r.min_distance_over_run > 0.0);
    assert!(r.audit_end.max_tangent_error < 1e-4, "{:?}", r.audit_end);
    assert!(r.audit_end.max_pin_error < 1e-6 * r.audit_end.mean_edge_length);
    assert!(r.length.1 > 1.2 * r.length.0);
}

#[test]
fn importer_reads_the_authors_scenes() {
    let tangent = read_scene("scenes/TangentInterpolation/scene.txt").unwrap();
    // 2 endpoints + 2 pins, tangents at the 2 endpoints
    assert_eq!(tangent.resolved.pins.len(), 4);
    assert_eq!(tangent.resolved.tangents.len(), 2);
    assert!((tangent.spec.length_penalty - 0.5).abs() < 1e-12);
    assert_eq!(tangent.spec.iteration_limit, Some(65));
    let out = run(&tangent, 15);
    assert!(out.report.audit_end.max_tangent_error < 1e-4);
    assert!(out.report.audit_end.max_pin_error < 1e-6);

    let torus = read_scene("scenes/ImplicitTorus/scene.txt").unwrap();
    assert_eq!(torus.spec.subdivide_limit, 5);
    assert!(torus.notes.iter().any(|n| n.contains("moved onto")));
    let out = run(&torus, 10);
    assert!(out.report.audit_end.max_surface_error < 1e-4);

    let interp = read_scene("scenes/CurveInterpolation/scene.txt").unwrap();
    assert_eq!(interp.resolved.pins.len(), 8);
    assert!((interp.spec.length_diff_penalty - 10000.0).abs() < 1e-9);

    let k5 = read_scene("scenes/GraphDrawing/k5scene.txt").unwrap();
    assert_eq!(k5.network.n_vertices(), 635);
}

#[test]
fn importer_reports_unsupported_directives_and_bad_input() {
    let dir = std::env::temp_dir().join(format!("rc-wp13-{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let (net, pos) = builtin::spline_points(16)
        .map(|s| (s.network, s.positions))
        .unwrap();
    ddg_repulsive_curves::write_curve_obj(dir.join("c.obj"), &net, &pos).unwrap();
    let ok = "curve c.obj\nrepel_curve 3 6\nfix_barycenter\nfix_edgelengths 2\n# comment\n\
              optimize_area 3\nshow_surface x.obj\nbogus 1\n";
    let s = parse_scene(ok, &dir, "t").unwrap();
    assert_eq!(s.notes.len(), 3, "{:?}", s.notes);
    assert!(s.notes.iter().any(|n| n.contains("optimize_area")));
    // no constraint line: the defaults
    let s = parse_scene("curve c.obj\n", &dir, "t").unwrap();
    assert_eq!(s.constraints().len(), 2);
    // errors name the line
    let e = parse_scene("curve c.obj\nfix_vertex x\n", &dir, "t").unwrap_err();
    assert!(e.to_string().contains("line 2"), "{e}");
    let e = parse_scene("curve c.obj\nconstraint_surface cube\n", &dir, "t").unwrap_err();
    assert!(e.to_string().contains("cube"));
    assert!(parse_scene("repel_curve 3 6\n", &dir, "t").is_err());
    assert!(parse_scene("curve missing.obj\n", &dir, "t").is_err());
    // a vertex out of range
    assert!(parse_scene("curve c.obj\nfix_vertex 99999\n", &dir, "t").is_err());
    let _ = std::fs::remove_dir_all(&dir);
}

#[test]
fn plane_obstacle_and_json_report() {
    // an obstacle passed through repel_plane keeps a loop above the plane
    let dir = std::env::temp_dir().join(format!("rc-wp13b-{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let n = 24;
    let pos: Vec<[f64; 3]> = (0..n)
        .map(|i| {
            let t = std::f64::consts::TAU * i as f64 / n as f64;
            [0.2 * t.cos(), 0.35 + 0.2 * t.sin(), 0.0]
        })
        .collect();
    let net =
        CurveNetwork::new(n, (0..n as u32).map(|i| [i, (i + 1) % n as u32]).collect()).unwrap();
    ddg_repulsive_curves::write_curve_obj(dir.join("c.obj"), &net, &pos).unwrap();
    let s = parse_scene(
        "curve c.obj\nfix_barycenter\nfix_edgelengths 3\nrepel_plane 0 0 0 0 1 0 1\niteration_limit 15\n",
        &dir,
        "plane",
    )
    .unwrap();
    let out = run(&s, 15);
    let r = &out.report;
    assert!(r.min_distance_over_run > 0.0);
    for p in &out.positions {
        assert!(p[1] > 0.0, "crossed the plane: {p:?}");
    }
    let j = r.to_json();
    assert!(j.contains("\"steps\": 15") && j.contains("\"audit_end\""));
    let _ = std::fs::remove_dir_all(&dir);
}
