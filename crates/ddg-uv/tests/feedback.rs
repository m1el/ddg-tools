//! The feedback loop of the cone solver (`ddg_optimal_cones::Feedback`), driven by the
//! pipeline (`UnwrapOptions::cone_feedback`) and by `solve_with_feedback` with Boundary First
//! Flattening: it finds the tips of thin spikes, changes nothing on smooth meshes, and refuses
//! invalid options.

mod common;

use common::{check, without_times};
use ddg_bff::{ConeFlattenInput, ConeFlattenOptions, CutSource, flatten_with_cones};
use ddg_core::generators::{bumpy_grid_disk, icosphere};
use ddg_core::{Geometry, TriMesh};
use ddg_optimal_cones::{
    BoundaryTreatment, ConeError, ConeProblem, FeedbackError, FeedbackOptions, FeedbackStop,
    SolveOptions, solve_with_feedback,
};
use ddg_uv::primitives::{capsule, uv_sphere};
use ddg_uv::{
    PolygonMesh, SpikeAction, SpikeCheck, UnwrapOptions, UvError, Value,
    realized_face_log_scale, unwrap,
};

/// A 24 x 24 grid with two vertices pulled 2.0 out of a sheet of side 1 (`h/ρ ≈ 48`), behind
/// three loose vertices (so that the numbers of the caller differ from those of the pipeline's
/// meshes). The same mesh as in `tests/spikes.rs`.
fn spiked_disk() -> (PolygonMesh, [u32; 2]) {
    let g = bumpy_grid_disk(24, 24, 0.02).unwrap();
    let at = |i: usize, j: usize| (i * 25 + j) as u32;
    let spikes = [at(12, 12), at(5, 18)];
    let mut positions = vec![[7.0, 7.0, 7.0], [8.0, 8.0, 8.0], [9.0, 9.0, 9.0]];
    positions.extend_from_slice(&g.positions);
    for s in spikes {
        positions[s as usize + 3][2] += 2.0;
    }
    let faces: Vec<[u32; 3]> = g.faces.iter().map(|f| f.map(|v| v + 3)).collect();
    let mesh = PolygonMesh::from_triangles(positions, &faces);
    let mut tips = spikes.map(|s| s + 3);
    tips.sort_unstable();
    (mesh, tips)
}

fn with_feedback(spike_check: Option<SpikeCheck>) -> UnwrapOptions {
    UnwrapOptions {
        spike_check,
        cone_feedback: Some(FeedbackOptions::default()),
        ..UnwrapOptions::default()
    }
}

#[test]
fn the_loop_finds_the_tips_of_thin_spikes() {
    let (mesh, tips) = spiked_disk();
    let input = mesh.input();
    // The spike check off: the loop alone has to find the tips.
    let options = with_feedback(None);
    let map = unwrap(&input, &options).unwrap();
    check(&mesh, &input, &map, &options);
    let feedback = map.report.charts[0].feedback.as_ref().expect("the loop ran");
    assert_eq!(feedback.reweighted, tips.to_vec(), "{:?}", feedback.report);
    assert_eq!(feedback.report.stop, FeedbackStop::Stable);
    let rounds = &feedback.report.rounds;
    assert!(rounds.len() >= 2, "{rounds:?}");
    assert!(rounds[0].max_abs_eligible > 10.0 * options.cone_feedback.unwrap().threshold);
    assert!(rounds.last().unwrap().max_abs_eligible <= options.cone_feedback.unwrap().threshold);
    // The final weights are zero at the tips (in the piece's numbering) and one elsewhere.
    let zero = feedback.report.weights.iter().filter(|&&w| w == 0.0).count();
    assert_eq!(zero, 2);

    // It ends where the spike check's exclusion of the same tips starts: the same cones.
    let excluded = UnwrapOptions {
        spike_check: Some(SpikeCheck {
            action: SpikeAction::ExcludeFromDistortion,
            ..SpikeCheck::default()
        }),
        ..UnwrapOptions::default()
    };
    let by_check = unwrap(&input, &excluded).unwrap();
    assert_eq!(by_check.report.charts[0].spikes.as_ref().unwrap().excluded, tips);
    assert_eq!(map.cones, by_check.cones);

    // The report and the warning.
    let summary = map.report.summary();
    let get = |key: &str| {
        summary
            .iter()
            .find(|(k, _)| *k == key)
            .map(|(_, v)| v.clone())
    };
    assert_eq!(get("n_cone_feedback_reweighted"), Some(Value::Int(2)));
    assert_eq!(
        get("cone_feedback_vertices"),
        Some(Value::Text(format!("{}, {}", tips[0], tips[1])))
    );
    assert!(matches!(get("cone_feedback_rounds"), Some(Value::Int(n)) if n >= 2));
    assert!(matches!(get("seconds_cone_feedback"), Some(Value::Float(s)) if s > 0.0));
    let warnings: Vec<&String> = map
        .report
        .warnings
        .iter()
        .filter(|w| w.contains("cone feedback"))
        .collect();
    assert_eq!(warnings.len(), 1, "{:?}", map.report.warnings);
    assert!(warnings[0].contains(&format!("vertices {}, {}", tips[0], tips[1])));
}

#[test]
fn with_the_spike_check_the_loop_has_nothing_left_to_do_at_the_tips() {
    let (mesh, tips) = spiked_disk();
    let options = with_feedback(Some(SpikeCheck::default()));
    let map = unwrap(&mesh.input(), &options).unwrap();
    let feedback = map.report.charts[0].feedback.as_ref().unwrap();
    assert!(feedback.reweighted.iter().all(|v| !tips.contains(v)));
    let plain = unwrap(&mesh.input(), &UnwrapOptions::default()).unwrap();
    if feedback.reweighted.is_empty() {
        assert_eq!(feedback.report.rounds.len(), 1);
        assert_eq!(map.uvs, plain.uvs);
        assert_eq!(map.cones, plain.cones);
    }
}

/// A sphere with one vertex pulled out to three times the radius: a closed piece, handled by
/// the default policy for closed surfaces; the loop finds the tip there too.
#[test]
fn the_loop_finds_a_spike_on_a_sphere() {
    let ico = icosphere(3).unwrap(); // 642 vertices
    let mut positions = ico.positions.clone();
    let tip = 17usize;
    positions[tip] = positions[tip].map(|x| 3.0 * x);
    let mesh = PolygonMesh::from_triangles(positions, &ico.faces);
    let input = mesh.input();
    let options = with_feedback(None);
    let map = unwrap(&input, &options).unwrap();
    check(&mesh, &input, &map, &options);
    let feedback = map.report.charts[0].feedback.as_ref().expect("the loop ran");
    assert!(
        feedback.reweighted.contains(&(tip as u32)),
        "{:?}",
        feedback.report
    );
    assert!(feedback.reweighted.len() <= 6, "{:?}", feedback.reweighted);
}

#[test]
fn smooth_meshes_are_not_changed() {
    let shapes = [
        (
            "bumpy disk",
            PolygonMesh::from_triangles(
                bumpy_grid_disk(24, 24, 0.3).unwrap().positions,
                &bumpy_grid_disk(24, 24, 0.3).unwrap().faces,
            ),
        ),
        ("uv sphere", uv_sphere(32, 16, 1.0).unwrap()),
        ("capsule", capsule(24, 8, 0.5, 1.0).unwrap()),
    ];
    for (name, mesh) in shapes {
        let input = mesh.input();
        let plain = unwrap(&input, &UnwrapOptions::default()).unwrap();
        let options = with_feedback(Some(SpikeCheck::default()));
        let map = unwrap(&input, &options).unwrap();
        check(&mesh, &input, &map, &options);
        for chart in &map.report.charts {
            let Some(feedback) = &chart.feedback else {
                continue;
            };
            assert!(feedback.reweighted.is_empty(), "{name}: {:?}", feedback.report);
            assert_eq!(feedback.report.rounds.len(), 1, "{name}");
            assert_eq!(feedback.report.stop, FeedbackStop::Stable, "{name}");
        }
        // One round compares what the pipeline made anyway: the map is the same.
        assert_eq!(map.uvs, plain.uvs, "{name}");
        assert_eq!(map.cones, plain.cones, "{name}");
        let mut a = without_times(map);
        let b = without_times(plain);
        for chart in &mut a.report.charts {
            chart.feedback = None;
        }
        assert_eq!(a.report.charts, b.report.charts, "{name}");
    }
}

#[test]
fn invalid_options_are_errors_not_panics() {
    let (mesh, _) = spiked_disk();
    let d = FeedbackOptions::default();
    for (bad, name) in [
        (FeedbackOptions { threshold: 0.0, ..d }, "cone_feedback.threshold"),
        (FeedbackOptions { threshold: f64::NAN, ..d }, "cone_feedback.threshold"),
        (FeedbackOptions { reduced_weight: 1.0, ..d }, "cone_feedback.reduced_weight"),
        (FeedbackOptions { max_rounds: 0, ..d }, "cone_feedback.max_rounds"),
        (
            FeedbackOptions {
                max_reweighted_fraction: -0.5,
                ..d
            },
            "cone_feedback.max_reweighted_fraction",
        ),
    ] {
        let options = UnwrapOptions {
            cone_feedback: Some(bad),
            ..UnwrapOptions::default()
        };
        match unwrap(&mesh.input(), &options) {
            Err(UvError::InvalidOption { name: got, .. }) => assert_eq!(got, name),
            other => panic!("{bad:?}: {other:?}"),
        }
    }
}

/// `solve_with_feedback` with Boundary First Flattening as the realization: a flat grid with one
/// lifted vertex. The loop finds the tip, and the second solve is the one with zero weight there.
#[test]
fn solve_with_feedback_and_boundary_first_flattening() {
    let g = bumpy_grid_disk(20, 20, 0.0).unwrap();
    let tip = 10 * 21 + 10;
    let mut positions = g.positions.clone();
    positions[tip][2] = 1.5;
    let mesh = TriMesh::new(positions.len(), g.faces.clone()).unwrap();
    let geometry = Geometry::from_positions(&mesh, &positions).unwrap();
    let problem = ConeProblem::new(&mesh, &geometry, 0.1);
    let mut flattenings = 0;
    let realize = |s: &ddg_optimal_cones::Solution| {
        flattenings += 1;
        let flat = flatten_with_cones(
            &ConeFlattenInput {
                mesh: &mesh,
                geometry: &geometry,
                cones: &s.cones,
                scale_factors: Some(&s.log_scale),
                cut: CutSource::Plan {
                    may_cut: None,
                    must_cut: &[],
                },
                anchors: &[],
            },
            &ConeFlattenOptions::default(),
        )?;
        let faces = realized_face_log_scale(&geometry, &flat.corner_uvs);
        Ok::<_, ddg_bff::BffError>((faces, flat.corner_uvs))
    };
    let result = solve_with_feedback(
        &problem,
        BoundaryTreatment::Fixed,
        &SolveOptions::default(),
        &FeedbackOptions::default(),
        realize,
    )
    .unwrap();
    let report = &result.feedback;
    assert_eq!(report.reweighted, vec![tip as u32], "{report:?}");
    assert_eq!(report.stop, FeedbackStop::Stable);
    assert_eq!(report.rounds.len(), 2);
    assert_eq!(flattenings, 2);
    assert_eq!(result.realized.len(), 3 * mesh.n_faces());
    // The kept solution is the optimum of the problem with zero weight at the tip.
    let mut weights = vec![1.0; mesh.n_vertices()];
    weights[tip] = 0.0;
    assert_eq!(report.weights, weights);
    let direct = ddg_optimal_cones::solve(
        &ConeProblem {
            distortion_weights: Some(weights),
            ..problem.clone()
        },
        &SolveOptions::default(),
    )
    .unwrap();
    assert_eq!(result.solution.cones, direct.cones);

    // Invalid options are an error of the set-up, before anything is solved.
    let bad = FeedbackOptions {
        max_rounds: 0,
        ..FeedbackOptions::default()
    };
    let error = solve_with_feedback(
        &problem,
        BoundaryTreatment::Fixed,
        &SolveOptions::default(),
        &bad,
        |_| Ok::<_, String>((Vec::new(), ())),
    )
    .unwrap_err();
    assert!(matches!(
        error,
        FeedbackError::Cones(ConeError::InvalidOptions { .. })
    ));
    // A realization of the wrong length is an error, not a panic.
    let error = solve_with_feedback(
        &problem,
        BoundaryTreatment::Fixed,
        &SolveOptions::default(),
        &FeedbackOptions::default(),
        |_| Ok::<_, String>((vec![0.0; 3], ())),
    )
    .unwrap_err();
    assert!(matches!(
        error,
        FeedbackError::Cones(ConeError::InvalidArgument { .. })
    ));
    // A failing first realization is the caller's error.
    let error = solve_with_feedback(
        &problem,
        BoundaryTreatment::Fixed,
        &SolveOptions::default(),
        &FeedbackOptions::default(),
        |_| Err::<(Vec<f64>, ()), _>("no flattener".to_owned()),
    )
    .unwrap_err();
    assert_eq!(error, FeedbackError::Realize("no flattener".to_owned()));
}
