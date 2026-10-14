//! The stop rule of a run (SPEC §4, D27): a run whose energy rises is not reported as stalled,
//! it is set back to its state of lowest energy, and the report says where that was.

mod common;

use common::{polar_angle, sphere};
use ddg_core::{Geometry, TriMesh};
use ddg_variational_cuts::{
    CutError, CutProblem, Energy, Flow, FlowOptions, RegionFlow, RegionWeights, RunReport,
    StepSize, StopReason, StopRule, energy_criterion, init,
};

/// The energies of the states of a run: before every step, and after the last one.
fn energies(report: &RunReport) -> Vec<f64> {
    report
        .steps
        .iter()
        .map(|r| r.energy.total)
        .chain([report.last_energy])
        .collect()
}

/// What every report must satisfy.
fn check_report(report: &RunReport) {
    let all = energies(report);
    let vanished = report.stop == StopReason::CutVanished;
    // A state without a cut has no energy to compare.
    let counted = if vanished {
        &all[..all.len() - 1]
    } else {
        &all[..]
    };
    let lowest = counted.iter().copied().fold(f64::INFINITY, f64::min);
    assert_eq!(report.lowest_energy, lowest);
    assert_eq!(
        report.lowest_energy_step,
        counted.iter().position(|e| *e == lowest).unwrap()
    );
    if report.stop == StopReason::EnergyRose {
        assert_eq!(report.returned_step, report.lowest_energy_step);
        assert_eq!(report.final_energy.total, report.lowest_energy);
    } else {
        assert_eq!(report.returned_step, report.steps.len());
        assert_eq!(report.final_energy.total, report.last_energy);
    }
}

/// A run whose energy rises: the stages of the paper as printed with a fixed step that is too
/// large for them (sphere, a cap with three lobes, Dirichlet energy, `α_L = 1.5`). The energy
/// falls for nine steps and then rises until the cut is lost. With a window of two steps the
/// run ends when it has risen, says so, and is back in the state of lowest energy.
#[test]
fn a_rising_energy_is_not_reported_as_stalled() {
    let s = sphere(3);
    let phi0 = s.sample(|p| polar_angle(p) - 1.2 - 0.05 * (3.0 * p[1].atan2(p[0])).cos());
    let problem = CutProblem::new(Energy::Dirichlet, 1.5);
    let tau = 0.005;
    let options = FlowOptions::paper(StepSize::Fixed(tau));
    let stop = StopRule {
        max_steps: 80,
        energy_window: 2,
        ..StopRule::default()
    };

    let mut flow = Flow::new(&s.mesh, &s.geometry, &phi0).unwrap();
    let report = flow.run(&problem, &options, &stop).unwrap();
    check_report(&report);
    assert_eq!(report.stop, StopReason::EnergyRose);
    let all = energies(&report);
    let n = report.steps.len();
    // It rose over the window by more than the tolerance, and the last state is not the best.
    assert!(all[n - 1] - all[n - 3] > stop.energy_tolerance * all[n - 1]);
    assert!(report.lowest_energy_step > 3 && report.lowest_energy_step < n);
    assert!(report.last_energy > report.lowest_energy * 1.01);
    assert!(report.lowest_energy < all[0]);
    // The flow is in the state of lowest energy: its energy, its time, its step count.
    let k = report.returned_step;
    let now = flow.evaluate(&problem, &options.eval).unwrap();
    assert_eq!(now.energy.total, report.lowest_energy);
    assert_eq!(now.energy, report.final_energy);
    assert!((flow.time() - k as f64 * tau).abs() < 1e-15);
    assert_eq!(flow.n_steps(), k);
    // It is the state that the same steps lead to.
    let mut again = Flow::new(&s.mesh, &s.geometry, &phi0).unwrap();
    for _ in 0..k {
        again.step(&problem, &options).unwrap();
    }
    assert_eq!(again.level_set(), flow.level_set());

    // Without the energy criterion the same run goes on and loses the cut; the report still
    // says where the energy was lowest, and the state without a cut does not count.
    let mut unchecked = Flow::new(&s.mesh, &s.geometry, &phi0).unwrap();
    let unstopped = StopRule {
        energy_tolerance: 0.0,
        ..stop
    };
    let lost = unchecked.run(&problem, &options, &unstopped).unwrap();
    check_report(&lost);
    assert_eq!(lost.stop, StopReason::CutVanished);
    assert_eq!(lost.lowest_energy, report.lowest_energy);
    assert_eq!(lost.lowest_energy_step, report.lowest_energy_step);
    assert_eq!(lost.returned_step, lost.steps.len());

    // The same through the flow for several regions (two functions, φ and −φ).
    let phis = vec![phi0.clone(), phi0.iter().map(|x| -x).collect()];
    let mut regions = RegionFlow::new(&s.mesh, &s.geometry, &phis).unwrap();
    let two = CutProblem::for_regions(Energy::Dirichlet, 1.5, 2);
    let report_regions = regions.run(&two, &options, &stop).unwrap();
    check_report(&report_regions);
    assert_eq!(report_regions.stop, StopReason::EnergyRose);
    assert_eq!(report_regions.returned_step, k);
    assert_eq!(regions.level_sets()[0], flow.level_set());
    assert!((regions.time() - k as f64 * tau).abs() < 1e-15);
}

/// A run that comes to rest is still reported as stalled, and stays in its last state.
#[test]
fn a_run_at_rest_is_reported_as_stalled() {
    let s = sphere(3);
    let phi0 = s.sample(|p| polar_angle(p) - 1.2 - 0.05 * (3.0 * p[1].atan2(p[0])).cos());
    let problem = CutProblem::new(Energy::Dirichlet, 1.5);
    let options = FlowOptions::new(StepSize::Fixed(0.005));
    let stop = StopRule {
        max_steps: 80,
        energy_window: 5,
        ..StopRule::default()
    };
    let mut flow = Flow::new(&s.mesh, &s.geometry, &phi0).unwrap();
    let report = flow.run(&problem, &options, &stop).unwrap();
    check_report(&report);
    assert_eq!(report.stop, StopReason::EnergyStalled);
    assert_eq!(report.stop.name(), "energy_stalled");
    assert_eq!(StopReason::EnergyRose.name(), "energy_rose");
    let all = energies(&report);
    let n = report.steps.len();
    assert!((all[n - 1] - all[n - 6]).abs() <= stop.energy_tolerance * all[n - 1]);
    assert_eq!(flow.n_steps(), n);
}

/// The run that showed the defect of SPEC D26 (it needs the model `spot` under `data/models`
/// and is skipped without it, and in a debug build): Hencky energy, `α_L = 7.4·10⁻⁴`, eight
/// clusters of normals, the stages that `FlowOptions::for_regions` selects. The energy fell
/// from `1.1955·10⁻²` to `1.104·10⁻²` in 8 steps and then rose in every step, to `1.82·10⁻²`
/// after 26 steps and `2.68·10⁻²` after 40, and the run was reported as stalled.
#[test]
fn the_run_on_spot_that_rose_falls() {
    let path = std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../data/models/spot/spot/spot_triangulated.obj");
    if !path.exists() {
        println!("{} not found: skipped", path.display());
        return;
    }
    if cfg!(debug_assertions) {
        println!("skipped in a debug build");
        return;
    }
    let obj = ddg_core::io::read_obj(path.to_str().unwrap()).unwrap();
    let mut faces = Vec::new();
    for f in 0..obj.n_faces() {
        let corners = &obj.corner_verts[obj.face_corners(f)];
        for k in 1..corners.len().saturating_sub(1) {
            faces.push([corners[0], corners[k], corners[k + 1]]);
        }
    }
    let mesh = TriMesh::new(obj.positions.len(), faces).unwrap();
    let geometry = Geometry::from_positions(&mesh, &obj.positions).unwrap();
    let labels = init::normal_clusters(&mesh, &obj.positions, 8, 100).unwrap();
    let phis = init::regions_from_vertex_labels(&mesh, &geometry, &labels, 8).unwrap();
    let problem = CutProblem::for_regions(Energy::Hencky, 7.4e-4, 8);
    let options = FlowOptions::for_regions(StepSize::default());

    let mut flow = RegionFlow::new(&mesh, &geometry, &phis).unwrap();
    let stop = StopRule {
        max_steps: 60,
        ..StopRule::default()
    };
    let report = flow.run(&problem, &options, &stop).unwrap();
    check_report(&report);
    let all = energies(&report);
    println!(
        "spot: {:?} after {} steps, energy {:.4e} -> {:.4e}, lowest {:.4e} after {} steps",
        report.stop,
        report.steps.len(),
        all[0],
        report.final_energy.total,
        report.lowest_energy,
        report.lowest_energy_step
    );
    assert!((all[0] - 1.1955e-2).abs() < 1e-5);
    // It was 1.82e-2 after 26 steps and 4.35e-2 after 60.
    assert!(all[26] < 1.0e-2 && report.final_energy.total < 0.95e-2);
    assert!(report.lowest_energy_step > 40);
    // No single step raises the energy by more than a few percent.
    let worst = all
        .windows(2)
        .map(|w| (w[1] - w[0]) / w[0])
        .fold(0.0, f64::max);
    assert!(worst < 0.05, "a step raised the energy by {worst:e}");
}

/// The three cases of the energy criterion (SPEC §4, D27): for a positive tolerance and finite
/// energies exactly one holds; a tolerance that is not positive turns the criterion off; an
/// energy that is not finite is never a reason to stop.
#[test]
fn the_three_cases_of_the_energy_criterion() {
    let values = [
        -3.0,
        -1.0,
        -1e-9,
        0.0,
        1e-9,
        0.5,
        1.0,
        1.0 + 1e-6,
        2.0,
        1e12,
    ];
    for tolerance in [1e-12, 1e-4, 0.3, 2.0] {
        for before in values {
            for now in values {
                let (fell, rose) = (before - now, now - before);
                let limit = tolerance * f64::abs(now);
                let cases = [fell > limit, fell.abs() <= limit, rose > limit];
                assert_eq!(cases.iter().filter(|c| **c).count(), 1);
                let expected = if cases[0] {
                    None
                } else if cases[1] {
                    Some(StopReason::EnergyStalled)
                } else {
                    Some(StopReason::EnergyRose)
                };
                assert_eq!(energy_criterion(tolerance, before, now), expected);
            }
        }
    }
    // Off: the energy has not fallen here, and with a tolerance of 0 that is no reason to stop.
    for tolerance in [0.0, -1.0, f64::NAN] {
        for (before, now) in [(1.0, 1.0), (1.0, 2.0), (2.0, 1.0)] {
            assert_eq!(energy_criterion(tolerance, before, now), None);
        }
    }
    for bad in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        assert_eq!(energy_criterion(1e-4, bad, 1.0), None);
        assert_eq!(energy_criterion(1e-4, 1.0, bad), None);
        assert_eq!(energy_criterion(1e-4, bad, bad), None);
    }
}

/// A number of a step that is not finite is an error that names the step and the number, and
/// the flow stays in the state before the step (SPEC §4): never a stop reason.
#[test]
fn a_non_finite_number_is_an_error_and_not_a_stop_reason() {
    let s = sphere(2);
    let phi0 = s.sample(|p| polar_angle(p) - 1.1);
    let good = CutProblem::new(Energy::Dirichlet, 1.0);
    let options = FlowOptions::new(StepSize::default());
    let stop = StopRule {
        max_steps: 40,
        ..StopRule::default()
    };
    // Weights that are valid one by one and overflow in the energy.
    let overflow = CutProblem {
        regions: vec![
            RegionWeights {
                distortion: 1e308,
                area: 0.0,
            };
            2
        ],
        ..good.clone()
    };
    let heavy = CutProblem::new(Energy::Dirichlet, f64::MAX);

    let mut flow = Flow::new(&s.mesh, &s.geometry, &phi0).unwrap();
    for bad in [&overflow, &heavy] {
        let error = flow.run(bad, &options, &stop).unwrap_err();
        println!("{error}");
        assert!(
            matches!(&error, CutError::NonFiniteStep { step: 0, value, .. } if !value.is_finite()),
            "{error:?}"
        );
        assert_eq!(flow.n_steps(), 0);
        assert_eq!(flow.time(), 0.0);
    }
    // After two good steps the error names step 2, and the state is the one after them.
    flow.step(&good, &options).unwrap();
    flow.step(&good, &options).unwrap();
    let state = flow.level_set();
    let error = flow.step(&overflow, &options).unwrap_err();
    assert!(
        matches!(error, CutError::NonFiniteStep { step: 2, .. }),
        "{error:?}"
    );
    assert!(error.to_string().contains("after 2 steps"));
    assert_eq!(flow.level_set(), state);
    assert_eq!(flow.n_steps(), 2);
    // The flow goes on with a valid problem.
    flow.step(&good, &options).unwrap();

    // A level-set function that leaves the range of `f64` in a step: a step size near the
    // limit, no length term.
    let free = CutProblem::new(Energy::Dirichlet, 0.0);
    let huge = FlowOptions::new(StepSize::Fixed(1e308));
    let mut flow = Flow::new(&s.mesh, &s.geometry, &phi0).unwrap();
    let state = flow.level_set();
    let error = flow.step(&free, &huge).unwrap_err();
    println!("{error}");
    assert!(
        matches!(&error, CutError::NonFiniteStep { step: 0, what, .. } if what.contains("level-set")),
        "{error:?}"
    );
    assert_eq!(flow.level_set(), state);

    // The same through the flow for several regions.
    let phis = vec![phi0.clone(), phi0.iter().map(|x| -x).collect()];
    let mut regions = RegionFlow::new(&s.mesh, &s.geometry, &phis).unwrap();
    let state = regions.level_sets();
    for (problem, options) in [(&overflow, &options), (&heavy, &options), (&free, &huge)] {
        let error = regions.run(problem, options, &stop).unwrap_err();
        assert!(
            matches!(error, CutError::NonFiniteStep { step: 0, .. }),
            "{error:?}"
        );
        assert_eq!(regions.level_sets(), state);
        assert_eq!(regions.n_steps(), 0);
    }
}

/// A tolerance of zero turns a criterion off: a run whose steps move nothing (no distortion
/// weight, no length weight) goes on to its last step with `motion_tolerance = 0`, and ends
/// as converged after one step with a positive tolerance.
#[test]
fn a_tolerance_of_zero_turns_a_criterion_off() {
    let s = sphere(2);
    let phi0 = s.sample(|p| polar_angle(p) - 1.1);
    let still = CutProblem {
        regions: vec![
            RegionWeights {
                distortion: 0.0,
                area: 0.0,
            };
            2
        ],
        ..CutProblem::new(Energy::Dirichlet, 0.0)
    };
    let options = FlowOptions::new(StepSize::default());
    let off = StopRule {
        max_steps: 3,
        motion_tolerance: 0.0,
        energy_tolerance: 0.0,
        energy_window: 1,
    };
    let mut flow = Flow::new(&s.mesh, &s.geometry, &phi0).unwrap();
    let report = flow.run(&still, &options, &off).unwrap();
    check_report(&report);
    assert!(report.steps.iter().all(|r| r.cut_motion == 0.0));
    assert_eq!(report.stop, StopReason::MaxSteps);
    assert_eq!(report.steps.len(), 3);

    let on = StopRule {
        motion_tolerance: 1e-9,
        ..off
    };
    let mut flow = Flow::new(&s.mesh, &s.geometry, &phi0).unwrap();
    let report = flow.run(&still, &options, &on).unwrap();
    assert_eq!(report.stop, StopReason::Converged);
    assert_eq!(report.steps.len(), 1);
    // The energy criterion: the energy stands still, which a tolerance of zero does not stop.
    let energy = StopRule {
        energy_tolerance: 1e-4,
        ..off
    };
    let mut flow = Flow::new(&s.mesh, &s.geometry, &phi0).unwrap();
    let report = flow.run(&still, &options, &energy).unwrap();
    assert_eq!(report.stop, StopReason::EnergyStalled);
    assert_eq!(report.steps.len(), 2);
}

/// The motion criterion (SPEC §4, D31): both the motion of the step and the motion that its
/// rate predicts must be within a positive tolerance.
#[test]
fn the_motion_criterion() {
    use ddg_variational_cuts::motion_criterion;
    assert!(motion_criterion(1e-5, 1e-6, 1e-6));
    assert!(motion_criterion(1e-5, 1e-5, 1e-5));
    assert!(motion_criterion(1e-5, 0.0, 0.0));
    // A short step of a cut that is not at rest.
    assert!(!motion_criterion(1e-5, 1e-10, 8e-3));
    // A step that moved the function further than its rate predicts.
    assert!(!motion_criterion(1e-5, 8e-3, 1e-10));
    // Off.
    for tolerance in [0.0, -1.0, f64::NAN] {
        assert!(!motion_criterion(tolerance, 0.0, 0.0));
    }
    assert!(!motion_criterion(1e-5, f64::NAN, 0.0));
    assert!(!motion_criterion(1e-5, 0.0, f64::NAN));
}

/// A mesh from a fixture: `mesh NV NF NE`, `centre V`, `f A B C`, `e A B LENGTH` with the
/// length as the bits of an `f64` in hexadecimal.
fn fixture_mesh(name: &str) -> (TriMesh, Geometry, usize) {
    let path = std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests/fixtures")
        .join(name);
    let text = std::fs::read_to_string(&path).unwrap_or_else(|e| panic!("{}: {e}", path.display()));
    let (mut n_vertices, mut centre) = (0, 0);
    let mut faces = Vec::new();
    let mut edges = Vec::new();
    for line in text.lines() {
        let words: Vec<&str> = line.split_whitespace().collect();
        let number = |k: usize| words[k].parse::<usize>().unwrap();
        match words.first().copied() {
            Some("mesh") => n_vertices = number(1),
            Some("centre") => centre = number(1),
            Some("f") => faces.push([number(1) as u32, number(2) as u32, number(3) as u32]),
            Some("e") => edges.push((
                number(1),
                number(2),
                f64::from_bits(u64::from_str_radix(words[3], 16).unwrap()),
            )),
            _ => {}
        }
    }
    let mesh = TriMesh::new(n_vertices, faces).unwrap();
    let mut lengths = vec![f64::NAN; mesh.n_edges()];
    for (a, b, length) in edges {
        lengths[mesh.find_edge(a, b).unwrap()] = length;
    }
    let geometry = Geometry::from_edge_lengths(&mesh, lengths).unwrap();
    (mesh, geometry, centre)
}

/// The run that was reported as converged after one step (SPEC D30, D31): the needles-and-caps
/// grid of `examples/bad_mesh_study.rs` with 16 cells, flipped to its intrinsic Delaunay
/// triangulation and exported by `ddg-intrinsic` (the fixture has its connectivity and edge
/// lengths), a geodesic circle around the middle vertex, Hencky energy, `α_L = 1`, default
/// options. Six vertices lie on the circle up to rounding, their crossings are moved by `δ`,
/// and the slope `|φ_a| / dist(a, pq)` of their faces was `2·10⁻¹²`: the mean slope of the
/// whole cut was `2·10⁻⁸`, the first step moved the cut by `1.5·10⁻¹⁰`, and the run ended.
#[test]
fn a_cut_through_vertices_moves_and_is_not_called_converged() {
    let (mesh, geometry, centre) = fixture_mesh("needles_and_caps_16_intrinsic.txt");
    let radius = 0.2 * geometry.total_area().sqrt();
    let phi0 = init::geodesic_circle(&mesh, &geometry, centre, radius).unwrap();
    let problem = CutProblem::new(Energy::Hencky, 1.0);
    let options = FlowOptions::new(StepSize::default());

    let mut flow = Flow::new(&mesh, &geometry, &phi0).unwrap();
    let start = flow.evaluate(&problem, &options.eval).unwrap();
    assert!((start.energy.total - 1.167016).abs() < 1e-6);
    // The crossings that were moved are flagged, and their faces do not count in the slope.
    assert_eq!(start.cut.n_clamped, 6);
    assert_eq!(start.cut.moved.iter().filter(|m| **m).count(), 6);
    let (first, stages) = flow.step_with_stages(&problem, &options).unwrap();
    assert!(stages.gradient_scale.iter().all(|g| (0.9..1.3).contains(g)));
    let slopes = first.gradient_norm.unwrap();
    assert!(slopes.min > 0.5, "{slopes:?}");
    assert!(
        first.cut_motion > 5e-3,
        "the first step moved the cut by {:e}",
        first.cut_motion
    );
    assert!(first.predicted_motion > 5e-3);

    let mut flow = Flow::new(&mesh, &geometry, &phi0).unwrap();
    let stop = StopRule {
        max_steps: 8,
        ..StopRule::default()
    };
    let report = flow.run(&problem, &options, &stop).unwrap();
    check_report(&report);
    assert_eq!(report.stop, StopReason::MaxSteps);
    assert_eq!(report.steps.len(), 8);
    // It was 1.167 after one step; with the criterion off, 0.894 after eight.
    assert!(
        report.final_energy.total < 0.87,
        "{}",
        report.final_energy.total
    );
}

/// A run that comes to rest still ends by the motion criterion, in as many steps as the motion
/// of its steps alone would have ended it (the rule before SPEC D31): a cap with three lobes on
/// the sphere under a length weight that makes it round. At its last step both the motion and
/// the predicted motion are within the tolerance, and at no earlier step.
#[test]
fn a_cut_at_rest_is_still_called_converged() {
    use ddg_variational_cuts::motion_criterion;
    let s = sphere(3);
    let cap = s.sample(|p| polar_angle(p) - 1.2 - 0.02 * (3.0 * p[1].atan2(p[0])).cos());
    let problem = CutProblem::new(Energy::Dirichlet, 6.0);
    let options = FlowOptions::new(StepSize::default());
    let stop = StopRule {
        max_steps: 400,
        energy_tolerance: 0.0,
        ..StopRule::default()
    };
    let mut flow = Flow::new(&s.mesh, &s.geometry, &cap).unwrap();
    let report = flow.run(&problem, &options, &stop).unwrap();
    check_report(&report);
    assert_eq!(report.stop, StopReason::Converged);
    let n = report.steps.len();
    for (k, r) in report.steps.iter().enumerate() {
        assert_eq!(
            motion_criterion(stop.motion_tolerance, r.cut_motion, r.predicted_motion),
            k + 1 == n
        );
        // The prediction is the motion of the step where nothing holds it back.
        assert!(r.predicted_motion.is_finite() && r.predicted_motion >= 0.0);
    }
    // The first step whose motion alone is within the tolerance.
    let by_motion_alone = report
        .steps
        .iter()
        .position(|r| r.cut_motion <= stop.motion_tolerance)
        .unwrap()
        + 1;
    println!(
        "sphere: converged after {n} steps; by the motion of the step alone after {by_motion_alone}"
    );
    assert!(n <= by_motion_alone + by_motion_alone / 10);
}
