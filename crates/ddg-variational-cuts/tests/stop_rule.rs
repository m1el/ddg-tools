//! The stop rule of a run (SPEC §4, D27): a run whose energy rises is not reported as stalled,
//! it is set back to its state of lowest energy, and the report says where that was.

mod common;

use common::{polar_angle, sphere};
use ddg_core::{Geometry, TriMesh};
use ddg_variational_cuts::{
    CutProblem, Energy, Flow, FlowOptions, RegionFlow, RunReport, StepSize, StopReason, StopRule,
    init,
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
