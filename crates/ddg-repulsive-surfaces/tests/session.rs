//! The session API (`SPEC.md` §7.1, §7.4, §9.1): determinism across thread counts, stepping by
//! hand equals [`flow`], cancellation, status codes, stuck handling, option and setup errors.

use std::sync::atomic::Ordering;
use std::time::Duration;

use ddg_repulsive_surfaces::constraints::Constraint;
use ddg_repulsive_surfaces::generate::primitives as generate;
use ddg_repulsive_surfaces::generate::Rng;
use ddg_repulsive_surfaces::solve::iterative::{Backend, HierarchicalOptions};
use ddg_repulsive_surfaces::{
    BhOptions, Control, EnergyEval, FlowOptions, LineSearch, Method, RemeshOptions, RsError,
    Schedule, Session, StartRule, Status, StepStatus, StopRule, Surface, flow,
};

fn bumpy_sphere(level: u32, a: f64, seed: u64) -> (Surface, Vec<[f64; 3]>) {
    let m = generate::sphere(level, 1.0, [0.0; 3]).unwrap();
    let mut rng = Rng::new(seed);
    let p: Vec<[f64; 3]> = m
        .positions
        .iter()
        .map(|q| {
            let s = 1.0 + a * rng.uniform(-1.0, 1.0);
            [q[0] * s, q[1] * s * 0.9, q[2] * s * 1.2]
        })
        .collect();
    (Surface::new(p.len(), m.faces.clone()).unwrap(), p)
}

fn base() -> FlowOptions {
    FlowOptions {
        constraints: vec![
            Constraint::Barycenter { component: None, target: None },
            Constraint::Area { component: None, target: None },
        ],
        stop: StopRule::Off,
        threads: Some(2),
        ..FlowOptions::default()
    }
}

/// A step report without its wall times.
fn untimed(mut r: ddg_repulsive_surfaces::StepReport) -> ddg_repulsive_surfaces::StepReport {
    r.seconds = Default::default();
    r.linear.factor_seconds = 0.0;
    r.linear.solve_seconds = 0.0;
    if let Some(m) = r.remesh.as_mut() {
        m.seconds = 0.0;
    }
    if let Some(p) = r.products.as_mut() {
        p.1 = 0.0;
    }
    if let Some(c) = r.ccd.as_mut() {
        c.seconds = 0.0;
    }
    r
}

fn bits(p: &[[f64; 3]]) -> Vec<u64> {
    p.iter().flatten().map(|x| x.to_bits()).collect()
}

#[test]
fn trajectories_do_not_depend_on_the_thread_count() {
    // Hierarchical back end (1280 faces > dense_max_faces), Barnes–Hut, remeshing every step.
    let (s, p) = bumpy_sphere(3, 0.08, 1);
    type Trajectory = (Vec<[u32; 3]>, Vec<u64>, Vec<u64>);
    let mut reference: Option<Trajectory> = None;
    for threads in [1, 3, 8] {
        let o = FlowOptions {
            threads: Some(threads),
            remesh: Some(RemeshOptions::default()),
            energy_eval: EnergyEval::BarnesHut(BhOptions::default()),
            backend: Backend::Auto { dense_max_faces: 1000, hierarchical: HierarchicalOptions::default() },
            max_iterations: 4,
            ..base()
        };
        let r = flow(s.clone(), p.clone(), o, None).unwrap();
        assert_eq!(r.report.status, Status::MaxIterations);
        assert!(r.report.steps[0].linear.solver.contains("gmres"));
        let energies: Vec<u64> = r.report.steps.iter().map(|s| s.energy_after.to_bits()).collect();
        let got = (r.surface.faces().to_vec(), bits(&r.positions), energies);
        match &reference {
            None => reference = Some(got),
            Some(want) => assert!(&got == want, "{threads} threads differ from 1 thread"),
        }
    }
}

#[test]
fn stepping_by_hand_equals_flow() {
    let (s, p) = bumpy_sphere(2, 0.1, 2);
    let o = FlowOptions {
        max_iterations: 6,
        remesh: Some(RemeshOptions::default()),
        schedules: vec![Schedule { constraint: 1, factor: 1.1, offset: 0.0, steps: 3 }],
        ..base()
    };
    let r = flow(s.clone(), p.clone(), o.clone(), None).unwrap();
    let mut sess = Session::new(s, p, o).unwrap();
    for k in 0..6 {
        let rep = sess.step().unwrap();
        assert_eq!(untimed(rep), untimed(r.report.steps[k].clone()), "step {k}");
    }
    assert_eq!(sess.status(), Status::MaxIterations);
    assert_eq!(sess.surface().faces(), r.surface.faces());
    assert_eq!(bits(sess.positions()), bits(&r.positions));
    assert_eq!(sess.energy().to_bits(), r.report.final_energy.to_bits());
}

#[test]
fn cancellation() {
    let (s, p) = bumpy_sphere(2, 0.1, 3);
    // Flag set before the run: no step.
    let mut sess = Session::new(s.clone(), p.clone(), base()).unwrap();
    sess.cancel_flag().store(true, Ordering::Relaxed);
    let r = sess.run(None).unwrap();
    assert_eq!(r.status, Status::Cancelled);
    assert!(r.steps.is_empty());
    assert_eq!(bits(sess.positions()), bits(&p));
    // step() with the flag set: the error, the positions of the last completed step.
    assert!(matches!(sess.step(), Err(RsError::Cancelled)));
    assert_eq!(sess.status(), Status::Cancelled);
    // The progress callback stops the run after two steps.
    let mut sess = Session::new(s.clone(), p.clone(), base()).unwrap();
    let mut n = 0;
    let mut cb = |_: &ddg_repulsive_surfaces::StepReport| {
        n += 1;
        if n == 2 { Control::Stop } else { Control::Continue }
    };
    let r = sess.run(Some(&mut cb)).unwrap();
    assert_eq!(r.status, Status::Cancelled);
    assert_eq!(r.steps.len(), 2);
    assert_eq!(sess.iteration(), 2);
    // The flag set from the callback (as another thread would): the run ends before the next
    // step, and the session can be resumed after clearing it.
    let mut sess = Session::new(s, p, base()).unwrap();
    let flag = sess.cancel_flag();
    let mut cb = |_: &ddg_repulsive_surfaces::StepReport| {
        flag.store(true, Ordering::Relaxed);
        Control::Continue
    };
    let r = sess.run(Some(&mut cb)).unwrap();
    assert_eq!(r.status, Status::Cancelled);
    assert_eq!(r.steps.len(), 1);
    sess.cancel_flag().store(false, Ordering::Relaxed);
    let before = sess.iteration();
    sess.step().unwrap();
    assert_eq!(sess.iteration(), before + 1);
}

#[test]
fn status_codes() {
    let (s, p) = bumpy_sphere(2, 0.1, 4);
    // MaxIterations.
    let r = flow(s.clone(), p.clone(), FlowOptions { max_iterations: 3, ..base() }, None).unwrap();
    assert_eq!((r.report.status, r.report.iterations), (Status::MaxIterations, 3));
    // TimeLimit (zero: before the first step).
    let o = FlowOptions { time_limit: Some(Duration::ZERO), ..base() };
    let r = flow(s.clone(), p.clone(), o, None).unwrap();
    assert_eq!((r.report.status, r.report.steps.len()), (Status::TimeLimit, 0));
    assert!(r.report.final_energy.is_finite());
    // Converged: a loose energy stall.
    let o = FlowOptions {
        stop: StopRule::EnergyStall { relative: 0.5, window: 2 },
        ..base()
    };
    let r = flow(s.clone(), p.clone(), o, None).unwrap();
    assert_eq!(r.report.status, Status::Converged);
    assert!(r.report.iterations >= 2 && r.report.iterations < 50);
    // Converged: the Sobolev norm rule.
    let o = FlowOptions { stop: StopRule::SobolevNorm { eps: 0.9 }, max_iterations: 200, ..base() };
    let r = flow(s.clone(), p.clone(), o, None).unwrap();
    assert_eq!(r.report.status, Status::Converged);
    // Not converged while a schedule is moving its target.
    let o = FlowOptions {
        stop: StopRule::EnergyStall { relative: 0.5, window: 1 },
        schedules: vec![Schedule { constraint: 1, factor: 1.05, offset: 0.0, steps: 6 }],
        ..base()
    };
    let r = flow(s.clone(), p.clone(), o, None).unwrap();
    assert_eq!(r.report.status, Status::Converged);
    assert!(r.report.iterations >= 6, "converged at {} while the schedule moved", r.report.iterations);
    // Error: a step error ends the run with `Status::Error` (GMRES with one iteration and a
    // tiny tolerance diverges, `SolverDiverged`), `run` itself succeeds.
    let mut h = HierarchicalOptions::default();
    h.gmres.max_iterations = 1;
    h.gmres.tolerance = 1e-12;
    let o = FlowOptions { backend: Backend::Hierarchical(h), ..base() };
    let r = flow(s, p, o, None).unwrap();
    match &r.report.status {
        Status::Error(RsError::SolverDiverged { .. }) => {}
        other => panic!("expected a solver error, got {other:?}"),
    }
    assert!(r.report.error().is_some());
    assert!(r.report.steps.is_empty());
}

#[test]
fn stuck_steps_end_the_run() {
    let (s, p) = bumpy_sphere(2, 0.1, 5);
    // A huge first step and at most two halvings: every line search fails.
    let o = FlowOptions {
        line_search: LineSearch {
            start: StartRule::Absolute { tau: 1e6 },
            max_halvings: Some(2),
            // Without CCD (it would cap the huge start).
            ccd: None,
            ..LineSearch::default()
        },
        remesh: None,
        stuck_limit: 3,
        max_iterations: 100,
        ..base()
    };
    let r = flow(s, p.clone(), o, None).unwrap();
    assert_eq!(r.report.status, Status::Stuck);
    assert_eq!(r.report.steps.len(), 3);
    for st in &r.report.steps {
        assert_eq!(st.status, StepStatus::Stuck);
        assert_eq!((st.tau, st.halvings, st.trials.len()), (0.0, 3, 3));
        assert_eq!(st.energy_trial, st.energy_before);
    }
    // The positions moved only by the projection (nothing to project: unchanged up to round-off).
    let d = r.positions.iter().zip(&p).flat_map(|(a, b)| (0..3).map(move |k| (a[k] - b[k]).abs()));
    assert!(d.fold(0.0f64, f64::max) < 1e-14);
}

#[test]
fn setup_errors() {
    let (s, p) = bumpy_sphere(1, 0.1, 6);
    let new = |o: FlowOptions| Session::new(s.clone(), p.clone(), o);
    // No translation fixed.
    let e = new(FlowOptions { constraints: vec![], ..base() }).unwrap_err();
    assert!(matches!(e, RsError::TranslationNotFixed { component: 0 }), "{e}");
    // A schedule on a barycenter.
    let o = FlowOptions { schedules: vec![Schedule { constraint: 0, factor: 2.0, offset: 0.0, steps: 3 }], ..base() };
    assert!(matches!(new(o).unwrap_err(), RsError::InvalidOption { name: "schedule.constraint", .. }));
    // A schedule out of range.
    let o = FlowOptions { schedules: vec![Schedule { constraint: 7, factor: 2.0, offset: 0.0, steps: 3 }], ..base() };
    assert!(matches!(new(o).unwrap_err(), RsError::ConstraintOutOfRange { constraint: 7, .. }));
    // Bad option values.
    for o in [
        FlowOptions { tp_weight: 0.0, ..base() },
        FlowOptions { stuck_limit: 0, ..base() },
        FlowOptions { threads: Some(0), ..base() },
        FlowOptions { stop: StopRule::EnergyStall { relative: 1e-3, window: 0 }, ..base() },
        FlowOptions {
            line_search: LineSearch { start: StartRule::MaxDisplacement { lengths: -1.0 }, ..LineSearch::default() },
            ..base()
        },
    ] {
        assert!(matches!(new(o).unwrap_err(), RsError::InvalidOption { .. }));
    }
    // The comparison optimizers need a soft problem.
    let o = FlowOptions { method: Method::Aqp { kappa: 100.0 }, ..base() };
    assert!(matches!(new(o).unwrap_err(), RsError::MethodNeedsSoftProblem { .. }));
    let o = FlowOptions { method: Method::Aqp { kappa: 100.0 }, constraints: vec![Constraint::Barycenter { component: None, target: None }], ..base() };
    assert!(new(o).is_ok());
    // Pins out of range.
    let o = FlowOptions { constraints: vec![Constraint::Pins { vertices: vec![10_000], targets: None }], ..base() };
    assert!(matches!(new(o).unwrap_err(), RsError::ConstraintOutOfRange { constraint: 0, index: 10_000 }));
}

#[test]
fn accessors() {
    let (s, p) = bumpy_sphere(2, 0.1, 8);
    let o = FlowOptions {
        constraints: vec![Constraint::Pins { vertices: vec![3, 1, 3], targets: None }],
        remesh: Some(RemeshOptions::default()),
        ..base()
    };
    let mut sess = Session::new(s.clone(), p.clone(), o).unwrap();
    assert_eq!(sess.pinned(), &[3, 1]);
    assert_eq!(sess.status(), Status::Running);
    assert!(sess.energy().is_finite());
    assert!(sess.mean_edge_length0() > 0.0);
    assert_eq!(sess.face_origin().len(), s.n_faces());
    assert_eq!(sess.options().constraints.len(), 1);
    sess.step().unwrap();
    assert_eq!(sess.iteration(), 1);
    assert_eq!(sess.pinned().len(), 2);
    for &v in sess.pinned() {
        assert!((v as usize) < sess.surface().n_vertices());
    }
}
