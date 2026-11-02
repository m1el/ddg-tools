//! The flow (`SPEC.md` §7, §13.4): energy decrease, constraints held (also under growth schedules
//! and moving pins), a torus of revolution stays one, a perturbed sphere becomes round (dense and
//! hierarchical back ends), two linked components do not intersect, remeshing every step keeps
//! manifoldness and pins, the comparison metrics, the other energies, penalties and obstacles.
//!
//! The oracle's `step` fixtures are replayed in `tests/oracle.rs` (`step_fixtures_are_reproduced`).

use std::sync::Arc;

use ddg_repulsive_curves::obstacle::Obstacle;
use ddg_repulsive_surfaces::ccd::{intersects, min_distance};
use ddg_repulsive_surfaces::constraints::{BarycenterMode, Constraint, ConstraintSet};
use ddg_repulsive_surfaces::generate::primitives as generate;
use ddg_repulsive_surfaces::generate::{Mesh, Rng};
use ddg_repulsive_surfaces::geometry::{FaceGeometry, vertex_areas};
use ddg_repulsive_surfaces::solve::iterative::{Backend, HierarchicalOptions};
use ddg_repulsive_surfaces::{
    BhOptions, CcdOptions, Energy, EnergyEval, FlowOptions, L2Kind, Metric, ObstacleKind,
    Penalty, PinSchedule, RemeshOptions, Schedule, ScheduleTiming, Session, StaticSet,
    StepReport, StepStatus, StopRule, Surface,
};

// ----- meshes -----

/// An icosphere (`20·4^level` faces) of radius 1 with seeded radial noise of amplitude `a`.
fn bumpy_sphere(level: u32, a: f64, seed: u64) -> (Surface, Vec<[f64; 3]>) {
    let m = generate::sphere(level, 1.0, [0.0; 3]).unwrap();
    let mut rng = Rng::new(seed);
    let p: Vec<[f64; 3]> = m
        .positions
        .iter()
        .map(|q| {
            let s = 1.0 + a * rng.uniform(-1.0, 1.0) + 0.5 * a * (3.0 * q[0]).sin();
            [q[0] * s, q[1] * s * 0.9, q[2] * s * 1.1]
        })
        .collect();
    (Surface::new(p.len(), m.faces.clone()).unwrap(), p)
}

fn torus(nu: usize, nv: usize, major: f64, minor: f64) -> Mesh {
    generate::torus(nu, nv, major, minor).unwrap()
}

fn surface_of(m: &Mesh) -> (Surface, Vec<[f64; 3]>) {
    (m.surface().unwrap(), m.positions.clone())
}

/// Two linked tori (a Hopf link): one around `z`, one around `y` through the first's tube.
fn hopf(nu: usize, nv: usize) -> (Surface, Vec<[f64; 3]>) {
    let a = torus(nu, nv, 1.0, 0.25);
    let mut b = torus(nu, nv, 1.0, 0.25);
    // Rotate b by 90° about x (its axis becomes y) and shift it by +1 in x.
    for p in &mut b.positions {
        *p = [p[0] + 1.0, -p[2], p[1]];
    }
    let m = a.merged(&b);
    surface_of(&m)
}

// ----- helpers -----

fn dense_options() -> FlowOptions {
    FlowOptions {
        energy_eval: EnergyEval::Exact,
        backend: Backend::Dense,
        remesh: None,
        stop: StopRule::Off,
        threads: Some(4),
        ..FlowOptions::default()
    }
}

fn area(s: &Surface, p: &[[f64; 3]]) -> f64 {
    FaceGeometry::new(s, p).unwrap().areas.iter().sum()
}

fn volume(s: &Surface, p: &[[f64; 3]]) -> f64 {
    s.faces()
        .iter()
        .map(|t| {
            let [a, b, c] = t.map(|i| p[i as usize]);
            (a[0] * (b[1] * c[2] - b[2] * c[1]) - a[1] * (b[0] * c[2] - b[2] * c[0])
                + a[2] * (b[0] * c[1] - b[1] * c[0]))
                / 6.0
        })
        .sum()
}

fn barycenter(s: &Surface, p: &[[f64; 3]]) -> [f64; 3] {
    let g = FaceGeometry::new(s, p).unwrap();
    let va = vertex_areas(s, &g.areas);
    let tot: f64 = va.iter().sum();
    let mut b = [0.0; 3];
    for (v, q) in p.iter().enumerate() {
        for k in 0..3 {
            b[k] += va[v] * q[k] / tot;
        }
    }
    b
}

/// `std/mean` of the distance to the barycenter.
fn roundness(s: &Surface, p: &[[f64; 3]]) -> f64 {
    let b = barycenter(s, p);
    let r: Vec<f64> = p
        .iter()
        .map(|q| ((q[0] - b[0]).powi(2) + (q[1] - b[1]).powi(2) + (q[2] - b[2]).powi(2)).sqrt())
        .collect();
    let m = r.iter().sum::<f64>() / r.len() as f64;
    (r.iter().map(|x| (x - m).powi(2)).sum::<f64>() / r.len() as f64).sqrt() / m
}

fn assert_accepted_steps_decrease(reps: &[StepReport]) {
    for r in reps {
        if r.status == StepStatus::Accepted {
            assert!(
                r.energy_trial < r.energy_before,
                "step {}: F {} -> {}",
                r.iteration,
                r.energy_before,
                r.energy_trial
            );
        }
    }
}

fn run(s: &mut Session, n: usize) -> Vec<StepReport> {
    (0..n).map(|_| s.step().unwrap()).collect()
}

// ----- tests -----

#[test]
fn energy_decreases_and_constraints_hold() {
    let (s, p) = bumpy_sphere(2, 0.15, 1);
    let a0 = area(&s, &p);
    let b0 = barycenter(&s, &p);
    let o = FlowOptions {
        constraints: vec![
            Constraint::Barycenter { component: None, target: None },
            Constraint::Area { component: None, target: None },
        ],
        ..dense_options()
    };
    let mut sess = Session::new(s, p, o).unwrap();
    let reps = run(&mut sess, 12);
    assert_accepted_steps_decrease(&reps);
    assert!(reps.iter().all(|r| r.status == StepStatus::Accepted));
    // Without remeshing and with constraints that keep F's value through the projection up to
    // the Newton correction, F decreases over the run.
    assert!(reps.last().unwrap().energy_after < 0.7 * reps[0].energy_before);
    let (s, p) = (sess.surface(), sess.positions());
    assert!((area(s, p) - a0).abs() <= 1e-8 * a0, "area {} vs {a0}", area(s, p));
    let b = barycenter(s, p);
    for k in 0..3 {
        assert!((b[k] - b0[k]).abs() < 1e-12, "barycenter {b:?} vs {b0:?}");
    }
    for r in &reps {
        assert!(r.max_violation <= 1e-8 * a0, "step {}: violation {}", r.iteration, r.max_violation);
    }
}

#[test]
fn area_and_volume_on_a_torus_hold() {
    let m = torus(24, 12, 1.0, 0.4);
    let (s, mut p) = surface_of(&m);
    for q in &mut p {
        q[2] *= 1.3;
    }
    let (a0, v0) = (area(&s, &p), volume(&s, &p));
    let o = FlowOptions {
        constraints: vec![
            Constraint::Barycenter { component: None, target: None },
            Constraint::Area { component: None, target: None },
            Constraint::Volume { component: None, target: None },
        ],
        ..dense_options()
    };
    let mut sess = Session::new(s, p, o).unwrap();
    let reps = run(&mut sess, 8);
    assert_accepted_steps_decrease(&reps);
    let (s, p) = (sess.surface(), sess.positions());
    assert!((area(s, p) / a0 - 1.0).abs() < 1e-8);
    assert!((volume(s, p) / v0 - 1.0).abs() < 1e-8);
}

#[test]
fn perturbed_sphere_becomes_round_dense_and_hierarchical() {
    for (name, backend, eval) in [
        ("dense", Backend::Dense, EnergyEval::Exact),
        (
            "hierarchical",
            Backend::Hierarchical(HierarchicalOptions::default()),
            EnergyEval::BarnesHut(BhOptions::default()),
        ),
    ] {
        let (s, p) = bumpy_sphere(2, 0.12, 7);
        let r0 = roundness(&s, &p);
        let o = FlowOptions {
            backend,
            energy_eval: eval,
            constraints: vec![
                Constraint::Barycenter { component: None, target: None },
                Constraint::Area { component: None, target: None },
            ],
            ..dense_options()
        };
        let mut sess = Session::new(s, p, o).unwrap();
        let reps = run(&mut sess, 25);
        assert_accepted_steps_decrease(&reps);
        let r1 = roundness(sess.surface(), sess.positions());
        println!("{name}: roundness {r0:.4} -> {r1:.4}, solver {}", reps[0].linear.solver);
        assert!(r1 < 0.25 * r0, "{name}: roundness {r0} -> {r1}");
        if name == "hierarchical" {
            assert!(reps[0].linear.solver.contains("gmres"));
            assert!(reps.iter().all(|r| r.linear.converged()));
        }
    }
}

#[test]
fn torus_of_revolution_stays_one() {
    // Exact energy, dense solve, no remeshing: the symmetry is kept to round-off. Barnes–Hut,
    // the default back end and remeshing every step: kept to the remeshing's noise.
    for (name, eval, backend, remesh, tol) in [
        ("exact", EnergyEval::Exact, Backend::Dense, None, 1e-10),
        (
            "bh+remesh",
            EnergyEval::BarnesHut(BhOptions::default()),
            Backend::default(),
            Some(RemeshOptions::default()),
            1e-2,
        ),
    ] {
        let m = torus(48, 16, 1.0, 0.35);
        let (s, p) = surface_of(&m);
        let o = FlowOptions {
            constraints: vec![
                Constraint::Barycenter { component: None, target: None },
                Constraint::Area { component: None, target: None },
            ],
            remesh,
            energy_eval: eval,
            backend,
            ..dense_options()
        };
        let mut sess = Session::new(s, p, o).unwrap();
        let reps = run(&mut sess, 10);
        assert_accepted_steps_decrease(&reps);
        assert!(reps.last().unwrap().energy_after < 0.7 * reps[0].energy_before, "{name}");
        let (s, p) = (sess.surface(), sess.positions());
        // Still a torus: one closed component of Euler characteristic 0, not self-intersecting.
        assert_eq!(s.n_components(), 1);
        assert!(s.boundary_loops().is_empty());
        assert_eq!(s.n_vertices() as i64 - s.edges().len() as i64 + s.n_faces() as i64, 0);
        assert!(!intersects(s, p).unwrap());
        // Still of revolution about z: the angular Fourier modes m = 1..4 of the area-weighted
        // distance from the axis vanish.
        let modes = angular_modes(s, p);
        println!("{name}: angular modes {modes:?}");
        assert!(modes.iter().all(|&c| c < tol), "{name}: not of revolution: {modes:?}");
    }
}

/// `|Σ a_S ρ_S e^{i m φ_S}| / Σ a_S ρ_S` for `m = 1..4` (face barycenters about the
/// barycenter, axis `z`).
fn angular_modes(s: &Surface, p: &[[f64; 3]]) -> [f64; 4] {
    let b = barycenter(s, p);
    let g = FaceGeometry::new(s, p).unwrap();
    let mut out = [0.0; 4];
    for m in 1..=4 {
        let (mut re, mut im, mut w) = (0.0, 0.0, 0.0);
        for (f, q) in g.barycenters.iter().enumerate() {
            let (x, y) = (q[0] - b[0], q[1] - b[1]);
            let (phi, rho) = (y.atan2(x), (x * x + y * y).sqrt());
            re += g.areas[f] * rho * (m as f64 * phi).cos();
            im += g.areas[f] * rho * (m as f64 * phi).sin();
            w += g.areas[f] * rho;
        }
        out[m - 1] = (re * re + im * im).sqrt() / w;
    }
    out
}

#[test]
fn linked_components_do_not_intersect() {
    let (s, p) = hopf(24, 8);
    assert!(!intersects(&s, &p).unwrap());
    assert_eq!(s.n_components(), 2);
    let o = FlowOptions {
        constraints: vec![
            Constraint::Barycenter { component: None, target: None },
            Constraint::Area { component: None, target: None },
        ],
        remesh: Some(RemeshOptions::default()),
        energy_eval: EnergyEval::BarnesHut(BhOptions::default()),
        backend: Backend::default(),
        ..dense_options()
    };
    let mut sess = Session::new(s, p, o).unwrap();
    for _ in 0..15 {
        let r = sess.step().unwrap();
        assert!(r.status == StepStatus::Accepted || r.iteration > 1);
        assert_eq!(sess.surface().n_components(), 2);
        assert!(
            !intersects(sess.surface(), sess.positions()).unwrap(),
            "step {}: the components intersect",
            r.iteration
        );
    }
    let d = min_distance(sess.surface(), sess.positions()).unwrap();
    assert!(d > 0.0);
}

#[test]
fn linked_components_with_ccd_keep_the_gap() {
    let (s, p) = hopf(24, 8);
    let ccd = CcdOptions::default();
    let o = FlowOptions {
        constraints: vec![Constraint::Barycenter { component: None, target: None }],
        remesh: Some(RemeshOptions::default()),
        energy_eval: EnergyEval::BarnesHut(BhOptions::default()),
        backend: Backend::default(),
        line_search: ddg_repulsive_surfaces::LineSearch {
            ccd: Some(ccd),
            // Large first steps so that the CCD cap matters.
            start: ddg_repulsive_surfaces::StartRule::MaxDisplacement { lengths: 4.0 },
            ..Default::default()
        },
        ..dense_options()
    };
    let mut sess = Session::new(s, p, o).unwrap();
    let reps = run(&mut sess, 6);
    assert!(reps.iter().all(|r| r.tau_ccd.is_some() && r.ccd.is_some()));
    assert!(!intersects(sess.surface(), sess.positions()).unwrap());
}

#[test]
fn remeshing_every_step_keeps_manifoldness_and_pins() {
    let m = torus(20, 10, 1.0, 0.35);
    let (s, mut p) = surface_of(&m);
    // Stretch it so that remeshing has work to do.
    for q in &mut p {
        q[0] *= 1.6;
    }
    let pins = vec![0u32, 37, 101];
    let targets: Vec<[f64; 3]> = pins.iter().map(|&v| p[v as usize]).collect();
    let o = FlowOptions {
        constraints: vec![Constraint::Pins { vertices: pins.clone(), targets: None }],
        remesh: Some(RemeshOptions::default()),
        energy_eval: EnergyEval::BarnesHut(BhOptions::default()),
        backend: Backend::default(),
        ..dense_options()
    };
    let mut sess = Session::new(s, p, o).unwrap();
    let mut ops = 0;
    for _ in 0..8 {
        let r = sess.step().unwrap();
        let rm = r.remesh.as_ref().expect("remeshed");
        ops += rm.splits + rm.collapses + rm.flips;
        // Valid manifold: rebuild from the faces and validate the positions.
        let s2 = Surface::new(sess.surface().n_vertices(), sess.surface().faces().to_vec()).unwrap();
        s2.validate_positions(sess.positions()).unwrap();
        assert_eq!(sess.pinned().len(), 3);
        for (k, &v) in sess.pinned().iter().enumerate() {
            assert_eq!(sess.positions()[v as usize], targets[k], "pin {k} moved");
        }
        assert_eq!(sess.face_origin().len(), sess.surface().n_faces());
    }
    assert!(ops > 0, "remeshing did nothing");
}

#[test]
fn growth_schedule_and_moving_pins() {
    // Area growth ×1.3 over 4 steps on a sphere; the targets are hit exactly at each step.
    for timing in [ScheduleTiming::BeforeProjection, ScheduleTiming::AfterProjection] {
        let (s, p) = bumpy_sphere(2, 0.05, 3);
        let a0 = area(&s, &p);
        let o = FlowOptions {
            constraints: vec![
                Constraint::Barycenter { component: None, target: None },
                Constraint::Area { component: None, target: None },
            ],
            schedules: vec![Schedule { constraint: 1, factor: 1.3, offset: 0.0, steps: 4 }],
            schedule_timing: timing,
            ..dense_options()
        };
        let mut sess = Session::new(s, p, o).unwrap();
        for k in 1..=6usize {
            let r = sess.step().unwrap();
            let kk = match timing {
                ScheduleTiming::BeforeProjection => k,
                ScheduleTiming::AfterProjection => k - 1,
            };
            let want = a0 * (1.0 + 0.3 * (kk.min(4) as f64) / 4.0);
            let got = area(sess.surface(), sess.positions());
            assert!(
                (got / want - 1.0).abs() < 1e-7,
                "{timing:?} step {k}: area {got} vs {want} ({:?})",
                r.projection_status
            );
        }
    }
    // Moving pins: a patch of a torus, three pins moved by (0, 0, 0.2) over 3 steps.
    let m = torus(20, 10, 1.0, 0.35);
    let (s, p) = surface_of(&m);
    let pins = vec![0u32, 1, 2];
    let p0: Vec<[f64; 3]> = pins.iter().map(|&v| p[v as usize]).collect();
    let o = FlowOptions {
        constraints: vec![Constraint::Pins { vertices: pins.clone(), targets: None }],
        pin_schedules: vec![PinSchedule { constraint: 0, displacement: vec![[0.0, 0.0, 0.2]; 3], steps: 3 }],
        ..dense_options()
    };
    let mut sess = Session::new(s, p, o).unwrap();
    for k in 1..=4usize {
        sess.step().unwrap();
        let t = 0.2 * (k.min(3) as f64) / 3.0;
        for (i, &v) in pins.iter().enumerate() {
            let q = sess.positions()[v as usize];
            assert_eq!([q[0], q[1]], [p0[i][0], p0[i][1]]);
            assert!((q[2] - (p0[i][2] + t)).abs() < 1e-15, "pin {i} at step {k}");
        }
    }
}

#[test]
fn comparison_metrics_descend_with_all_constraints() {
    for metric in [Metric::L2(L2Kind::Mass), Metric::L2(L2Kind::Identity), Metric::H1, Metric::H2] {
        let (s, p) = bumpy_sphere(2, 0.1, 11);
        let a0 = area(&s, &p);
        let o = FlowOptions {
            metric,
            constraints: vec![
                Constraint::Barycenter { component: None, target: None },
                Constraint::Area { component: None, target: None },
            ],
            ..dense_options()
        };
        let mut sess = Session::new(s, p, o).unwrap();
        let reps = run(&mut sess, 5);
        assert_accepted_steps_decrease(&reps);
        assert!(reps.iter().any(|r| r.status == StepStatus::Accepted), "{metric:?}: no step");
        assert!(reps[0].linear.solver.contains("sparse"));
        assert!(reps.iter().all(|r| r.linear.converged()), "{metric:?}: {:?}", reps[0].linear);
        let got = area(sess.surface(), sess.positions());
        assert!((got / a0 - 1.0).abs() < 1e-8, "{metric:?}");
    }
}

#[test]
fn other_energies_and_penalties() {
    // Coulomb energy (dense) and the Willmore energy with the H² metric.
    for (energy, metric) in [
        (Energy::Coulomb { alpha: 1.0 }, Metric::default()),
        (Energy::Willmore, Metric::H2),
    ] {
        let (s, p) = bumpy_sphere(1, 0.1, 5);
        let o = FlowOptions {
            energy,
            metric,
            constraints: vec![
                Constraint::Barycenter { component: None, target: None },
                Constraint::Area { component: None, target: None },
            ],
            ..dense_options()
        };
        let mut sess = Session::new(s, p, o).unwrap();
        let reps = run(&mut sess, 4);
        assert_accepted_steps_decrease(&reps);
        assert!(reps.iter().any(|r| r.status == StepStatus::Accepted), "{energy:?}");
    }
    // Penalties with remeshing: squared error (targets follow the vertex map), a Willmore
    // penalty (its term in the metric), a tangent-point mesh obstacle (diagonal in the metric,
    // frozen cross lists in the line search).
    let (s, p) = bumpy_sphere(2, 0.1, 9);
    let ob = generate::sphere(1, 0.3, [2.0, 0.0, 0.0]).unwrap();
    let obstacle = Obstacle::new(ob.positions.clone(), ob.faces.clone()).unwrap();
    let set = Arc::new(StaticSet::from_obstacle(&obstacle, 8).unwrap());
    let o = FlowOptions {
        constraints: vec![Constraint::Barycenter { component: None, target: None }],
        penalties: vec![
            Penalty::SquaredError { weight: 0.1, targets: p.clone() },
            Penalty::Willmore { weight: 0.01 },
            Penalty::MeshObstacle { weight: 1.0, obstacle: set, kind: ObstacleKind::TangentPoint },
            Penalty::AreaDeviation { weight: 10.0, target: area(&s, &p) },
        ],
        remesh: Some(RemeshOptions { split: 1.2, collapse: 0.55, ..RemeshOptions::default() }),
        energy_eval: EnergyEval::BarnesHut(BhOptions::default()),
        ccd_obstacles: vec![Arc::new(obstacle)],
        ..dense_options()
    };
    let mut sess = Session::new(s, p, o).unwrap();
    let reps = run(&mut sess, 4);
    assert_accepted_steps_decrease(&reps);
    for r in &reps {
        assert_eq!(r.penalties_after.len(), 4);
        assert!(r.penalties_after.iter().all(|v| v.is_finite()));
    }
    match &sess.penalties()[0] {
        Penalty::SquaredError { targets, .. } => assert_eq!(targets.len(), sess.surface().n_vertices()),
        _ => unreachable!(),
    }
}

#[test]
fn barycenter_shift_and_reference_compat_run() {
    let (s, p) = bumpy_sphere(2, 0.1, 13);
    let o = FlowOptions {
        constraints: vec![
            Constraint::Barycenter { component: None, target: None },
            Constraint::Area { component: None, target: None },
        ],
        barycenter_shift: true,
        ..dense_options()
    };
    let mut sess = Session::new(s.clone(), p.clone(), o).unwrap();
    let reps = run(&mut sess, 4);
    assert_accepted_steps_decrease(&reps);
    let o = FlowOptions {
        constraints: vec![
            Constraint::Barycenter { component: None, target: None },
            Constraint::Area { component: None, target: None },
        ],
        reference_compat: true,
        ..dense_options()
    };
    let mut sess = Session::new(s, p, o).unwrap();
    let e = sess.effective_options().clone();
    assert_eq!(e.barycenter_mode, Some(BarycenterMode::Global));
    assert_eq!(e.newton_steps, ddg_repulsive_surfaces::NewtonSteps::Fixed(1));
    let reps = run(&mut sess, 3);
    assert!(reps.iter().any(|r| r.status == StepStatus::Accepted));
    assert!(reps.iter().all(|r| r.newton_steps <= 1));
}

#[test]
fn global_and_per_component_barycenters() {
    // Auto = Global with B⁰ on (one global barycenter for both tori), PerComponent otherwise.
    let (s, p) = hopf(16, 6);
    for (metric, mode) in [
        (Metric::Sobolev { low_order: true }, BarycenterMode::Global),
        (Metric::Sobolev { low_order: false }, BarycenterMode::PerComponent),
    ] {
        let o = FlowOptions { metric, ..dense_options() };
        let mut sess = Session::new(s.clone(), p.clone(), o).unwrap();
        assert_eq!(sess.barycenter_mode(), mode);
        let r = sess.step().unwrap();
        let set = ConstraintSet::new(sess.surface(), sess.positions(), &sess.constraints(), mode).unwrap();
        let phi = set.values(sess.positions()).unwrap();
        assert!(phi.iter().all(|x| x.abs() < 1e-12), "{mode:?}: {phi:?}");
        assert_eq!(set.n_rows(), if mode == BarycenterMode::Global { 3 } else { 6 });
        assert!(r.energy_trial < r.energy_before);
    }
}

// ----- CCD on inputs that already intersect (SPEC 7.3, the intersection band) -----

const KLEIN_SCENE: &str =
    "/var/data/discrete-geometry/RepulsiveSurfaces/repulsive-surfaces/scenes/KleinBottle/scene.txt";

/// Two tori that cross each other (an immersion-like input).
fn crossing_tori() -> (Surface, Vec<[f64; 3]>) {
    let a = torus(24, 8, 1.0, 0.3);
    let mut b = torus(24, 8, 1.0, 0.3);
    // b is a copy of a rotated by 90 degrees about x and shifted by 0.5 in x: its tube passes
    // through a's tube.
    for p in &mut b.positions {
        *p = [p[0] + 0.5, -p[2], p[1]];
    }
    surface_of(&a.merged(&b))
}

fn ccd_flow_options(threads: usize) -> FlowOptions {
    FlowOptions {
        constraints: vec![Constraint::Barycenter { component: None, target: None }],
        remesh: Some(RemeshOptions::default()),
        energy_eval: EnergyEval::BarnesHut(BhOptions::default()),
        backend: Backend::default(),
        threads: Some(threads),
        line_search: ddg_repulsive_surfaces::LineSearch {
            ccd: Some(CcdOptions::default()),
            start: ddg_repulsive_surfaces::StartRule::MaxDisplacement { lengths: 2.0 },
            ..Default::default()
        },
        ..dense_options()
    }
}

fn bits(p: &[[f64; 3]]) -> Vec<u64> {
    p.iter().flatten().map(|x| x.to_bits()).collect()
}

#[test]
fn ccd_on_an_intersecting_input_exempts_and_is_deterministic() {
    let (s, p) = crossing_tori();
    assert!(intersects(&s, &p).unwrap());
    let mut ends = Vec::new();
    for threads in [1, 3] {
        let mut sess = Session::new(s.clone(), p.clone(), ccd_flow_options(threads)).unwrap();
        let reps = run(&mut sess, 5);
        let exempt: Vec<u64> = reps.iter().map(|r| r.ccd.as_ref().unwrap().exempt_pairs).collect();
        assert!(exempt[0] > 0, "{exempt:?}");
        assert!(reps.iter().all(|r| r.status == StepStatus::Accepted));
        // The flow moves: the energy goes down.
        assert!(reps.last().unwrap().energy_after < reps[0].energy_before);
        ends.push((bits(sess.positions()), exempt));
    }
    assert_eq!(ends[0], ends[1], "trajectory differs between 1 and 3 threads");
}

#[test]
fn ccd_on_an_embedded_input_has_an_empty_exemption_and_is_deterministic() {
    let (s, p) = hopf(24, 8);
    assert!(!intersects(&s, &p).unwrap());
    let mut ends = Vec::new();
    for threads in [1, 3] {
        let mut sess = Session::new(s.clone(), p.clone(), ccd_flow_options(threads)).unwrap();
        let reps = run(&mut sess, 6);
        // Embedded at every step start: no pair is ever exempt.
        assert!(reps.iter().all(|r| r.ccd.as_ref().unwrap().exempt_pairs == 0));
        assert!(!intersects(sess.surface(), sess.positions()).unwrap());
        ends.push(bits(sess.positions()));
    }
    assert_eq!(ends[0], ends[1]);
}

/// The Klein bottle scene (an immersion, exponents (3.5, 7)): the CCD stays on, and after 200
/// steps the objective is comparable to the run without CCD (SPEC 16: the forced CCD without the
/// rule reached 631 against 156). Release build: `cargo test --release -p ddg-repulsive-surfaces
/// --test flow -- --ignored klein`.
#[test]
#[ignore = "200 steps of the Klein bottle; run in release"]
fn klein_bottle_with_ccd_reaches_the_objective_of_the_run_without() {
    let spec = ddg_repulsive_surfaces::SceneSpec::load(KLEIN_SCENE).unwrap();
    let base = ddg_repulsive_surfaces::scene_options(&spec).unwrap();
    assert!(base.line_search.ccd.is_some(), "CCD is on for subcritical exponents too");
    let mut finals = Vec::new();
    for ccd_on in [true, false] {
        let mut o = base.clone();
        o.max_iterations = 200;
        o.stop = StopRule::Off;
        o.threads = Some(8);
        if !ccd_on {
            o.line_search.ccd = None;
        }
        let mut sess = Session::new(spec.surface.clone(), spec.positions.clone(), o).unwrap();
        let reps = run(&mut sess, 200);
        let f = reps.last().unwrap().energy_after;
        let exempt: u64 = reps.iter().map(|r| r.ccd.as_ref().map_or(0, |c| c.exempt_pairs)).sum();
        eprintln!("klein ccd_on={ccd_on}: F(200) = {f:.3}, exempt pairs summed = {exempt}");
        finals.push(f);
    }
    // The immersion is free to move with the rule (631 vs 156 without it).
    assert!(finals[0] < 2.0 * finals[1], "{finals:?}");
}

/// Regression (WP18b defect 1): the blocky Klein bottle in `reference_compat` with splits and
/// collapses off, state after step 299 (the flow had shrunk an edge to 7e-4 and the final flip
/// of the step made two needle faces: min angle 0.1 degrees, `F` 24.7 -> 4.4e9). The
/// printed-preconditioner left GMRES ends at relative residual 1.466e3 there (above 1e-1), which
/// used to abort the run with `SolverDiverged`. The step now recovers (right GMRES with the
/// augmented preconditioner, else a dense LU), says so in `linear.fallbacks`, and the flow goes
/// on and brings `F` down. Release build: `cargo test --release -p ddg-repulsive-surfaces
/// --test flow -- solver_fallback`.
#[test]
fn solver_fallback_recovers_the_klein_bottle_needle_state() {
    let dir = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures/solver_fallback");
    let spec = ddg_repulsive_surfaces::SceneSpec::load(dir.join("scene.txt")).unwrap();
    let mut o = ddg_repulsive_surfaces::scene_options(&spec).unwrap();
    o.max_iterations = 4;
    o.stop = StopRule::Off;
    o.threads = Some(2);
    o.reference_compat = true;
    o.line_search.ccd = None;
    let r = o.remesh.as_mut().expect("the scene remeshes");
    r.enabled_ops.split = false;
    r.enabled_ops.collapse = false;
    let mut sess = Session::new(spec.surface.clone(), spec.positions.clone(), o).unwrap();
    let f0 = sess.energy();
    assert!(f0 > 1e9, "the needle state: F = {f0:e}");
    let reps = run(&mut sess, 4);
    assert!(!reps[0].linear.fallbacks.is_empty(), "step 1 used a fallback: {:?}", reps[0].linear);
    assert!(reps[0].linear.fallbacks[0].contains("diverged"));
    assert!(reps[0].linear.max_relative_residual() < 1e-1);
    let f = reps.last().unwrap().energy_after;
    assert!(f < 1e3, "F after the recovery steps: {f}");
}
