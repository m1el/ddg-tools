//! Regressions of the solver and robustness fixes of `NOTES-SOLVERFIX.md`: the scale of the
//! augmented preconditioner (hard pins), the coarse correction on per-component translations
//! (many components with a global barycenter), the projection guard (obstacles), the remesher's
//! growth budget and the dense back end's size limit.

#[path = "../examples/apps/apps.rs"]
#[allow(unused_imports)]
mod apps;

use std::sync::Arc;

use ddg_core::generators::icosphere;
use ddg_repulsive_surfaces::constraints::{BarycenterMode, Constraint, ConstraintSet};
use ddg_repulsive_surfaces::energy::BhOptions;
use ddg_repulsive_surfaces::energy::barnes_hut::BhTree;
use ddg_repulsive_surfaces::energy::exact::make_pool;
use ddg_repulsive_surfaces::generate::primitives;
use ddg_repulsive_surfaces::metric::SobolevWeights;
use ddg_repulsive_surfaces::metric::dense::{DENSE_MAX_VERTICES, dense_metric};
use ddg_repulsive_surfaces::remesh::remesh;
use ddg_repulsive_surfaces::solve::iterative::{HierarchicalOptions, MetricExtras, hierarchical_saddle};
use ddg_repulsive_surfaces::solve::precond::{FractionalScale, fractional_symbol_ratio};
use ddg_repulsive_surfaces::solve::SchurComplement;
use ddg_repulsive_surfaces::step::{FlowOptions, Schedule, Session, StopRule};
use ddg_repulsive_surfaces::{
    Backend, Exponents, ProjectionStatus, RemeshOptions, RsError, Surface, TrialOutcome,
};

fn flat(v: &[[f64; 3]]) -> Vec<f64> {
    v.iter().flatten().copied().collect()
}

/// `(iterations of the Schur column, iterations of the direction)` of one constrained direction
/// with the hierarchical back end, or the error.
fn direction_iterations(
    s: &Surface,
    p: &[[f64; 3]],
    cons: &[Constraint],
    mode: BarycenterMode,
    o: &HierarchicalOptions,
) -> Result<Vec<usize>, RsError> {
    let pool = Arc::new(make_pool(Some(4)).unwrap());
    let e = Exponents::default();
    let set = ConstraintSet::new(s, p, cons, mode).unwrap();
    let cq = set.jacobian(p, &set.positional_rows()).unwrap();
    let nl = set.nonlinear_rows();
    let cn = set.jacobian(p, &nl).unwrap();
    let bh = BhTree::build(&pool, s, p, BhOptions::default()).unwrap();
    let g = flat(&bh.energy_and_differential(&pool, &e).unwrap().differential);
    let mut solver = hierarchical_saddle(
        pool.clone(),
        s,
        p,
        &e,
        &SobolevWeights::default(),
        &MetricExtras::default(),
        &cq,
        o,
        None,
    )?;
    let schur = SchurComplement::new(&mut solver, nl, cn)?;
    let d = schur.direction(&mut solver, &g)?;
    let mut its: Vec<usize> = schur.records().iter().map(|r| r.iterations).collect();
    its.push(d.record.iterations);
    assert!(d.record.relative_residual <= 1e-4, "residual {}", d.record.relative_residual);
    Ok(its)
}

#[test]
fn symbol_ratio_matches_the_closed_form() {
    // 4 I(2/3) I(1/3) at s = 5/3, I(t) = π Γ(1−t) / (t 4^t Γ(1+t)).
    let r = fractional_symbol_ratio(&Exponents::default());
    assert!((r - 200.0).abs() < 2.0, "{r}");
    let r2 = fractional_symbol_ratio(&Exponents::new(3.5, 7.0).unwrap());
    assert!(r2 > 100.0 && r2 < 250.0, "{r2}");
}

/// Failure 1 of WP15: GMRES diverged at step 0 on the H of Fig. 22 with its 24 corner pins
/// (3840 faces). The augmented preconditioner's constraint part has eigenvalue 1, its
/// fractional part ≈ 200: the unscaled form diverges, the scaled one converges.
#[test]
fn hard_pins_on_the_h_converge() {
    let (m, _) = apps::h_mesh(3);
    let (s, p) = apps::mesh_pair(&m).unwrap();
    let cons = vec![
        Constraint::Pins { vertices: (0..24).collect(), targets: None },
        Constraint::Volume { component: None, target: None },
    ];
    let o = HierarchicalOptions::default();
    let its = direction_iterations(&s, &p, &cons, BarycenterMode::PerComponent, &o).unwrap();
    eprintln!("scaled: iterations {its:?}");
    assert!(its.iter().all(|&i| i <= 120), "{its:?}");
    let mut unit = o;
    unit.preconditioner.scale = FractionalScale::Unit;
    let r = direction_iterations(&s, &p, &cons, BarycenterMode::PerComponent, &unit);
    eprintln!("unscaled: {r:?}");
    assert!(matches!(r, Err(RsError::SolverDiverged { .. })), "{r:?}");
}

/// Failure 2 of WP15: 14 spheres with one global barycenter and an area row each needed ≈ 200
/// GMRES iterations per solve; the coarse correction on the per-component translations brings
/// it to the single-sphere count.
#[test]
fn many_components_with_a_global_barycenter() {
    let c = primitives::cobblestones(14, 0.5, 0.15, [3.2, 3.2, 0.9], 1, 11).unwrap();
    let (s, p) = apps::mesh_pair(&c.spheres).unwrap();
    let mut cons = vec![Constraint::Barycenter { component: None, target: None }];
    for k in 0..s.n_components() {
        cons.push(Constraint::Area { component: Some(k as u32), target: None });
    }
    let o = HierarchicalOptions::default();
    let its = direction_iterations(&s, &p, &cons, BarycenterMode::Global, &o).unwrap();
    eprintln!("coarse: {its:?}");
    assert!(its.iter().all(|&i| i <= 15), "{its:?}");
    let mut plain = o;
    plain.preconditioner.coarse = false;
    let its0 = direction_iterations(&s, &p, &cons, BarycenterMode::Global, &plain).unwrap();
    eprintln!("no coarse space: {its0:?}");
    assert!(its0.iter().sum::<usize>() > 3 * its.iter().sum::<usize>(), "{its0:?} vs {its:?}");
}

/// A sphere in a slab whose area target doubles in one step (failure 3 of WP15, cobblestones):
/// without the guard the Newton projection pushes it through the slab; with it the surface
/// stays inside (the area lags behind its target, reported).
#[test]
fn projection_respects_obstacles() {
    let sphere = primitives::sphere(2, 0.5, [0.0; 3]).unwrap();
    // A weak slab 0.02 above and below the sphere; the area target doubles in one step.
    let (half, wslab, fac, nst) = (0.52, 0.001, 2.0, 1);
    let slab = primitives::box_mesh([2.0, 2.0, half], 6).unwrap();
    let (s, p) = apps::mesh_pair(&sphere).unwrap();
    let run = |guard: bool| -> (f64, Vec<ddg_repulsive_surfaces::StepReport>) {
        let (pen, ob) = apps::mesh_obstacle(&slab, wslab).unwrap();
        let mut o = FlowOptions {
            max_iterations: 3,
            stop: StopRule::Off,
            backend: Backend::Dense,
            remesh: None,
            threads: Some(2),
            constraints: vec![
                Constraint::Barycenter { component: None, target: None },
                Constraint::Area { component: None, target: None },
            ],
            ..FlowOptions::default()
        };
        o.schedules.push(Schedule { constraint: 1, factor: fac, offset: 0.0, steps: nst });
        o.penalties.push(pen);
        o.ccd_obstacles.push(ob);
        o.line_search.guard_projection = guard;
        let mut sess = Session::new(s.clone(), p.clone(), o).unwrap();
        let mut reps = Vec::new();
        for _ in 0..3 {
            reps.push(sess.step().unwrap());
        }
        let zmax = sess.positions().iter().map(|q| q[2].abs()).fold(0.0, f64::max);
        (zmax, reps)
    };
    let (z_off, _) = run(false);
    let (z_on, reps) = run(true);
    eprintln!("max |z|: guard off {z_off:.4}, on {z_on:.4} (slab half {half})");
    assert!(z_off > half, "the unguarded projection should cross the slab here: {z_off}");
    assert!(z_on < half, "{z_on}");
    assert!(reps.iter().any(|r| r.projection_blocked > 0));
    assert!(
        reps.iter().any(|r| r.projection_status == ProjectionStatus::Blocked
            || r.trials.iter().any(|t| t.outcome == TrialOutcome::ProjectionBlocked))
    );
}

/// The remesher's growth budget: a pass may at most double the face count.
#[test]
fn remesh_growth_is_bounded() {
    let m = icosphere(2).unwrap();
    let s = Surface::new(m.positions.len(), m.faces.clone()).unwrap();
    // Every edge is "long" for a tiny L₀: an unbounded pass would split all of them (×4).
    let l0 = 0.01;
    let free = RemeshOptions { max_growth: f64::INFINITY, rounds: 0, ..RemeshOptions::default() };
    let (out, rep) = remesh(&s, &m.positions, &[], l0, &free, None).unwrap();
    assert_eq!(out.surface.n_faces(), 4 * s.n_faces());
    assert_eq!(rep.splits_deferred, 0);
    let bounded = RemeshOptions { rounds: 0, ..RemeshOptions::default() };
    let (out, rep) = remesh(&s, &m.positions, &[], l0, &bounded, None).unwrap();
    eprintln!("bounded: {} -> {} faces, deferred {}", s.n_faces(), out.surface.n_faces(), rep.splits_deferred);
    assert!(out.surface.n_faces() <= 2 * s.n_faces());
    assert!(rep.splits_deferred > 0);
    assert!(RemeshOptions { max_growth: 0.5, ..RemeshOptions::default() }.validate().is_err());
}

/// The dense back end refuses meshes whose `|V| × |V|` matrices would not fit, with a clear
/// error, before allocating.
#[test]
fn dense_metric_refuses_huge_meshes() {
    let m = icosphere(6).unwrap();
    assert!(m.positions.len() > DENSE_MAX_VERTICES);
    let s = Surface::new(m.positions.len(), m.faces.clone()).unwrap();
    let pool = make_pool(Some(1)).unwrap();
    let r = dense_metric(&pool, &s, &m.positions, &Exponents::default(), &SobolevWeights::default(), None);
    match r {
        Err(RsError::InvalidOption { name, value, .. }) => {
            assert_eq!(name, "dense metric");
            assert!(value.contains("GB"), "{value}");
        }
        other => panic!("expected a size error, got {:?}", other.map(|_| ())),
    }
}
