//! The back ends of the session (WP3c, `SPEC.md` §8.5, §17h): the thread pools reach the
//! solvers, `Solver::Multigrid` with live subdivision and constraints, the dense fallback,
//! `Backend::Auto`.

use std::sync::Arc;

use ddg_repulsive_curves::{
    Backend, BhOptions, CurveNetwork, EnergyEval, FlowOptions, MgOptions, Session, Solver,
};

type P = [f64; 3];

fn trefoil(n: usize) -> (CurveNetwork, Vec<P>) {
    let pos = (0..n)
        .map(|k| {
            let t = std::f64::consts::TAU * k as f64 / n as f64;
            [
                t.sin() + 2.0 * (2.0 * t).sin(),
                t.cos() - 2.0 * (2.0 * t).cos(),
                -(3.0 * t).sin(),
            ]
        })
        .collect();
    let edges = (0..n as u32).map(|k| [k, (k + 1) % n as u32]).collect();
    (CurveNetwork::new(n, edges).unwrap(), pos)
}

fn small_mg() -> MgOptions {
    MgOptions {
        coarsest_vertices: 40,
        ..MgOptions::default()
    }
}

fn options(energy: EnergyEval, solver: Solver, threads: usize) -> FlowOptions {
    FlowOptions {
        backend: Backend::Explicit { energy, solver },
        threads: Some(threads),
        ..FlowOptions::default()
    }
}

/// Three steps with every energy / solver pairing give bitwise identical curves at 1, 3 and 8
/// threads: the solvers run on the session's pools only.
#[test]
fn every_backend_is_bitwise_independent_of_the_thread_count() {
    let (net, pos) = trefoil(200);
    let energies = [
        ("exact", EnergyEval::Exact),
        ("bh", EnergyEval::BarnesHut(BhOptions::default())),
    ];
    for (en, energy) in energies {
        for sn in ["dense", "multigrid"] {
            let mut curves = Vec::new();
            for t in [1, 3, 8] {
                let solver = match sn {
                    "dense" => Solver::Dense,
                    _ => Solver::Multigrid(small_mg()),
                };
                let mut s =
                    Session::new(net.clone(), pos.clone(), options(energy, solver, t)).unwrap();
                for _ in 0..3 {
                    s.step().unwrap();
                }
                curves.push(s.positions().to_vec());
            }
            assert!(
                curves[1] == curves[0] && curves[2] == curves[0],
                "{en} + {sn}: curves differ between 1, 3 and 8 threads"
            );
        }
    }
}

// ------------------------------------------------------------------ Solver::Multigrid

use ddg_repulsive_curves::constraints::{Constraint, TangentSite};
use ddg_repulsive_curves::{LengthSchedule, RcError, SaddleSolution, SaddleSolver, SolverFactory};
use ddg_repulsive_curves::{DirectSolvers, SolverRequest};

fn mg_options(energy: EnergyEval, threads: usize) -> FlowOptions {
    options(energy, Solver::Multigrid(small_mg()), threads)
}

fn max_violation(s: &Session) -> f64 {
    let phi = s.constraints().values(s.network(), s.positions()).unwrap();
    phi.iter().fold(0.0, |m, p| m.max(p.abs()))
}

/// `Solver::Multigrid` follows the dense solver closely with both energy models, reports the
/// multigrid solve and never falls back on a healthy problem.
#[test]
fn multigrid_session_follows_the_dense_session() {
    let (net, pos) = trefoil(160);
    for energy in [EnergyEval::Exact, EnergyEval::BarnesHut(BhOptions::default())] {
        let mut mg = Session::new(net.clone(), pos.clone(), mg_options(energy, 4)).unwrap();
        let mut dense = Session::new(
            net.clone(),
            pos.clone(),
            options(energy, Solver::Dense, 4),
        )
        .unwrap();
        for _ in 0..6 {
            let rm = mg.step().unwrap();
            let rd = dense.step().unwrap();
            assert_eq!(rm.solver_backend, "multigrid");
            assert_eq!(rm.linear.solver, "multigrid");
            assert!(rm.linear.levels >= 3, "{:?}", rm.linear);
            assert!(rm.solver_fallback.is_none());
            assert!(rd.solver_fallback.is_none());
            let rel = (rm.energy_after - rd.energy_after).abs() / rd.energy_after.abs();
            assert!(rel < 1e-3, "energies {} {}", rm.energy_after, rd.energy_after);
        }
        assert!(max_violation(&mg) < 1e-3);
    }
}

/// Live subdivision: the hierarchy is rebuilt for the refined network (more levels), the
/// constraints follow, and no step needs the fallback. Both energy models.
#[test]
fn multigrid_session_survives_live_subdivision() {
    let (net, pos) = trefoil(60);
    for energy in [EnergyEval::Exact, EnergyEval::BarnesHut(BhOptions::default())] {
        let opts = FlowOptions {
            schedule: Some(LengthSchedule {
                factor: 2.5,
                rate: 0.1,
            }),
            max_subdivisions: 1,
            max_iterations: 40,
            ..mg_options(energy, 4)
        };
        let mut s = Session::new(net.clone(), pos.clone(), opts).unwrap();
        let (mut before, mut after) = (0, 0);
        for _ in 0..40 {
            let r = s.step().unwrap();
            assert!(r.solver_fallback.is_none(), "{:?}", r.solver_fallback);
            assert_eq!(s.positions().len(), s.network().n_vertices());
            assert_eq!(r.solver_backend, "multigrid");
            if r.n_vertices == 60 {
                before = before.max(r.linear.levels);
            } else {
                after = after.max(r.linear.levels);
            }
        }
        assert_eq!(s.network().n_vertices(), 120);
        assert!(before >= 2 && after > before, "levels {before} -> {after}");
        assert!(max_violation(&s) < 1e-3);
    }
}

/// Constraints on an open curve: both ends pinned, one end tangent fixed, edge lengths fixed.
/// The coarse levels keep the pins and the tangent site (SPEC §8.4) and the flow agrees with the
/// dense one.
#[test]
fn multigrid_session_handles_pins_and_tangents() {
    let n = 120usize;
    let pos: Vec<P> = (0..n)
        .map(|k| {
            let t = k as f64 / (n - 1) as f64;
            [
                3.0 * t,
                0.6 * (6.0 * t).sin(),
                0.4 * (9.0 * t).cos() * (1.0 - t),
            ]
        })
        .collect();
    let edges = (0..n as u32 - 1).map(|k| [k, k + 1]).collect();
    let net = CurveNetwork::new(n, edges).unwrap();
    let constraints = vec![
        Constraint::EdgeLengths {
            edges: None,
            targets: None,
        },
        Constraint::Point {
            vertex: 0,
            target: None,
        },
        Constraint::Point {
            vertex: n as u32 - 1,
            target: None,
        },
        Constraint::Tangent {
            site: TangentSite::Vertex(0),
            direction: None,
        },
    ];
    let mk = |solver| FlowOptions {
        constraints: constraints.clone(),
        ..options(EnergyEval::Exact, solver, 4)
    };
    let mut mg = Session::new(net.clone(), pos.clone(), mk(Solver::Multigrid(small_mg()))).unwrap();
    let mut dense = Session::new(net, pos.clone(), mk(Solver::Dense)).unwrap();
    // The block cluster tree's product is an approximation of the dense metric (θ = 0.5): the
    // projected gradients of the first step agree to 2 %, later steps are different flows.
    let (mut previous, e0) = (mg.energy(), mg.energy());
    for k in 0..8 {
        let rm = mg.step().unwrap();
        let rd = dense.step().unwrap();
        assert!(rm.solver_fallback.is_none(), "{:?}", rm.solver_fallback);
        assert!(rm.linear.levels >= 2, "{:?}", rm.linear);
        if k == 0 {
            let rel = (rm.gradient_norm - rd.gradient_norm).abs() / rd.gradient_norm;
            assert!(rel < 0.02, "gradient norms {} {}", rm.gradient_norm, rd.gradient_norm);
        }
        assert!(rm.energy_after <= previous);
        previous = rm.energy_after;
    }
    assert!(previous < 0.97 * e0, "{previous} {e0}");
    assert!(max_violation(&mg) < 1e-3);
    assert!(mg.positions()[0] == pos[0] || {
        let d = mg.positions()[0]
            .iter()
            .zip(&pos[0])
            .map(|(a, b)| (a - b).abs())
            .fold(0.0, f64::max);
        d < 1e-3
    });
}

// ------------------------------------------------------------------ fallback

/// Multigrid options that cannot converge in one cycle at a 1e-14 tolerance.
fn failing_mg() -> MgOptions {
    MgOptions {
        tolerance: 1e-14,
        failure_residual: 1e-13,
        max_cycles: 1,
        ..small_mg()
    }
}

/// A multigrid solve that stagnates above `failure_residual` ends in the dense solver of the
/// same step: the step's result equals the dense session's bit for bit, and the report says so.
#[test]
fn stagnating_multigrid_falls_back_to_the_dense_solver() {
    let (net, pos) = trefoil(100);
    let mut mg = Session::new(
        net.clone(),
        pos.clone(),
        options(EnergyEval::Exact, Solver::Multigrid(failing_mg()), 4),
    )
    .unwrap();
    let mut dense = Session::new(net, pos, options(EnergyEval::Exact, Solver::Dense, 4)).unwrap();
    for _ in 0..3 {
        let rm = mg.step().unwrap();
        let rd = dense.step().unwrap();
        let f = rm.solver_fallback.as_ref().expect("the step fell back");
        assert_eq!((f.from, f.to), ("multigrid", "direct"));
        assert!(matches!(f.cause, RcError::MultigridNotConverged { .. }), "{:?}", f.cause);
        assert!(f.dense_solves >= 1);
        assert!(f.cause.to_string().contains("did not converge"));
        assert_eq!(rm.linear.solver, rd.linear.solver);
        assert!(rd.solver_fallback.is_none());
    }
    assert!(mg.positions() == dense.positions());
}

/// Above `dense_fallback_max_vertices` the error goes up and the status is `Error`.
#[test]
fn fallback_is_refused_above_the_size_limit() {
    let (net, pos) = trefoil(100);
    let opts = FlowOptions {
        dense_fallback_max_vertices: 99,
        ..options(EnergyEval::Exact, Solver::Multigrid(failing_mg()), 4)
    };
    let mut s = Session::new(net, pos, opts).unwrap();
    let e = s.step().unwrap_err();
    assert!(matches!(e, RcError::MultigridNotConverged { .. }), "{e:?}");
    assert!(matches!(s.status(), ddg_repulsive_curves::Status::Error(_)));
}

/// A saddle solver that fails with `MultigridDiverged` on the first solve.
struct Diverging(Box<dyn SaddleSolver>, bool);

impl SaddleSolver for Diverging {
    fn n_dofs(&self) -> usize {
        self.0.n_dofs()
    }
    fn n_constraints(&self) -> usize {
        self.0.n_constraints()
    }
    fn solve(&mut self, a: &[f64], b: &[f64]) -> Result<SaddleSolution, RcError> {
        if self.1 {
            return Err(RcError::MultigridDiverged {
                residuals: vec![1e-2, 1.0, 40.0],
            });
        }
        self.0.solve(a, b)
    }
}

struct DivergingFactory;

impl SolverFactory for DivergingFactory {
    fn name(&self) -> &'static str {
        "diverging"
    }
    fn build(&self, r: &SolverRequest<'_>) -> Result<Box<dyn SaddleSolver>, RcError> {
        Ok(Box::new(Diverging(DirectSolvers::new().build(r)?, true)))
    }
}

/// `MultigridDiverged` from any iterative factory also falls back, and the report names the
/// factory.
#[test]
fn diverging_solver_falls_back_and_names_the_cause() {
    let (net, pos) = trefoil(60);
    let mut s = Session::new(
        net,
        pos,
        options(EnergyEval::Exact, Solver::Custom(Arc::new(DivergingFactory)), 2),
    )
    .unwrap();
    let r = s.step().unwrap();
    let f = r.solver_fallback.expect("fell back");
    assert_eq!(f.from, "diverging");
    assert!(matches!(f.cause, RcError::MultigridDiverged { .. }));
    assert!(f.cause.to_string().contains("diverged"));
    assert!(f.dense_solves >= 1);
}
