//! The line-search options of WP14b (`SPEC.md` §17l): the second-order correction of the trials,
//! the Newton polish of the projection, the back-projection switch, and the multigrid thread cap.

use ddg_repulsive_curves::ccd::path_contact;
use ddg_repulsive_curves::{
    Backend, BhOptions, CcdOptions, CcdReport, CurveNetwork, EnergyEval, FlowOptions, MgOptions,
    Session, Solver, StepStatus,
};

type P = [f64; 3];

fn trefoil(n: usize) -> (CurveNetwork, Vec<P>) {
    let pos = (0..n)
        .map(|k| {
            let t = std::f64::consts::TAU * k as f64 / n as f64;
            [
                t.sin() + 2.0 * (2.0 * t).sin() + 0.05 * (7.0 * t).sin(),
                t.cos() - 2.0 * (2.0 * t).cos(),
                -(3.0 * t).sin() + 0.05 * (5.0 * t).cos(),
            ]
        })
        .collect();
    let edges = (0..n as u32).map(|k| [k, (k + 1) % n as u32]).collect();
    (CurveNetwork::new(n, edges).unwrap(), pos)
}

fn options(solver: Solver, threads: usize) -> FlowOptions {
    FlowOptions {
        backend: Backend::Explicit {
            energy: EnergyEval::BarnesHut(BhOptions::default()),
            solver,
        },
        threads: Some(threads),
        record_path: true,
        dense_fallback_max_vertices: 0,
        ..FlowOptions::default()
    }
}

fn small_mg() -> MgOptions {
    MgOptions {
        coarsest_vertices: 40,
        ..MgOptions::default()
    }
}

/// Runs `steps` steps; checks every accepted step: within the tolerance (when projecting), lower
/// energy, an embedded path (CCD audit). Returns the session and the summed projection solves.
fn run(opts: FlowOptions, steps: usize) -> (Session, usize) {
    let (net, pos) = trefoil(150);
    let project = opts.back_projection;
    let mut s = Session::new(net, pos, opts).unwrap();
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(2)
        .build()
        .unwrap();
    let audit = CcdOptions {
        gap: 1e-6,
        max_advances: 100_000,
    };
    let mut solves = 0;
    for _ in 0..steps {
        let r = s.step().unwrap();
        solves += r.projection_solves;
        if r.status != StepStatus::Accepted {
            continue;
        }
        assert!(r.energy_after < r.energy_before, "step {}", r.iteration);
        if project {
            let tol = 1e-4 * s.mean_edge_length0();
            assert!(
                r.max_violation <= tol * (1.0 + 1e-12),
                "step {}: violation {:.3e} > {tol:.3e}",
                r.iteration,
                r.max_violation
            );
        }
        let mut rep = CcdReport::default();
        let hit = path_contact(
            &pool,
            s.network(),
            r.path.as_ref().unwrap(),
            &[],
            &audit,
            &mut rep,
        )
        .unwrap();
        assert!(
            hit.is_none(),
            "step {}: the path has a contact",
            r.iteration
        );
    }
    (s, solves)
}

#[test]
fn second_order_correction_keeps_the_constraints_and_the_isotopy() {
    for solver in [Solver::Dense, Solver::Multigrid(small_mg())] {
        let o = FlowOptions {
            second_order_correction: true,
            newton_polish: None,
            ..options(solver.clone(), 4)
        };
        let (s, _) = run(o, 8);
        assert!(s.iteration() == 8);
        // The correction's solve is reported.
        let mut s2 = Session::new(
            s.network().clone(),
            s.positions().to_vec(),
            FlowOptions {
                second_order_correction: true,
                ..options(solver, 4)
            },
        )
        .unwrap();
        let r = s2.step().unwrap();
        assert!(r.second_order.is_some());
    }
}

#[test]
fn newton_polish_keeps_the_constraints_and_the_isotopy() {
    for soc in [false, true] {
        for solver in [Solver::Dense, Solver::Multigrid(small_mg())] {
            let o = FlowOptions {
                second_order_correction: soc,
                newton_polish: Some(100.0),
                ..options(solver, 4)
            };
            run(o, 8);
        }
    }
}

#[test]
fn without_back_projection_no_projection_solve_runs_and_the_lengths_drift() {
    let o = FlowOptions {
        back_projection: false,
        second_order_correction: true, // ignored without back-projection
        ..options(Solver::Multigrid(small_mg()), 4)
    };
    let (s, solves) = run(o, 6);
    assert_eq!(solves, 0);
    // The edge lengths are not held exactly any more (O(τ²) per step) ...
    let phi = s.constraints().values(s.network(), s.positions()).unwrap();
    let worst = phi.iter().fold(0.0f64, |m, p| m.max(p.abs()));
    assert!(
        worst > 1e-4 * s.mean_edge_length0(),
        "violation {worst:.3e}"
    );
    // ... but they drift only to second order: far below an edge length.
    assert!(worst < 0.2 * s.mean_edge_length0(), "violation {worst:.3e}");
}

#[test]
fn new_trial_rules_are_bitwise_independent_of_the_thread_count() {
    for (soc, polish) in [(true, None), (true, Some(100.0)), (false, Some(100.0))] {
        let mut first: Option<Vec<P>> = None;
        for t in [1, 3, 8] {
            let o = FlowOptions {
                second_order_correction: soc,
                newton_polish: polish,
                record_path: false,
                ..options(Solver::Multigrid(small_mg()), t)
            };
            let (net, pos) = trefoil(150);
            let mut s = Session::new(net, pos, o).unwrap();
            for _ in 0..3 {
                s.step().unwrap();
            }
            match &first {
                None => first = Some(s.positions().to_vec()),
                Some(p) => assert!(
                    p == s.positions(),
                    "soc {soc} polish {polish:?}: {t} threads"
                ),
            }
        }
    }
}

#[test]
fn multigrid_thread_cap_does_not_change_the_result() {
    let mut first: Option<Vec<P>> = None;
    for cap in [None, Some(2), Some(5)] {
        let mg = MgOptions {
            max_threads: cap,
            ..small_mg()
        };
        let o = FlowOptions {
            record_path: false,
            ..options(Solver::Multigrid(mg), 8)
        };
        let (net, pos) = trefoil(150);
        let mut s = Session::new(net, pos, o).unwrap();
        for _ in 0..3 {
            s.step().unwrap();
        }
        match &first {
            None => first = Some(s.positions().to_vec()),
            Some(p) => assert!(p == s.positions(), "cap {cap:?}"),
        }
    }
    assert!(
        MgOptions {
            max_threads: Some(0),
            ..MgOptions::default()
        }
        .validate()
        .is_err()
    );
    assert!(
        FlowOptions {
            newton_polish: Some(0.5),
            ..FlowOptions::default()
        }
        .validate()
        .is_err()
    );
}
