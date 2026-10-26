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

/// Pulls vertices `a` and `b` together: `w |x_a − x_b|`.
#[derive(Debug)]
struct Spring {
    a: usize,
    b: usize,
    w: f64,
}

impl ddg_repulsive_curves::PenaltyTerm for Spring {
    fn name(&self) -> &str {
        "spring"
    }
    fn evaluate(
        &self,
        _network: &CurveNetwork,
        positions: &[[f64; 3]],
        gradient: Option<&mut [[f64; 3]]>,
    ) -> Result<f64, ddg_repulsive_curves::RcError> {
        let d: Vec<f64> = (0..3)
            .map(|c| positions[self.a][c] - positions[self.b][c])
            .collect();
        let n = (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt();
        if let Some(g) = gradient
            && n > 0.0
        {
            for c in 0..3 {
                g[self.a][c] += self.w * d[c] / n;
                g[self.b][c] -= self.w * d[c] / n;
            }
        }
        Ok(self.w * n)
    }
    fn subdivided(
        &self,
        _map: &ddg_repulsive_curves::SubdivisionMap,
    ) -> Result<std::sync::Arc<dyn ddg_repulsive_curves::PenaltyTerm>, ddg_repulsive_curves::RcError>
    {
        Err(ddg_repulsive_curves::RcError::Unsupported {
            what: "spring subdivision",
        })
    }
}

/// A metric blind to the repulsion (L²) and a penalty that overpowers the tangent-point energy
/// drive two strands together; each guarded step can only close a fraction of the gap, so before
/// the contact floor (`SPEC.md` §17l) the distance shrank geometrically to rounding (with the
/// floor set to 0 the curve fails `validate_positions` after step 25: edges 9 and 29 at 1.8e-13).
/// Now the flow ends `Stuck` at the floor (step 20, distance 5.3e-7 ℓ̄₀) with a valid curve and an
/// embedded path at every step.
#[test]
fn contact_floor_stops_a_flow_that_drives_two_edges_together() {
    let n = 40;
    let pos: Vec<P> = (0..n)
        .map(|k| {
            let t = std::f64::consts::TAU * k as f64 / n as f64;
            [3.0 * t.cos(), 0.6 * t.sin(), 0.05 * (2.0 * t).sin()]
        })
        .collect();
    let edges = (0..n as u32).map(|k| [k, (k + 1) % n as u32]).collect();
    let net = CurveNetwork::new(n, edges).unwrap();
    let o = FlowOptions {
        metric: ddg_repulsive_curves::Metric::L2,
        backend: Backend::Explicit {
            energy: EnergyEval::Exact,
            solver: Solver::Dense,
        },
        penalties: vec![std::sync::Arc::new(Spring {
            a: 10,
            b: 30,
            w: 1e60,
        })],
        threads: Some(2),
        record_path: true,
        max_iterations: 400,
        ..FlowOptions::default()
    };
    let mut s = Session::new(net, pos, o).unwrap();
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(2)
        .build()
        .unwrap();
    let audit = CcdOptions {
        gap: 1e-6,
        max_advances: 100_000,
    };
    let mut floor_hit = false;
    for _ in 0..400 {
        let r = s.step().unwrap();
        if let Some(p) = &r.path {
            let mut rep = CcdReport::default();
            let hit = path_contact(&pool, s.network(), p, &[], &audit, &mut rep).unwrap();
            assert!(
                hit.is_none(),
                "step {}: the path has a contact",
                r.iteration
            );
        }
        s.network().validate_positions(s.positions()).unwrap();
        if r.iteration.is_multiple_of(25) || r.contact_floor {
            eprintln!(
                "step {}: {:?}, tau {:.3e}, min distance / l0 {:.3e}",
                r.iteration,
                r.status,
                r.tau,
                ddg_repulsive_curves::ccd::min_pair_distance(s.network(), s.positions())
                    / s.mean_edge_length0()
            );
        }
        if r.contact_floor {
            floor_hit = true;
            assert_eq!(r.status, StepStatus::Stuck);
            break;
        }
        if s.status().is_terminal() {
            break;
        }
    }
    assert!(
        floor_hit,
        "the flow ended {:?} before the floor",
        s.status()
    );
    assert_eq!(s.status(), ddg_repulsive_curves::Status::Stuck);
    let d = ddg_repulsive_curves::ccd::min_pair_distance(s.network(), s.positions());
    let l0 = s.mean_edge_length0();
    assert!(
        d < 1e-6 * l0 && d > 1e-8 * l0,
        "distance {d:.3e}, ℓ̄₀ {l0:.3e}"
    );
}

/// A Sobolev flow never comes near the floor: no step reports it.
#[test]
fn contact_floor_is_never_hit_by_a_sobolev_flow() {
    let (_, solves) = run(options(Solver::Dense, 2), 10);
    let _ = solves;
    let (net, pos) = trefoil(150);
    let mut s = Session::new(net, pos, options(Solver::Dense, 2)).unwrap();
    for _ in 0..10 {
        assert!(!s.step().unwrap().contact_floor);
    }
}
