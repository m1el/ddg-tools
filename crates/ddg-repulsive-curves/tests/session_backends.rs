//! The back ends of the session (WP3c, `SPEC.md` §8.5, §17h): the thread pools reach the
//! solvers, `Solver::Multigrid` with live subdivision and constraints, the dense fallback,
//! `Backend::Auto`.

use std::sync::Arc;

use ddg_repulsive_curves::{
    Backend, BhOptions, CurveNetwork, EnergyEval, FlowOptions, MgOptions, MultigridSolvers,
    Session, Solver,
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
        for sn in ["dense", "custom-mg"] {
            let mut curves = Vec::new();
            for t in [1, 3, 8] {
                let solver = match sn {
                    "dense" => Solver::Dense,
                    _ => Solver::Custom(Arc::new(MultigridSolvers::new(small_mg()).unwrap())),
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
