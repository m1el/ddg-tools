//! The factor-once solver returns the reference solver's optimum on generated meshes.

use ddg_optimal_cones::{ConeProblem, SolveOptions, solve};
use ddg_solver_lab::lasso::{self, Options};
use ddg_solver_lab::meshes::load;
use ddg_solver_lab::prepare::{Constraints, Prepared};

fn check(spec: &str, lambda: f64, cons: &Constraints, threads: usize) {
    let case = load(spec).unwrap();
    let nv = case.mesh.n_vertices();
    let prep = Prepared::with_constraints(&case, cons).unwrap();
    let s = lasso::solve(&prep, lambda, &Options::for_threads(threads));
    assert!(s.stats.converged, "{spec}: not converged");
    let mut p = ConeProblem::new(&case.mesh, &case.geometry, lambda);
    p.cone_costs = cons.cost.clone();
    if cons.nonneg.is_some() {
        p.lower_bounds = Some(vec![0.0; nv]);
    }
    let cert = ddg_optimal_cones::certify(&p, &s.u, &s.phi, &s.mu).unwrap();
    assert!(cert.holds(1e-9), "{spec}: certificate {cert:?}");
    let r = solve(&p, &SolveOptions::default()).unwrap();
    let dmu = s.mu.iter().zip(&r.measure).map(|(a, b)| (a - b).abs()).fold(0.0, f64::max);
    let du = s.u.iter().zip(&r.log_scale_optimal).map(|(a, b)| (a - b).abs()).fold(0.0, f64::max);
    assert!(dmu < 1e-7 && du < 1e-7, "{spec} lambda {lambda}: dmu {dmu:e} du {du:e}");
    assert!((s.objective - r.objective).abs() < 1e-10 * r.objective.abs().max(1.0));
}

#[test]
fn bounded_disk() {
    for lambda in [0.5, 0.02, 0.004, 0.001] {
        check("bump:40", lambda, &Constraints::default(), 1);
    }
}

#[test]
fn closed_sphere() {
    check("sphere:3", 0.05, &Constraints::default(), 1);
    check("sphere:4", 0.05, &Constraints::default(), 4);
}

#[test]
fn steering() {
    let case = load("sphere:4").unwrap();
    let nv = case.mesh.n_vertices();
    let mut cost = vec![f64::INFINITY; nv];
    for i in 0..8 {
        cost[i * nv / 8 + nv / 16] = 0.0;
    }
    check("sphere:4", 0.05, &Constraints { cost: Some(cost), nonneg: Some(vec![true; nv]), nonpos: None }, 1);
    let case = load("bump:40").unwrap();
    let nv = case.mesh.n_vertices();
    let mut cost = vec![1.0; nv];
    for v in (0..nv).step_by(37) {
        cost[v] = 0.0;
    }
    for c in cost.iter_mut().take(nv / 4).skip(1).step_by(2) {
        *c = f64::INFINITY;
    }
    check("bump:40", 0.002, &Constraints { cost: Some(cost), nonneg: None, nonpos: None }, 1);
}
