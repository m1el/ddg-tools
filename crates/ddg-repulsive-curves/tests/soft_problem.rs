//! The soft problem of the Sec. 7 evaluation (`SPEC.md` §6.4, §14.2, §17m).
//!
//! At a minimizer of the hard problem (barycenter + edge lengths) Eq. 22 gives the edge
//! multipliers `λ_I` with `dÊ = −Σ_I λ_I ∇ℓ_I` (the rows are `ℓ⁰_I − ℓ_I`). The soft problem
//! `Ê + w Σ_I ℓ⁰_I ln²(ℓ_I/ℓ⁰_I)` is stationary where `dÊ + Σ_I 2w (ℓ⁰_I/ℓ_I) ln(ℓ_I/ℓ⁰_I) ∇ℓ_I
//! = 0`, so to first order in `1/w` its minimizer near the hard one has the strains
//! `ln(ℓ_I/ℓ⁰_I) ≈ λ_I / (2w)`: the Hencky weight sets the equilibrium strain through the edge
//! tensions, which is how §17m chooses `w_H`. Both flows here run to their stop test.

use std::sync::Arc;

use ddg_repulsive_curves::energy::exact::{energy_and_differential, make_pool};
use ddg_repulsive_curves::metric::{GramQuadrature, sobolev_gram};
use ddg_repulsive_curves::optimize::{Hencky, Method};
use ddg_repulsive_curves::solve::{ProjectionOptions, Tolerance};
use ddg_repulsive_curves::step::stop::{GradNorm, Status, StopRule};
use ddg_repulsive_curves::{
    Backend, ConstraintSet, CurveNetwork, DenseSaddle, EnergyEval, Exponents, FlowOptions, Session,
    Solver, default_constraints, projected_gradient,
};

type P = [f64; 3];

fn rough_trefoil(n: usize) -> (CurveNetwork, Vec<P>) {
    let pos = (0..n)
        .map(|i| {
            let t = std::f64::consts::TAU * i as f64 / n as f64;
            let w = 0.35 * (4.0 * t).sin();
            [
                t.sin() + 2.0 * (2.0 * t).sin() + w,
                t.cos() - 2.0 * (2.0 * t).cos() + 0.5 * w,
                -(3.0 * t).sin() + 0.3 * (7.0 * t).cos(),
            ]
        })
        .collect();
    let edges = (0..n as u32).map(|i| [i, (i + 1) % n as u32]).collect();
    (CurveNetwork::new(n, edges).unwrap(), pos)
}

fn exact_dense() -> Backend {
    Backend::Explicit {
        energy: EnergyEval::Exact,
        solver: Solver::Dense,
    }
}

fn lengths(net: &CurveNetwork, pos: &[P]) -> Vec<f64> {
    net.edges()
        .iter()
        .map(|&[a, b]| {
            let (p, q) = (pos[a as usize], pos[b as usize]);
            ((p[0] - q[0]).powi(2) + (p[1] - q[1]).powi(2) + (p[2] - q[2]).powi(2)).sqrt()
        })
        .collect()
}

fn run(net: &CurveNetwork, pos: &[P], opts: FlowOptions) -> (Vec<P>, Status, usize) {
    let mut s = Session::new(net.clone(), pos.to_vec(), opts).unwrap();
    loop {
        let r = s.step().unwrap();
        if std::env::var("RC_TRACE").is_ok() && (r.iteration % 100 == 1 || r.iteration < 20) {
            eprintln!(
                "TRACE {} F {:.12} g {:.3e} dot {:.3e} tau {:.3e} h {} viol {:.1e}",
                r.iteration,
                r.energy_after,
                r.gradient_norm,
                r.sobolev_norm_sq,
                r.tau,
                r.halvings,
                r.max_violation
            );
        }
        if s.status().is_terminal() {
            break;
        }
    }
    (s.positions().to_vec(), s.status(), s.iteration())
}

#[test]
fn soft_equilibrium_strain_is_tension_over_two_w() {
    let exps = Exponents::new(2.0, 4.5).unwrap();
    let (net, pos0) = rough_trefoil(40);
    let stop = |eps| StopRule::GradientNorm {
        norm: GradNorm::Euclidean,
        eps,
    };
    // 1. The hard minimizer.
    let hard = FlowOptions {
        exponents: exps,
        backend: exact_dense(),
        // A tight projection tolerance: with the default one the flow stalls at the tolerance
        // boundary (`SPEC.md` §17m, Finding B) before `‖g̃‖` is small.
        projection: ProjectionOptions {
            tolerance: Tolerance::MeanEdgeLength(1e-10),
            max_iterations: 3,
        },
        stop: stop(3e-5),
        max_iterations: 3000,
        threads: Some(2),
        dense_threads: Some(1),
        ..FlowOptions::default()
    };
    let (xh, st, it) = run(&net, &pos0, hard);
    assert_eq!(st, Status::Converged, "hard flow: {it} iterations");

    // 2. Its edge multipliers.
    let pool = Arc::new(make_pool(Some(2)).unwrap());
    let (_, diff) = energy_and_differential(&pool, &net, &xh, exps);
    let set = ConstraintSet::new(&net, &xh, &default_constraints(), None).unwrap();
    set.analyze(&net, &xh).unwrap();
    let c = set.jacobian(&net, &xh).unwrap();
    let g = sobolev_gram(&pool, &net, &xh, exps, GramQuadrature::Midpoint);
    let mut solver = DenseSaddle::new(pool.clone(), &g.a, &c, set.layout()).unwrap();
    let pg = projected_gradient(&mut solver, set.layout(), &diff).unwrap();
    assert!(pg.gradient_norm < 3e-5);
    let lambda = &pg.lambda[3..]; // three barycenter rows first
    assert_eq!(lambda.len(), net.n_edges());
    // Scaling identity: Σ λ_I ℓ_I = −dÊ·γ = −(2 + α − β) Ê = ½ Ê > 0 (the curve wants to grow).
    let l = lengths(&net, &xh);
    let e_hard = energy_and_differential(&pool, &net, &xh, exps).0;
    let sum: f64 = lambda.iter().zip(&l).map(|(a, b)| a * b).sum();
    assert!(
        (sum - 0.5 * e_hard).abs() < 1e-6 * e_hard,
        "{sum} vs {}",
        0.5 * e_hard
    );

    // 3. The soft minimizer (Hs NCG) from the hard one, rest lengths = the hard lengths.
    // w = 10 and 30 converge in a few hundred iterations; w = 100 is slow (2702 iterations with
    // the dense LU, more than 4000 with the Schur solver: the slow tail is path sensitive), so
    // it is only required to be much slower, not to converge.
    let cap = 6000;
    let mut iterations = Vec::new();
    for w in [10.0, 30.0, 100.0] {
        let h = Hencky::at_rest(w, &net, &xh).unwrap();
        let soft = FlowOptions {
            exponents: exps,
            backend: exact_dense(),
            method: Method::Ncg,
            constraints: Vec::new(),
            penalties: vec![Arc::new(h)],
            stop: stop(1e-6),
            max_iterations: cap,
            threads: Some(2),
            dense_threads: Some(1),
            ..FlowOptions::default()
        };
        let (xs, st, it) = run(&net, &xh, soft);
        iterations.push(it);
        let ls = lengths(&net, &xs);
        let mut worst = 0.0f64;
        let lmax = lambda.iter().fold(0.0f64, |m, x| m.max(x.abs()));
        for ((a, b), lam) in ls.iter().zip(&l).zip(lambda) {
            let strain = (a / b).ln();
            let predicted = lam / (2.0 * w);
            worst = worst.max((strain - predicted).abs());
        }
        // first order in 1/w: the error is O(1/w²) relative to λ/(2w)
        let scale = lmax / (2.0 * w);
        eprintln!(
            "w {w}: {st:?} after {it} NCG iterations, max predicted strain {scale:.3e}, \
             worst error {worst:.3e}"
        );
        if w < 50.0 {
            assert_eq!(st, Status::Converged, "soft w = {w}: {it} iterations");
            assert!(worst < 0.15 * scale, "w {w}: {worst} vs {scale}");
        }
    }
    // The Hencky stiffness (a second-order operator of size 2w along the tangents) is not in the
    // Hs preconditioner: ten times the weight costs many times the iterations (measured 180
    // against 2702 or more), the conditioning effect of §17m Finding A.
    assert!(
        iterations[2] > 4 * iterations[0],
        "Hs NCG iterations at w = 10, 30, 100: {iterations:?}"
    );
}
