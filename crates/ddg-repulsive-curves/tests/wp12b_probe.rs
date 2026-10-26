//! WP12b probes (`SPEC.md` §17m): the soft problem of the Sec. 7 evaluation and the stop test
//! of the hard one. Env-driven and `#[ignore]`d; run one probe per process, e.g.
//!
//! ```text
//! RC_START=final.obj RC_BACKEND=exact-dense RC_ITERS=50 cargo test --release \
//!   -p ddg-repulsive-curves --test wp12b_probe flow_probe -- --ignored --nocapture
//! ```
//!
//! * `flow_probe`: a flow from `RC_START` (an OBJ; the network is taken from it), hard
//!   (barycenter + edge lengths at their start values) unless `RC_SOFT=ncg|lbfgs` (then `Ê` +
//!   Hencky with weight `RC_W` and rest lengths of `RC_REST`, default the start curve), metric
//!   `RC_METRIC=hs|h1`, back end `RC_BACKEND=exact-dense|bh-dense|exact-mg|bh-mg|auto`, multigrid
//!   tolerance `RC_MG_TOL`, `RC_ITERS`, `RC_SECONDS`, `RC_THREADS`; one `STEP` line per step;
//!   `RC_OUT` writes the final curve.
//! * `multipliers_probe`: at `RC_START`, exact energy + dense solve of Eq. 22 with barycenter +
//!   edge lengths: `‖g̃‖`, `dE·g̃` and the edge multipliers `λ_I` (the edge tensions; the soft
//!   equilibrium strain is `≈ −λ_I / (2 w_H)`).

use std::sync::Arc;
use std::time::{Duration, Instant};

use ddg_repulsive_curves::energy::exact::{energy_and_differential, make_pool};
use ddg_repulsive_curves::metric::{GramQuadrature, sobolev_gram};
use ddg_repulsive_curves::optimize::{Hencky, Method};
use ddg_repulsive_curves::solve::{ProjectionOptions, Tolerance};
use ddg_repulsive_curves::step::stop::StepStatus;
use ddg_repulsive_curves::{
    Backend, ConstraintSet, DenseSaddle, EnergyEval, Exponents, FlowOptions, Metric, MgOptions,
    PenaltyTerm, Session, Solver, default_constraints, projected_gradient, read_curve_obj,
    write_curve_obj,
};

fn env(name: &str) -> Option<String> {
    std::env::var(name).ok().filter(|s| !s.is_empty())
}

fn envf(name: &str, default: f64) -> f64 {
    env(name).and_then(|s| s.parse().ok()).unwrap_or(default)
}

fn exps() -> Exponents {
    Exponents::new(envf("RC_ALPHA", 2.0), envf("RC_BETA", 4.5)).unwrap()
}

#[test]
#[ignore]
fn flow_probe() {
    let Some(start) = env("RC_START") else {
        eprintln!("RC_START not set; skipped");
        return;
    };
    let (net, pos) = read_curve_obj(&start).unwrap();
    let threads = envf("RC_THREADS", 1.0) as usize;
    let mut mg = MgOptions::default();
    if let Some(t) = env("RC_MG_TOL").and_then(|s| s.parse().ok()) {
        mg.tolerance = t;
        mg.failure_residual = mg.failure_residual.max(t);
    }
    let backend = match env("RC_BACKEND").as_deref().unwrap_or("auto") {
        "exact-dense" => Backend::Explicit {
            energy: EnergyEval::Exact,
            solver: Solver::Dense,
        },
        "bh-dense" => Backend::Explicit {
            energy: EnergyEval::BarnesHut(Default::default()),
            solver: Solver::Dense,
        },
        "exact-mg" => Backend::Explicit {
            energy: EnergyEval::Exact,
            solver: Solver::Multigrid(mg),
        },
        "bh-mg" => Backend::Explicit {
            energy: EnergyEval::BarnesHut(Default::default()),
            solver: Solver::Multigrid(mg),
        },
        _ => Backend::default(),
    };
    let metric = match env("RC_METRIC").as_deref() {
        Some("h1") => Metric::H1,
        _ => Metric::default(),
    };
    let method = match env("RC_SOFT").as_deref() {
        Some("ncg") => Method::Ncg,
        Some("lbfgs") => Method::Lbfgs { memory: 30 },
        _ => Method::ProjectedGradient,
    };
    let mut opts = FlowOptions {
        exponents: exps(),
        metric,
        backend,
        method,
        projection: ProjectionOptions {
            tolerance: Tolerance::Absolute(1e-4),
            ..ProjectionOptions::default()
        },
        max_iterations: envf("RC_ITERS", 50.0) as usize,
        time_limit: env("RC_SECONDS")
            .and_then(|s| s.parse::<f64>().ok())
            .map(Duration::from_secs_f64),
        threads: Some(threads),
        dense_fallback_max_vertices: 0,
        ..FlowOptions::default()
    };
    let mut hencky = None;
    if method.is_soft() {
        let rest = match env("RC_REST") {
            Some(p) => read_curve_obj(&p).unwrap().1,
            None => pos.clone(),
        };
        let h = Hencky::at_rest(envf("RC_W", 1e4), &net, &rest).unwrap();
        opts.constraints = Vec::new();
        opts.penalties = vec![Arc::new(h.clone())];
        hencky = Some(h);
    }
    let mut s = Session::new(net, pos, opts).unwrap();
    let t0 = Instant::now();
    eprintln!(
        "HEAD it t F E g gmass soboDot tau tau_start halvings lin_it lin_res max_viol status strain"
    );
    loop {
        let r = match s.step() {
            Ok(r) => r,
            Err(e) => {
                eprintln!("ERROR {e}");
                break;
            }
        };
        let strain = hencky
            .as_ref()
            .map_or(0.0, |h| h.max_log_strain(s.network(), s.positions()));
        eprintln!(
            "STEP {} {:.2} {:.10} {:.10} {:.4e} {:.4e} {:.4e} {:.3e} {:.3e} {} {} {:.2e} {:.2e} {:?} {:.3e}",
            r.iteration,
            t0.elapsed().as_secs_f64(),
            r.energy_after,
            r.tangent_point_after,
            r.gradient_norm,
            r.gradient_mass_norm,
            r.sobolev_norm_sq,
            r.tau,
            r.tau_start,
            r.halvings,
            r.linear.iterations,
            r.linear.relative_residual,
            r.max_violation,
            r.status,
            strain
        );
        // `Session::step` does not apply `time_limit` (only `run` does): the wall cap is here.
        let over = s.options().time_limit.is_some_and(|l| t0.elapsed() >= l);
        if over || s.status().is_terminal() || r.status == StepStatus::Converged {
            break;
        }
    }
    eprintln!("END status {:?} iterations {}", s.status(), s.iteration());
    if let Some(out) = env("RC_OUT") {
        write_curve_obj(&out, s.network(), s.positions()).unwrap();
    }
}

#[test]
#[ignore]
fn multipliers_probe() {
    let Some(start) = env("RC_START") else {
        eprintln!("RC_START not set; skipped");
        return;
    };
    let (net, pos) = read_curve_obj(&start).unwrap();
    let threads = envf("RC_THREADS", 1.0) as usize;
    let pool = Arc::new(make_pool(Some(threads)).unwrap());
    let e = exps();
    let (en, diff) = energy_and_differential(&pool, &net, &pos, e);
    let set = ConstraintSet::new(&net, &pos, &default_constraints(), None).unwrap();
    set.analyze(&net, &pos).unwrap();
    let c = set.jacobian(&net, &pos).unwrap();
    let g = sobolev_gram(&pool, &net, &pos, e, GramQuadrature::Midpoint);
    let mut solver = DenseSaddle::new(pool.clone(), &g.a, &c, set.layout()).unwrap();
    let pg = projected_gradient(&mut solver, set.layout(), &diff).unwrap();
    let dnorm = diff.iter().flatten().map(|x| x * x).sum::<f64>().sqrt();
    // rows: 3 barycenter rows, then one per edge
    let mut lam: Vec<f64> = pg.lambda[3..].to_vec();
    lam.sort_by(f64::total_cmp);
    let pct = |q: f64| lam[((lam.len() - 1) as f64 * q) as usize];
    let amax = lam.iter().fold(0.0f64, |m, x| m.max(x.abs()));
    eprintln!(
        "MULT E {en:.8} |dE| {dnorm:.4e} |g| {:.4e} dE.g {:.4e} lambda min {:.4e} p5 {:.4e} \
         median {:.4e} p95 {:.4e} max {:.4e} max|.| {amax:.4e} residual {:.1e}",
        pg.gradient_norm,
        pg.sobolev_norm_sq,
        lam[0],
        pct(0.05),
        pct(0.5),
        pct(0.95),
        lam[lam.len() - 1],
        pg.linear.relative_residual
    );
    for w in [1.0, 10.0, 100.0, 1e3, 1e4] {
        eprintln!(
            "MULT w {w:e}: predicted soft equilibrium max |ln(l/l0)| ~ {:.3e}",
            amax / (2.0 * w)
        );
    }
    // the Hencky differential on this curve against central differences of its value
    let h = Hencky::at_rest(
        1e4,
        &net,
        &read_curve_obj(env("RC_REST").as_deref().unwrap_or(&start))
            .unwrap()
            .1,
    )
    .unwrap();
    let mut grad = vec![[0.0; 3]; pos.len()];
    h.evaluate(&net, &pos, Some(&mut grad)).unwrap();
    let mut worst = 0.0f64;
    for (k, v) in [3usize, 500, 1000].iter().enumerate() {
        let v = (*v).min(pos.len() - 1);
        for a in 0..3 {
            let step = 1e-6;
            let mut p = pos.clone();
            p[v][a] += step;
            let fp = h.evaluate(&net, &p, None).unwrap();
            p[v][a] -= 2.0 * step;
            let fm = h.evaluate(&net, &p, None).unwrap();
            let fd = (fp - fm) / (2.0 * step);
            let err = (fd - grad[v][a]).abs() / grad[v][a].abs().max(1e-3);
            worst = worst.max(err);
            if k == 0 {
                eprintln!(
                    "MULT hencky fd v{v} a{a}: fd {fd:.6e} grad {:.6e}",
                    grad[v][a]
                );
            }
        }
    }
    eprintln!("MULT hencky fd worst relative error {worst:.2e}");
}

/// The tail of the hard Hs flow on a 40-vertex trefoil (exact + dense) for several projection
/// tolerances `RC_TOLS` (comma separated, `MeanEdgeLength`): `‖g̃‖`, `dE·g̃`, `τ`, the largest
/// violation at the base and the share of accepted steps whose trial needed no projection.
#[test]
#[ignore]
fn small_hard_tail() {
    let n = 40;
    let pos: Vec<[f64; 3]> = (0..n)
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
    let net = ddg_repulsive_curves::CurveNetwork::new(n, edges).unwrap();
    let tols: Vec<f64> = env("RC_TOLS")
        .unwrap_or_else(|| "1e-4,1e-6,1e-8".into())
        .split(',')
        .filter_map(|s| s.trim().parse().ok())
        .collect();
    let iters = envf("RC_ITERS", 1500.0) as usize;
    for tol in tols {
        let opts = FlowOptions {
            exponents: exps(),
            backend: Backend::Explicit {
                energy: EnergyEval::Exact,
                solver: Solver::Dense,
            },
            projection: ProjectionOptions {
                tolerance: Tolerance::MeanEdgeLength(tol),
                max_iterations: envf("RC_PROJ_ITERS", 3.0) as u32,
            },
            stop: ddg_repulsive_curves::StopRule::GradientNorm {
                norm: ddg_repulsive_curves::GradNorm::Euclidean,
                eps: 1e-9,
            },
            max_iterations: iters,
            threads: Some(2),
            dense_threads: Some(1),
            ..FlowOptions::default()
        };
        let mut s = Session::new(net.clone(), pos.clone(), opts).unwrap();
        let (mut unprojected, mut accepted) = (0usize, 0usize);
        loop {
            let r = s.step().unwrap();
            if r.status == StepStatus::Accepted {
                accepted += 1;
                if r.projection_iterations == 0 {
                    unprojected += 1;
                }
            }
            if [10, 30, 100, 300, 1000, 3000].contains(&r.iteration) || s.status().is_terminal() {
                eprintln!(
                    "TAIL tol {tol:.0e} it {} F {:.12} g {:.3e} dot {:.3e} tau {:.3e} viol {:.2e} \
                     unprojected {unprojected}/{accepted} {:?}",
                    r.iteration,
                    r.energy_after,
                    r.gradient_norm,
                    r.sobolev_norm_sq,
                    r.tau,
                    r.max_violation,
                    s.status()
                );
            }
            if s.status().is_terminal() {
                break;
            }
        }
    }
}
