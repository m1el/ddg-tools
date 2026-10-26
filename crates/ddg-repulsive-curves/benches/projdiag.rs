//! Diagnostic of the constraint projection of the line search (`SPEC.md` §17l): after `--warm`
//! steps of a multigrid flow, the chord projection (Eq. 23) of the trials `γ − τ p` for
//! `τ = τ₀ / 2^j` with several solvers, printing the violation history of every chord
//! iteration, the cycles and the seconds per solve.
//!
//! ```text
//! cargo bench --bench projdiag -- --curve fhw-r --size 2000 [--warm 3] [--threads 1]
//!     [--solvers mg,mg-1e-1,dense] [--max-it 8]
//! ```

mod common;

use std::time::Instant;

use common::components::pool;
use common::*;
use ddg_repulsive_curves::metric::sobolev_gram;
use ddg_repulsive_curves::multigrid::HierarchyInput;
use ddg_repulsive_curves::solve::SaddleSolver;
use ddg_repulsive_curves::{
    Backend, BhOptions, DenseSaddle, EnergyEval, Exponents, FlowOptions, GramQuadrature, MgOptions,
    MultigridSaddle, ProjectionOptions, Session, Solver, project_constraints, projected_gradient,
};

/// The projection with `chord` Ā-chord steps first, then Newton steps with the Jacobian at the
/// iterate in the lumped-mass metric (`x = M⁻¹ Cᵀ (C M⁻¹ Cᵀ)⁻¹ (−Φ)`). Linear reports only for
/// the chord steps.
fn polished(
    solver: &mut dyn SaddleSolver,
    cset: &ddg_repulsive_curves::ConstraintSet,
    net: &ddg_repulsive_curves::CurveNetwork,
    trial: &[[f64; 3]],
    opts: &ProjectionOptions,
    mel0: f64,
    chord: u32,
) -> ddg_repulsive_curves::Projection {
    let layout = cset.layout();
    let d = layout.dim();
    let mut pos = trial.to_vec();
    let mut history = Vec::new();
    let mut linear = Vec::new();
    let mut it = 0u32;
    loop {
        let phi = cset.values(net, &pos).unwrap();
        let tol = cset
            .row_tolerances(net, &pos, opts.tolerance, mel0)
            .unwrap();
        let rel = phi
            .iter()
            .zip(&tol)
            .map(|(p, t)| p.abs() / t)
            .fold(0.0f64, f64::max);
        history.push(rel);
        let status = if rel <= 1.0 {
            Some(ddg_repulsive_curves::ProjectionStatus::Converged)
        } else if it >= opts.max_iterations + 3 {
            Some(ddg_repulsive_curves::ProjectionStatus::MaxIterations)
        } else {
            None
        };
        if let Some(status) = status {
            return ddg_repulsive_curves::Projection {
                positions: pos,
                status,
                iterations: it,
                max_violation: phi.iter().fold(0.0f64, |m, p| m.max(p.abs())),
                history,
                linear,
            };
        }
        let b: Vec<f64> = phi.iter().map(|p| -p).collect();
        let x = if it < chord {
            let sol = solver.solve(&vec![0.0; solver.n_dofs()], &b).unwrap();
            linear.push(sol.report);
            sol.x
        } else {
            let mut c = cset.jacobian(net, &pos).unwrap();
            let m = net.dual_lengths(&pos);
            let cp = c.col_ptr().to_vec();
            let vals = c.values_mut();
            for j in 0..cp.len() - 1 {
                let w = 1.0 / m[j / d].sqrt();
                for v in &mut vals[cp[j]..cp[j + 1]] {
                    *v *= w;
                }
            }
            let pr = ddg_repulsive_curves::Projector::new(&c, None).unwrap();
            let mut z = pr.pseudo_inverse(&b).unwrap();
            for (j, v) in z.iter_mut().enumerate() {
                *v /= m[j / d].sqrt();
            }
            z
        };
        layout.add_scaled(&mut pos, 1.0, &x);
        it += 1;
    }
}

fn main() {
    let a = args();
    let curve_name = a.get("curve").cloned().unwrap_or("fhw-r".into());
    let size: usize = get(&a, "size", 2000);
    let warm: usize = get(&a, "warm", 3);
    let threads: usize = get(&a, "threads", 1);
    let max_it: u32 = get(&a, "max-it", 8);
    let ntrials: usize = get(&a, "trials", 4);
    let solvers = a
        .get("solvers")
        .cloned()
        .unwrap_or("mg,mg-1e-1,dense".into());
    let paper = curve_name.starts_with("fhw");
    let alpha: f64 = get(&a, "alpha", if paper { 2.0 } else { 3.0 });
    let beta: f64 = get(&a, "beta", if paper { 4.5 } else { 6.0 });
    let exps = Exponents::new(alpha, beta).unwrap();
    let (net, pos) = curve(&curve_name, size).unwrap();
    let opts = FlowOptions {
        exponents: exps,
        backend: Backend::Explicit {
            energy: EnergyEval::BarnesHut(BhOptions::default()),
            solver: Solver::Multigrid(MgOptions::default()),
        },
        threads: Some(threads),
        dense_threads: Some(threads.min(16)),
        dense_fallback_max_vertices: 0,
        ..FlowOptions::default()
    };
    let mut s = Session::new(net.clone(), pos, opts.clone()).unwrap();
    let mut tau0 = 0.0;
    for _ in 0..warm {
        let r = s.step().unwrap();
        tau0 = r.tau_start.max(r.tau);
        eprintln!(
            "warm step {}: tau {:.3e} tau0 {:.3e} trials {}",
            r.iteration,
            r.tau,
            r.tau_start,
            r.trials.len()
        );
    }
    let pos = s.positions().to_vec();
    let cset = s.constraints().clone();
    let jac = cset.jacobian(&net, &pos).unwrap();
    let layout = cset.layout();
    let p = pool(threads);
    let dp = pool(threads.min(16));
    let ev = ddg_repulsive_curves::energy(
        &net,
        &pos,
        exps,
        &EnergyEval::BarnesHut(BhOptions::default()),
        Some(threads),
    )
    .unwrap();
    let input = HierarchyInput {
        network: &net,
        positions: &pos,
        exponents: exps,
        quadrature: GramQuadrature::Midpoint,
        layout,
        jacobian: &jac,
        constraints: Some(&cset),
    };
    let mut mg = MultigridSaddle::new(&p, &dp, &input, MgOptions::default()).unwrap();
    for (l, lv) in mg.hierarchy().levels().iter().enumerate() {
        let n = lv.network().n_vertices();
        match lv.hierarchical_metric() {
            Some(h) => {
                let st = h.stats();
                println!(
                    "level {l}: {n} vertices, near entries {} ({:.1}/row), far blocks {}, near pairs corrected {}",
                    st.near_entries,
                    st.near_entries as f64 / n as f64,
                    st.admissible_blocks,
                    lv.near_correction().pairs().len()
                );
            }
            None => println!("level {l}: {n} vertices, dense"),
        }
    }
    let t = Instant::now();
    let pg = projected_gradient(&mut mg, layout, &ev.differential).unwrap();
    println!(
        "direction: {:.3} s, {} cycles, |g| {:.3e}, tau0 {:.3e}",
        t.elapsed().as_secs_f64(),
        pg.linear.iterations,
        pg.gradient_norm,
        tau0
    );
    let dense_dir = if a.contains_key("compare") {
        let g = sobolev_gram(&p, &net, &pos, exps, GramQuadrature::Midpoint);
        let mut d = DenseSaddle::new(dp.clone(), &g.a, &jac, layout).unwrap();
        let r = projected_gradient(&mut d, layout, &ev.differential).unwrap();
        Some(r.direction)
    } else {
        None
    };
    let mel0 = s.mean_edge_length0();
    let popts = ProjectionOptions {
        max_iterations: max_it,
        ..s.options().projection
    };
    let mut plain: Vec<Vec<[f64; 3]>> = Vec::new();
    for full in solvers.split(',') {
        // suffixes: +polish (first chord step with Ā, then Newton steps with the current
        // Jacobian in the lumped-mass metric), +polish0 (Newton steps only), then +soc
        let (rest, polish) = if let Some(n) = full.strip_suffix("+polish") {
            (n, Some(1u32))
        } else if let Some(n) = full.strip_suffix("+polish0") {
            (n, Some(0u32))
        } else {
            (full, None)
        };
        let (name, soc) = match rest.strip_suffix("+soc") {
            Some(n) => (n, true),
            None => (rest, false),
        };
        let t = Instant::now();
        let mut solver: Box<dyn SaddleSolver> = match name {
            "dense" => {
                let g = sobolev_gram(&p, &net, &pos, exps, GramQuadrature::Midpoint);
                Box::new(DenseSaddle::new(dp.clone(), &g.a, &jac, layout).unwrap())
            }
            _ => {
                // mg[:tol=…][:tt=…][:tx=…]
                let mut o = MgOptions::default();
                for kv in name.split(':').skip(1) {
                    let (k, v) = kv.split_once('=').unwrap();
                    let v: f64 = v.parse().unwrap();
                    match k {
                        "tol" => {
                            o.tolerance = v;
                            o.failure_residual = o.failure_residual.max(v);
                        }
                        "tt" => o.bct.theta_t = v,
                        "tx" => o.bct.theta_x = v,
                        _ => panic!("unknown multigrid key {k}"),
                    }
                }
                let mut m = MultigridSaddle::new(&p, &dp, &input, o).unwrap();
                let st = m.hierarchy().levels()[0]
                    .hierarchical_metric()
                    .map(|h| h.stats().near_entries)
                    .unwrap_or(0);
                let t = Instant::now();
                let g = projected_gradient(&mut m, layout, &ev.differential).unwrap();
                let secs = t.elapsed().as_secs_f64();
                let rel = |r: &[[f64; 3]]| -> f64 {
                    let (mut num, mut den) = (0.0, 0.0);
                    for (a, b) in g.direction.iter().zip(r) {
                        for c in 0..3 {
                            num += (a[c] - b[c]).powi(2);
                            den += b[c] * b[c];
                        }
                    }
                    (num / den).sqrt()
                };
                println!(
                    "  {name}: level-0 near entries {st}, direction {secs:.3} s, {} cycles, vs default mg {:.2e}{}",
                    g.linear.iterations,
                    rel(&pg.direction),
                    match &dense_dir {
                        Some(d) => format!(", vs dense {:.2e}", rel(d)),
                        None => String::new(),
                    }
                );
                Box::new(m)
            }
        };
        println!("solver {full}: build {:.3} s", t.elapsed().as_secs_f64());
        // Second-order correction: q with C q = −Φ₂, Φ₂ = ½ Φ''[p, p] by central differences.
        let q: Vec<[f64; 3]> = if soc {
            let t = Instant::now();
            let pmax = pg
                .direction
                .iter()
                .map(|d| (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt())
                .fold(0.0f64, f64::max);
            let h = 1e-3 * mel0 / pmax;
            let shift = |sgn: f64| -> Vec<[f64; 3]> {
                pos.iter()
                    .zip(&pg.direction)
                    .map(|(b, d)| {
                        [
                            b[0] + sgn * h * d[0],
                            b[1] + sgn * h * d[1],
                            b[2] + sgn * h * d[2],
                        ]
                    })
                    .collect()
            };
            let fp = cset.values(&net, &shift(1.0)).unwrap();
            let fm = cset.values(&net, &shift(-1.0)).unwrap();
            let f0 = cset.values(&net, &pos).unwrap();
            let b: Vec<f64> = (0..f0.len())
                .map(|r| -(fp[r] + fm[r] - 2.0 * f0[r]) / (2.0 * h * h))
                .collect();
            let sol = solver.solve(&vec![0.0; solver.n_dofs()], &b).unwrap();
            println!(
                "  soc solve {:.3} s, {} cycles",
                t.elapsed().as_secs_f64(),
                sol.report.iterations
            );
            layout.scatter(&sol.x)
        } else {
            vec![[0.0; 3]; pos.len()]
        };
        for j in 0..ntrials {
            let tau = tau0 / f64::from(1u32 << j);
            let trial: Vec<[f64; 3]> = pos
                .iter()
                .zip(&pg.direction)
                .zip(&q)
                .map(|((b, d), q)| {
                    let t2 = tau * tau;
                    [
                        b[0] - tau * d[0] + t2 * q[0],
                        b[1] - tau * d[1] + t2 * q[1],
                        b[2] - tau * d[2] + t2 * q[2],
                    ]
                })
                .collect();
            let t = Instant::now();
            let pr = match polish {
                None => {
                    project_constraints(solver.as_mut(), &cset, &net, &trial, &popts, mel0).unwrap()
                }
                Some(chord) => polished(solver.as_mut(), &cset, &net, &trial, &popts, mel0, chord),
            };
            let secs = t.elapsed().as_secs_f64();
            let cycles: Vec<usize> = pr.linear.iter().map(|l| l.iterations).collect();
            let hist: Vec<String> = pr.history.iter().map(|h| format!("{h:.2e}")).collect();
            let moved: f64 = pr
                .positions
                .iter()
                .zip(&trial)
                .map(|(a, b)| (0..3).map(|c| (a[c] - b[c]).powi(2)).sum::<f64>())
                .sum::<f64>()
                .sqrt();
            let diff = if plain.len() > j {
                pr.positions
                    .iter()
                    .zip(&plain[j])
                    .map(|(a, b)| (0..3).map(|c| (a[c] - b[c]).powi(2)).sum::<f64>())
                    .sum::<f64>()
                    .sqrt()
            } else {
                plain.push(pr.positions.clone());
                0.0
            };
            println!(
                "  tau0/{:<2} {:?} it {} {:.3} s cycles {:?} |x| {:.3e} |vs first| {:.2e} history [{}]",
                1u32 << j,
                pr.status,
                pr.iterations,
                secs,
                cycles,
                moved,
                diff,
                hist.join(" ")
            );
        }
    }
}
