//! Fig. 20 scaling (`SPEC.md` §14.1, §17k): mean time per iteration of a running flow, one
//! configuration per invocation, one CSV row on stdout (header with `--header 1`).
//!
//! ```text
//! cargo bench --bench fig20 -- --curve fhw|fhw-r|knot --size N --config exact|bh|mg|mg-noproj|mg-bary
//!     [--threads T] [--dense-threads D] [--steps 20] [--min-steps 5] [--cap-s 240] [--trace 1]
//!     [--soc 0|1] [--proj-it 3] [--polish F] [--abort 0|1] [--bct-tt 2] [--bct-tx 0.25] [--mg-threads 16]
//! ```
//!
//! `--trace 1` prints one line per step on stderr (phases, cycles, every trial of the line
//! search as `τ/τ₀:chord iterations:outcome`).
//!
//! `exact`: exact energy + dense; `bh`: Barnes–Hut energy + dense; `mg`: BH + block cluster tree +
//! multigrid; `mg-noproj`: the same without the back-projection of the line search's trials
//! (`FlowOptions::back_projection = false`, the paper's "Multigrid without projection");
//! `mg-bary`: multigrid with the barycenter as the only constraint (the proxy of WP14, `SPEC.md`
//! §17k). `--soc 0` turns the second-order correction of the trials off (`SPEC.md` §17l).
//! Exponents (2, 4.5) on `fhw`/`fhw-r`,
//! (3, 6) otherwise (`--alpha/--beta` override). Barycenter + edge lengths otherwise.

mod common;

use std::time::Instant;

use common::*;
use ddg_repulsive_curves::{
    Backend, BhOptions, Constraint, EnergyEval, Exponents, FlowOptions, MgOptions, Session, Solver,
    TrialOutcome,
};

fn main() {
    let a = args();
    if a.contains_key("header") {
        println!(
            "curve,n,config,threads,dense_threads,steps,mean_total,median_total,energy,metric,solve,\
             line_search,ls_energy,projection,ccd,mean_mg_cycles,fallbacks,peak_rss_gb,load_start,load_end,\
             soc,second_order,proj_solves,mean_tau,max_violation,energy_first,energy_last,note"
        );
        return;
    }
    let curve_name = a.get("curve").cloned().unwrap_or("fhw".into());
    let size: usize = get(&a, "size", 600);
    let config = a.get("config").cloned().unwrap_or("mg".into());
    let threads: usize = get(&a, "threads", 1);
    let dense_threads: usize = get(&a, "dense-threads", threads.min(16));
    let max_steps: usize = get(&a, "steps", 20);
    let min_steps: usize = get(&a, "min-steps", 5);
    let cap: f64 = get(&a, "cap-s", 240.0);
    let trace = a.contains_key("trace");
    let paper = curve_name.starts_with("fhw");
    let alpha: f64 = get(&a, "alpha", if paper { 2.0 } else { 3.0 });
    let beta: f64 = get(&a, "beta", if paper { 4.5 } else { 6.0 });

    let (net, pos) = match curve(&curve_name, size) {
        Ok(c) => c,
        Err(e) => {
            eprintln!("{e}");
            std::process::exit(2);
        }
    };
    let (energy, solver) = match config.as_str() {
        "exact" => (EnergyEval::Exact, Solver::Dense),
        "bh" => (EnergyEval::BarnesHut(BhOptions::default()), Solver::Dense),
        _ => {
            // --bct-tt / --bct-tx: admissibility of the multigrid's block cluster trees
            let mut mg = MgOptions::default();
            mg.bct.theta_t = get(&a, "bct-tt", mg.bct.theta_t);
            mg.bct.theta_x = get(&a, "bct-tx", mg.bct.theta_x);
            // --mg-threads N: MgOptions::max_threads (0: no cap)
            if let Some(t) = a.get("mg-threads").and_then(|s| s.parse::<usize>().ok()) {
                mg.max_threads = (t > 0).then_some(t);
            }
            (
                EnergyEval::BarnesHut(BhOptions::default()),
                Solver::Multigrid(mg),
            )
        }
    };
    let mut opts = FlowOptions {
        exponents: Exponents::new(alpha, beta).unwrap(),
        backend: Backend::Explicit { energy, solver },
        threads: Some(threads),
        dense_threads: Some(dense_threads),
        // the measurement must not silently switch solver
        dense_fallback_max_vertices: 0,
        ..FlowOptions::default()
    };
    match config.as_str() {
        // the paper's "Multigrid without projection"
        "mg-noproj" => opts.back_projection = false,
        // the WP14 proxy: barycenter as the only constraint
        "mg-bary" => opts.constraints = vec![Constraint::Barycenter { target: None }],
        _ => {}
    }
    // --soc 0|1: straight trials (before WP14b) or the second-order correction (default: the
    // library's)
    if let Some(v) = a.get("soc") {
        opts.second_order_correction = v != "0";
    }
    // --polish F: Newton polish below F times the tolerance (0: off; default: the library's)
    if let Some(p) = a.get("polish").and_then(|s| s.parse::<f64>().ok()) {
        opts.newton_polish = (p > 0.0).then_some(p);
    }
    // --abort 0|1: the early abort of hopeless chord projections (default: the library's)
    if let Some(v) = a.get("abort") {
        opts.projection_early_abort = v != "0";
    }
    // --proj-it N: chord iterations of the projection
    opts.projection.max_iterations = get(&a, "proj-it", opts.projection.max_iterations);
    let opts_soc = opts.second_order_correction;
    let mut note = match opts.newton_polish {
        Some(p) => format!("polish={p}"),
        None => String::new(),
    };
    let n = net.n_vertices();
    let load0 = load_avg();
    let t_all = Instant::now();
    let mut s = Session::new(net, pos, opts).expect("session");
    let mut rows = Vec::new();
    let mut cycles = Vec::new();
    let mut fallbacks = 0;
    let (mut solves, mut taus, mut viol) = (0usize, 0.0f64, 0.0f64);
    let (mut e_first, mut e_last) = (f64::NAN, f64::NAN);
    while rows.len() < max_steps {
        let k = rows.len();
        let el = t_all.elapsed().as_secs_f64();
        if k >= 1 {
            let mean = el / k as f64;
            if (k >= min_steps && el + mean > cap) || (k >= 2 && el + mean > 1.5 * cap) {
                break;
            }
        }
        match s.step() {
            Ok(r) => {
                if trace {
                    let trials: Vec<String> = r
                        .trials
                        .iter()
                        .map(|t| {
                            format!(
                                "{:.3}:{}:{}",
                                t.tau / r.tau_start.max(f64::MIN_POSITIVE),
                                t.projection_iterations,
                                match t.outcome {
                                    TrialOutcome::Accepted => "ok".to_string(),
                                    TrialOutcome::ArmijoFailed => "armijo".to_string(),
                                    TrialOutcome::NonFinite => "nan".to_string(),
                                    TrialOutcome::PathCollision { .. } => "path".to_string(),
                                    TrialOutcome::ProjectionFailed(st) => format!("{st:?}"),
                                }
                            )
                        })
                        .collect();
                    eprintln!(
                        "step {} total {:.3} energy {:.3} metric {:.3} solve {:.3} ({} cyc) ls {:.3} soc {:.3} proj {:.3} \
                         ({} solves, {} cyc) ccd {:.3} f {:.3} tau {:.3e} tau0 {:.3e} viol {:.2e} E {:.6e} trials [{}]",
                        r.iteration,
                        r.seconds.total,
                        r.seconds.energy,
                        r.seconds.metric,
                        r.seconds.solve,
                        r.linear.iterations,
                        r.seconds.line_search,
                        r.seconds.second_order,
                        r.seconds.projection,
                        r.projection_solves,
                        r.projection_cycles,
                        r.seconds.ccd,
                        r.seconds.line_search_energy,
                        r.tau,
                        r.tau_start,
                        r.max_violation,
                        r.energy_after,
                        trials.join(" ")
                    );
                }
                fallbacks += usize::from(r.solver_fallback.is_some());
                solves += r.projection_solves + usize::from(r.second_order.is_some());
                taus += r.tau;
                viol = viol.max(r.max_violation);
                if e_first.is_nan() {
                    e_first = r.energy_before;
                }
                e_last = r.energy_after;
                cycles.push(r.linear.iterations as f64);
                rows.push(r.seconds);
            }
            Err(e) => {
                note = format!("{note} step {k}: {e}").replace(',', ";");
                break;
            }
        }
    }
    let m = |f: fn(&ddg_repulsive_curves::Timings) -> f64| -> f64 {
        rows.iter().map(f).sum::<f64>() / rows.len().max(1) as f64
    };
    let tot: Vec<f64> = rows.iter().map(|t| t.total).collect();
    println!(
        "{curve_name},{n},{config},{threads},{dense_threads},{},{:.4},{:.4},{:.4},{:.4},{:.4},{:.4},{:.4},{:.4},{:.4},{:.2},{fallbacks},{:.2},{load0:.1},{:.1},{},{:.4},{:.2},{:.4e},{:.3e},{:.8e},{:.8e},{note}",
        rows.len(),
        m(|t| t.total),
        median(&tot),
        m(|t| t.energy),
        m(|t| t.metric),
        m(|t| t.solve),
        m(|t| t.line_search),
        m(|t| t.line_search_energy),
        m(|t| t.projection),
        m(|t| t.ccd),
        cycles.iter().sum::<f64>() / cycles.len().max(1) as f64,
        peak_rss_gb(),
        load_avg(),
        u8::from(opts_soc),
        m(|t| t.second_order),
        solves as f64 / rows.len().max(1) as f64,
        taus / rows.len().max(1) as f64,
        viol,
        e_first,
        e_last,
    );
}
