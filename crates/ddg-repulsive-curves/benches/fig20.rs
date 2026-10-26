//! Fig. 20 scaling (`SPEC.md` §14.1, §17k): mean time per iteration of a running flow, one
//! configuration per invocation, one CSV row on stdout (header with `--header 1`).
//!
//! ```text
//! cargo bench --bench fig20 -- --curve fhw|fhw-r|knot --size N --config exact|bh|mg|mg-noproj
//!     [--threads T] [--dense-threads D] [--steps 20] [--min-steps 5] [--cap-s 240]
//! ```
//!
//! `exact`: exact energy + dense; `bh`: Barnes–Hut energy + dense; `mg`: BH + block cluster tree +
//! multigrid; `mg-noproj`: the same with the barycenter as the only constraint (no edge-length
//! rows, so the line search's projections need no solve). Exponents (2, 4.5) on `fhw`/`fhw-r`,
//! (3, 6) otherwise (`--alpha/--beta` override). Barycenter + edge lengths otherwise.

mod common;

use std::time::Instant;

use common::*;
use ddg_repulsive_curves::{
    Backend, BhOptions, Constraint, EnergyEval, Exponents, FlowOptions, MgOptions, Session, Solver,
};

fn main() {
    let a = args();
    if a.contains_key("header") {
        println!(
            "curve,n,config,threads,dense_threads,steps,mean_total,median_total,energy,metric,solve,\
             line_search,ls_energy,projection,ccd,mean_mg_cycles,fallbacks,peak_rss_gb,load_start,load_end,note"
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
        _ => (
            EnergyEval::BarnesHut(BhOptions::default()),
            Solver::Multigrid(MgOptions::default()),
        ),
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
    if config == "mg-noproj" {
        opts.constraints = vec![Constraint::Barycenter { target: None }];
    }
    let n = net.n_vertices();
    let load0 = load_avg();
    let mut note = String::new();
    let t_all = Instant::now();
    let mut s = Session::new(net, pos, opts).expect("session");
    let mut rows = Vec::new();
    let mut cycles = Vec::new();
    let mut fallbacks = 0;
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
                fallbacks += usize::from(r.solver_fallback.is_some());
                cycles.push(r.linear.iterations as f64);
                rows.push(r.seconds);
            }
            Err(e) => {
                note = format!("step {k}: {e}").replace(',', ";");
                break;
            }
        }
    }
    let m = |f: fn(&ddg_repulsive_curves::Timings) -> f64| -> f64 {
        rows.iter().map(f).sum::<f64>() / rows.len().max(1) as f64
    };
    let tot: Vec<f64> = rows.iter().map(|t| t.total).collect();
    println!(
        "{curve_name},{n},{config},{threads},{dense_threads},{},{:.4},{:.4},{:.4},{:.4},{:.4},{:.4},{:.4},{:.4},{:.4},{:.2},{fallbacks},{:.2},{load0:.1},{:.1},{note}",
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
    );
}
