//! Cost of the near-pair quadrature (WP16, `SPEC.md` §3.4, §17o): seconds per iteration of a
//! running flow and the part spent in the near-pair passes (`Timings::near_pair`: the base curve
//! and every trial energy of the line search), one configuration per invocation, one CSV row on
//! stdout (header with `--header 1`). Run alone (single-thread timings).
//!
//! ```text
//! cargo bench --bench nearpair -- (--curve fhw|fhw-r|knot --size N | --obj FILE.obj)
//!     --config exact|bh|mg [--near 0|1] [--threads 1] [--steps 10] [--alpha A --beta B]
//! ```
//!
//! Exponents (2, 4.5) unless given; barycenter + edge lengths; the library's defaults otherwise.

mod common;

use std::time::Instant;

use common::*;
use ddg_repulsive_curves::energy::NearPairOptions;
use ddg_repulsive_curves::{
    Backend, BhOptions, EnergyEval, Exponents, FlowOptions, MgOptions, Session, Solver,
    read_curve_obj,
};

fn main() {
    let a = args();
    if a.contains_key("header") {
        println!(
            "curve,n,config,near,threads,steps,mean_total,mean_energy,mean_ls_energy,mean_near_pair,\
             near_share,mean_near_pairs,max_near_pairs,energy_first,energy_last,load_start,load_end,note"
        );
        return;
    }
    let (name, (net, pos)) = if let Some(f) = a.get("obj") {
        let c = read_curve_obj(f).unwrap_or_else(|e| {
            eprintln!("{f}: {e}");
            std::process::exit(2)
        });
        (
            std::path::Path::new(f)
                .file_stem()
                .map(|s| s.to_string_lossy().to_string())
                .unwrap_or_default(),
            c,
        )
    } else {
        let curve_name = a.get("curve").cloned().unwrap_or("fhw".into());
        let size: usize = get(&a, "size", 600);
        match curve(&curve_name, size) {
            Ok(c) => (curve_name, c),
            Err(e) => {
                eprintln!("{e}");
                std::process::exit(2);
            }
        }
    };
    let config = a.get("config").cloned().unwrap_or("exact".into());
    let near = get(&a, "near", 1usize) != 0;
    let threads: usize = get(&a, "threads", 1);
    let steps: usize = get(&a, "steps", 10);
    let alpha: f64 = get(&a, "alpha", 2.0);
    let beta: f64 = get(&a, "beta", 4.5);
    let (energy, solver) = match config.as_str() {
        "exact" => (EnergyEval::Exact, Solver::Dense),
        "bh" => (EnergyEval::BarnesHut(BhOptions::default()), Solver::Dense),
        _ => (
            EnergyEval::BarnesHut(BhOptions::default()),
            Solver::Multigrid(MgOptions::default()),
        ),
    };
    let opts = FlowOptions {
        exponents: Exponents::new(alpha, beta).unwrap(),
        backend: Backend::Explicit { energy, solver },
        threads: Some(threads),
        dense_threads: Some(threads.min(16)),
        dense_fallback_max_vertices: 0,
        near_pair: near.then(NearPairOptions::default),
        ..FlowOptions::default()
    };
    let n = net.n_vertices();
    let load0 = load_avg();
    let mut s = Session::new(net, pos, opts).expect("session");
    let t0 = Instant::now();
    let mut rows = Vec::new();
    let mut pairs = Vec::new();
    let (mut e_first, mut e_last) = (f64::NAN, f64::NAN);
    let mut note = String::new();
    while rows.len() < steps {
        match s.step() {
            Ok(r) => {
                if e_first.is_nan() {
                    e_first = r.energy_before;
                }
                e_last = r.energy_after;
                pairs.push(r.near_pairs as f64);
                rows.push(r.seconds);
            }
            Err(e) => {
                note = format!("step {}: {e}", rows.len()).replace(',', ";");
                break;
            }
        }
    }
    let _ = t0;
    let m = |f: fn(&ddg_repulsive_curves::Timings) -> f64| -> f64 {
        rows.iter().map(f).sum::<f64>() / rows.len().max(1) as f64
    };
    let total = m(|t| t.total);
    let np = m(|t| t.near_pair);
    println!(
        "{name},{n},{config},{},{threads},{},{:.4},{:.4},{:.4},{:.5},{:.4},{:.1},{},{:.8e},{:.8e},{load0:.1},{:.1},{note}",
        u8::from(near),
        rows.len(),
        total,
        m(|t| t.energy),
        m(|t| t.line_search_energy),
        np,
        np / total.max(f64::MIN_POSITIVE),
        pairs.iter().sum::<f64>() / pairs.len().max(1) as f64,
        pairs.iter().cloned().fold(0.0, f64::max),
        e_first,
        e_last,
        load_avg(),
    );
}
