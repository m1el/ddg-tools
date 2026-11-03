//! The trial rule of the line search on one knot (`SPEC.md` §17l): straight trials or the
//! second-order correction, with or without the Newton polish, as the start-rule study of §17f
//! (`tests/flow.rs`, `ls_compare_one_knot`): `(α, β) = (3, 6)`, `Absolute(1e-4)`, 3 chord
//! iterations, default back end, a CCD audit of every accepted path (`path_contact`, gap 1e-6).
//!
//! ```text
//! cargo bench --bench trialrule -- --knot PATH [--soc 0|1] [--polish F] [--abort 0|1] [--steps 100] [--threads 4]
//!     [--library-defaults 1]
//! ```
//!
//! Prints `TRAJ iteration seconds energy |g| trials solves backend` per step (seconds without the
//! audit) and a final `RESULT` line.

mod common;

use std::time::Instant;

use common::*;
use ddg_repulsive_curves::{
    CcdOptions, CcdReport, Exponents, FlowOptions, ProjectionOptions, Session, StepStatus,
    Tolerance, read_curve_obj,
};

fn main() {
    let a = args();
    if let Some(p) = a.get("eval") {
        // --eval PATH: exact and Barnes-Hut energies, length, bounding box of a curve file
        let (net, pos) = read_curve_obj(p).expect("curve file");
        let ex = Exponents::new(3.0, 6.0).unwrap();
        let e = ddg_repulsive_curves::energy(
            &net,
            &pos,
            ex,
            &ddg_repulsive_curves::EnergyEval::Exact,
            Some(16),
        )
        .unwrap();
        let b =
            ddg_repulsive_curves::energy(
                &net,
                &pos,
                ex,
                &ddg_repulsive_curves::EnergyEval::BarnesHut(
                    ddg_repulsive_curves::BhOptions::default(),
                ),
                Some(16),
            )
            .unwrap();
        let g = net.edge_geometry(&pos);
        let (lo, hi) = pos
            .iter()
            .fold(([f64::MAX; 3], [f64::MIN; 3]), |(lo, hi), x| {
                (
                    [0, 1, 2].map(|c| lo[c].min(x[c])),
                    [0, 1, 2].map(|c| hi[c].max(x[c])),
                )
            });
        println!(
            "EVAL exact {:.6} bh {:.6} length {:.4} min edge {:.4e} max edge {:.4e} bbox {:?} {:?}",
            e.energy,
            b.energy,
            g.length.iter().sum::<f64>(),
            g.length.iter().cloned().fold(f64::MAX, f64::min),
            g.length.iter().cloned().fold(0.0, f64::max),
            lo,
            hi
        );
        return;
    }
    let path = a.get("knot").cloned().expect("--knot PATH");
    let soc = get(&a, "soc", 0usize) != 0;
    let polish: f64 = get(&a, "polish", 0.0);
    let steps: usize = get(&a, "steps", 100);
    let abort = get(&a, "abort", 0usize) != 0;
    let threads: usize = get(&a, "threads", 4);
    let (net, pos) = read_curve_obj(&path).expect("knot file");
    // --library-defaults 1: FlowOptions::default() (only the step cap, the threads and the
    // path for the audit set), the check of SPEC §17t; the flags above are then ignored.
    let library = get(&a, "library-defaults", 0usize) != 0;
    let opts = if library {
        FlowOptions {
            max_iterations: steps,
            record_path: true,
            threads: Some(threads),
            dense_threads: Some(threads),
            ..FlowOptions::default()
        }
    } else {
        FlowOptions {
            exponents: Exponents::new(3.0, 6.0).unwrap(),
            projection: ProjectionOptions {
                tolerance: Tolerance::Absolute(1e-4),
                max_iterations: 3,
            },
            second_order_correction: soc,
            newton_polish: (polish > 0.0).then_some(polish),
            projection_early_abort: abort,
            max_iterations: steps,
            record_path: true,
            threads: Some(threads),
            dense_threads: Some(threads),
            ..FlowOptions::default()
        }
    };
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(threads)
        .build()
        .unwrap();
    let audit = CcdOptions {
        gap: 1e-6,
        max_advances: 100_000,
    };
    let t0 = Instant::now();
    let mut s = Session::new(net.clone(), pos, opts).expect("session");
    let e0 = s.energy();
    let (mut audit_s, mut trials, mut solves, mut stuck, mut violations) = (0.0, 0, 0, 0, 0);
    let mut error = String::new();
    loop {
        let r = match s.step() {
            Ok(r) => r,
            Err(e) => {
                error = format!(" error {e}").replace('\n', " ");
                break;
            }
        };
        trials += r.trials.len();
        let sv = r.projection_solves + usize::from(r.second_order.is_some());
        solves += sv;
        if r.status == StepStatus::Stuck {
            stuck += 1;
        }
        println!(
            "TRAJ {} {:.4} {:.10} {:.3e} {} {} {} {:.3e} {:.3e}",
            r.iteration,
            t0.elapsed().as_secs_f64() - audit_s,
            r.energy_after,
            r.gradient_norm,
            r.trials.len(),
            sv,
            r.solver_backend,
            r.max_violation,
            r.tau
        );
        if let Some(p) = &r.path {
            let ta = Instant::now();
            let mut rep = CcdReport::default();
            let hit = ddg_repulsive_curves::ccd::path_contact(
                &pool,
                s.network(),
                p,
                &[],
                &audit,
                &mut rep,
            )
            .unwrap();
            audit_s += ta.elapsed().as_secs_f64();
            if hit.is_some() {
                violations += 1;
            }
        }
        if s.status().is_terminal() {
            break;
        }
    }
    if let Some(out) = a.get("out") {
        ddg_repulsive_curves::write_curve_obj(out, s.network(), s.positions()).unwrap();
    }
    println!(
        "RESULT status {:?} iterations {} seconds {:.3} energy0 {:.6} energy {:.6} trials {} solves {} \
         stuck_steps {} violations {} vertices {}{}",
        s.status(),
        s.iteration(),
        t0.elapsed().as_secs_f64() - audit_s,
        e0,
        s.energy(),
        trials,
        solves,
        stuck,
        violations,
        net.n_vertices(),
        error
    );
}
