//! Runs a reference scene file (`scene.txt` of `repulsive-surfaces/scenes`, fetched by
//! `scripts/rs-fetch-scenes.sh` into `data/repulsive-surfaces/scenes/`).
//!
//! ```text
//! cargo run --release -p ddg-repulsive-surfaces --example run_scene -- SCENE [options]
//!   --steps N        overrides the scene's iteration_limit
//!   --threads T      thread count (default: all cores)
//!   --no-remesh      keep the topology
//!   --ccd [GAP]      continuous collision detection in the line search (gap in mean edge lengths, 0.05)
//!   --no-stop        no energy-stall stop rule
//!   --csv PATH       per-step CSV log (default: the scene's `log`)
//!   --obj PATH       write the final mesh
//! ```
//!
//! Prints one line per step and a summary line
//! `summary scene=… faces0=… steps=… status=… energy0=… energy=… ms_per_step=… gmres_per_step=…`.

use std::path::PathBuf;

use ddg_repulsive_surfaces::step::{RunSceneOptions, csv_line, run_scene};
use ddg_repulsive_surfaces::CcdOptions;

fn main() {
    let mut args = std::env::args().skip(1);
    let Some(scene) = args.next() else {
        eprintln!("usage: run_scene SCENE [--steps N] [--threads T] [--no-remesh] [--ccd [GAP]] [--no-stop] [--csv PATH] [--obj PATH]");
        std::process::exit(2);
    };
    let mut o = RunSceneOptions::default();
    let rest: Vec<String> = args.collect();
    let mut i = 0;
    while i < rest.len() {
        let a = rest[i].as_str();
        let next = |i: usize| rest.get(i + 1).cloned().unwrap_or_else(|| panic!("{a} needs a value"));
        match a {
            "--steps" => {
                o.max_steps = Some(next(i).parse().expect("--steps N"));
                i += 1;
            }
            "--threads" => {
                o.threads = Some(next(i).parse().expect("--threads T"));
                i += 1;
            }
            "--no-remesh" => o.no_remesh = true,
            "--no-stop" => o.no_stop_rule = true,
            "--ccd" => {
                let mut c = CcdOptions::default();
                if let Some(g) = rest.get(i + 1).and_then(|s| s.parse::<f64>().ok()) {
                    c.gap = g;
                    i += 1;
                }
                o.ccd = Some(c);
            }
            "--csv" => {
                o.csv = Some(PathBuf::from(next(i)));
                i += 1;
            }
            "--obj" => {
                o.output_obj = Some(PathBuf::from(next(i)));
                i += 1;
            }
            other => panic!("unknown option {other}"),
        }
        i += 1;
    }
    let r = match run_scene(&scene, &o) {
        Ok(r) => r,
        Err(e) => {
            eprintln!("error: {e}");
            std::process::exit(1);
        }
    };
    for s in &r.flow.steps {
        println!("{}", csv_line(s));
    }
    let n = r.flow.steps.len().max(1) as f64;
    let ms: f64 = r.flow.steps.iter().map(|s| s.seconds.total).sum::<f64>() * 1e3 / n;
    let gm: f64 = r.flow.steps.iter().map(|s| s.linear.total_iterations() as f64).sum::<f64>() / n;
    let remeshes = r.flow.steps.iter().filter(|s| s.remesh.is_some()).count();
    let splits: u64 = r.flow.steps.iter().filter_map(|s| s.remesh.as_ref().map(|m| m.splits)).sum();
    let collapses: u64 = r.flow.steps.iter().filter_map(|s| s.remesh.as_ref().map(|m| m.collapses)).sum();
    let flips: u64 = r.flow.steps.iter().filter_map(|s| s.remesh.as_ref().map(|m| m.flips)).sum();
    let ccd_frac = {
        let mut v: Vec<f64> = r
            .flow
            .steps
            .iter()
            .filter(|s| s.seconds.total > 0.0)
            .map(|s| s.seconds.ccd / s.seconds.total)
            .collect();
        v.sort_by(f64::total_cmp);
        v.get(v.len() / 2).copied().unwrap_or(0.0)
    };
    let e0 = r.flow.steps.first().map_or(f64::NAN, |s| s.energy_before);
    let stuck = r
        .flow
        .steps
        .iter()
        .filter(|s| s.status == ddg_repulsive_surfaces::StepStatus::Stuck)
        .count();
    println!(
        "summary scene={} faces0={} faces={} steps={} stuck={} status={:?} energy0={:.6e} energy={:.6e} \
         ms_per_step={:.1} gmres_per_step={:.1} remeshes={} splits={} collapses={} flips={} ccd_median_fraction={:.4} seconds={:.1}",
        scene,
        r.initial_size.1,
        r.surface.n_faces(),
        r.flow.steps.len(),
        stuck,
        r.flow.status,
        e0,
        r.flow.final_energy,
        ms,
        gm,
        remeshes,
        splits,
        collapses,
        flips,
        ccd_frac,
        r.flow.seconds
    );
}
