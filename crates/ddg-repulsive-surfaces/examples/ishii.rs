//! Fig. 16 of "Repulsive Surfaces" (Yu, Brakensiek, Schumacher, Crane 2021): the genus-two
//! handlebody-knots of Ishii, Kishimoto, Moriuchi, Suzuki 2012 as tubes around trivalent graphs,
//! flowed to minimizers.
//!
//! ```text
//! ishii list
//! ishii check                                   # transcription report for every entry
//! ishii run NAME [--faces N] [--threads T] [--max-steps S] [--window W] [--rel R] [--out DIR]
//! ```
//! `run`: area + barycenter constraints, default options (CCD on, remeshing), no attractor; stops when the mean of
//! the energy over the last `W` steps (default 200) differs from the mean over the `W` before by less than `R`
//! (default 1e-5) relative, or after `S` steps (default 3000). Prints one `RESULT` line.
#[path = "apps/apps.rs"]
mod apps;

use std::time::Instant;

use ddg_repulsive_surfaces::generate::ishii;
use ddg_repulsive_surfaces::step::{FlowOptions, Session, StopRule};
use ddg_repulsive_surfaces::{Constraint, ccd, write_obj};

fn main() {
    let args: Vec<String> = std::env::args().skip(1).collect();
    match args.first().map(String::as_str) {
        None | Some("list") => {
            for n in ishii::NAMES {
                println!("{n}");
            }
        }
        Some("check") => check(),
        Some("run") => run(&args[1..]),
        Some(o) => {
            eprintln!("unknown command {o}");
            std::process::exit(2);
        }
    }
}

fn check() {
    for e in ishii::all_entries() {
        let status = match ishii::check_projection(&e) {
            Ok(()) => "OK".to_string(),
            Err(m) => format!("MISMATCH {m}"),
        };
        println!("{} | projection {} | clearance {:.3} radius {:.3}", ishii::describe(&e), status, ishii::clearance(&e), ishii::tube_radius(&e));
    }
}

fn run(args: &[String]) {
    let name = args[0].clone();
    let (mut faces, mut threads, mut max_steps, mut window, mut rel, mut out) = (6000usize, 4usize, 3000usize, 200usize, 1e-5, None);
    let mut i = 1;
    while i < args.len() {
        let v = args.get(i + 1).cloned().expect("option needs a value");
        match args[i].as_str() {
            "--faces" => faces = v.parse().unwrap(),
            "--threads" => threads = v.parse().unwrap(),
            "--max-steps" => max_steps = v.parse().unwrap(),
            "--window" => window = v.parse().unwrap(),
            "--rel" => rel = v.parse().unwrap(),
            "--out" => out = Some(std::path::PathBuf::from(v)),
            o => panic!("unknown option {o}"),
        }
        i += 2;
    }
    let fail = |m: String| -> ! {
        println!("RESULT name={name} status=ERROR error={m}");
        std::process::exit(1)
    };
    let e = ishii::entry(&name).unwrap_or_else(|x| fail(x.to_string()));
    let (mesh, radius) = ishii::tube(&e, faces).unwrap_or_else(|x| fail(x.to_string()));
    let surface = mesh.surface().unwrap_or_else(|x| fail(x.to_string()));
    let faces0 = surface.n_faces();
    let genus0 = 1 - surface.euler_characteristic() / 2;
    let mut o = FlowOptions { max_iterations: max_steps, stop: StopRule::Off, ..FlowOptions::default() };
    o.threads = Some(threads);
    o.constraints = vec![Constraint::Barycenter { component: None, target: None }, Constraint::Area { component: None, target: None }];
    let t0 = Instant::now();
    let mut s = Session::new(surface, mesh.positions.clone(), o).unwrap_or_else(|x| fail(x.to_string()));
    let (mut e0, mut e1) = (f64::NAN, f64::NAN);
    let mut hist: Vec<f64> = Vec::new();
    let mut status = "MaxSteps".to_string();
    while s.iteration() < max_steps {
        match s.step() {
            Ok(r) => {
                if e0.is_nan() {
                    e0 = r.energy_before;
                }
                e1 = r.energy_after;
                hist.push(e1);
                if s.status().is_terminal() {
                    status = format!("{:?}", s.status());
                    break;
                }
                let n = hist.len();
                if n >= 2 * window {
                    let m = |a: &[f64]| a.iter().sum::<f64>() / a.len() as f64;
                    let (last, prev) = (m(&hist[n - window..]), m(&hist[n - 2 * window..n - window]));
                    if (prev - last).abs() <= rel * last.abs() {
                        status = "Converged".into();
                        break;
                    }
                }
            }
            Err(x) => {
                status = format!("Error({x})");
                break;
            }
        }
    }
    let secs = t0.elapsed().as_secs_f64();
    let (_, _, desc) = apps::inertia_indicators(s.surface().faces(), s.positions());
    let inter = ccd::intersects(s.surface(), s.positions()).map_or("n/a".to_string(), |b| b.to_string());
    let genus1 = 1 - s.surface().euler_characteristic() / 2;
    if let Some(d) = out {
        std::fs::create_dir_all(&d).unwrap();
        write_obj(d.join(format!("{name}_start.obj")), &mesh.surface().unwrap(), &mesh.positions).unwrap();
        write_obj(d.join(format!("{name}_final.obj")), s.surface(), s.positions()).unwrap();
    }
    println!(
        "RESULT name={name} kind={:?} crossings={} radius={radius:.3} faces={faces0}->{} genus={genus0}->{genus1} steps={} status={status} seconds={secs:.1} energy={e0:.5e}->{e1:.5e} intersects={inter} {desc}",
        e.kind,
        e.crossings,
        s.surface().n_faces(),
        s.iteration()
    );
}
