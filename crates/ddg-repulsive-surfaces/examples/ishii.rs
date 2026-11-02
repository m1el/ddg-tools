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
//! (default 1e-5) relative, or after `S` steps (default 3000), or when the file named by `ISHII_STOP_FILE` exists.
//! With `--out DIR` it also writes `DIR/steps.csv` (one line per step, `step::CSV_HEADER`).
//! `genus G [--threads T] [--max-steps S] [--window W] [--rel R] [--out DIR]`: the WP15 genus-G run (start
//! `shapes::genus_surface(G, 3000 + 600 G)`, area + barycenter, default options) with the same stop rule
//! (default cap 8000), OBJ frames every 50 steps and `final.obj` in DIR.
//! Prints `PROGRESS` lines every 25 steps and one `RESULT` line.
#[path = "apps/apps.rs"]
mod apps;

use std::time::Instant;

use ddg_repulsive_surfaces::generate::{ishii, shapes};
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
        Some("genus") => genus(&args[1..]),
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
    let (status, e0, e1) = drive(&mut s, max_steps, window, rel, out.as_deref(), 0, t0);
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

/// The step loop shared by `run` and `genus`: windowed stop rule, `PROGRESS` lines every 25 steps, per-step
/// CSV and (if `frame_every > 0`) OBJ frames in `out`. Returns (status, first energy, last energy).
fn drive(s: &mut Session, max_steps: usize, window: usize, rel: f64, out: Option<&std::path::Path>, frame_every: usize, t0: Instant) -> (String, f64, f64) {
    use std::io::Write;
    let mut csv = out.map(|d| {
        std::fs::create_dir_all(d).unwrap();
        let mut f = std::io::BufWriter::new(std::fs::File::create(d.join("steps.csv")).unwrap());
        writeln!(f, "{}", ddg_repulsive_surfaces::step::CSV_HEADER).unwrap();
        f
    });
    if let (Some(d), true) = (out, frame_every > 0) {
        write_obj(d.join("frame_00000.obj"), s.surface(), s.positions()).unwrap();
    }
    let (mut e0, mut e1) = (f64::NAN, f64::NAN);
    let mut hist: Vec<f64> = Vec::new();
    let mut status = "MaxSteps".to_string();
    let stop_file = std::env::var("ISHII_STOP_FILE").ok();
    while s.iteration() < max_steps {
        if stop_file.as_ref().is_some_and(|f| std::path::Path::new(f).exists()) {
            status = "StoppedByFile".into();
            break;
        }
        match s.step() {
            Ok(r) => {
                if e0.is_nan() {
                    e0 = r.energy_before;
                }
                e1 = r.energy_after;
                hist.push(e1);
                if let Some(f) = csv.as_mut() {
                    writeln!(f, "{}", ddg_repulsive_surfaces::step::csv_line(&r)).unwrap();
                    f.flush().unwrap();
                }
                let it = s.iteration();
                if it.is_multiple_of(25) {
                    println!("PROGRESS step={} faces={} energy={e1:.6e} t={:.0}s", it, s.surface().n_faces(), t0.elapsed().as_secs_f64());
                }
                if let (Some(d), true) = (out, frame_every > 0 && it.is_multiple_of(frame_every.max(1))) {
                    write_obj(d.join(format!("frame_{it:05}.obj")), s.surface(), s.positions()).unwrap();
                }
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
    (status, e0, e1)
}

fn genus(args: &[String]) {
    let g: u32 = args[0].parse().unwrap();
    let (mut threads, mut max_steps, mut window, mut rel, mut out) = (4usize, 8000usize, 200usize, 1e-5, None);
    let mut i = 1;
    while i < args.len() {
        let v = args.get(i + 1).cloned().expect("option needs a value");
        match args[i].as_str() {
            "--threads" => threads = v.parse().unwrap(),
            "--max-steps" => max_steps = v.parse().unwrap(),
            "--window" => window = v.parse().unwrap(),
            "--rel" => rel = v.parse().unwrap(),
            "--out" => out = Some(std::path::PathBuf::from(v)),
            o => panic!("unknown option {o}"),
        }
        i += 2;
    }
    let mesh = shapes::genus_surface(g, 3000 + 600 * g as usize).unwrap();
    let surface = mesh.surface().unwrap();
    let faces0 = surface.n_faces();
    let mut o = FlowOptions { max_iterations: max_steps, stop: StopRule::Off, ..FlowOptions::default() };
    o.threads = Some(threads);
    o.constraints = vec![Constraint::Barycenter { component: None, target: None }, Constraint::Area { component: None, target: None }];
    let t0 = Instant::now();
    let mut s = Session::new(surface, mesh.positions.clone(), o).unwrap();
    let (status, e0, e1) = drive(&mut s, max_steps, window, rel, out.as_deref(), 50, t0);
    let (_, _, desc) = apps::inertia_indicators(s.surface().faces(), s.positions());
    let inter = ccd::intersects(s.surface(), s.positions()).map_or("n/a".to_string(), |b| b.to_string());
    if let Some(d) = out {
        write_obj(d.join("final.obj"), s.surface(), s.positions()).unwrap();
    }
    println!(
        "RESULT genus={g} faces={faces0}->{} genus_out={} steps={} status={status} seconds={:.1} energy={e0:.5e}->{e1:.5e} intersects={inter} {desc}",
        s.surface().n_faces(),
        1 - s.surface().euler_characteristic() / 2,
        s.iteration(),
        t0.elapsed().as_secs_f64()
    );
}
