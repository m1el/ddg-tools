//! Comparison of the optimizers (`SPEC.md` §14.2, WP9 sanity run, WP16B evaluation): one method on
//! one OBJ, energy against wall time.
//!
//! ```text
//! compare_methods OBJ METHOD [STEPS] [THREADS] [OPTIONS]
//! compare_methods gen SHAPE OUT.obj
//! METHOD = hs | l2 | h1 | h2 | aqp | h1-lbfgs | h1-lbfgs-ref | bcqn | bcqn-ref
//! SHAPE  = fig12:N (N = 1..5, the five Fig. 12 analogues) | genus:G:FACES | handcuffs:K:FACES
//! OPTIONS (WP16B):
//!   --budget SECS         wall budget of the stepping (checkpoints not counted); STEPS = 0: no step limit
//!   --csv PATH            one row per step (default: space-separated rows on stdout)
//!   --exact-time T1,T2,.. exact energy at the first step at/after each wall time
//!   --exact-steps S1,S2,. exact energy after those steps
//!   --exact-end           exact energy at the end (summary line)
//! ```
//!
//! `hs`, `l2`, `h1`, `h2` use hard barycenter + area constraints; the others the soft problem
//! (`soft_problem`: area penalty, weight 100). CSV columns: `step wall E_bh F area faces stuck
//! gmres_per_solve solves tau E_exact` (`E_bh` = the unweighted Barnes-Hut energy after the step's
//! projection, `E_exact` only at checkpoints). The exact evaluations (1 thread) are excluded from `wall`.
//! The last line is `summary ...` with the ending (`budget`, `steps`, `Stuck` after 5 stuck steps in a
//! row, or `Error: ...`).

use std::io::Write;
use std::time::Instant;

use ddg_repulsive_surfaces::constraints::Constraint;
use ddg_repulsive_surfaces::optimize::soft_problem;
use ddg_repulsive_surfaces::{
    BcqnVariant, EnergyEval, Exponents, FaceGeometry, FlowOptions, L2Kind, LbfgsVariant, Method, Metric, Session,
    StepStatus, StopRule, generate, read_obj, write_obj,
};

fn flag<'a>(a: &'a [String], name: &str) -> Option<&'a str> {
    a.iter().position(|x| x == name).and_then(|i| a.get(i + 1)).map(String::as_str)
}

fn list(a: &[String], name: &str) -> Vec<f64> {
    flag(a, name).map_or(vec![], |s| s.split(',').map(|x| x.parse().expect(name)).collect())
}

fn exact_energy(ses: &Session) -> f64 {
    ddg_repulsive_surfaces::energy(ses.surface(), ses.positions(), Exponents::default(), &EnergyEval::Exact, &[], Some(1))
        .expect("exact energy")
        .energy
}

fn area_of(ses: &Session) -> f64 {
    FaceGeometry::new(ses.surface(), ses.positions()).unwrap().areas.iter().sum()
}

fn generate_obj(a: &[String]) {
    let spec: Vec<&str> = a[2].split(':').collect();
    let mesh = match spec[0] {
        "fig12" => generate::shapes::fig12(spec[1].parse::<usize>().expect("N") - 1).expect("fig12").mesh,
        "genus" => generate::shapes::genus_surface(spec[1].parse().expect("G"), spec[2].parse().expect("FACES")).expect("genus"),
        "handcuffs" => {
            generate::shapes::handcuffs(spec[1].parse().expect("K"), spec[2].parse().expect("FACES")).expect("handcuffs").mesh
        }
        _ => panic!("unknown shape {}", a[2]),
    };
    let surface = mesh.surface().expect("surface");
    write_obj(&a[3], &surface, &mesh.positions).expect("write");
    println!("{} faces {}", a[3], surface.n_faces());
}

fn main() {
    let a: Vec<String> = std::env::args().collect();
    if a[1] == "gen" {
        return generate_obj(&a);
    }
    let (obj, name) = (&a[1], a[2].as_str());
    let num = |i: usize| a.get(i).filter(|s| !s.starts_with("--"));
    let steps: usize = num(3).map_or(50, |s| s.parse().expect("STEPS"));
    let threads: usize = num(4).map_or(4, |s| s.parse().expect("THREADS"));
    let budget: f64 = flag(&a, "--budget").map_or(f64::INFINITY, |s| s.parse().expect("budget"));
    let ex_time = list(&a, "--exact-time");
    let ex_steps = list(&a, "--exact-steps");
    let mut csv = flag(&a, "--csv").map(|p| std::fs::File::create(p).expect("csv"));

    let (surface, positions) = read_obj(obj).expect("read OBJ");
    let hard = FlowOptions {
        constraints: vec![
            Constraint::Barycenter { component: None, target: None },
            Constraint::Area { component: None, target: None },
        ],
        stop: StopRule::Off,
        threads: Some(threads),
        ..FlowOptions::default()
    };
    let mut o = hard.clone();
    match name {
        "hs" => {}
        "l2" => o.metric = Metric::L2(L2Kind::Mass),
        "h1" => o.metric = Metric::H1,
        "h2" => o.metric = Metric::H2,
        _ => {
            o = soft_problem(&surface, &positions, &hard).expect("soft problem");
            o.method = match name {
                "aqp" => Method::Aqp { kappa: 100.0 },
                "h1-lbfgs" => Method::H1Lbfgs { memory: 20, variant: LbfgsVariant::Standard },
                "h1-lbfgs-ref" => Method::H1Lbfgs { memory: 20, variant: LbfgsVariant::Reference },
                "bcqn" => Method::Bcqn { memory: 20, variant: BcqnVariant::Standard },
                "bcqn-ref" => Method::Bcqn { memory: 20, variant: BcqnVariant::Reference },
                _ => panic!("unknown method {name}"),
            };
        }
    }
    let faces0 = surface.n_faces();
    let mut ses = Session::new(surface, positions, o).expect("session");
    if let Some(f) = csv.as_mut() {
        writeln!(f, "step,wall,E_bh,F,area,faces,stuck,gmres_per_solve,solves,tau,E_exact").unwrap();
    }
    let area0 = area_of(&ses);
    let (mut clock, mut e0, mut e_last) = (0.0f64, f64::NAN, f64::NAN);
    let (mut next_t, mut done, mut stuck_run, mut max_stuck, mut n_stuck) = (0usize, 0usize, 0usize, 0usize, 0usize);
    let (mut git, mut gsolves) = (0usize, 0usize);
    let mut ending = "budget".to_string();
    loop {
        if steps > 0 && done >= steps {
            ending = "steps".into();
            break;
        }
        if clock >= budget {
            break;
        }
        let t = Instant::now();
        let r = match ses.step() {
            Ok(r) => r,
            Err(e) => {
                ending = format!("Error: {e}");
                break;
            }
        };
        clock += t.elapsed().as_secs_f64();
        done += 1;
        let area = area_of(&ses);
        if e0.is_nan() {
            e0 = r.energy_before;
        }
        e_last = r.energy_term_after;
        let st = r.status == StepStatus::Stuck;
        n_stuck += usize::from(st);
        stuck_run = if st { stuck_run + 1 } else { 0 };
        max_stuck = max_stuck.max(stuck_run);
        git += r.linear.total_iterations();
        gsolves += r.linear.solves.len();
        let mut due = ex_steps.contains(&(done as f64));
        while next_t < ex_time.len() && clock >= ex_time[next_t] {
            next_t += 1;
            due = true;
        }
        let ex = if due { format!("{:.8e}", exact_energy(&ses)) } else { String::new() };
        let per = r.linear.total_iterations() as f64 / r.linear.solves.len().max(1) as f64;
        let line = format!(
            "{},{:.3},{:.6e},{:.6e},{:.6e},{},{},{:.2},{},{:.4e},{}",
            r.iteration, clock, e_last, r.energy_after, area, r.faces, u8::from(st), per, r.linear.solves.len(), r.tau, ex
        );
        match csv.as_mut() {
            Some(f) => writeln!(f, "{line}").unwrap(),
            None => println!("{}", line.replace(',', " ")),
        }
        if stuck_run >= 5 {
            ending = "Stuck".into();
            break;
        }
    }
    let exact_end = if a.iter().any(|x| x == "--exact-end") && !ending.starts_with("Error") { exact_energy(&ses) } else { f64::NAN };
    println!(
        "summary method={name} faces0={faces0} faces={} steps={done} seconds={clock:.2} ending={ending} F0={e0:.5e} E_bh={e_last:.5e} E_exact_end={exact_end:.8e} area0={area0:.6e} area_ratio={:.4} stuck_steps={n_stuck} max_stuck_run={max_stuck} gmres_iters_per_solve={:.2} diag={:?}",
        ses.surface().n_faces(),
        area_of(&ses) / area0,
        git as f64 / gsolves.max(1) as f64,
        ses.method_diagnostics()
    );
}
