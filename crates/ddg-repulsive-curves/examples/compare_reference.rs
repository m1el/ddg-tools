//! Final-shape comparison with the authors' program (`SPEC.md` §14.7, README "Comparison with the
//! authors' program").
//!
//! ```text
//! compare_reference flow   IN.obj OUT.obj [--alpha A] [--beta B] [--threads T] [--max-steps N]
//!                          [--time-limit-s S] [--log steps.csv]
//! compare_reference energy FILE.obj [--alpha A] [--beta B] [--threads T]
//! ```
//!
//! `flow` runs our flow with the settings of the reference's batch mode as far as they can be
//! matched: the exponents, the constraints (all edge lengths at their initial values and the
//! barycenter), the exact energy and the dense solver, the absolute projection tolerance 1e-4;
//! the stop rule (gradient norm below 1e-4) and the line search are ours. It writes the final
//! curve and a CSV with one line per step. `energy` prints the exact energy (Eq. 17, our edge-pair
//! discretization) of any curve file, e.g. of an output of the reference.

use std::io::Write;

use ddg_repulsive_curves::solve::{ProjectionOptions, Tolerance};
use ddg_repulsive_curves::{
    Backend, Control, EnergyEval, Exponents, FlowOptions, Session, Solver, StepReport, energy,
    read_curve_obj, write_curve_obj,
};

struct Args {
    alpha: f64,
    beta: f64,
    threads: Option<usize>,
    max_steps: usize,
    time_limit_s: Option<f64>,
    log: Option<String>,
    positional: Vec<String>,
}

fn parse(args: impl Iterator<Item = String>) -> Result<Args, String> {
    let mut a = Args {
        alpha: 3.0,
        beta: 6.0,
        threads: None,
        max_steps: 1000,
        time_limit_s: None,
        log: None,
        positional: Vec::new(),
    };
    let mut it = args.peekable();
    while let Some(s) = it.next() {
        let mut value = |name: &str| it.next().ok_or(format!("{name} needs a value"));
        match s.as_str() {
            "--alpha" => a.alpha = value("--alpha")?.parse().map_err(|e| format!("{e}"))?,
            "--beta" => a.beta = value("--beta")?.parse().map_err(|e| format!("{e}"))?,
            "--threads" => {
                a.threads = Some(value("--threads")?.parse().map_err(|e| format!("{e}"))?)
            }
            "--max-steps" => {
                a.max_steps = value("--max-steps")?.parse().map_err(|e| format!("{e}"))?
            }
            "--time-limit-s" => {
                a.time_limit_s = Some(
                    value("--time-limit-s")?
                        .parse()
                        .map_err(|e| format!("{e}"))?,
                )
            }
            "--log" => a.log = Some(value("--log")?),
            _ => a.positional.push(s),
        }
    }
    Ok(a)
}

fn main() {
    if let Err(e) = run() {
        eprintln!("error: {e}");
        std::process::exit(1);
    }
}

fn run() -> Result<(), Box<dyn std::error::Error>> {
    let mut argv = std::env::args().skip(1);
    let cmd = argv
        .next()
        .ok_or("usage: compare_reference flow|energy ...")?;
    let a = parse(argv)?;
    let exponents = Exponents::new(a.alpha, a.beta)?;
    match cmd.as_str() {
        "energy" => {
            let file = a.positional.first().ok_or("energy needs a file")?;
            let (net, pos) = read_curve_obj(file)?;
            let r = energy(&net, &pos, exponents, &EnergyEval::Exact, a.threads)?;
            println!("{:.17e}", r.energy);
        }
        "flow" => {
            let [input, output] = [0, 1].map(|i| a.positional.get(i).cloned());
            let (input, output) = (
                input.ok_or("flow needs IN.obj OUT.obj")?,
                output.ok_or("flow needs IN.obj OUT.obj")?,
            );
            let (net, pos) = read_curve_obj(&input)?;
            let options = FlowOptions {
                exponents,
                backend: Backend::Explicit {
                    energy: EnergyEval::Exact,
                    solver: Solver::Dense,
                },
                projection: ProjectionOptions {
                    tolerance: Tolerance::Absolute(1e-4),
                    max_iterations: 3,
                },
                max_iterations: a.max_steps,
                time_limit: a.time_limit_s.map(std::time::Duration::from_secs_f64),
                threads: a.threads,
                dense_threads: a.threads,
                ..FlowOptions::default()
            };
            let mut log = match &a.log {
                Some(p) => {
                    let mut f = std::io::BufWriter::new(std::fs::File::create(p)?);
                    writeln!(f, "step,energy,gradient_norm,tau,halvings,max_violation")?;
                    Some(f)
                }
                None => None,
            };
            let mut session = Session::new(net, pos, options)?;
            let mut progress = |s: &StepReport| {
                if let Some(f) = log.as_mut() {
                    let _ = writeln!(
                        f,
                        "{},{:.17e},{:.6e},{:.6e},{},{:.3e}",
                        s.iteration,
                        s.energy_after,
                        s.gradient_norm,
                        s.tau,
                        s.halvings,
                        s.max_violation
                    );
                    let _ = f.flush();
                }
                Control::Continue
            };
            let report = session.run(Some(&mut progress))?;
            write_curve_obj(&output, session.network(), session.positions())?;
            println!(
                "status {:?} steps {} seconds {:.1} energy {:.17e}",
                report.status, report.iterations, report.seconds, report.final_energy
            );
        }
        other => return Err(format!("unknown command {other}").into()),
    }
    Ok(())
}
