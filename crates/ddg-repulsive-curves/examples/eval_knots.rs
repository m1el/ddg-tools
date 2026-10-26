//! Sec. 7 comparison driver (`SPEC.md` §14.2): one knot, one method, one process.
//!
//! ```text
//! eval_knots --knot FILE.obj --method NAME [--threads T] [--max-iters N] [--max-seconds S]
//!            [--alpha A] [--beta B] [--hencky-w W] [--backend auto|exact|bh-mg]
//!            [--no-audit] [--out RESULT.jsonl] [--final-obj FINAL.obj]
//!            [--frames-every-seconds S --frames-out FRAMES.rcf] [--tag TEXT]
//! ```
//!
//! Methods (names as in the tables): hard constraints (barycenter + edge lengths), collision
//! guarded line search: `hs-pg` (reference), `l2-pg`, `h1-pg`, `h2-pg`; `hs-armijo` (Hs projected
//! gradient with the normalized Armijo search); soft (`Ê` + Hencky, weight `--hencky-w`, default
//! 1e3): `soft-hs-ncg` (reference), `soft-hs-lbfgs` (memory 30), `soft-h1-ncg`, `soft-h1-lbfgs`.
//! Soft methods use the guarded search too (isotopy). Exponents (2, 4.5), projection tolerance
//! `Absolute(1e-4)`, stop `‖g̃‖₂ < 1e-4`.
//!
//! One JSON line is appended to `--out` (stdout when absent) with: knot, method, status
//! (`converged` | `stuck` | `nonconvergent` | `error`), `end_reason`, iterations, wall seconds
//! (without the audit and the frame writing), the history `[iteration, seconds, F, Ê, ‖g̃‖]`
//! (`F` the objective, `Ê` = `F` minus the penalties), exact initial and final energy (Eq. 17), the
//! isotopy audit (CCD on every accepted step, as `path_contact` in §17f), the constraint residual
//! (largest `|ℓ/ℓ⁰ − 1|`, the last step's `max_violation`, the largest `|ln(ℓ/ℓ⁰)|`), threads and
//! back ends. Progress lines `TRAJ it seconds F gnorm` go to stderr.
//!
//! Frames (`--frames-every-seconds`): little-endian file `RCF1`, `u32` vertices `n`, `u32` edges
//! `m`, `f32` interval, `u32` 0, then `m` pairs of `u32`, then records `f32` seconds, `u32`
//! iteration, `u32` flag (1 = final curve), `3n` `f32`. A record is written after the first step
//! that ends at least one interval after the previous record (plus the start and the final curve),
//! so a renderer picks, for a video time `t`, the last record with seconds `<= t`.

use std::io::Write;
use std::path::Path;
use std::time::Instant;

use ddg_repulsive_curves::ccd::{CcdOptions, CcdReport, path_contact};
use ddg_repulsive_curves::optimize::{Hencky, Method};
use ddg_repulsive_curves::solve::{ProjectionOptions, Tolerance};
use ddg_repulsive_curves::step::stop::{Status, StepStatus};
use ddg_repulsive_curves::{
    Backend, CurveNetwork, EnergyEval, Exponents, FlowOptions, LineSearch, Metric, RcError,
    Session, Solver, energy, read_curve_obj, write_curve_obj,
};

struct Args {
    knot: String,
    method: String,
    threads: usize,
    max_iters: usize,
    max_seconds: Option<f64>,
    alpha: f64,
    beta: f64,
    hencky_w: f64,
    backend: String,
    audit: bool,
    out: Option<String>,
    final_obj: Option<String>,
    frames_every: Option<f64>,
    frames_out: Option<String>,
    tag: String,
}

fn parse(args: impl Iterator<Item = String>) -> Result<Args, String> {
    let mut a = Args {
        knot: String::new(),
        method: String::new(),
        threads: 1,
        max_iters: 1_000_000,
        max_seconds: None,
        alpha: 2.0,
        beta: 4.5,
        hencky_w: 1e3,
        backend: "auto".into(),
        audit: true,
        out: None,
        final_obj: None,
        frames_every: None,
        frames_out: None,
        tag: String::new(),
    };
    let mut it = args;
    while let Some(s) = it.next() {
        let mut value = |name: &str| it.next().ok_or(format!("{name} needs a value"));
        let num = |s: String| s.parse::<f64>().map_err(|e| format!("{s}: {e}"));
        match s.as_str() {
            "--knot" => a.knot = value("--knot")?,
            "--method" => a.method = value("--method")?,
            "--threads" => {
                a.threads = value("--threads")?
                    .parse()
                    .map_err(|e| format!("--threads: {e}"))?
            }
            "--max-iters" => {
                a.max_iters = value("--max-iters")?
                    .parse()
                    .map_err(|e| format!("--max-iters: {e}"))?
            }
            "--max-seconds" => a.max_seconds = Some(num(value("--max-seconds")?)?),
            "--alpha" => a.alpha = num(value("--alpha")?)?,
            "--beta" => a.beta = num(value("--beta")?)?,
            "--hencky-w" => a.hencky_w = num(value("--hencky-w")?)?,
            "--backend" => a.backend = value("--backend")?,
            "--no-audit" => a.audit = false,
            "--out" => a.out = Some(value("--out")?),
            "--final-obj" => a.final_obj = Some(value("--final-obj")?),
            "--frames-every-seconds" => {
                a.frames_every = Some(num(value("--frames-every-seconds")?)?)
            }
            "--frames-out" => a.frames_out = Some(value("--frames-out")?),
            "--tag" => a.tag = value("--tag")?,
            other => return Err(format!("unknown argument {other}")),
        }
    }
    if a.knot.is_empty() || a.method.is_empty() {
        return Err("--knot and --method are required".into());
    }
    if a.threads == 0 {
        return Err("--threads must be positive".into());
    }
    if a.frames_every.is_some() != a.frames_out.is_some() {
        return Err("--frames-every-seconds and --frames-out go together".into());
    }
    if a.frames_every.is_some_and(|s| !(s.is_finite() && s > 0.0)) {
        return Err("--frames-every-seconds must be positive".into());
    }
    Ok(a)
}

struct MethodSpec {
    soft: bool,
    metric: Metric,
    method: Method,
    armijo: bool,
}

fn method_spec(name: &str) -> Option<MethodSpec> {
    let sob = Metric::default();
    let (soft, metric, method, armijo) = match name {
        "hs-pg" => (false, sob, Method::ProjectedGradient, false),
        "l2-pg" => (false, Metric::L2, Method::ProjectedGradient, false),
        "h1-pg" => (false, Metric::H1, Method::ProjectedGradient, false),
        "h2-pg" => (false, Metric::H2, Method::ProjectedGradient, false),
        "hs-armijo" => (false, sob, Method::ProjectedGradient, true),
        "soft-hs-ncg" => (true, sob, Method::Ncg, false),
        "soft-hs-lbfgs" => (true, sob, Method::Lbfgs { memory: 30 }, false),
        "soft-h1-ncg" => (true, Metric::H1, Method::Ncg, false),
        "soft-h1-lbfgs" => (true, Metric::H1, Method::Lbfgs { memory: 30 }, false),
        _ => return None,
    };
    Some(MethodSpec {
        soft,
        metric,
        method,
        armijo,
    })
}

fn json_str(s: &str) -> String {
    let mut o = String::with_capacity(s.len() + 2);
    o.push('"');
    for c in s.chars() {
        match c {
            '"' => o.push_str("\\\""),
            '\\' => o.push_str("\\\\"),
            '\n' => o.push_str("\\n"),
            c if (c as u32) < 0x20 => o.push(' '),
            c => o.push(c),
        }
    }
    o.push('"');
    o
}

fn jnum(x: f64) -> String {
    if x.is_finite() {
        format!("{x:.10e}")
    } else {
        "null".into()
    }
}

struct Frames {
    file: std::io::BufWriter<std::fs::File>,
    every: f64,
    next: f64,
}

impl Frames {
    fn create(
        path: &str,
        every: f64,
        net: &CurveNetwork,
        n: usize,
    ) -> Result<Frames, std::io::Error> {
        let mut file = std::io::BufWriter::new(std::fs::File::create(path)?);
        file.write_all(b"RCF1")?;
        file.write_all(&(n as u32).to_le_bytes())?;
        file.write_all(&(net.edges().len() as u32).to_le_bytes())?;
        file.write_all(&(every as f32).to_le_bytes())?;
        file.write_all(&0u32.to_le_bytes())?;
        for e in net.edges() {
            file.write_all(&e[0].to_le_bytes())?;
            file.write_all(&e[1].to_le_bytes())?;
        }
        Ok(Frames {
            file,
            every,
            next: every,
        })
    }

    fn record(
        &mut self,
        t: f64,
        iteration: usize,
        flag: u32,
        pos: &[[f64; 3]],
    ) -> Result<(), std::io::Error> {
        self.file.write_all(&(t as f32).to_le_bytes())?;
        self.file.write_all(&(iteration as u32).to_le_bytes())?;
        self.file.write_all(&flag.to_le_bytes())?;
        for p in pos {
            for c in p {
                self.file.write_all(&(*c as f32).to_le_bytes())?;
            }
        }
        Ok(())
    }

    fn maybe(&mut self, t: f64, iteration: usize, pos: &[[f64; 3]]) -> Result<(), std::io::Error> {
        if t >= self.next {
            self.record(t, iteration, 0, pos)?;
            self.next = t + self.every;
        }
        Ok(())
    }
}

fn run(a: &Args) -> Result<(), String> {
    let spec = method_spec(&a.method).ok_or(format!("unknown method {}", a.method))?;
    let (net, pos0) = read_curve_obj(&a.knot).map_err(|e| format!("{}: {e}", a.knot))?;
    let n = pos0.len();
    let exps = Exponents::new(a.alpha, a.beta).map_err(|e| e.to_string())?;
    let backend = match a.backend.as_str() {
        "auto" => Backend::default(),
        "exact" => Backend::Explicit {
            energy: EnergyEval::Exact,
            solver: Solver::Dense,
        },
        "bh-mg" => Backend::Explicit {
            energy: EnergyEval::BarnesHut(Default::default()),
            solver: Solver::Multigrid(Default::default()),
        },
        other => return Err(format!("unknown backend {other}")),
    };
    let mut opts = FlowOptions {
        exponents: exps,
        metric: spec.metric,
        backend,
        method: spec.method,
        line_search: if spec.armijo {
            LineSearch::normalized_armijo()
        } else {
            LineSearch::default()
        },
        projection: ProjectionOptions {
            tolerance: Tolerance::Absolute(1e-4),
            ..ProjectionOptions::default()
        },
        max_iterations: a.max_iters,
        record_path: a.audit,
        threads: Some(a.threads),
        ..FlowOptions::default()
    };
    let hencky = if spec.soft {
        let h = Hencky::at_rest(a.hencky_w, &net, &pos0).map_err(|e| e.to_string())?;
        opts.constraints = Vec::new();
        opts.penalties = vec![std::sync::Arc::new(h.clone())];
        Some(h)
    } else {
        None
    };
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(a.threads)
        .build()
        .map_err(|e| e.to_string())?;
    let audit_opts = CcdOptions {
        gap: 1e-6,
        max_advances: 100_000,
    };
    let exact_energy = |net: &CurveNetwork, p: &[[f64; 3]]| -> Result<f64, RcError> {
        Ok(energy(net, p, exps, &EnergyEval::Exact, Some(a.threads))?.energy)
    };
    let e0_exact = exact_energy(&net, &pos0).map_err(|e| e.to_string())?;
    let rest: Vec<f64> = net
        .edges()
        .iter()
        .map(|&[i, j]| {
            let (p, q) = (pos0[i as usize], pos0[j as usize]);
            ((p[0] - q[0]).powi(2) + (p[1] - q[1]).powi(2) + (p[2] - q[2]).powi(2)).sqrt()
        })
        .collect();

    let mut frames = match (&a.frames_out, a.frames_every) {
        (Some(p), Some(s)) => Some(Frames::create(p, s, &net, n).map_err(|e| format!("{p}: {e}"))?),
        _ => None,
    };
    let mut session = Session::new(net.clone(), pos0.clone(), opts).map_err(|e| e.to_string())?;
    let f0 = session.energy();

    let (mut history, mut audit_s) = (Vec::<[f64; 5]>::new(), 0.0f64);
    let (mut audited, mut violations, mut first_violation) = (0usize, 0usize, None);
    let (mut audit_errors, mut stuck_steps, mut trials) = (0usize, 0usize, 0usize);
    let mut last_violation = 0.0;
    let mut error: Option<String> = None;
    let (mut eb, mut sb) = ("", "");
    let mut end_reason = "";
    let t0 = Instant::now();
    let mut frame_s = 0.0f64;
    let wall = |audit_s: f64, frame_s: f64| t0.elapsed().as_secs_f64() - audit_s - frame_s;
    if let Some(f) = frames.as_mut() {
        let tf = Instant::now();
        f.record(0.0, 0, 0, session.positions())
            .map_err(|e| e.to_string())?;
        frame_s += tf.elapsed().as_secs_f64();
    }
    history.push([0.0, 0.0, f0, f0, f64::NAN]);
    loop {
        if let Some(limit) = a.max_seconds
            && wall(audit_s, frame_s) >= limit
        {
            end_reason = "time";
            break;
        }
        let r = match session.step() {
            Ok(r) => r,
            Err(e) => {
                error = Some(e.to_string());
                end_reason = "error";
                break;
            }
        };
        let t = wall(audit_s, frame_s);
        trials += r.trials.len();
        (eb, sb) = (r.energy_backend, r.solver_backend);
        last_violation = r.max_violation;
        if r.status == StepStatus::Stuck {
            stuck_steps += 1;
        }
        let pen: f64 = r.penalties_after.iter().sum();
        history.push([
            r.iteration as f64,
            t,
            r.energy_after,
            r.energy_after - pen,
            r.gradient_norm,
        ]);
        eprintln!(
            "TRAJ {} {:.3} {:.8} {:.3e}",
            r.iteration, t, r.energy_after, r.gradient_norm
        );
        if let Some(p) = &r.path {
            let ta = Instant::now();
            let mut rep = CcdReport::default();
            match path_contact(&pool, session.network(), p, &[], &audit_opts, &mut rep) {
                Ok(Some((piece, c))) => {
                    violations += 1;
                    first_violation.get_or_insert((
                        r.iteration,
                        piece,
                        format!("{:?}", c.primitive),
                    ));
                }
                Ok(None) => {}
                Err(_) => audit_errors += 1,
            }
            audited += 1;
            audit_s += ta.elapsed().as_secs_f64();
        }
        if let Some(f) = frames.as_mut() {
            let tf = Instant::now();
            f.maybe(t, r.iteration, session.positions())
                .map_err(|e| e.to_string())?;
            frame_s += tf.elapsed().as_secs_f64();
        }
        if session.status().is_terminal() {
            break;
        }
    }
    let seconds = wall(audit_s, frame_s);
    let status_dbg = format!("{:?}", session.status());
    let (status, reason) = match session.status() {
        Status::Converged => ("converged", "converged"),
        Status::Stuck => ("stuck", "stuck"),
        Status::MaxIterations => ("nonconvergent", "max_iters"),
        Status::TimeLimit => ("nonconvergent", "time"),
        Status::Error(e) => {
            error.get_or_insert_with(|| e.to_string());
            ("error", "error")
        }
        Status::Cancelled => ("error", "cancelled"),
        Status::Running => {
            if error.is_some() {
                ("error", "error")
            } else {
                (
                    "nonconvergent",
                    if end_reason.is_empty() {
                        "time"
                    } else {
                        end_reason
                    },
                )
            }
        }
    };
    let pos = session.positions().to_vec();
    if let Some(f) = frames.as_mut() {
        f.record(seconds, session.iteration(), 1, &pos)
            .map_err(|e| e.to_string())?;
        f.file.flush().map_err(|e| e.to_string())?;
    }
    if let Some(p) = &a.final_obj {
        write_curve_obj(p, session.network(), &pos).map_err(|e| e.to_string())?;
    }
    let ef_exact = exact_energy(session.network(), &pos).map_err(|e| e.to_string())?;
    let (mut len_err, mut log_strain) = (0.0f64, 0.0f64);
    for (&[i, j], r) in session.network().edges().iter().zip(&rest) {
        let (p, q) = (pos[i as usize], pos[j as usize]);
        let l = ((p[0] - q[0]).powi(2) + (p[1] - q[1]).powi(2) + (p[2] - q[2]).powi(2)).sqrt();
        len_err = len_err.max((l / r - 1.0).abs());
        log_strain = log_strain.max((l / r).ln().abs());
    }
    let log_strain_lib = hencky
        .as_ref()
        .map(|h| h.max_log_strain(session.network(), &pos));
    let (f_last, e_last) = history.last().map(|h| (h[2], h[3])).unwrap_or((f0, f0));
    let knot_name = Path::new(&a.knot)
        .file_stem()
        .map(|s| s.to_string_lossy().into_owned())
        .unwrap_or_default();
    let backends = match &session.options().backend {
        Backend::Auto { .. } => "auto".to_string(),
        other => format!("{other:?}").chars().take(60).collect(),
    };
    let mut o = String::new();
    o.push_str(&format!(
        "{{\"knot\":{},\"method\":{},\"tag\":{},\"status\":{},\"end_reason\":{},\"status_debug\":{},\
\"error\":{},\"iterations\":{},\"seconds\":{},\"audit_seconds\":{},\"threads\":{},\"backend_option\":{},\
\"vertices\":{},\"alpha\":{},\"beta\":{},\"hencky_w\":{},\"max_iters\":{},\"max_seconds\":{},\
\"energy0_exact\":{},\"final_energy_exact\":{},\"final_objective\":{},\"final_ehat\":{},\
\"audit\":{{\"enabled\":{},\"audited_steps\":{},\"violations\":{},\"first_violation\":{},\"errors\":{}}},\
\"constraint\":{{\"max_len_rel_err\":{},\"max_abs_log_strain\":{},\"last_step_max_violation\":{}}},\
\"stuck_steps\":{},\"trials\":{},\"energy_backend\":{},\"solver_backend\":{},",
        json_str(&knot_name),
        json_str(&a.method),
        json_str(&a.tag),
        json_str(status),
        json_str(reason),
        json_str(&status_dbg),
        error.as_deref().map_or("null".into(), json_str),
        session.iteration(),
        jnum(seconds),
        jnum(audit_s),
        a.threads,
        json_str(&backends),
        n,
        jnum(a.alpha),
        jnum(a.beta),
        if spec.soft { jnum(a.hencky_w) } else { "null".into() },
        if a.max_iters >= 1_000_000 { "null".into() } else { a.max_iters.to_string() },
        a.max_seconds.map_or("null".into(), jnum),
        jnum(e0_exact),
        jnum(ef_exact),
        jnum(f_last),
        jnum(e_last),
        a.audit,
        audited,
        violations,
        first_violation
            .as_ref()
            .map_or("null".into(), |(i, p, c)| json_str(&format!("iter {i} piece {p} {c}"))),
        audit_errors,
        jnum(len_err),
        jnum(log_strain_lib.unwrap_or(log_strain)),
        jnum(last_violation),
        stuck_steps,
        trials,
        json_str(eb),
        json_str(sb),
    ));
    o.push_str("\"history\":[");
    for (k, h) in history.iter().enumerate() {
        if k > 0 {
            o.push(',');
        }
        o.push_str(&format!(
            "[{},{},{},{},{}]",
            h[0] as usize,
            jnum(h[1]),
            jnum(h[2]),
            jnum(h[3]),
            jnum(h[4])
        ));
    }
    o.push_str("]}");
    match &a.out {
        Some(p) => {
            let mut f = std::fs::OpenOptions::new()
                .create(true)
                .append(true)
                .open(p)
                .map_err(|e| format!("{p}: {e}"))?;
            writeln!(f, "{o}").map_err(|e| e.to_string())?;
        }
        None => println!("{o}"),
    }
    eprintln!(
        "RESULT knot {knot_name} method {} status {status} ({reason}) iterations {} seconds {seconds:.1} \
         E0 {e0_exact:.6} E {ef_exact:.6} violations {violations}",
        a.method,
        session.iteration()
    );
    Ok(())
}

fn main() {
    let a = match parse(std::env::args().skip(1)) {
        Ok(a) => a,
        Err(e) => {
            eprintln!("eval_knots: {e}");
            std::process::exit(2);
        }
    };
    if let Err(e) = run(&a) {
        eprintln!("eval_knots: {e}");
        std::process::exit(1);
    }
}
