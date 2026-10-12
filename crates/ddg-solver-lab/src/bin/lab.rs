use ddg_solver_lab::{lasso, meshes, prepare::{Constraints, Prepared}};
use ddg_optimal_cones::{ConeProblem, SolveOptions, solve};
use std::time::Instant;

fn main() {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let cmd = args.first().map(String::as_str).unwrap_or("");
    match cmd {
        "compare" => compare(&args[1..]),
        "table" => table(&args[1..]),
        "scale" => scale(&args[1..]),
        "refopts" => refopts(&args[1..]),
        _ => eprintln!("usage: lab compare MESH [--lambda X] [--threads N] [--no-ref] ..."),
    }
}

fn flag(args: &[String], name: &str) -> Option<String> {
    args.iter().position(|a| a == name).map(|i| args[i + 1].clone())
}

/// Steering from flags: `--nonneg`, `--free-every K` (cost 0 at every K-th vertex),
/// `--sites K` (K free sites spread over the mesh, every other vertex forbidden), `--forbid-first F`
/// (the first fraction F of the vertices forbidden).
fn constraints_from(args: &[String], nv: usize) -> Constraints {
    let mut c = Constraints::default();
    if args.iter().any(|a| a == "--nonneg") { c.nonneg = Some(vec![true; nv]); }
    if let Some(k) = flag(args, "--free-every") {
        let k: usize = k.parse().unwrap();
        let mut cost = c.cost.take().unwrap_or(vec![1.0; nv]);
        for v in (0..nv).step_by(k) { cost[v] = 0.0; }
        c.cost = Some(cost);
    }
    if let Some(k) = flag(args, "--sites") {
        let k: usize = k.parse().unwrap();
        let mut cost = vec![f64::INFINITY; nv];
        for i in 0..k { cost[(i * nv) / k + nv / (2 * k)] = 0.0; }
        c.cost = Some(cost);
    }
    if let Some(f) = flag(args, "--forbid-first") {
        let f: f64 = f.parse().unwrap();
        let mut cost = c.cost.take().unwrap_or(vec![1.0; nv]);
        for v in 0..(f * nv as f64) as usize { cost[v] = f64::INFINITY; }
        c.cost = Some(cost);
    }
    c
}

fn problem_with<'a>(case: &'a meshes::Case, lambda: f64, c: &Constraints) -> ConeProblem<'a> {
    let mut p = ConeProblem::new(&case.mesh, &case.geometry, lambda);
    p.cone_costs = c.cost.clone();
    if c.nonneg.is_some() { p.lower_bounds = Some(vec![0.0; case.mesh.n_vertices()]); }
    p
}

fn opts_from(args: &[String]) -> lasso::Options {
    let mut o = lasso::Options::default();
    if let Some(v) = flag(args, "--threads") { o.threads = v.parse().unwrap(); }
    if let Some(v) = flag(args, "--chunk") { o.chunk = v.parse().unwrap(); }
    if let Some(v) = flag(args, "--steps") { o.cont_steps = v.parse().unwrap(); }
    if let Some(v) = flag(args, "--factor") { o.cont_factor = v.parse().unwrap(); }
    if let Some(v) = flag(args, "--batch-min") { o.batch_min = v.parse().unwrap(); }
    if let Some(v) = flag(args, "--batch-frac") { o.batch_frac = v.parse().unwrap(); }
    if let Some(v) = flag(args, "--nms") { o.nms_radius = v.parse().unwrap(); }
    if let Some(v) = flag(args, "--tol") { o.tol = v.parse().unwrap(); }
    if args.iter().any(|a| a == "--no-gemv") { o.gemv_u = false; }
    if let Some(v) = flag(args, "--multi-add") { o.multi_add = v.parse().unwrap(); }
    o
}

/// One line per mesh: fastest of `--repeat` runs of the lab solver and of the reference.
fn table(args: &[String]) {
    let lambda: f64 = flag(args, "--lambda").map(|v| v.parse().unwrap()).unwrap_or(0.5);
    let repeat: usize = flag(args, "--repeat").map(|v| v.parse().unwrap()).unwrap_or(3);
    let ref_repeat: usize = flag(args, "--ref-repeat").map(|v| v.parse().unwrap()).unwrap_or(repeat);
    let o = opts_from(args);
    let names: Vec<String> = args.iter().take_while(|a| !a.starts_with("--")).cloned().collect();
    println!("lambda {lambda}, lab threads {}, lab best of {repeat}, ref best of {ref_repeat} (1 thread)", o.threads);
    println!("{:<18} {:>6} {:>7} {:>5} | {:>7} {:>7} {:>5} {:>5} {:>5} {:>5} | {:>7} {:>4} {:>5} | {:>6} | {:>9} {:>9} {:>9} {:>5}",
        "mesh", "V", "F", "cones", "lab s", "setup", "cols", "grad", "peakA", "W", "ref s", "fact", "slv", "x", "max dmu", "max du", "obj rel", "cert");
    for spec in names {
        let case = meshes::load(&spec).unwrap_or_else(|e| panic!("{e}"));
        let prep = Prepared::new(&case).unwrap();
        let mut best: Option<lasso::Solution> = None;
        for _ in 0..repeat {
            let s = lasso::solve(&prep, lambda, &o);
            if best.as_ref().is_none_or(|b| s.stats.t_total < b.stats.t_total) { best = Some(s); }
        }
        let s = best.unwrap();
        let problem = ConeProblem::new(&case.mesh, &case.geometry, lambda);
        let cert = ddg_optimal_cones::certify(&problem, &s.u, &s.phi, &s.mu).map(|c| c.holds(1e-9)).unwrap_or(false);
        let mut rbest: Option<ddg_optimal_cones::Solution> = None;
        for _ in 0..ref_repeat {
            let r = solve(&problem, &SolveOptions::default()).unwrap();
            if rbest.as_ref().is_none_or(|b| r.report.timings.total < b.report.timings.total) { rbest = Some(r); }
        }
        let r = rbest.unwrap();
        let st = &s.stats;
        let dmu = s.mu.iter().zip(&r.measure).map(|(a, b)| (a - b).abs()).fold(0.0, f64::max);
        let du = s.u.iter().zip(&r.log_scale_optimal).map(|(a, b)| (a - b).abs()).fold(0.0, f64::max);
        let rl = &r.report.linear_solver;
        println!("{:<18} {:>6} {:>7} {:>5} | {:>7.3} {:>7.3} {:>5} {:>5} {:>5} {:>5} | {:>7.3} {:>4} {:>5} | {:>5.1}x | {:>9.1e} {:>9.1e} {:>9.1e} {:>5}",
            spec, case.mesh.n_vertices(), case.mesh.n_faces(), s.mu.iter().filter(|m| m.abs() > 1e-12).count(),
            st.t_total, st.t_setup, st.column_solves, st.gradient_solves, st.peak_active, st.working_set,
            r.report.timings.total, rl.ldlt_factorizations + rl.lu_factorizations, rl.solves,
            r.report.timings.total / st.t_total, dmu, du, (s.objective - r.objective).abs() / r.objective.abs(), cert);
    }
}

fn compare(args: &[String]) {
    let spec = &args[0];
    let lambda: f64 = flag(args, "--lambda").map(|v| v.parse().unwrap()).unwrap_or(0.5);
    let o = opts_from(args);
    let no_ref = args.iter().any(|a| a == "--no-ref");
    let case = meshes::load(spec).unwrap_or_else(|e| panic!("{e}"));
    let cons = constraints_from(args, case.mesh.n_vertices());
    let prep = Prepared::with_constraints(&case, &cons).unwrap();
    println!("{spec}: V={} F={} closed={} setup {:.3}s (factor {:.3}s)", case.mesh.n_vertices(), case.mesh.n_faces(), prep.closed, prep.setup_seconds, prep.factor_seconds);
    let t = Instant::now();
    let s = lasso::solve(&prep, lambda, &o);
    let wall = t.elapsed().as_secs_f64();
    let st = &s.stats;
    let nz = s.mu.iter().filter(|m| m.abs() > 1e-12).count();
    println!("lab: wall {wall:.3}s (+setup {:.3}) converged={} cones={nz} peakA={} W={} rounds={} stages={} solves: cols {} grad {} total {} excess {:.1e}",
        prep.setup_seconds, st.converged, st.peak_active, st.working_set, st.rounds, st.stages, st.column_solves, st.gradient_solves, st.total_solves, st.max_excess);
    println!("  times: columns {:.3} gram {:.3} dense {:.3} gradient {:.3}; dense iters {} adds {} drops {} failed {}", st.t_columns, st.t_gram, st.t_dense, st.t_gradient, st.dense.iterations, st.dense.adds, st.dense.drops, st.dense.failed_adds);
    let problem = problem_with(&case, lambda, &cons);
    match ddg_optimal_cones::certify(&problem, &s.u, &s.phi, &s.mu) {
        Ok(c) => println!("  certificate: state {:.1e} adjoint {:.1e} cone {:.1e} dualinf {:.1e} gap {:.1e} compl {:.1e} holds(1e-9)={}", c.state_residual, c.adjoint_residual, c.cone_residual, c.dual_infeasibility, c.duality_gap, c.complementarity_gap, c.holds(1e-9)),
        Err(e) => println!("  certify error {e}"),
    }
    println!("  objective {:.12e}", s.objective);
    if !no_ref {
        let t = Instant::now();
        let r = solve(&problem, &SolveOptions::default()).unwrap();
        let rw = t.elapsed().as_secs_f64();
        let rl = &r.report.linear_solver;
        println!("ref: wall {rw:.3}s status {:?} cones {} factorizations {} solves {} objective {:.12e}", r.report.status, r.measure.iter().filter(|m| m.abs() > 1e-12).count(), rl.ldlt_factorizations + rl.lu_factorizations, rl.solves, r.objective);
        let dmu = s.mu.iter().zip(&r.measure).map(|(a, b)| (a - b).abs()).fold(0.0, f64::max);
        let du = s.u.iter().zip(&r.log_scale_optimal).map(|(a, b)| (a - b).abs()).fold(0.0, f64::max);
        let umax = r.log_scale_optimal.iter().fold(0.0f64, |a, b| a.max(b.abs()));
        let supp_ref: Vec<bool> = r.measure.iter().map(|m| m.abs() > 1e-12).collect();
        let supp: Vec<bool> = s.mu.iter().map(|m| m.abs() > 1e-12).collect();
        let sdiff = supp.iter().zip(&supp_ref).filter(|(a, b)| a != b).count();
        println!("  diff: max|dmu| {dmu:.2e} max|du| {du:.2e} (max|u| {umax:.2e}) objective rel {:.2e} support differs at {sdiff}", (s.objective - r.objective).abs() / r.objective.abs());
        println!("  speedup (lab incl. setup / ref): {:.2}x", rw / (wall + prep.setup_seconds));
    }
}

/// Thread scaling of the pieces the factor-once method is made of.
fn scale(args: &[String]) {
    let spec = &args[0];
    let case = meshes::load(spec).unwrap_or_else(|e| panic!("{e}"));
    let prep = Prepared::new(&case).unwrap();
    let ns = prep.n_slots();
    println!("{spec}: V={} setup {:.3}s (factor {:.3}s)", case.mesh.n_vertices(), prep.setup_seconds, prep.factor_seconds);
    let reps = 3;
    // (a) multi-rhs solves on one thread
    let slots: Vec<usize> = (0..ns).step_by((ns / 256).max(1)).take(256).collect();
    println!("(a) one thread, time per right-hand side vs rhs per call (best of {reps})");
    for nrhs in [1usize, 2, 4, 8, 16, 32] {
        let mut best = f64::MAX;
        for _ in 0..reps {
            let t = Instant::now();
            let cols = &slots[..(nrhs * (32 / nrhs).max(1)).min(slots.len())];
            let _ = prep.columns(cols, None, nrhs);
            best = best.min(t.elapsed().as_secs_f64() / cols.len() as f64);
        }
        println!("  nrhs {nrhs:>2}: {:.2} ms per rhs", best * 1e3);
    }
    // (b) independent solves on a pool
    println!("(b) 256 right-hand sides, chunks of 4 per call, one call per task (best of {reps})");
    let mut t1 = 0.0;
    for threads in [1usize, 2, 4, 8, 16, 32, 64] {
        let pool = rayon::ThreadPoolBuilder::new().num_threads(threads).build().unwrap();
        let mut best = f64::MAX;
        for _ in 0..reps {
            let t = Instant::now();
            let _ = prep.columns(&slots, Some(&pool), 4);
            best = best.min(t.elapsed().as_secs_f64());
        }
        if threads == 1 { t1 = best; }
        println!("  threads {threads:>2}: {:.3} s, speedup {:.2}x, {:.2} ms per rhs", best, t1 / best, best / slots.len() as f64 * 1e3);
    }
    // (c) Gram cross product: 192 old columns, 64 new
    println!("(c) Gram cross terms, 64 new columns against 192 old (own row-chunked product; best of {reps})");
    let pool1 = rayon::ThreadPoolBuilder::new().num_threads(1).build().unwrap();
    let old = prep.columns(&slots[..192], Some(&pool1), 4);
    let new = prep.columns(&slots[192..256], Some(&pool1), 4);
    let blocks = vec![(old, 192usize)];
    let mut t1 = 0.0;
    for threads in [1usize, 2, 4, 8, 16, 32, 64] {
        let pool = rayon::ThreadPoolBuilder::new().num_threads(threads).build().unwrap();
        let mut best = f64::MAX;
        for _ in 0..reps {
            let t = Instant::now();
            let _ = prep.gram_cross(&blocks, &new, 64, Some(&pool), threads);
            best = best.min(t.elapsed().as_secs_f64());
        }
        if threads == 1 { t1 = best; }
        println!("  threads {threads:>2}: {:.4} s, speedup {:.2}x", best, t1 / best);
    }
    // (d) numeric Cholesky of the Laplacian with the library's pool
    println!("(d) ddg-core numeric Cholesky of the restricted Laplacian with a ThreadPool of n threads (best of {reps})");
    {
        use ddg_core::sparse::{Method, Parallelism, Symbolic, ThreadPool};
        let l = ddg_core::operators::cotan_laplacian(&case.mesh, &case.geometry).unwrap();
        let unk: Vec<u32> = if prep.closed { (1..case.mesh.n_vertices() as u32).collect() } else { case.mesh.interior_vertices() };
        let a = l.submatrix(&unk, &unk).unwrap();
        let sym = Symbolic::analyze(&a, Method::Llt).unwrap();
        println!("  factor values {:?}", sym.info().factor_values);
        for threads in [1usize, 2, 4, 8, 16] {
            let pool = ThreadPool::new(threads).unwrap();
            let par = if threads == 1 { Parallelism::Sequential } else { Parallelism::Pool(&pool) };
            let mut best = f64::MAX;
            for _ in 0..reps {
                let t = Instant::now();
                let _ = sym.factor(&a, par).unwrap();
                best = best.min(t.elapsed().as_secs_f64());
            }
            println!("  threads {threads:>2}: {:.3} s", best);
        }
    }
}

/// Reference solver with other regularization schedules / tolerances; compares with the default.
/// `lab refopts MESH [--lambda X] [--repeat R] SCHEDULE...` where SCHEDULE is `default`, `none`
/// or comma separated gammas, optionally `@tol` for the Newton tolerance (e.g. `1e-3,1e-10@1e-9`).
fn refopts(args: &[String]) {
    let spec = &args[0];
    let lambda: f64 = flag(args, "--lambda").map(|v| v.parse().unwrap()).unwrap_or(0.5);
    let repeat: usize = flag(args, "--repeat").map(|v| v.parse().unwrap()).unwrap_or(2);
    let case = meshes::load(spec).unwrap_or_else(|e| panic!("{e}"));
    let problem = ConeProblem::new(&case.mesh, &case.geometry, lambda);
    let scheds: Vec<&String> = {
        let mut v = Vec::new();
        let mut i = 1;
        while i < args.len() {
            if args[i].starts_with("--") { i += 2; continue; }
            v.push(&args[i]);
            i += 1;
        }
        v
    };
    let base = solve(&problem, &SolveOptions::default()).unwrap();
    println!("{spec} lambda {lambda}: {:<34} {:>7} {:>5} {:>5} {:>6} {:>9} {:>9} {:>9}", "schedule", "time s", "fact", "slv", "status", "max dmu", "obj rel", "stages");
    for sch in scheds {
        let mut o = SolveOptions::default();
        let (g, tol) = match sch.split_once('@') { Some((g, t)) => (g, Some(t)), None => (sch.as_str(), None) };
        match g {
            "default" => {}
            "none" => o.regularization = vec![],
            _ => o.regularization = g.split(',').map(|x| x.parse().unwrap()).collect(),
        }
        if let Some(t) = tol { o.newton_tolerance = t.parse().unwrap(); }
        let mut best: Option<ddg_optimal_cones::Solution> = None;
        for _ in 0..repeat {
            match solve(&problem, &o) {
                Ok(r) => if best.as_ref().is_none_or(|b| r.report.timings.total < b.report.timings.total) { best = Some(r); },
                Err(e) => { println!("  {sch}: error {e}"); break; }
            }
        }
        let Some(r) = best else { continue };
        let rl = &r.report.linear_solver;
        let dmu = r.measure.iter().zip(&base.measure).map(|(a, b)| (a - b).abs()).fold(0.0, f64::max);
        println!("{:>width$} {:<34} {:>7.3} {:>5} {:>5} {:>6} {:>9.1e} {:>9.1e} {:>9}", "", sch, r.report.timings.total, rl.ldlt_factorizations + rl.lu_factorizations, rl.solves,
            format!("{:?}", r.report.status).chars().take(6).collect::<String>(), dmu, (r.objective - base.objective).abs() / base.objective.abs(), r.report.stages.len(), width = spec.len() + 14);
    }
}
