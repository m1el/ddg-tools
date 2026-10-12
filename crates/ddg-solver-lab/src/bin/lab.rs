use ddg_solver_lab::{lasso, meshes, prepare::Prepared};
use ddg_optimal_cones::{ConeProblem, SolveOptions, solve};
use std::time::Instant;

fn main() {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let cmd = args.first().map(String::as_str).unwrap_or("");
    match cmd {
        "compare" => compare(&args[1..]),
        _ => eprintln!("usage: lab compare MESH [--lambda X] [--threads N] [--no-ref] ..."),
    }
}

fn flag(args: &[String], name: &str) -> Option<String> {
    args.iter().position(|a| a == name).map(|i| args[i + 1].clone())
}

fn compare(args: &[String]) {
    let spec = &args[0];
    let lambda: f64 = flag(args, "--lambda").map(|v| v.parse().unwrap()).unwrap_or(0.5);
    let mut o = lasso::Options::default();
    if let Some(v) = flag(args, "--threads") { o.threads = v.parse().unwrap(); }
    if let Some(v) = flag(args, "--chunk") { o.chunk = v.parse().unwrap(); }
    if let Some(v) = flag(args, "--steps") { o.cont_steps = v.parse().unwrap(); }
    if let Some(v) = flag(args, "--factor") { o.cont_factor = v.parse().unwrap(); }
    if let Some(v) = flag(args, "--batch-min") { o.batch_min = v.parse().unwrap(); }
    if let Some(v) = flag(args, "--batch-frac") { o.batch_frac = v.parse().unwrap(); }
    if let Some(v) = flag(args, "--nms") { o.nms_radius = v.parse().unwrap(); }
    if let Some(v) = flag(args, "--tol") { o.tol = v.parse().unwrap(); }
    let no_ref = args.iter().any(|a| a == "--no-ref");
    let case = meshes::load(spec).unwrap_or_else(|e| panic!("{e}"));
    let prep = Prepared::new(&case).unwrap();
    println!("{spec}: V={} F={} closed={} setup {:.3}s (factor {:.3}s)", case.mesh.n_vertices(), case.mesh.n_faces(), prep.closed, prep.setup_seconds, prep.factor_seconds);
    let t = Instant::now();
    let s = lasso::solve(&prep, lambda, &o);
    let wall = t.elapsed().as_secs_f64();
    let st = &s.stats;
    let nz = s.mu.iter().filter(|m| m.abs() > 1e-12).count();
    println!("lab: wall {wall:.3}s (+setup {:.3}) converged={} cones={nz} peakA={} W={} rounds={} stages={} solves: cols {} grad {} total {} excess {:.1e}",
        prep.setup_seconds, st.converged, st.peak_active, st.working_set, st.rounds, st.stages, st.column_solves, st.gradient_solves, st.total_solves, st.max_excess);
    println!("  times: columns {:.3} gram {:.3} dense {:.3} gradient {:.3}; dense iters {} adds {} drops {} failed {}", st.t_columns, st.t_gram, st.t_dense, st.t_gradient, st.dense.iterations, st.dense.adds, st.dense.drops, st.dense.failed_adds);
    let problem = ConeProblem::new(&case.mesh, &case.geometry, lambda);
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
