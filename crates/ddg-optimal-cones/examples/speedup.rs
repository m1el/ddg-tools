//! The comparison table of the factor-once method against the reference method.
//!
//! ```text
//! cargo run --release -p ddg-optimal-cones --example speedup -- [options] [mesh names...]
//!
//!   --lambda X[,Y,...]  the parameters λ (default 0.5); several give a sweep
//!   --data DIR          dataset directory (default data/optimal-cones/optimal-cones-data/data)
//!   --threads N         size of the thread pool for the factor-once and automatic runs
//!                       (default 1: sequential); the reference method always runs sequentially
//!                       (threads do not help it)
//!   --repeat N          solve N times and keep the fastest (default 3)
//!   --ref-repeat N      the same for the reference method (default: --repeat)
//!   --skip-reference    do not run the reference method (no agreement columns)
//!   --old-schedule      the reference method with the old default schedule 1, 1e-1, ..., 1e-10
//! ```
//!
//! For every mesh and λ one line: the reference method (`SolveMethod::Reference`), the
//! factor-once method forced with its limits lifted (`SolveMethod::FactorOnce`), and the
//! automatic mode (`SolveMethod::Auto`, default options) with the method that produced its
//! result; times are `report.timings.total` in seconds. The last columns compare the automatic
//! result with the reference result (largest difference of the cone angles and of the scale
//! factor, relative difference of the objective) and tell whether the certificate holds.

use std::path::PathBuf;

use ddg_core::io::read_obj;
use ddg_core::{Geometry, TriMesh};
use ddg_optimal_cones::{
    ConeProblem, MethodUsed, Parallelism, Solution, SolveMethod, SolveOptions, ThreadPool, solve,
};

fn best_of(
    n: usize,
    problem: &ConeProblem<'_>,
    options: &SolveOptions<'_>,
) -> Solution {
    let mut best: Option<Solution> = None;
    for _ in 0..n.max(1) {
        let s = solve(problem, options).unwrap_or_else(|e| panic!("{e}"));
        if best
            .as_ref()
            .is_none_or(|b| s.report.timings.total < b.report.timings.total)
        {
            best = Some(s);
        }
    }
    best.expect("at least one run")
}

fn max_diff(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b).map(|(x, y)| (x - y).abs()).fold(0.0, f64::max)
}

fn main() {
    let mut lambdas = vec![0.5];
    let mut dir = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../data/optimal-cones/optimal-cones-data/data");
    let mut names: Vec<String> = Vec::new();
    let mut threads = 1usize;
    let mut repeat = 3usize;
    let mut ref_repeat = None;
    let mut skip_reference = false;
    let mut old_schedule = false;
    let mut args = std::env::args().skip(1);
    while let Some(arg) = args.next() {
        let mut value = |what: &str| args.next().unwrap_or_else(|| panic!("{what} needs a value"));
        match arg.as_str() {
            "--lambda" => {
                lambdas = value("--lambda")
                    .split(',')
                    .map(|s| s.parse().expect("--lambda: numbers"))
                    .collect()
            }
            "--data" => dir = PathBuf::from(value("--data")),
            "--threads" => threads = value("--threads").parse().expect("--threads: an integer"),
            "--repeat" => repeat = value("--repeat").parse().expect("--repeat: an integer"),
            "--ref-repeat" => {
                ref_repeat = Some(value("--ref-repeat").parse().expect("--ref-repeat: an integer"))
            }
            "--skip-reference" => skip_reference = true,
            "--old-schedule" => old_schedule = true,
            name => names.push(name.to_string()),
        }
    }
    let ref_repeat = ref_repeat.unwrap_or(repeat);
    if !dir.is_dir() {
        eprintln!("dataset directory {} not found", dir.display());
        std::process::exit(2);
    }
    if names.is_empty() {
        names = std::fs::read_dir(&dir)
            .expect("readable dataset directory")
            .filter_map(|e| e.ok())
            .filter_map(|e| {
                e.file_name()
                    .to_str()
                    .and_then(|n| n.strip_suffix(".obj"))
                    .map(str::to_string)
            })
            .collect();
        names.sort();
    }
    let pool = (threads > 1).then(|| ThreadPool::new(threads).expect("thread pool"));
    let parallelism = pool.as_ref().map_or(Parallelism::Sequential, Parallelism::Pool);

    println!(
        "threads = {threads}, best of {repeat} (reference: best of {ref_repeat}, sequential), times in seconds"
    );
    println!(
        "{:<17} {:>6} {:>6} {:>6} | {:>7} {:>5} | {:>7} {:>5} {:>5} | {:>7} {:>9} {:>5} {:>7} | {:>6} | {:>9} {:>9} {:>9} {:>4}",
        "mesh", "V", "lambda", "cones", "ref", "fact", "fast", "ws", "round", "auto", "used", "ws", "factor", "x", "max dmu", "max du", "obj rel", "cert"
    );
    for name in &names {
        let obj = read_obj(dir.join(format!("{name}.obj"))).unwrap_or_else(|e| panic!("{e}"));
        let faces = obj.triangles().unwrap_or_else(|| panic!("{name}: not a triangle mesh"));
        let mesh = TriMesh::new(obj.positions.len(), faces).unwrap_or_else(|e| panic!("{name}: {e}"));
        let geometry =
            Geometry::from_positions(&mesh, &obj.positions).unwrap_or_else(|e| panic!("{name}: {e}"));
        for &lambda in &lambdas {
            let problem = ConeProblem::new(&mesh, &geometry, lambda);
            let mut reference_options = SolveOptions {
                method: SolveMethod::Reference,
                ..SolveOptions::default()
            };
            if old_schedule {
                reference_options.regularization = (0..=10).map(|k| 10f64.powi(-k)).collect();
            }
            let reference = (!skip_reference).then(|| best_of(ref_repeat, &problem, &reference_options));
            let mut fast_options = SolveOptions {
                method: SolveMethod::FactorOnce,
                parallelism,
                ..SolveOptions::default()
            };
            fast_options.factor_once.max_working_set_fraction = 1e9;
            fast_options.factor_once.min_working_set = usize::MAX / 2;
            fast_options.factor_once.memory_limit_bytes = usize::MAX;
            let fast = best_of(repeat, &problem, &fast_options);
            let auto = best_of(
                repeat,
                &problem,
                &SolveOptions {
                    parallelism,
                    ..SolveOptions::default()
                },
            );
            let fast_report = fast.report.method.factor_once.as_ref();
            let auto_report = auto.report.method.factor_once.as_ref();
            let used = match (auto.report.method.used, &auto.report.method.fallback) {
                (MethodUsed::FactorOnce, _) => "fast".to_string(),
                (MethodUsed::Reference, None) => "reference".to_string(),
                (MethodUsed::Reference, Some(_)) => "ref(fb)".to_string(),
            };
            let (rt, rf, dmu, du, obj_rel) = match &reference {
                Some(r) => {
                    let l = &r.report.linear_solver;
                    (
                        format!("{:.3}", r.report.timings.total),
                        format!("{}", l.ldlt_factorizations + l.lu_factorizations),
                        max_diff(&auto.measure, &r.measure),
                        max_diff(&auto.log_scale_optimal, &r.log_scale_optimal),
                        (auto.objective - r.objective).abs() / r.objective.abs().max(1e-300),
                    )
                }
                None => ("-".to_string(), "-".to_string(), f64::NAN, f64::NAN, f64::NAN),
            };
            let ratio = reference
                .as_ref()
                .map_or(f64::NAN, |r| r.report.timings.total / auto.report.timings.total);
            println!(
                "{:<17} {:>6} {:>6} {:>6} | {:>7} {:>5} | {:>7.3} {:>5} {:>5} | {:>7.3} {:>9} {:>5} {:>7.3} | {:>5.1}x | {:>9.1e} {:>9.1e} {:>9.1e} {:>4}",
                name,
                mesh.n_vertices(),
                lambda,
                auto.measure.iter().filter(|m| **m != 0.0).count(),
                rt,
                rf,
                fast.report.timings.total,
                fast_report.map_or(0, |r| r.working_set),
                fast_report.map_or(0, |r| r.rounds),
                auto.report.timings.total,
                used,
                auto_report.map_or(0, |r| r.working_set),
                auto.report.timings.factor_once,
                ratio,
                dmu,
                du,
                obj_rel,
                if auto.report.certificate.holds(1e-9) { "ok" } else { "FAIL" }
            );
        }
    }
}
