//! Baseline performance of the reference solver on the optimal-cones dataset.
//!
//! ```text
//! cargo run --release -p ddg-optimal-cones --example baseline -- [options] [mesh names...]
//!
//!   --lambda X          the parameter λ (default 0.5)
//!   --data DIR          dataset directory (default data/optimal-cones/optimal-cones-data/data)
//!   --metric mass|unweighted
//!   --gammas a,b,c      regularization schedule (default: the solver's default)
//!   --threads N         size of the thread pool (default 1: sequential)
//!   --lu                force LU for the Newton systems
//!   --merge             merge cone clusters (without a limit on the cluster size)
//!   --repeat N          solve N times and report the fastest (default 1)
//!   --stages            print one line per stage
//!   --method M          reference (default: this is the baseline of the reference method), auto, fast
//! ```
//!
//! Without mesh names every `*.obj` of the directory is run. One line per mesh is printed:
//! size, status, wall times (total, numeric factorization, triangular solves, symbolic
//! analysis, set-up, rest), Newton iterations of Stage I and Stage II, safeguard rounds,
//! factorizations, solves, cones before and after merging, and the certificate.

use std::path::PathBuf;

use ddg_core::io::read_obj;
use ddg_core::{Geometry, TriMesh};
use ddg_optimal_cones::{
    ConeProblem, LinearSolver, Parallelism, RegularizationMetric, SolveMethod, SolveOptions,
    ThreadPool, solve,
};

fn main() {
    let mut lambda = 0.5;
    let mut dir = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../data/optimal-cones/optimal-cones-data/data");
    let mut names: Vec<String> = Vec::new();
    let mut options = SolveOptions {
        method: SolveMethod::Reference,
        ..SolveOptions::default()
    };
    let mut threads = 1usize;
    let mut print_stages = false;
    let mut repeat = 1usize;
    let mut args = std::env::args().skip(1);
    while let Some(arg) = args.next() {
        let mut value = |what: &str| {
            args.next()
                .unwrap_or_else(|| panic!("{what} needs a value"))
        };
        match arg.as_str() {
            "--lambda" => lambda = value("--lambda").parse().expect("--lambda: a number"),
            "--data" => dir = PathBuf::from(value("--data")),
            "--metric" => {
                options.regularization_metric = match value("--metric").as_str() {
                    "mass" => RegularizationMetric::Mass,
                    "unweighted" => RegularizationMetric::Unweighted,
                    other => panic!("--metric: expected mass or unweighted, got {other}"),
                }
            }
            "--gammas" => {
                options.regularization = value("--gammas")
                    .split(',')
                    .filter(|s| !s.is_empty())
                    .map(|s| s.parse().expect("--gammas: comma separated numbers"))
                    .collect()
            }
            "--threads" => threads = value("--threads").parse().expect("--threads: an integer"),
            "--lu" => options.linear_solver = LinearSolver::SaddleLu,
            "--merge" => {
                options.merge_cone_clusters = true;
                options.max_merge_cluster_size = usize::MAX;
            }
            "--repeat" => repeat = value("--repeat").parse().expect("--repeat: an integer"),
            "--stages" => print_stages = true,
            "--method" => {
                options.method = match value("--method").as_str() {
                    "reference" => SolveMethod::Reference,
                    "auto" => SolveMethod::Auto,
                    "fast" => SolveMethod::FactorOnce,
                    other => panic!("--method: expected reference, auto or fast, got {other}"),
                }
            }
            name => names.push(name.to_string()),
        }
    }
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
    if let Some(pool) = &pool {
        options.parallelism = Parallelism::Pool(pool);
    }

    println!(
        "lambda = {lambda}, metric = {:?}, gammas = {:?}, threads = {threads}, best of {repeat}, \
         merging {}",
        options.regularization_metric,
        options.regularization,
        if options.merge_cone_clusters {
            "on (no size limit)"
        } else {
            "off"
        }
    );
    println!(
        "{:<19} {:>6} {:>6} {:>9} | {:>7} {:>7} {:>7} {:>7} {:>7} {:>6} | {:>4} {:>4} {:>2} {:>4} {:>6} {:>4} | {:>5} {:>4} {:>6} {:>6} {:>10} {:>10} | {:>7} {:>7}  Newton steps per regularized problem",
        "mesh",
        "V",
        "F",
        "status",
        "total",
        "factor",
        "solves",
        "symbol",
        "setup",
        "rest",
        "itI",
        "itII",
        "sg",
        "fact",
        "solves",
        "refi",
        "cones",
        "clus",
        "merge8",
        "merged",
        "E",
        "E merged",
        "state",
        "gap"
    );
    for name in &names {
        let obj = read_obj(dir.join(format!("{name}.obj"))).unwrap_or_else(|e| panic!("{e}"));
        let faces = obj
            .triangles()
            .unwrap_or_else(|| panic!("{name}: not a triangle mesh"));
        let mesh =
            TriMesh::new(obj.positions.len(), faces).unwrap_or_else(|e| panic!("{name}: {e}"));
        let geometry = Geometry::from_positions(&mesh, &obj.positions)
            .unwrap_or_else(|e| panic!("{name}: {e}"));
        let problem = ConeProblem::new(&mesh, &geometry, lambda);
        // The fastest of the runs; everything but the times is the same in each.
        let mut best: Option<ddg_optimal_cones::Solution> = None;
        let mut failed = false;
        for _ in 0..repeat {
            match solve(&problem, &options) {
                Ok(s) => {
                    if best
                        .as_ref()
                        .is_none_or(|b| s.report.timings.total < b.report.timings.total)
                    {
                        best = Some(s);
                    }
                }
                Err(e) => {
                    println!("{name:<19} ERROR {e}");
                    failed = true;
                    break;
                }
            }
        }
        if failed {
            continue;
        }
        let solution = best.expect("at least one run");
        let r = &solution.report;
        let c = &r.certificate;
        let l = &r.linear_solver;
        let (factor, solves, symbolic) = (
            l.factorization_time.as_secs_f64(),
            l.solve_time.as_secs_f64(),
            l.symbolic_time.as_secs_f64(),
        );
        let raw = solution.measure.iter().filter(|m| **m != 0.0).count();
        // What the default limit of 8 cones per cluster would have merged.
        let merged_small = raw
            - solution
                .extraction
                .merged
                .iter()
                .filter(|cluster| cluster.members.len() <= 8)
                .map(|cluster| cluster.members.len() - 1)
                .sum::<usize>();
        let per_stage: Vec<String> = r
            .stages
            .iter()
            .filter(|s| !s.safeguard && s.kind != ddg_optimal_cones::StageKind::Exact)
            .map(|s| s.newton_iterations.to_string())
            .collect();
        println!(
            "{:<19} {:>6} {:>6} {:>9} | {:>7.3} {:>7.3} {:>7.3} {:>7.3} {:>7.3} {:>6.3} | {:>4} {:>4} {:>2} {:>4} {:>6} {:>4} | {:>5} {:>4} {:>6} {:>6} {:>10.4e} {:>10.4e} | {:>7.0e} {:>7.0e}  {}",
            name,
            mesh.n_vertices(),
            mesh.n_faces(),
            format!("{:?}", r.status),
            r.timings.total,
            factor,
            solves,
            symbolic,
            r.timings.setup,
            r.timings.total - factor - solves - symbolic - r.timings.setup,
            r.regularized_iterations,
            r.exact_iterations,
            r.safeguard_rounds,
            l.ldlt_factorizations + l.lu_factorizations,
            l.solves,
            l.refinement_steps,
            raw,
            solution.extraction.n_clusters,
            if options.merge_cone_clusters {
                merged_small.to_string()
            } else {
                "-".to_string()
            },
            if options.merge_cone_clusters {
                solution.cones.len().to_string()
            } else {
                "-".to_string()
            },
            solution.distortion_optimal,
            solution.distortion,
            c.state_residual,
            c.duality_gap,
            per_stage.join(",")
        );
        if print_stages {
            for s in &r.stages {
                println!(
                    "    {:<40} newton {:>3}  damped {:>3}  shifts {:>2}  {:?}  {:.3}s",
                    format!("{:?}", s.kind),
                    s.newton_iterations,
                    s.damped_steps,
                    s.constant_shifts,
                    s.stop,
                    s.seconds
                );
            }
            if let Some(m) = &r.message {
                println!("    message: {m}");
            }
            println!(
                "    timings: {:?}\n    certificate: {c:?}\n    extraction: clusters {} dropped {}",
                r.timings, solution.extraction.n_clusters, solution.extraction.n_dropped
            );
        }
    }
    println!(
        "\ntimes in seconds, wall clock: total = everything inside solve(); factor = numeric \
         factorizations of the Newton matrix; solves = triangular solves and refinement;\nsymbol = \
         symbolic analysis and ordering; setup = validation and assembly of the operators; rest = \
         Newton bookkeeping, line searches, certificate, extraction.\nitI / itII: Newton steps of \
         Stage I (all regularized problems) and Stage II; sg: safeguard rounds; fact: numeric \
         factorizations; solves: triangular solve pairs,\nof which refi are refinement steps. \
         cones: cones of the optimal measure; clus: clusters (edge-connected cones of one sign); \
         merge8 / merged: cones after merging the clusters of at most 8 cones / all clusters; \
         E / E merged: distortion of the optimum / of the cones as returned."
    );
}
