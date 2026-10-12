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
//!   --merge             merge cone clusters
//!   --stages            print one line per stage
//! ```
//!
//! Without mesh names every `*.obj` of the directory is run. One line per mesh is printed:
//! size, status, Newton iterations of Stage I and Stage II, safeguard rounds, factorizations,
//! solves, wall times, and the certificate.

use std::path::PathBuf;
use std::time::Instant;

use ddg_core::io::read_obj;
use ddg_core::{Geometry, TriMesh};
use ddg_optimal_cones::{
    ConeProblem, LinearSolver, Parallelism, RegularizationMetric, SolveOptions, ThreadPool, solve,
};

fn main() {
    let mut lambda = 0.5;
    let mut dir = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../data/optimal-cones/optimal-cones-data/data");
    let mut names: Vec<String> = Vec::new();
    let mut options = SolveOptions::default();
    let mut threads = 1usize;
    let mut print_stages = false;
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
            "--merge" => options.merge_cone_clusters = true,
            "--stages" => print_stages = true,
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
        "lambda = {lambda}, metric = {:?}, gammas = {:?}, threads = {threads}",
        options.regularization_metric, options.regularization
    );
    println!(
        "{:<20} {:>7} {:>7} {:>13} {:>4} {:>4} {:>3} {:>5} {:>3} {:>6} {:>7} {:>7} {:>7} {:>6} {:>9} {:>10} {:>8} {:>8} {:>8}",
        "mesh",
        "V",
        "F",
        "status",
        "itI",
        "itII",
        "sg",
        "ldlt",
        "lu",
        "solves",
        "load[s]",
        "solve[s]",
        "fact[s]",
        "cones",
        "angle/pi",
        "distortion",
        "state",
        "cone",
        "gap"
    );
    for name in &names {
        let start = Instant::now();
        let obj = read_obj(dir.join(format!("{name}.obj"))).unwrap_or_else(|e| panic!("{e}"));
        let faces = obj
            .triangles()
            .unwrap_or_else(|| panic!("{name}: not a triangle mesh"));
        let mesh =
            TriMesh::new(obj.positions.len(), faces).unwrap_or_else(|e| panic!("{name}: {e}"));
        let geometry = Geometry::from_positions(&mesh, &obj.positions)
            .unwrap_or_else(|e| panic!("{name}: {e}"));
        let load = start.elapsed().as_secs_f64();
        let problem = ConeProblem::new(&mesh, &geometry, lambda);
        let solution = match solve(&problem, &options) {
            Ok(s) => s,
            Err(e) => {
                println!("{name:<20} ERROR {e}");
                continue;
            }
        };
        let r = &solution.report;
        let c = &r.certificate;
        println!(
            "{:<20} {:>7} {:>7} {:>13} {:>4} {:>4} {:>3} {:>5} {:>3} {:>6} {:>7.3} {:>7.3} {:>7.3} {:>6} {:>9.4} {:>10.4e} {:>8.1e} {:>8.1e} {:>8.1e}",
            name,
            mesh.n_vertices(),
            mesh.n_faces(),
            format!("{:?}", r.status),
            r.regularized_iterations,
            r.exact_iterations,
            r.safeguard_rounds,
            r.linear_solver.ldlt_factorizations,
            r.linear_solver.lu_factorizations,
            r.linear_solver.solves,
            load,
            r.timings.total,
            r.linear_solver.factorization_time.as_secs_f64(),
            solution.cones.len(),
            solution.total_cone_angle / std::f64::consts::PI,
            solution.distortion,
            c.state_residual,
            c.cone_residual,
            c.duality_gap,
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
                "    certificate: {c:?}\n    extraction: clusters {} dropped {}",
                solution.extraction.n_clusters, solution.extraction.n_dropped
            );
        }
    }
}
