//! Solves one dataset mesh for a given `λ` and compares the cones with the authors' file.
//!
//! ```text
//! cargo run --release -p ddg-optimal-cones --example compare_cones -- NAME LAMBDA [options]
//!
//!   --data DIR        dataset directory
//!   --raw-area        divide LAMBDA by the total area of the mesh (i.e. LAMBDA is meant for
//!                     masses that are not area normalized)
//!   --merge N         merge cone clusters of up to N cones
//!   --list            list every cone of both sets
//! ```

use std::collections::BTreeMap;
use std::path::PathBuf;

use ddg_core::io::{read_cones, read_obj};
use ddg_core::{Geometry, TriMesh};
use ddg_optimal_cones::{ConeProblem, SolveOptions, solve};

fn main() {
    let mut dir = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../data/optimal-cones/optimal-cones-data/data");
    let mut positional: Vec<String> = Vec::new();
    let mut raw_area = false;
    let mut list = false;
    let mut options = SolveOptions::default();
    let mut args = std::env::args().skip(1);
    while let Some(arg) = args.next() {
        match arg.as_str() {
            "--data" => dir = PathBuf::from(args.next().expect("--data needs a value")),
            "--raw-area" => raw_area = true,
            "--list" => list = true,
            "--merge" => {
                options.merge_cone_clusters = true;
                options.max_merge_cluster_size = args
                    .next()
                    .expect("--merge needs a size")
                    .parse()
                    .expect("--merge: an integer");
            }
            other => positional.push(other.to_string()),
        }
    }
    let [name, lambda] = positional.as_slice() else {
        eprintln!("usage: compare_cones NAME LAMBDA [--raw-area] [--merge N] [--list]");
        std::process::exit(2);
    };
    let mut lambda: f64 = lambda.parse().expect("LAMBDA: a number");
    let obj = read_obj(dir.join(format!("{name}.obj"))).unwrap_or_else(|e| panic!("{e}"));
    let faces = obj.triangles().expect("triangle mesh");
    let mesh = TriMesh::new(obj.positions.len(), faces).unwrap_or_else(|e| panic!("{e}"));
    let geometry =
        Geometry::from_positions(&mesh, &obj.positions).unwrap_or_else(|e| panic!("{e}"));
    let theirs = read_cones(dir.join(format!("{name}-cones.txt")), mesh.n_vertices())
        .unwrap_or_else(|e| panic!("{e}"));
    if raw_area {
        lambda /= geometry.total_area();
    }
    let problem = ConeProblem::new(&mesh, &geometry, lambda);
    let solution = solve(&problem, &options).unwrap_or_else(|e| panic!("{e}"));
    let ours: BTreeMap<u32, f64> = solution.cones.iter().map(|c| (c.vertex, c.angle)).collect();
    let file: BTreeMap<u32, f64> = theirs.iter().map(|c| (c.vertex, c.angle)).collect();
    let common: Vec<u32> = ours
        .keys()
        .filter(|v| file.contains_key(v))
        .copied()
        .collect();
    let worst = common
        .iter()
        .map(|v| (ours[v] - file[v]).abs())
        .fold(0.0, f64::max);
    let total = |m: &BTreeMap<u32, f64>| m.values().map(|a| a.abs()).sum::<f64>();
    println!(
        "{name}: area {:.6}, lambda (area normalized) {lambda:.6e}, status {:?}, clusters {}, merged {}",
        geometry.total_area(),
        solution.report.status,
        solution.extraction.n_clusters,
        solution.extraction.merged.len(),
    );
    println!(
        "  ours: {} cones, total angle {:.6}   file: {} cones, total angle {:.6}",
        ours.len(),
        total(&ours),
        file.len(),
        total(&file)
    );
    println!(
        "  common vertices: {}, largest angle difference there: {:.3e}",
        common.len(),
        worst
    );
    let only_ours: Vec<_> = ours.iter().filter(|(v, _)| !file.contains_key(v)).collect();
    let only_file: Vec<_> = file.iter().filter(|(v, _)| !ours.contains_key(v)).collect();
    let sum = |x: &[(&u32, &f64)]| x.iter().map(|(_, a)| a.abs()).sum::<f64>();
    println!(
        "  only ours: {} (total |angle| {:.3e})   only file: {} (total |angle| {:.3e})",
        only_ours.len(),
        sum(&only_ours),
        only_file.len(),
        sum(&only_file)
    );
    if list {
        for (v, a) in &ours {
            match file.get(v) {
                Some(b) => println!("    {:>7} ours {:>12.6} file {:>12.6}", v + 1, a, b),
                None => println!("    {:>7} ours {:>12.6}", v + 1, a),
            }
        }
        for (v, b) in &only_file {
            println!("    {:>7} {:>17} file {:>12.6}", *v + 1, "", b);
        }
    }
}
