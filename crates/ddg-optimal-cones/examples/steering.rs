//! What the solver offers where the plain problem gives nothing usable: closed convex meshes,
//! and a mesh with needle-like spikes.
//!
//! ```text
//! cargo run --release -p ddg-optimal-cones --example steering -- closed
//! cargo run --release -p ddg-optimal-cones --example steering -- spikes [--data DIR]
//! cargo run --release -p ddg-optimal-cones --example steering -- closed-dataset [--data DIR] [meshes]
//! ```
//!
//! `closed`: an icosphere, an ellipsoid and a capped cylinder (all closed, convex, genus 0),
//! each solved for a range of `λ` under eight settings: plain; positive cones only; one free
//! cone site, without and with an upper bound of `π` on every angle, without and with positive
//! cones only; an upper bound alone; and a handful of free sites with cones forbidden
//! elsewhere. One line per solve: cones, clusters, the smallest and largest angle, total angle,
//! distortion, and the same after merging the clusters.
//!
//! `closed-dataset`: the closed meshes of the dataset (`axes`, `bunny`, `chair`; `brain` on
//! request) with positive cones only: optimum, clusters merged, angles re-optimized.
//!
//! `spikes`: the dataset's `camelshotnoise` (22 needles of length 0.5 and radius `10⁻⁴`),
//! solved as it is, with merging (this crate's rule, and the rule that also merges cones of
//! opposite sign), with positive cones only, and with zero distortion weight at the spike
//! tips.

use std::f64::consts::PI;
use std::path::PathBuf;

use ddg_core::generators::{icosphere, open_cylinder};
use ddg_core::io::{Cone, read_obj};
use ddg_core::operators::angle_defects;
use ddg_core::{Geometry, TriMesh};
use ddg_optimal_cones::{
    ConeProblem, Parallelism, Solution, SolveOptions, Status, evaluate_cones, solve,
};

struct Shape {
    name: &'static str,
    positions: Vec<[f64; 3]>,
    mesh: TriMesh,
    geometry: Geometry,
}

impl Shape {
    fn new(name: &'static str, positions: Vec<[f64; 3]>, faces: Vec<[u32; 3]>) -> Shape {
        let mesh = TriMesh::new(positions.len(), faces).expect("valid mesh");
        let geometry = Geometry::from_positions(&mesh, &positions).expect("valid geometry");
        Shape {
            name,
            positions,
            mesh,
            geometry,
        }
    }
}

/// A cylinder closed by a fan over each end.
fn capped_cylinder(n_around: usize, n_height: usize, radius: f64, height: f64) -> Shape {
    let open = open_cylinder(n_around, n_height, radius, height).expect("a cylinder");
    let mesh = TriMesh::new(open.positions.len(), open.faces.clone()).expect("valid mesh");
    let mut positions = open.positions.clone();
    let mut faces = open.faces.clone();
    let bottom = positions.len() as u32;
    positions.push([0.0, 0.0, 0.0]);
    let top = positions.len() as u32;
    positions.push([0.0, 0.0, height]);
    for h in 0..3 * mesh.n_faces() {
        if mesh.twin(h).is_none() {
            let (a, b) = (mesh.from_vertex(h), mesh.to_vertex(h));
            let centre = if positions[a][2] < 0.5 * height {
                bottom
            } else {
                top
            };
            faces.push([b as u32, a as u32, centre]);
        }
    }
    Shape::new("capped cylinder", positions, faces)
}

fn shapes() -> Vec<Shape> {
    let sphere = icosphere(3).expect("an icosphere");
    let ellipsoid: Vec<[f64; 3]> = sphere
        .positions
        .iter()
        .map(|p| [1.6 * p[0], p[1], 0.6 * p[2]])
        .collect();
    vec![
        Shape::new("icosphere", sphere.positions.clone(), sphere.faces.clone()),
        Shape::new("ellipsoid", ellipsoid, sphere.faces.clone()),
        capped_cylinder(24, 12, 0.5, 2.0),
    ]
}

/// `k` vertices spread over the surface: the first is the one of largest angle defect, each
/// further one is the vertex farthest (in a straight line) from those chosen so far.
fn spread_sites(shape: &Shape, k: usize) -> Vec<usize> {
    let omega = angle_defects(&shape.mesh, &shape.geometry).unwrap();
    let nv = shape.mesh.n_vertices();
    let first = (0..nv)
        .max_by(|&a, &b| omega[a].total_cmp(&omega[b]))
        .expect("a vertex");
    let distance = |a: usize, b: usize| -> f64 {
        (0..3)
            .map(|c| (shape.positions[a][c] - shape.positions[b][c]).powi(2))
            .sum::<f64>()
            .sqrt()
    };
    let mut sites = vec![first];
    let mut nearest: Vec<f64> = (0..nv).map(|v| distance(v, first)).collect();
    while sites.len() < k {
        let next = (0..nv)
            .max_by(|&a, &b| nearest[a].total_cmp(&nearest[b]))
            .expect("a vertex");
        sites.push(next);
        for (v, d) in nearest.iter_mut().enumerate() {
            *d = d.min(distance(v, next));
        }
    }
    sites
}

struct Setting {
    name: &'static str,
    free_sites: usize,
    forbid_elsewhere: bool,
    positive: bool,
    upper: Option<f64>,
}

fn describe(solution: &Solution) -> String {
    let mu: Vec<f64> = solution
        .measure
        .iter()
        .cloned()
        .filter(|m| *m != 0.0)
        .collect();
    let min = mu.iter().cloned().fold(f64::INFINITY, f64::min);
    let max = mu.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    let beyond = mu.iter().filter(|m| m.abs() >= 2.0 * PI).count();
    let e = &solution.extraction;
    let merged_range = solution
        .cones
        .iter()
        .fold((f64::INFINITY, f64::NEG_INFINITY), |(lo, hi), c| {
            (lo.min(c.angle), hi.max(c.angle))
        });
    format!(
        "{:>5} {:>5} {:>8.3} {:>8.3} {:>8.3} {:>3} {:>10.4e} | {:>5} {:>8.3} {:>8.3} {:>3} {:>10.4e}",
        mu.len(),
        e.n_clusters,
        min / PI,
        max / PI,
        mu.iter().map(|m| m.abs()).sum::<f64>() / PI,
        beyond,
        solution.distortion_optimal,
        solution.cones.len(),
        merged_range.0 / PI,
        merged_range.1 / PI,
        e.n_beyond_two_pi,
        solution.distortion
    )
}

fn closed() {
    let settings = [
        Setting {
            name: "plain",
            free_sites: 0,
            forbid_elsewhere: false,
            positive: false,
            upper: None,
        },
        Setting {
            name: "mu >= 0",
            free_sites: 0,
            forbid_elsewhere: false,
            positive: true,
            upper: None,
        },
        Setting {
            name: "mu <= pi",
            free_sites: 0,
            forbid_elsewhere: false,
            positive: false,
            upper: Some(PI),
        },
        Setting {
            name: "1 free site",
            free_sites: 1,
            forbid_elsewhere: false,
            positive: false,
            upper: None,
        },
        Setting {
            name: "1 free site, mu >= 0",
            free_sites: 1,
            forbid_elsewhere: false,
            positive: true,
            upper: None,
        },
        Setting {
            name: "1 free site, mu <= pi",
            free_sites: 1,
            forbid_elsewhere: false,
            positive: false,
            upper: Some(PI),
        },
        Setting {
            name: "1 free site, 0 <= mu <= pi",
            free_sites: 1,
            forbid_elsewhere: false,
            positive: true,
            upper: Some(PI),
        },
        Setting {
            name: "8 free sites, none elsewhere, 0 <= mu <= pi",
            free_sites: 8,
            forbid_elsewhere: true,
            positive: true,
            upper: Some(PI),
        },
    ];
    let options = SolveOptions {
        merge_cone_clusters: true,
        max_merge_cluster_size: usize::MAX,
        ..SolveOptions::default()
    };
    for shape in shapes() {
        let nv = shape.mesh.n_vertices();
        println!(
            "{}: {} vertices, {} faces, area {:.4}",
            shape.name,
            nv,
            shape.mesh.n_faces(),
            shape.geometry.total_area()
        );
        println!(
            "  {:<44} {:>7} {:>9} | {:>5} {:>5} {:>8} {:>8} {:>8} {:>3} {:>10} | {:>5} {:>8} {:>8} {:>3} {:>10}",
            "setting",
            "lambda",
            "status",
            "cones",
            "clus",
            "min/pi",
            "max/pi",
            "sum/pi",
            "x",
            "E",
            "merge",
            "min/pi",
            "max/pi",
            "x",
            "E merged"
        );
        for setting in &settings {
            let sites = spread_sites(&shape, setting.free_sites.max(1));
            let lambdas: &[f64] = if setting.forbid_elsewhere {
                &[0.1]
            } else {
                &[1.0, 0.3, 0.1, 0.03, 0.01, 0.003]
            };
            for &lambda in lambdas {
                let mut costs = vec![
                    if setting.forbid_elsewhere {
                        f64::INFINITY
                    } else {
                        1.0
                    };
                    nv
                ];
                for &v in sites.iter().take(setting.free_sites) {
                    costs[v] = 0.0;
                }
                let problem = ConeProblem {
                    cone_costs: Some(costs),
                    lower_bounds: setting.positive.then(|| vec![0.0; nv]),
                    upper_bounds: setting.upper.map(|hi| vec![hi; nv]),
                    ..ConeProblem::new(&shape.mesh, &shape.geometry, lambda)
                };
                match solve(&problem, &options) {
                    Ok(solution) => println!(
                        "  {:<44} {:>7} {:>9} | {}",
                        setting.name,
                        lambda,
                        format!("{:?}", solution.report.status),
                        describe(&solution)
                    ),
                    Err(e) => println!("  {:<44} {:>7} error: {e}", setting.name, lambda),
                }
            }
        }
    }
    println!(
        "\ncones, clus, min, max, sum, E: the optimal measure (number of cones, same-sign \
         clusters, smallest and largest angle and total |angle| in units of pi, distortion).\n\
         x: number of angles at or beyond 2 pi in magnitude (not realizable as a cone).\n\
         merge ...: the same after merging every cluster (no limit on the cluster size)."
    );
}

/// Merges edge-connected clusters of at most `limit` cones, of one sign (`by_sign`) or of any
/// sign, into the member of largest magnitude. (The crate merges by sign at the Fréchet mean;
/// the position within a cluster does not matter for what this example shows.)
fn merge(mesh: &TriMesh, mu: &[f64], by_sign: bool, limit: usize) -> Vec<Cone> {
    let nv = mesh.n_vertices();
    let mut visited = vec![false; nv];
    let mut cones = Vec::new();
    for seed in 0..nv {
        if visited[seed] || mu[seed] == 0.0 {
            continue;
        }
        visited[seed] = true;
        let mut members = vec![seed];
        let mut head = 0;
        while head < members.len() {
            let v = members[head];
            head += 1;
            for w in mesh.vertex_neighbors(v) {
                if !visited[w] && mu[w] != 0.0 && (!by_sign || (mu[w] > 0.0) == (mu[seed] > 0.0)) {
                    visited[w] = true;
                    members.push(w);
                }
            }
        }
        if members.len() > limit {
            cones.extend(members.iter().map(|&v| Cone {
                vertex: v as u32,
                angle: mu[v],
            }));
            continue;
        }
        let at = *members
            .iter()
            .max_by(|&&a, &&b| mu[a].abs().total_cmp(&mu[b].abs()))
            .expect("a member");
        let angle: f64 = members.iter().map(|&v| mu[v]).sum();
        if angle.abs() > 1e-12 {
            cones.push(Cone {
                vertex: at as u32,
                angle,
            });
        }
    }
    cones.sort_by_key(|c| c.vertex);
    cones
}

fn summary(cones: &[Cone]) -> String {
    let negative = cones.iter().filter(|c| c.angle < 0.0).count();
    let min = cones.iter().map(|c| c.angle).fold(f64::INFINITY, f64::min);
    let max = cones
        .iter()
        .map(|c| c.angle)
        .fold(f64::NEG_INFINITY, f64::max);
    let beyond = cones.iter().filter(|c| c.angle.abs() >= 2.0 * PI).count();
    format!(
        "{:>5} cones ({:>4} negative), angles in [{:>7.3}, {:>6.3}], {:>3} at or beyond 2 pi",
        cones.len(),
        negative,
        min,
        max,
        beyond
    )
}

fn spikes(dir: &std::path::Path) {
    let obj = read_obj(dir.join("camelshotnoise.obj")).expect("readable OBJ");
    let mesh = TriMesh::new(obj.positions.len(), obj.triangles().expect("triangles")).unwrap();
    let geometry = Geometry::from_positions(&mesh, &obj.positions).unwrap();
    let nv = mesh.n_vertices();
    let omega = angle_defects(&mesh, &geometry).unwrap();
    let tips: Vec<usize> = (0..nv)
        .filter(|&v| !mesh.is_boundary_vertex(v) && omega[v] > 1.5 * PI)
        .collect();
    println!(
        "camelshotnoise: {} vertices, {} spike tips (angle defect above 3 pi / 2)",
        nv,
        tips.len()
    );
    let plain = ConeProblem::new(&mesh, &geometry, 0.0);
    let distortion_of = |cones: &[Cone]| -> f64 {
        evaluate_cones(&plain, cones, Parallelism::Sequential)
            .expect("valid cones")
            .distortion
    };
    for lambda in [0.5, 0.1, 0.02] {
        println!("lambda = {lambda}");
        let problem = ConeProblem::new(&mesh, &geometry, lambda);
        let solution = solve(&problem, &SolveOptions::default()).expect("valid problem");
        assert_eq!(solution.report.status, Status::Converged);
        println!(
            "  optimum                               {}, distortion {:.4e}",
            summary(&solution.cones),
            solution.distortion
        );
        for (label, by_sign, limit) in [
            ("merged, any sign, clusters <= 8  ", false, 8),
            ("merged, any sign, all clusters   ", false, usize::MAX),
            ("merged, one sign, clusters <= 8  ", true, 8),
            ("merged, one sign, all clusters   ", true, usize::MAX),
        ] {
            let cones = merge(&mesh, &solution.measure, by_sign, limit);
            println!(
                "  {label}     {}, distortion {:.4e}",
                summary(&cones),
                distortion_of(&cones)
            );
        }
        // What the crate does (one sign, Fréchet mean), with its default limit and without.
        for (label, limit) in [
            ("crate, default limit 8", 8),
            ("crate, no limit      ", usize::MAX),
        ] {
            let merged = solve(
                &problem,
                &SolveOptions {
                    merge_cone_clusters: true,
                    max_merge_cluster_size: limit,
                    ..SolveOptions::default()
                },
            )
            .expect("valid problem");
            println!(
                "  {label}                {}, distortion {:.4e} (optimum {:.4e})",
                summary(&merged.cones),
                merged.distortion,
                merged.distortion_optimal
            );
        }
        let positive = solve(
            &ConeProblem {
                lower_bounds: Some(vec![0.0; nv]),
                ..problem.clone()
            },
            &SolveOptions::default(),
        )
        .expect("valid problem");
        println!(
            "  mu >= 0                               {}, distortion {:.4e}, status {:?}",
            summary(&positive.cones),
            positive.distortion,
            positive.report.status
        );
        // The spikes excluded from the distortion measure: wE = 0 at the tips.
        let mut weights = vec![1.0; nv];
        for &v in &tips {
            weights[v] = 0.0;
        }
        for (label, lower) in [
            ("wE = 0 at the tips          ", None),
            ("wE = 0 at the tips, mu >= 0 ", Some(vec![0.0; nv])),
        ] {
            let ignoring = solve(
                &ConeProblem {
                    distortion_weights: Some(weights.clone()),
                    lower_bounds: lower,
                    ..problem.clone()
                },
                &SolveOptions::default(),
            )
            .expect("valid problem");
            println!(
                "  {label}          {}, weighted distortion {:.4e}, status {:?}",
                summary(&ignoring.cones),
                ignoring.distortion,
                ignoring.report.status
            );
        }
    }
}

/// Closed meshes of the dataset with positive cones only: the optimum (which does not depend
/// on `λ` and comes in clusters), the clusters merged, and the angles re-optimized at the
/// merged positions within `0 ≤ µ ≤ hi`.
fn closed_dataset(dir: &std::path::Path, names: &[String]) {
    let hi = 1.9 * PI;
    println!(
        "{:<8} {:>7} | {:>5} {:>5} {:>10} | {:>6} {:>8} {:>3} {:>10} | {:>6} {:>8} {:>8} {:>10}",
        "mesh",
        "V",
        "cones",
        "clus",
        "E",
        "merged",
        "max/pi",
        "x",
        "E merged",
        "refit",
        "min/pi",
        "max/pi",
        "E refit"
    );
    for name in names {
        let obj = read_obj(dir.join(format!("{name}.obj"))).expect("readable OBJ");
        let mesh = TriMesh::new(obj.positions.len(), obj.triangles().expect("triangles")).unwrap();
        let geometry = Geometry::from_positions(&mesh, &obj.positions).unwrap();
        let nv = mesh.n_vertices();
        let positive = ConeProblem {
            lower_bounds: Some(vec![0.0; nv]),
            ..ConeProblem::new(&mesh, &geometry, 0.5)
        };
        let merged = solve(
            &positive,
            &SolveOptions {
                merge_cone_clusters: true,
                max_merge_cluster_size: usize::MAX,
                ..SolveOptions::default()
            },
        )
        .expect("valid problem");
        assert_eq!(merged.report.status, Status::Converged);
        // Free sites at the merged positions, no cone elsewhere, angles in [0, hi].
        let mut costs = vec![f64::INFINITY; nv];
        for c in &merged.cones {
            costs[c.vertex as usize] = 0.0;
        }
        let refit = solve(
            &ConeProblem {
                cone_costs: Some(costs),
                lower_bounds: Some(vec![0.0; nv]),
                upper_bounds: Some(vec![hi; nv]),
                ..ConeProblem::new(&mesh, &geometry, 0.5)
            },
            &SolveOptions::default(),
        );
        let raw = merged.measure.iter().filter(|m| **m != 0.0).count();
        let max_merged = merged.cones.iter().map(|c| c.angle).fold(0.0, f64::max);
        let refit_text = match &refit {
            Ok(r) => format!(
                "{:>6} {:>8.3} {:>8.3} {:>10.4e}{}",
                r.cones.len(),
                r.cones
                    .iter()
                    .map(|c| c.angle)
                    .fold(f64::INFINITY, f64::min)
                    / PI,
                r.cones.iter().map(|c| c.angle).fold(0.0, f64::max) / PI,
                r.distortion,
                if r.report.status == Status::Converged {
                    ""
                } else {
                    "  NOT CONVERGED"
                }
            ),
            Err(e) => format!("error: {e}"),
        };
        println!(
            "{:<8} {:>7} | {:>5} {:>5} {:>10.4e} | {:>6} {:>8.3} {:>3} {:>10.4e} | {}",
            name,
            nv,
            raw,
            merged.extraction.n_clusters,
            merged.distortion_optimal,
            merged.cones.len(),
            max_merged / PI,
            merged.extraction.n_beyond_two_pi,
            merged.distortion,
            refit_text
        );
    }
    println!(
        "\ncones, clus, E: the optimum with mu >= 0. merged, max, x, E merged: after merging \
         every cluster (x: angles at or beyond 2 pi).\nrefit: angles re-optimized at the merged \
         positions within [0, 1.9 pi], no cone elsewhere."
    );
}

fn main() {
    let mut dir = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../data/optimal-cones/optimal-cones-data/data");
    let mut args = std::env::args().skip(1);
    let command = args.next().unwrap_or_default();
    let mut names: Vec<String> = Vec::new();
    while let Some(arg) = args.next() {
        match arg.as_str() {
            "--data" => dir = PathBuf::from(args.next().expect("--data needs a value")),
            name => names.push(name.to_string()),
        }
    }
    match command.as_str() {
        "closed" => closed(),
        "spikes" => spikes(&dir),
        "closed-dataset" => {
            if names.is_empty() {
                names = ["axes", "bunny", "chair"].map(str::to_string).to_vec();
            }
            closed_dataset(&dir, &names)
        }
        _ => {
            eprintln!(
                "usage: steering closed | spikes [--data DIR] | closed-dataset [--data DIR] [meshes]"
            );
            std::process::exit(2);
        }
    }
}
