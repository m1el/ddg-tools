//! Verification layer 4: the authors' published results
//! (`data/optimal-cones/optimal-cones-data/data`: `<name>.obj` and `<name>-cones.txt`).
//!
//! The dataset is not part of the repository; without it every test here prints a notice and
//! passes. The solves are too slow for an unoptimized build, so in a debug build the tests
//! that solve print a notice and pass as well. Run everything with
//!
//! ```text
//! cargo test --release -p ddg-optimal-cones --test dataset -- --nocapture
//! ```
//!
//! What the tests establish is summarized in the crate README ("The published dataset").

use std::collections::BTreeMap;
use std::path::PathBuf;

use ddg_core::io::{read_cones, read_obj};
use ddg_core::{Geometry, TriMesh};
use ddg_optimal_cones::{
    Cone, ConeProblem, Parallelism, Solution, SolveOptions, Status, evaluate_cones, solve,
};

fn dataset_dir() -> Option<PathBuf> {
    let dir = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../data/optimal-cones/optimal-cones-data/data");
    if dir.is_dir() {
        Some(dir)
    } else {
        eprintln!(
            "SKIPPED: dataset directory {} not found; nothing was checked",
            dir.display()
        );
        None
    }
}

/// Whether the solves should run. They are skipped in unoptimized builds.
fn solves_enabled(test: &str) -> bool {
    if cfg!(debug_assertions) {
        eprintln!(
            "SKIPPED: {test} solves dataset meshes, which is too slow in a debug build; run \
             with --release"
        );
        false
    } else {
        true
    }
}

struct Dataset {
    name: String,
    mesh: TriMesh,
    geometry: Geometry,
    cones: Vec<Cone>,
}

fn load(dir: &std::path::Path, name: &str) -> Dataset {
    let obj = read_obj(dir.join(format!("{name}.obj"))).unwrap_or_else(|e| panic!("{e}"));
    let faces = obj
        .triangles()
        .unwrap_or_else(|| panic!("{name}: not a triangle mesh"));
    let mesh = TriMesh::new(obj.positions.len(), faces).unwrap_or_else(|e| panic!("{name}: {e}"));
    let geometry =
        Geometry::from_positions(&mesh, &obj.positions).unwrap_or_else(|e| panic!("{name}: {e}"));
    let cones = read_cones(dir.join(format!("{name}-cones.txt")), mesh.n_vertices())
        .unwrap_or_else(|e| panic!("{e}"));
    Dataset {
        name: name.to_string(),
        mesh,
        geometry,
        cones,
    }
}

fn names(dir: &std::path::Path) -> Vec<String> {
    let mut names: Vec<String> = std::fs::read_dir(dir)
        .unwrap()
        .filter_map(|e| e.ok())
        .filter_map(|e| {
            e.file_name()
                .to_str()
                .and_then(|n| n.strip_suffix(".obj"))
                .map(str::to_string)
        })
        .collect();
    names.sort();
    names
}

/// How well a cone file satisfies the optimality conditions of the canonical problem with unit
/// weights and no bounds.
///
/// If the cones are an optimum for some `λ`, the adjoint state `φ` computed from them has
/// `φ_i = λ·sign(µ_i)` at every cone and `|φ_i| ≤ λ` elsewhere; on a closed mesh the same holds
/// for `φ + c` with one constant `c`.
#[derive(Debug, Clone, Copy)]
struct LambdaFit {
    closed: bool,
    n_positive: usize,
    n_negative: usize,
    /// Least-squares fit of `sign(µ_i) (φ_i + c) = λ` over the cones.
    lambda: f64,
    shift: f64,
    /// Whether `λ` is determined by the cones. On a closed mesh whose cones all have the same
    /// sign it is not (only `λ ∓ c` is); then `lambda` is the smallest value for which
    /// `|φ + c| ≤ λ` can hold everywhere.
    identifiable: bool,
    /// `(max − min)` of `sign(µ_i) (φ_i + c)` over the cones, divided by `λ`.
    spread: f64,
    /// `max (|φ_i + c| − λ) / λ` over the interior vertices without a cone.
    excess: f64,
    gauss_bonnet_defect: f64,
}

fn fit_lambda(data: &Dataset) -> LambdaFit {
    let problem = ConeProblem::new(&data.mesh, &data.geometry, 0.0);
    let evaluation = evaluate_cones(&problem, &data.cones, Parallelism::Sequential)
        .unwrap_or_else(|e| panic!("{}: {e}", data.name));
    let phi = &evaluation.adjoint;
    let closed = data.mesh.boundary_loops().is_empty();
    let at = |positive: bool| -> Vec<f64> {
        data.cones
            .iter()
            .filter(|c| (c.angle > 0.0) == positive)
            .map(|c| phi[c.vertex as usize])
            .collect()
    };
    let (positive, negative) = (at(true), at(false));
    let mean = |x: &[f64]| x.iter().sum::<f64>() / x.len() as f64;
    let mut identifiable = true;
    let (lambda, shift) = if !closed {
        let signed: Vec<f64> = data
            .cones
            .iter()
            .map(|c| c.angle.signum() * phi[c.vertex as usize])
            .collect();
        (mean(&signed), 0.0)
    } else if !positive.is_empty() && !negative.is_empty() {
        (
            0.5 * (mean(&positive) - mean(&negative)),
            -0.5 * (mean(&positive) + mean(&negative)),
        )
    } else {
        identifiable = false;
        let (lo, hi) = phi
            .iter()
            .fold((f64::INFINITY, f64::NEG_INFINITY), |(lo, hi), &p| {
                (lo.min(p), hi.max(p))
            });
        if negative.is_empty() {
            let top = mean(&positive);
            (0.5 * (top - lo), -0.5 * (top + lo))
        } else {
            let bottom = mean(&negative);
            (0.5 * (hi - bottom), -0.5 * (hi + bottom))
        }
    };
    let signed: Vec<f64> = data
        .cones
        .iter()
        .map(|c| c.angle.signum() * (phi[c.vertex as usize] + shift))
        .collect();
    let min = signed.iter().cloned().fold(f64::INFINITY, f64::min);
    let max = signed.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    let mut is_cone = vec![false; data.mesh.n_vertices()];
    for c in &data.cones {
        is_cone[c.vertex as usize] = true;
    }
    let excess = (0..data.mesh.n_vertices())
        .filter(|&v| !is_cone[v] && !data.mesh.is_boundary_vertex(v))
        .map(|v| (phi[v] + shift).abs() - lambda)
        .fold(f64::NEG_INFINITY, f64::max);
    LambdaFit {
        closed,
        n_positive: positive.len(),
        n_negative: negative.len(),
        lambda,
        shift,
        identifiable,
        spread: (max - min) / lambda,
        excess: excess / lambda,
        gauss_bonnet_defect: evaluation.gauss_bonnet_defect,
    }
}

/// For each of the 23 cone files: the `λ` its cones imply, how consistently, and how far the
/// rest of the mesh is from satisfying `|φ| ≤ λ`.
///
/// Asserted: the two files that are optima of the canonical problem with `u = 0` on the
/// boundary pass the check (`cow01`, and the closed mesh `axes`); the check itself recovers
/// `λ` from our own solutions (see `properties.rs`). The other files do not satisfy the
/// conditions for any `λ`, which the table shows. The reason is in the tests below: their
/// angles were re-optimized without the cone penalty, so the adjoint *vanishes* at their cones
/// (the fitted `λ` of this table is then about zero), and half of them used optimized
/// boundary values.
#[test]
fn published_cones_and_the_optimality_conditions() {
    let Some(dir) = dataset_dir() else { return };
    let all = names(&dir);
    assert_eq!(
        all.len(),
        23,
        "expected the 23 meshes of the dataset: {all:?}"
    );
    println!(
        "{:<20} {:>7} {:>6} {:>6} {:>5} {:>5} {:>11} {:>11} {:>11} {:>10} {:>10} {:>9}",
        "mesh",
        "V",
        "closed",
        "cones",
        "pos",
        "neg",
        "lambda",
        "lambda*area",
        "shift c",
        "spread/lam",
        "excess/lam",
        "GB defect"
    );
    for name in &all {
        // In a debug build only the small meshes (the factorization is slow unoptimized).
        let size = std::fs::metadata(dir.join(format!("{name}.obj")))
            .unwrap()
            .len();
        if cfg!(debug_assertions) && size > 600_000 {
            println!("{name:<20} skipped in a debug build");
            continue;
        }
        let data = load(&dir, name);
        let fit = fit_lambda(&data);
        println!(
            "{:<20} {:>7} {:>6} {:>6} {:>5} {:>5} {:>11.6} {:>11.6} {:>11.3e} {:>10.2e} {:>10.2e} {:>9.1e}{}",
            name,
            data.mesh.n_vertices(),
            fit.closed,
            data.cones.len(),
            fit.n_positive,
            fit.n_negative,
            fit.lambda,
            fit.lambda * data.geometry.total_area(),
            fit.shift,
            fit.spread,
            fit.excess,
            fit.gauss_bonnet_defect,
            if fit.identifiable {
                ""
            } else {
                "   (one sign only: lambda is a lower bound)"
            }
        );
        match name.as_str() {
            // 439 cones of both signs: an optimum for λ·area = 0.01, i.e. for λ = 0.01 with
            // masses that are not area normalized.
            "cow01" => {
                assert!(fit.spread.abs() < 0.01, "{fit:?}");
                assert!(fit.excess < 0.01, "{fit:?}");
                let raw = fit.lambda * data.geometry.total_area();
                assert!((raw - 0.01).abs() < 2e-5, "lambda * area = {raw}");
            }
            // Four positive cones on a closed mesh: equal adjoint values at the cones.
            "axes" => {
                assert!(!fit.identifiable);
                assert!(fit.spread.abs() < 2e-3, "{fit:?}");
                assert!(fit.excess < 0.01, "{fit:?}");
                assert!(fit.gauss_bonnet_defect.abs() < 2e-4);
            }
            _ => {}
        }
    }
}

struct Comparison {
    ours: usize,
    file: usize,
    common: usize,
    /// Largest angle difference at the common vertices.
    worst: f64,
    /// Sum of `|angle|` of our cones that the file does not have, and vice versa.
    only_ours: f64,
    only_file: f64,
    total_ours: f64,
    total_file: f64,
}

fn compare_cones(ours: &[Cone], file: &[Cone]) -> Comparison {
    let a: BTreeMap<u32, f64> = ours.iter().map(|c| (c.vertex, c.angle)).collect();
    let b: BTreeMap<u32, f64> = file.iter().map(|c| (c.vertex, c.angle)).collect();
    let common: Vec<u32> = a.keys().filter(|v| b.contains_key(v)).copied().collect();
    let total = |m: &BTreeMap<u32, f64>| m.values().map(|x| x.abs()).sum::<f64>();
    let only = |m: &BTreeMap<u32, f64>, other: &BTreeMap<u32, f64>| {
        m.iter()
            .filter(|(v, _)| !other.contains_key(v))
            .map(|(_, x)| x.abs())
            .sum::<f64>()
    };
    Comparison {
        ours: a.len(),
        file: b.len(),
        common: common.len(),
        worst: common
            .iter()
            .map(|v| (a[v] - b[v]).abs())
            .fold(0.0, f64::max),
        only_ours: only(&a, &b),
        only_file: only(&b, &a),
        total_ours: total(&a),
        total_file: total(&b),
    }
}

fn solve_plain(data: &Dataset, lambda: f64, options: &SolveOptions<'_>) -> Solution {
    let problem = ConeProblem::new(&data.mesh, &data.geometry, lambda);
    let solution = solve(&problem, options).unwrap_or_else(|e| panic!("{}: {e}", data.name));
    assert_eq!(
        solution.report.status,
        Status::Converged,
        "{} at lambda {lambda}: {:?}",
        data.name,
        solution.report
    );
    solution
}

/// The hypothesis "the cow file names encode `λ`" (`cow1` = 1, `cow75e-2` = 0.75,
/// `cow5e-1` = 0.5, `cow01` = 0.1), tested by solving and comparing the cones with the files,
/// for the area-normalized `λ` of this crate and for `λ / area` (masses not normalized).
///
/// What holds: `cow01` is reproduced, cone for cone, by `λ = 0.01` with masses that are **not**
/// area normalized (`λ = 0.01 / area` here). Nothing else matches: with the names read as
/// 1, 0.75 and 0.5 the solver places cones at the same vertices as the files (240, 588, …) but
/// with different angles and counts, under either normalization.
#[test]
fn cow_file_names_and_lambda() {
    let Some(dir) = dataset_dir() else { return };
    if !solves_enabled("cow_file_names_and_lambda") {
        return;
    }
    let options = SolveOptions::default();
    println!(
        "{:<10} {:>8} {:<16} {:>11} {:>6} {:>6} {:>7} {:>10} {:>10} {:>10} {:>10} {:>10}",
        "file",
        "name as",
        "normalization",
        "lambda",
        "ours",
        "file",
        "common",
        "max diff",
        "only ours",
        "only file",
        "sum ours",
        "sum file"
    );
    let mut cow01_match = None;
    for (name, values) in [
        ("cow1", vec![1.0]),
        ("cow75e-2", vec![0.75]),
        ("cow5e-1", vec![0.5]),
        ("cow01", vec![0.1, 0.01]),
    ] {
        let data = load(&dir, name);
        let area = data.geometry.total_area();
        for value in values {
            for (normalization, lambda) in
                [("area normalized", value), ("lambda / area", value / area)]
            {
                let solution = solve_plain(&data, lambda, &options);
                let c = compare_cones(&solution.cones, &data.cones);
                println!(
                    "{:<10} {:>8} {:<16} {:>11.4e} {:>6} {:>6} {:>7} {:>10.2e} {:>10.2e} {:>10.2e} {:>10.4} {:>10.4}",
                    name,
                    value,
                    normalization,
                    lambda,
                    c.ours,
                    c.file,
                    c.common,
                    c.worst,
                    c.only_ours,
                    c.only_file,
                    c.total_ours,
                    c.total_file
                );
                if name == "cow01" && value == 0.01 && normalization == "lambda / area" {
                    cow01_match = Some(c);
                } else {
                    // None of the other readings reproduces its file.
                    assert!(
                        c.common < c.file || c.worst > 0.1,
                        "{name} read as {value} ({normalization}) unexpectedly matches its file"
                    );
                }
            }
        }
    }
    // cow01 = 0.01 with raw masses: every cone of the file is found, with its angle to the six
    // digits the file has; we have one cone more, of 1e-4 rad.
    let c = cow01_match.expect("cow01 was solved");
    assert_eq!(c.file, 439);
    assert_eq!(c.common, 439, "every cone of the file is a cone of ours");
    assert!(c.worst < 2e-3, "largest angle difference {}", c.worst);
    assert!(
        c.ours - c.common <= 2 && c.only_ours < 1e-3,
        "extra cones: {}",
        c.only_ours
    );
}

/// `axes` is closed and its optimum has only positive cones, so it does not depend on `λ`
/// (the cone cost is `4πλ` whatever the configuration): a comparison without any free
/// parameter. Our optimum has 14 cones in 4 edge-connected clusters; merged (paper Sec. 6.1)
/// they are the 4 cones of the file, with the same angles to `2·10⁻⁴`. One of the four sits at
/// a different vertex of its cluster: the authors' "Fréchet mean" is not ours.
#[test]
fn axes_is_reproduced_after_merging_clusters() {
    let Some(dir) = dataset_dir() else { return };
    if !solves_enabled("axes_is_reproduced_after_merging_clusters") {
        return;
    }
    let data = load(&dir, "axes");
    let raw = solve_plain(&data, 0.5, &SolveOptions::default());
    let merged_options = SolveOptions {
        merge_cone_clusters: true,
        ..SolveOptions::default()
    };
    let merged = solve_plain(&data, 0.5, &merged_options);
    println!(
        "axes: {} cones in the optimal measure, {} clusters, {} cones after merging; \
         distortion {:.6} -> {:.6}",
        raw.cones.len(),
        raw.extraction.n_clusters,
        merged.cones.len(),
        merged.distortion_optimal,
        merged.distortion
    );
    assert_eq!(raw.extraction.n_clusters, 4);
    assert_eq!(merged.cones.len(), 4);
    assert_eq!(data.cones.len(), 4);
    // λ does not matter here.
    let other = solve_plain(&data, 0.05, &SolveOptions::default());
    for (a, b) in raw.measure.iter().zip(&other.measure) {
        assert!((a - b).abs() < 1e-9);
    }
    // The angles agree, cluster by cluster.
    let mut ours: Vec<f64> = merged.cones.iter().map(|c| c.angle).collect();
    let mut file: Vec<f64> = data.cones.iter().map(|c| c.angle).collect();
    ours.sort_by(f64::total_cmp);
    file.sort_by(f64::total_cmp);
    for (a, b) in ours.iter().zip(&file) {
        println!(
            "  merged angle {a:.6}   file {b:.6}   difference {:.1e}",
            a - b
        );
        assert!((a - b).abs() < 2e-4);
    }
    // Every cone of the file lies in one of our clusters, at the vertex we chose or not.
    let same_vertex = data
        .cones
        .iter()
        .filter(|c| merged.cones.iter().any(|m| m.vertex == c.vertex))
        .count();
    println!("  {same_vertex} of 4 merged cones sit at the vertex of the file");
    assert!(same_vertex >= 3);
    for cone in &data.cones {
        assert!(
            merged
                .extraction
                .merged
                .iter()
                .any(|cluster| cluster.members.iter().any(|m| m.vertex == cone.vertex)),
            "the file's cone at vertex {} is in none of our clusters",
            cone.vertex + 1
        );
    }
    // Merging keeps the total angle and costs distortion, as reported.
    let sum: f64 = merged.cones.iter().map(|c| c.angle).sum();
    assert!((sum - 4.0 * std::f64::consts::PI).abs() < 1e-9);
    assert!(merged.extraction.distortion_change > 0.0);
    assert!(
        (merged.distortion - merged.distortion_optimal - merged.extraction.distortion_change).abs()
            < 1e-15
    );
}

/// Every mesh of the dataset is solved at the paper's default `λ = 0.5`, single-threaded, to
/// a certified optimum, without any safeguard round.
#[test]
fn every_dataset_mesh_is_solved_at_the_default_lambda() {
    let Some(dir) = dataset_dir() else { return };
    if !solves_enabled("every_dataset_mesh_is_solved_at_the_default_lambda") {
        return;
    }
    println!(
        "{:<20} {:>7} {:>7} {:>5} {:>5} {:>3} {:>5} {:>6} {:>6} {:>9} {:>11} {:>9} {:>9}",
        "mesh",
        "V",
        "F",
        "itI",
        "itII",
        "sg",
        "fact",
        "solves",
        "cones",
        "angle/pi",
        "distortion",
        "state",
        "gap"
    );
    for name in names(&dir) {
        let data = load(&dir, &name);
        let solution = solve_plain(&data, 0.5, &SolveOptions::default());
        let r = &solution.report;
        println!(
            "{:<20} {:>7} {:>7} {:>5} {:>5} {:>3} {:>5} {:>6} {:>6} {:>9.4} {:>11.4e} {:>9.1e} {:>9.1e}",
            name,
            data.mesh.n_vertices(),
            data.mesh.n_faces(),
            r.regularized_iterations,
            r.exact_iterations,
            r.safeguard_rounds,
            r.linear_solver.ldlt_factorizations + r.linear_solver.lu_factorizations,
            r.linear_solver.solves,
            solution.cones.len(),
            solution.total_cone_angle / std::f64::consts::PI,
            solution.distortion,
            r.certificate.state_residual,
            r.certificate.duality_gap
        );
        assert_eq!(r.safeguard_rounds, 0, "{name}");
        assert_eq!(r.linear_solver.lu_factorizations, 0, "{name}");
        assert!(r.certificate.holds(1e-9), "{name}: {:?}", r.certificate);
        // A certified optimum cannot be beaten: the published cones, evaluated in the same
        // problem with the same λ, have an objective that is not smaller. (Meshes with
        // boundary only: on a closed mesh the file's rounded angles do not sum to 2πχ, so
        // they are not feasible. cow00001 with its 13 836 cones is left out.)
        if !data.mesh.boundary_loops().is_empty() && data.cones.len() < 1000 {
            let file = evaluate_cones(
                &ConeProblem::new(&data.mesh, &data.geometry, 0.5),
                &data.cones,
                Parallelism::Sequential,
            )
            .unwrap();
            let file_objective = file.distortion + 0.5 * file.total_cone_angle;
            assert!(
                file_objective >= solution.objective - 1e-9,
                "{name}: the file's cones have objective {file_objective}, ours {}",
                solution.objective
            );
        }
        if data.mesh.boundary_loops().is_empty() {
            let chi = data.mesh.euler_characteristic() as f64;
            let sum: f64 = solution.cones.iter().map(|c| c.angle).sum();
            assert!(
                (sum - 2.0 * std::f64::consts::PI * chi).abs() < 1e-9,
                "{name}"
            );
        }
    }
}

// ----- the published angles are re-optimized angles -----

/// The angles that minimize the distortion for cones at the given vertices, with no cone
/// penalty: the canonical problem with these vertices as free cone sites (`wR = 0`) and a cone
/// forbidden everywhere else (`wR = +∞`). `λ` plays no role then.
fn least_distortion_angles(mesh: &TriMesh, geometry: &Geometry, vertices: &[u32]) -> Solution {
    let mut costs = vec![f64::INFINITY; mesh.n_vertices()];
    for &v in vertices {
        costs[v as usize] = 0.0;
    }
    let problem = ConeProblem {
        cone_costs: Some(costs),
        ..ConeProblem::new(mesh, geometry, 1.0)
    };
    let solution = solve(&problem, &SolveOptions::default()).unwrap_or_else(|e| panic!("{e}"));
    assert_eq!(
        solution.report.status,
        Status::Converged,
        "{:?}",
        solution.report
    );
    solution
}

fn worst_difference(solution: &Solution, cones: &[Cone]) -> f64 {
    cones
        .iter()
        .map(|c| (solution.measure[c.vertex as usize] - c.angle).abs())
        .fold(0.0, f64::max)
}

/// The published angles are not those of the penalized optimum: they are the angles that
/// minimize the distortion for the published positions when cones cost nothing there
/// (paper Sec. 8.1, "fixed cone points with free cone angles", applied to every cone).
///
/// For the files below (boundary treated as in the canonical problem, `u = 0`, or closed) the
/// re-optimized angles agree with the file to the accuracy that the six-digit positions of
/// the OBJ files allow. The bounds are the measured differences, rounded up.
#[test]
fn published_angles_minimize_the_distortion_at_the_published_positions() {
    let Some(dir) = dataset_dir() else { return };
    if !solves_enabled("published_angles_minimize_the_distortion_at_the_published_positions") {
        return;
    }
    println!(
        "{:<12} {:>5} {:>12} {:>14} {:>14}",
        "mesh", "cones", "max |d|", "max |phi| S", "max |phi| N"
    );
    for (name, bound) in [
        ("bump", 5e-4),
        ("face", 2e-4),
        ("mannequin", 1e-3),
        ("tree", 5e-3),
        ("axes", 5e-4),
        ("chair", 2e-3),
    ] {
        let data = load(&dir, name);
        let vertices: Vec<u32> = data.cones.iter().map(|c| c.vertex).collect();
        let refit = least_distortion_angles(&data.mesh, &data.geometry, &vertices);
        let worst = worst_difference(&refit, &data.cones);
        // The adjoint of the file's own cones: about zero at the cones (on a closed mesh:
        // about one common constant), far from it elsewhere.
        let evaluation = evaluate_cones(
            &ConeProblem::new(&data.mesh, &data.geometry, 0.0),
            &data.cones,
            Parallelism::Sequential,
        )
        .unwrap();
        let closed = data.mesh.boundary_loops().is_empty();
        let at_cones: Vec<f64> = vertices
            .iter()
            .map(|&v| evaluation.adjoint[v as usize])
            .collect();
        let offset = if closed {
            at_cones.iter().sum::<f64>() / at_cones.len() as f64
        } else {
            0.0
        };
        let on_cones = at_cones
            .iter()
            .fold(0.0f64, |m, p| m.max((p - offset).abs()));
        let mut is_cone = vec![false; data.mesh.n_vertices()];
        for &v in &vertices {
            is_cone[v as usize] = true;
        }
        let elsewhere = (0..data.mesh.n_vertices())
            .filter(|&v| !is_cone[v] && !data.mesh.is_boundary_vertex(v))
            .map(|v| (evaluation.adjoint[v] - offset).abs())
            .fold(0.0, f64::max);
        println!(
            "{:<12} {:>5} {:>12.2e} {:>14.3e} {:>14.3e}",
            name,
            data.cones.len(),
            worst,
            on_cones,
            elsewhere
        );
        assert!(worst < bound, "{name}: {worst}");
        assert!(
            on_cones < 0.05 * elsewhere,
            "{name}: {on_cones} vs {elsewhere}"
        );
    }
    // body: the same holds for six of its ten cones; the other four are off by up to 0.09.
    let data = load(&dir, "body");
    let vertices: Vec<u32> = data.cones.iter().map(|c| c.vertex).collect();
    let refit = least_distortion_angles(&data.mesh, &data.geometry, &vertices);
    let close = data
        .cones
        .iter()
        .filter(|c| (refit.measure[c.vertex as usize] - c.angle).abs() < 6e-3)
        .count();
    let worst = worst_difference(&refit, &data.cones);
    println!(
        "{:<12} {:>5} {:>12.2e}   ({close} of 10 within 6e-3)",
        "body", 10, worst
    );
    assert_eq!(close, 6);
    assert!(worst < 0.1);
}

/// The mesh glued to a mirror copy of itself along its boundary: a closed mesh. Vertex `v` of
/// the original keeps its index; `copy[v]` is its mirror image (`v` itself on the boundary).
///
/// The canonical problem on this double, with the seam vertices as free cone sites and `λ/2`,
/// is the problem on the original mesh in which **every boundary value of `u` is optimized**
/// (paper Sec. 8.3, Fig. 10 c) with `λ`: the solution is mirror symmetric, and each seam row of
/// the Laplacian, seam mass and seam curvature is twice that of the boundary vertex.
fn double(data: &Dataset) -> (TriMesh, Geometry, Vec<u32>) {
    let mesh = &data.mesh;
    let nv = mesh.n_vertices();
    let mut copy: Vec<u32> = (0..nv as u32).collect();
    let mut next = nv as u32;
    for (v, slot) in copy.iter_mut().enumerate() {
        if !mesh.is_boundary_vertex(v) {
            *slot = next;
            next += 1;
        }
    }
    let mut faces = mesh.faces().to_vec();
    faces.extend(mesh.faces().iter().map(|f| {
        [
            copy[f[0] as usize],
            copy[f[2] as usize],
            copy[f[1] as usize],
        ]
    }));
    let doubled =
        TriMesh::new(next as usize, faces).unwrap_or_else(|e| panic!("{}: {e}", data.name));
    let mut original = vec![0u32; next as usize];
    for v in 0..nv {
        original[v] = v as u32;
        original[copy[v] as usize] = v as u32;
    }
    let lengths: Vec<f64> = (0..doubled.n_edges())
        .map(|e| {
            let [a, b] = doubled.edge_vertices(e);
            let edge = mesh
                .find_edge(original[a] as usize, original[b] as usize)
                .expect("the copy of an edge");
            data.geometry.edge_lengths()[edge]
        })
        .collect();
    let geometry = Geometry::from_edge_lengths(&doubled, lengths).unwrap();
    (doubled, geometry, copy)
}

/// The second family of files: their angles are the least-distortion angles for their
/// positions when, in addition, the boundary values of `u` are optimized instead of fixed to
/// zero. This crate does not implement that boundary model; it is evaluated here on the
/// double of the mesh. The agreement is at the level of `10⁻²` rad, not at the level of the
/// first family, so some detail of the authors' boundary treatment is not captured; with
/// `u = 0` on the boundary the same angles are off by `0.3` to `1`.
#[test]
fn half_of_the_published_files_used_optimized_boundary_values() {
    let Some(dir) = dataset_dir() else { return };
    if !solves_enabled("half_of_the_published_files_used_optimized_boundary_values") {
        return;
    }
    println!(
        "{:<10} {:>5} {:>16} {:>16}",
        "mesh", "cones", "free boundary", "u = 0"
    );
    for name in [
        "cow1",
        "cow75e-2",
        "cow5e-1",
        "camel",
        "camel_noise1",
        "camelbad",
        "rhino",
    ] {
        let data = load(&dir, name);
        let nv = data.mesh.n_vertices();
        let (doubled, geometry, copy) = double(&data);
        let mut sites: Vec<u32> = (0..nv as u32)
            .filter(|&v| data.mesh.is_boundary_vertex(v as usize))
            .collect();
        for c in &data.cones {
            sites.push(c.vertex);
            sites.push(copy[c.vertex as usize]);
        }
        let free = worst_difference(
            &least_distortion_angles(&doubled, &geometry, &sites),
            &data.cones,
        );
        let vertices: Vec<u32> = data.cones.iter().map(|c| c.vertex).collect();
        let fixed = worst_difference(
            &least_distortion_angles(&data.mesh, &data.geometry, &vertices),
            &data.cones,
        );
        println!(
            "{:<10} {:>5} {:>16.2e} {:>16.2e}",
            name,
            data.cones.len(),
            free,
            fixed
        );
        assert!(free < 1.2e-2, "{name}: {free}");
        assert!(fixed > 0.2, "{name}: {fixed}");
    }
}

/// One file is reproduced from start to finish: `bump`. Solve the canonical problem with
/// `λ = 0.01` for the masses of the OBJ as it is (`λ = 0.01 / area` here), merge the cone
/// clusters, and re-optimize the angles at the merged positions: the five cones of the file,
/// at the same vertices, with the same angles to `3·10⁻⁴`.
#[test]
fn bump_is_reproduced_by_penalty_merging_and_reoptimized_angles() {
    let Some(dir) = dataset_dir() else { return };
    if !solves_enabled("bump_is_reproduced_by_penalty_merging_and_reoptimized_angles") {
        return;
    }
    let data = load(&dir, "bump");
    let lambda = 0.01 / data.geometry.total_area();
    let options = SolveOptions {
        merge_cone_clusters: true,
        max_merge_cluster_size: usize::MAX,
        ..SolveOptions::default()
    };
    let penalized = solve_plain(&data, lambda, &options);
    let raw = penalized.measure.iter().filter(|m| **m != 0.0).count();
    let mut ours: Vec<u32> = penalized.cones.iter().map(|c| c.vertex).collect();
    let mut file: Vec<u32> = data.cones.iter().map(|c| c.vertex).collect();
    ours.sort_unstable();
    file.sort_unstable();
    println!(
        "bump: {raw} cones in the optimal measure, {} clusters, {} cones after merging; \
         file {} cones",
        penalized.extraction.n_clusters,
        ours.len(),
        file.len()
    );
    assert!(raw > 20, "the optimum itself is not sparse here");
    assert_eq!(ours, file, "merged positions");
    // The merged cones keep the penalized (shrunk) angles: far from the file.
    let shrunk: f64 = penalized.cones.iter().map(|c| c.angle).sum();
    let published: f64 = data.cones.iter().map(|c| c.angle).sum();
    assert!(published - shrunk > 0.3, "{shrunk} vs {published}");
    // Re-optimized angles at the merged positions: the file.
    let refit = least_distortion_angles(&data.mesh, &data.geometry, &ours);
    let worst = worst_difference(&refit, &data.cones);
    println!(
        "  total angle: penalized and merged {shrunk:.4}, re-optimized {:.4}, file {published:.4}; largest angle difference {worst:.2e}",
        refit.measure.iter().sum::<f64>()
    );
    assert!(worst < 5e-4, "{worst}");
}

/// Where the published positions are exactly the support of the canonical optimum: the
/// three sparse cow files and `body`. The `λ` are for the masses of the OBJ as it is
/// (divide by the area for this crate's `λ`); each lies inside the interval for which the
/// support is optimal (`examples/published.rs refit` computes the intervals).
#[test]
fn published_positions_that_are_supports_of_the_canonical_optimum() {
    let Some(dir) = dataset_dir() else { return };
    if !solves_enabled("published_positions_that_are_supports_of_the_canonical_optimum") {
        return;
    }
    for (name, lambda_raw) in [
        ("cow1", 2.0),
        ("cow75e-2", 1.25),
        ("cow5e-1", 1.03),
        ("body", 0.0875),
    ] {
        let data = load(&dir, name);
        let lambda = lambda_raw / data.geometry.total_area();
        let solution = solve_plain(&data, lambda, &SolveOptions::default());
        let mut ours: Vec<u32> = solution.cones.iter().map(|c| c.vertex).collect();
        let mut file: Vec<u32> = data.cones.iter().map(|c| c.vertex).collect();
        ours.sort_unstable();
        file.sort_unstable();
        println!(
            "{name}: lambda = {lambda_raw} / area gives {} cones, the file has {}",
            ours.len(),
            file.len()
        );
        assert_eq!(ours, file, "{name}");
        // Same positions, different angles: ours carry the penalty, the file's do not.
        let total_ours: f64 = solution.cones.iter().map(|c| c.angle.abs()).sum();
        let total_file: f64 = data.cones.iter().map(|c| c.angle.abs()).sum();
        assert!(
            total_file > total_ours + 1.0,
            "{name}: {total_ours} vs {total_file}"
        );
    }
}
