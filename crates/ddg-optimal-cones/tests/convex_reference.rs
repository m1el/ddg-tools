//! Verification layer 2: the solver against a generic convex solver.
//!
//! The canonical problem is handed, as written in `CLAUDE.md`, to Clarabel (an interior-point
//! solver for conic quadratic programs) on small random meshes of every topology and with
//! every feature: weights, zero distortion weights, free sites, forbidden vertices, angle
//! bounds, prescribed cones. This checks the derivation (optimality system, signs, closed-mesh
//! treatment, bounds), not only the code.
//!
//! `cargo test -p ddg-optimal-cones --test convex_reference -- --nocapture` prints the table.

mod common;

use common::{
    FeatureMask, Features, N_KINDS, Reference, TestMesh, max_abs_diff, random_features,
    random_mesh, reference_solve, summary,
};
use ddg_optimal_cones::{ConeError, SolveOptions, Status, solve};
use proptest::prelude::*;

/// Agreement required of the objective, relative to `max(1, |objective|)`.
const OBJECTIVE_TOLERANCE: f64 = 1e-7;
/// Agreement required of `µ` (radians) and of `u` (relative to its largest value), where the
/// minimizer is unique. The interior-point solution is the limiting factor, see
/// [`GROWTH_SLACK`]: where the objective is flat, as on a closed convex mesh where `u = 0` is
/// optimal, its cone angles are off by up to about `10⁻⁴`. An error of convention (a sign, a
/// factor) would show as a difference of order one.
const SOLUTION_TOLERANCE: f64 = 2e-3;
/// The sharp comparison of the minimizers is quadratic growth (`optimality-system-sufficient`
/// in `lean/proofs`): every feasible `(v, ν)` has
/// `objective(v, ν) ≥ optimum + ½ Σ wE m (v − u)²`. With the reference solution as `(v, ν)`
/// this bounds its weighted distance from ours by the difference of the objectives, up to this
/// slack (relative to `max(1, |objective|)`) for the reference's own infeasibility.
const GROWTH_SLACK: f64 = 1e-8;

#[derive(Debug, Default, Clone, Copy)]
struct Tally {
    cases: usize,
    solved: usize,
    infeasible: usize,
    unreliable: usize,
    unique: usize,
    worst_objective: f64,
    worst_measure: f64,
    worst_log_scale: f64,
    safeguard_rounds: usize,
    newton_iterations: usize,
}

/// Compares one case; returns an error text on disagreement.
fn compare(mesh: &TestMesh, f: &Features, tally: &mut Tally) -> Result<(), String> {
    let problem = f.problem(mesh);
    let reference = reference_solve(mesh, f);
    let result = solve(&problem, &SolveOptions::default());
    tally.cases += 1;
    let describe = || format!("{} with {f:.3?}", mesh.name);
    match (result, reference) {
        (Err(ConeError::Infeasible { .. }), Reference::Infeasible) => {
            tally.infeasible += 1;
            Ok(())
        }
        (Err(e), reference) => Err(format!(
            "solver returned an error ({e}) but the reference says {reference:?}: {}",
            describe()
        )),
        (Ok(solution), Reference::Infeasible) => Err(format!(
            "reference says infeasible, solver returned {:?} with objective {}: {}",
            solution.report.status,
            solution.objective,
            describe()
        )),
        (Ok(solution), Reference::Unreliable(status)) => {
            // The generic solver gave up; the solver under test must still certify its result.
            tally.unreliable += 1;
            if solution.report.status != Status::Converged {
                return Err(format!(
                    "not converged (reference unreliable: {status}): {}\n  {}",
                    summary(&solution),
                    describe()
                ));
            }
            Ok(())
        }
        (
            Ok(solution),
            Reference::Solved {
                objective,
                log_scale,
                measure,
            },
        ) => {
            if solution.report.status != Status::Converged {
                return Err(format!(
                    "not converged (reference objective {objective:.12e}): {}\n  {}",
                    summary(&solution),
                    describe()
                ));
            }
            tally.solved += 1;
            tally.safeguard_rounds += solution.report.safeguard_rounds;
            tally.newton_iterations +=
                solution.report.regularized_iterations + solution.report.exact_iterations;
            let scale = objective.abs().max(1.0);
            let objective_error = (solution.objective - objective).abs() / scale;
            tally.worst_objective = tally.worst_objective.max(objective_error);
            if objective_error > OBJECTIVE_TOLERANCE {
                return Err(format!(
                    "objective {} vs reference {objective}: {}",
                    solution.objective,
                    describe()
                ));
            }
            // Quadratic growth around our (certified) minimizer.
            let mass =
                ddg_core::operators::lumped_mass_normalized(&mesh.mesh, &mesh.geometry).unwrap();
            let weighted_distance: f64 = (0..mesh.mesh.n_vertices())
                .map(|v| {
                    let w = f.distortion_weights.as_ref().map_or(1.0, |w| w[v]) * mass[v];
                    0.5 * w * (solution.log_scale_optimal[v] - log_scale[v]).powi(2)
                })
                .sum();
            if weighted_distance > (objective - solution.objective).max(0.0) + GROWTH_SLACK * scale
            {
                return Err(format!(
                    "quadratic growth violated: weighted distance {weighted_distance:.3e} of the \
                     reference u, objective difference {:.3e}: {}",
                    objective - solution.objective,
                    describe()
                ));
            }
            if f.strictly_weighted(mesh) {
                tally.unique += 1;
                let measure_error = max_abs_diff(&solution.measure, &measure);
                let u_scale = log_scale.iter().fold(1.0f64, |m, u| m.max(u.abs()));
                let log_scale_error =
                    max_abs_diff(&solution.log_scale_optimal, &log_scale) / u_scale;
                tally.worst_measure = tally.worst_measure.max(measure_error);
                tally.worst_log_scale = tally.worst_log_scale.max(log_scale_error);
                if measure_error > SOLUTION_TOLERANCE || log_scale_error > SOLUTION_TOLERANCE {
                    return Err(format!(
                        "measure differs by {measure_error}, u by {log_scale_error} (relative): \
                         {}",
                        describe()
                    ));
                }
            }
            Ok(())
        }
    }
}

/// A deterministic sweep over all topologies and feature sets, with a table of the worst
/// disagreements.
#[test]
fn sweep_against_generic_convex_solver() {
    let kinds = ["disk", "annulus", "cylinder", "sphere", "torus"];
    let masks: [(&str, FeatureMask); 3] = [
        (
            "plain",
            FeatureMask {
                weights: false,
                zero_distortion_weights: false,
                free_sites: false,
                forbidden: false,
                bounds: false,
                prescribed: false,
            },
        ),
        ("all features, wE > 0", FeatureMask::UNIQUE),
        ("all features, some wE = 0", FeatureMask::ALL),
    ];
    println!(
        "{:<10} {:<26} {:>5} {:>6} {:>6} {:>6} {:>10} {:>10} {:>10} {:>7} {:>6}",
        "topology",
        "features",
        "cases",
        "solved",
        "infeas",
        "unrel",
        "objective",
        "measure",
        "u (rel)",
        "newton",
        "safeg"
    );
    let mut failures = Vec::new();
    for (kind, kind_name) in kinds.iter().enumerate() {
        for (mask_name, mask) in &masks {
            let mut tally = Tally::default();
            for seed in 0..60u64 {
                let mesh = random_mesh(kind, seed);
                let features = random_features(&mesh, 1000 * kind as u64 + seed, *mask);
                if let Err(e) = compare(&mesh, &features, &mut tally) {
                    failures.push(e);
                }
            }
            println!(
                "{:<10} {:<26} {:>5} {:>6} {:>6} {:>6} {:>10.2e} {:>10.2e} {:>10.2e} {:>7.1} {:>6}",
                kind_name,
                mask_name,
                tally.cases,
                tally.solved,
                tally.infeasible,
                tally.unreliable,
                tally.worst_objective,
                tally.worst_measure,
                tally.worst_log_scale,
                tally.newton_iterations as f64 / tally.solved.max(1) as f64,
                tally.safeguard_rounds,
            );
        }
    }
    assert!(
        failures.is_empty(),
        "{} disagreements, first ones:\n{}",
        failures.len(),
        failures
            .iter()
            .take(5)
            .cloned()
            .collect::<Vec<_>>()
            .join("\n")
    );
}

proptest! {
    #![proptest_config(ProptestConfig { cases: 200, failure_persistence: None, ..ProptestConfig::default() })]

    /// Random meshes, random features: same objective as the generic solver, same minimizer
    /// where it is unique, and the same verdict on feasibility.
    #[test]
    fn solver_agrees_with_generic_convex_solver(
        kind in 0..N_KINDS,
        mesh_seed in 0..100_000u64,
        feature_seed in 0..100_000u64,
        unique in any::<bool>(),
    ) {
        let mesh = random_mesh(kind, mesh_seed);
        let mask = if unique { FeatureMask::UNIQUE } else { FeatureMask::ALL };
        let features = random_features(&mesh, feature_seed, mask);
        let mut tally = Tally::default();
        if let Err(e) = compare(&mesh, &features, &mut tally) {
            return Err(TestCaseError::fail(e));
        }
    }
}

/// Reruns one case of the property test and prints everything about it:
/// `CASE="kind mesh_seed feature_seed unique" cargo test ... -- --ignored debug_case --nocapture`.
#[test]
#[ignore = "diagnostic: set CASE=\"kind mesh_seed feature_seed unique\""]
fn debug_case() {
    let case = std::env::var("CASE").expect("CASE=\"kind mesh_seed feature_seed unique\"");
    let parts: Vec<&str> = case.split_whitespace().collect();
    let kind: usize = parts[0].parse().unwrap();
    let mesh_seed: u64 = parts[1].parse().unwrap();
    let feature_seed: u64 = parts[2].parse().unwrap();
    let unique: bool = parts[3].parse().unwrap();
    let mesh = random_mesh(kind, mesh_seed);
    let mask = if unique {
        FeatureMask::UNIQUE
    } else {
        FeatureMask::ALL
    };
    let features = random_features(&mesh, feature_seed, mask);
    println!(
        "{} ({} vertices): {features:.4?}",
        mesh.name,
        mesh.mesh.n_vertices()
    );
    let solution = solve(&features.problem(&mesh), &SolveOptions::default());
    let reference = reference_solve(&mesh, &features);
    match (&solution, &reference) {
        (
            Ok(solution),
            Reference::Solved {
                objective,
                log_scale,
                measure,
            },
        ) => {
            println!("{}", summary(solution));
            println!(
                "objective: ours {:.15e}, reference {:.15e}, ours - reference {:.3e}",
                solution.objective,
                objective,
                solution.objective - objective
            );
            println!(
                "{:>4} {:>14} {:>14} {:>10} {:>14} {:>14} {:>14}",
                "v", "mu ours", "mu ref", "diff", "u ours", "u ref", "phi ours"
            );
            for v in 0..mesh.mesh.n_vertices() {
                println!(
                    "{:>4} {:>14.9} {:>14.9} {:>10.2e} {:>14.9} {:>14.9} {:>14.9}",
                    v,
                    solution.measure[v],
                    measure[v],
                    solution.measure[v] - measure[v],
                    solution.log_scale_optimal[v],
                    log_scale[v],
                    solution.adjoint[v]
                );
            }
        }
        (solution, reference) => println!(
            "ours: {:?}\nreference: {reference:?}",
            solution.as_ref().map(summary)
        ),
    }
}
