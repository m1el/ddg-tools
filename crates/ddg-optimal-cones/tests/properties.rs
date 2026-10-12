//! Verification layer 5: properties the solution must have, on random meshes of every
//! topology and with random features.

// Per-vertex loops that index several arrays read best with the vertex index.
#![allow(clippy::needless_range_loop)]

mod common;

use std::f64::consts::PI;

use common::{
    FeatureMask, Features, N_KINDS, TestMesh, max_abs_diff, random_features, random_mesh, summary,
};
use ddg_core::generators::{SplitMix64, bumpy_grid_disk, icosphere};
use ddg_core::operators::angle_defects;
use ddg_optimal_cones::{
    BoundaryTreatment, Cone, ConeError, Parallelism, Solution, SolveOptions, Status, certify,
    evaluate_cones,
};
use proptest::prelude::*;

fn solve_ok(mesh: &TestMesh, f: &Features) -> Option<Solution> {
    match common::solve_t(&f.problem(mesh), BoundaryTreatment::Fixed, &SolveOptions::default()) {
        Ok(solution) => {
            assert_eq!(
                solution.report.status,
                Status::Converged,
                "{}: {}",
                mesh.name,
                summary(&solution)
            );
            Some(solution)
        }
        Err(ConeError::Infeasible { .. }) => None,
        Err(e) => panic!("{}: {e}", mesh.name),
    }
}

fn config(cases: u32) -> ProptestConfig {
    ProptestConfig {
        cases,
        failure_persistence: None,
        ..ProptestConfig::default()
    }
}

proptest! {
    #![proptest_config(config(150))]

    /// Scaling all positions by a common factor changes nothing: masses are area normalized
    /// and `L`, `Ω` are scale invariant.
    #[test]
    fn invariant_under_uniform_scaling(
        kind in 0..N_KINDS,
        mesh_seed in 0..100_000u64,
        feature_seed in 0..100_000u64,
        exponent in -4.0..4.0f64,
    ) {
        let mesh = random_mesh(kind, mesh_seed);
        let features = random_features(&mesh, feature_seed, FeatureMask::UNIQUE);
        let scale = 10f64.powf(exponent);
        let scaled = TestMesh::from_positions(
            "scaled",
            mesh.positions.iter().map(|p| p.map(|x| scale * x)).collect(),
            mesh.mesh.faces().to_vec(),
        );
        let (Some(a), Some(b)) = (solve_ok(&mesh, &features), solve_ok(&scaled, &features)) else {
            // Infeasibility does not depend on the scale either.
            prop_assert!(common::solve_t(&features.problem(&mesh), BoundaryTreatment::Fixed, &SolveOptions::default()).is_err());
            prop_assert!(common::solve_t(&features.problem(&scaled), BoundaryTreatment::Fixed, &SolveOptions::default()).is_err());
            return Ok(());
        };
        let tolerance = 1e-8;
        prop_assert!((a.objective - b.objective).abs() <= tolerance * a.objective.abs().max(1.0));
        prop_assert!(max_abs_diff(&a.measure, &b.measure) <= tolerance,
            "measure differs by {}", max_abs_diff(&a.measure, &b.measure));
        let u_scale = a.log_scale_optimal.iter().fold(1.0f64, |m, u| m.max(u.abs()));
        prop_assert!(max_abs_diff(&a.log_scale_optimal, &b.log_scale_optimal) <= tolerance * u_scale);
    }

    /// Renumbering the vertices renumbers the solution.
    #[test]
    fn invariant_under_vertex_relabeling(
        kind in 0..N_KINDS,
        mesh_seed in 0..100_000u64,
        feature_seed in 0..100_000u64,
        shuffle_seed in any::<u64>(),
    ) {
        let mesh = random_mesh(kind, mesh_seed);
        let features = random_features(&mesh, feature_seed, FeatureMask::UNIQUE);
        let nv = mesh.mesh.n_vertices();
        // new_of[old] = new index.
        let mut new_of: Vec<usize> = (0..nv).collect();
        let mut rng = SplitMix64::new(shuffle_seed);
        for i in (1..nv).rev() {
            new_of.swap(i, (rng.next_u64() % (i as u64 + 1)) as usize);
        }
        let mut positions = vec![[0.0; 3]; nv];
        for old in 0..nv {
            positions[new_of[old]] = mesh.positions[old];
        }
        let faces = mesh.mesh.faces().iter().map(|f| f.map(|v| new_of[v as usize] as u32)).collect();
        let relabeled = TestMesh::from_positions("relabeled", positions, faces);
        let permute = |x: &Option<Vec<f64>>| x.as_ref().map(|x| {
            let mut y = vec![0.0; nv];
            for old in 0..nv {
                y[new_of[old]] = x[old];
            }
            y
        });
        let permuted = Features {
            lambda: features.lambda,
            distortion_weights: permute(&features.distortion_weights),
            cone_costs: permute(&features.cone_costs),
            lower_bounds: permute(&features.lower_bounds),
            upper_bounds: permute(&features.upper_bounds),
            prescribed_cones: features.prescribed_cones.iter()
                .map(|c| Cone { vertex: new_of[c.vertex as usize] as u32, angle: c.angle })
                .collect(),
        };
        let (Some(a), Some(b)) = (solve_ok(&mesh, &features), solve_ok(&relabeled, &permuted)) else {
            return Ok(());
        };
        let tolerance = 1e-8;
        prop_assert!((a.objective - b.objective).abs() <= tolerance * a.objective.abs().max(1.0));
        for old in 0..nv {
            prop_assert!((a.measure[old] - b.measure[new_of[old]]).abs() <= tolerance,
                "vertex {old}: {} vs {}", a.measure[old], b.measure[new_of[old]]);
            prop_assert!((a.log_scale_optimal[old] - b.log_scale_optimal[new_of[old]]).abs()
                <= tolerance * (1.0 + a.log_scale_optimal[old].abs()));
        }
    }

    /// A larger `λ` never buys more (weighted) cone angle and never less distortion.
    #[test]
    fn cone_angle_decreases_and_distortion_increases_with_lambda(
        kind in 0..N_KINDS,
        mesh_seed in 0..100_000u64,
        feature_seed in 0..100_000u64,
        low in -3.5..0.5f64,
        ratio in 0.05..1.5f64,
    ) {
        let mesh = random_mesh(kind, mesh_seed);
        let mut features = random_features(&mesh, feature_seed, FeatureMask::UNIQUE);
        let lambdas = [10f64.powf(low), 10f64.powf(low + ratio)];
        let mut results = Vec::new();
        for lambda in lambdas {
            features.lambda = lambda;
            let Some(solution) = solve_ok(&mesh, &features) else { return Ok(()); };
            // `cone_cost / λ` is the weighted cone norm Σ wR |µ| of the canonical problem.
            results.push((solution.cone_cost / lambda, solution.distortion_optimal));
        }
        let slack = |x: f64| 1e-8 * x.abs().max(1.0);
        prop_assert!(results[1].0 <= results[0].0 + slack(results[0].0),
            "cone norm rose from {} to {} as lambda rose from {} to {}",
            results[0].0, results[1].0, lambdas[0], lambdas[1]);
        prop_assert!(results[1].1 >= results[0].1 - slack(results[0].1),
            "distortion fell from {} to {} as lambda rose from {} to {}",
            results[0].1, results[1].1, lambdas[0], lambdas[1]);
    }

    /// Paper Sec. 4.5: no cone configuration of equal or smaller total angle has less
    /// distortion; and in general the objective of any competitor is no smaller.
    #[test]
    fn no_competitor_of_no_larger_norm_has_less_distortion(
        kind in 0..N_KINDS,
        mesh_seed in 0..100_000u64,
        lambda_exponent in -3.0..0.0f64,
        competitor_seed in any::<u64>(),
    ) {
        let mesh = random_mesh(kind, mesh_seed);
        let features = Features::plain(10f64.powf(lambda_exponent));
        let solution = solve_ok(&mesh, &features).expect("the plain problem is feasible");
        let problem = features.problem(&mesh);
        let interior = mesh.interior();
        let mut rng = SplitMix64::new(competitor_seed);
        let norm: f64 = solution.measure.iter().map(|m| m.abs()).sum();
        for attempt in 0..6 {
            // Competitors: the optimal cones moved to neighbouring vertices (same norm), the
            // optimal measure shrunk (smaller norm; meshes with boundary only), or random.
            let mut nu = vec![0.0; mesh.mesh.n_vertices()];
            match attempt % 3 {
                0 => {
                    for cone in &solution.cones {
                        let v = cone.vertex as usize;
                        let neighbors: Vec<usize> = mesh.mesh.vertex_neighbors(v)
                            .filter(|&w| !mesh.mesh.is_boundary_vertex(w))
                            .collect();
                        let target = if neighbors.is_empty() || rng.next_f64() < 0.5 {
                            v
                        } else {
                            neighbors[(rng.next_u64() as usize) % neighbors.len()]
                        };
                        nu[target] += cone.angle;
                    }
                }
                1 if !mesh.is_closed() => {
                    let shrink = rng.uniform(0.0, 1.0);
                    for v in 0..nu.len() {
                        nu[v] = shrink * solution.measure[v];
                    }
                }
                _ => {
                    for &v in &interior {
                        if rng.next_f64() < 0.3 {
                            nu[v] = rng.uniform(-1.0, 1.0);
                        }
                    }
                }
            }
            if mesh.is_closed() {
                // Gauss–Bonnet: put the missing angle on one vertex.
                let omega = angle_defects(&mesh.mesh, &mesh.geometry).unwrap();
                let missing = omega.iter().sum::<f64>() - nu.iter().sum::<f64>();
                nu[interior[(rng.next_u64() as usize) % interior.len()]] += missing;
            }
            let cones: Vec<Cone> = (0..nu.len())
                .filter(|&v| nu[v] != 0.0)
                .map(|v| Cone { vertex: v as u32, angle: nu[v] })
                .collect();
            let competitor = evaluate_cones(&problem, &cones, Parallelism::Sequential).unwrap();
            prop_assert!(competitor.gauss_bonnet_defect.abs() < 1e-9);
            let competitor_norm = competitor.total_cone_angle;
            let slack = 1e-9 * solution.objective.abs().max(1.0);
            // Optimality of the objective.
            prop_assert!(
                solution.objective
                    <= competitor.distortion + features.lambda * competitor_norm + slack,
                "competitor has a smaller objective: {} < {}",
                competitor.distortion + features.lambda * competitor_norm, solution.objective);
            // The guarantee of Sec. 4.5.
            if competitor_norm <= norm {
                prop_assert!(solution.distortion_optimal <= competitor.distortion + slack,
                    "competitor with norm {competitor_norm} <= {norm} has distortion {} < {}",
                    competitor.distortion, solution.distortion_optimal);
            }
        }
    }

    /// Prescribed cones, forbidden vertices, bounds and the boundary are respected exactly, the
    /// certificate holds, and the certificate in the report is the one `certify` computes.
    #[test]
    fn constraints_are_respected_exactly(
        kind in 0..N_KINDS,
        mesh_seed in 0..100_000u64,
        feature_seed in 0..100_000u64,
    ) {
        let mesh = random_mesh(kind, mesh_seed);
        let features = random_features(&mesh, feature_seed, FeatureMask::ALL);
        let Some(solution) = solve_ok(&mesh, &features) else { return Ok(()); };
        let nv = mesh.mesh.n_vertices();
        let mut prescribed = vec![None; nv];
        for cone in &features.prescribed_cones {
            prescribed[cone.vertex as usize] = Some(cone.angle);
        }
        for v in 0..nv {
            let mu = solution.measure[v];
            if mesh.mesh.is_boundary_vertex(v) {
                prop_assert_eq!(mu, 0.0);
                prop_assert_eq!(solution.log_scale_optimal[v], 0.0);
                prop_assert_eq!(solution.adjoint[v], 0.0);
            } else if let Some(angle) = prescribed[v] {
                prop_assert_eq!(mu, angle, "prescribed cone at {}", v);
            } else if features.cone_costs.as_ref().is_some_and(|c| c[v] == f64::INFINITY) {
                prop_assert_eq!(mu, 0.0, "forbidden vertex {}", v);
            } else {
                let lo = features.lower_bounds.as_ref().map_or(f64::NEG_INFINITY, |b| b[v]);
                let hi = features.upper_bounds.as_ref().map_or(f64::INFINITY, |b| b[v]);
                prop_assert!(lo <= mu && mu <= hi, "vertex {}: {} not in [{}, {}]", v, mu, lo, hi);
            }
        }
        // Without merging, the cones are the measure.
        prop_assert!(!solution.extraction.cones_differ_from_measure || solution.extraction.n_dropped > 0);
        let mut from_cones = vec![0.0; nv];
        for cone in &solution.cones {
            from_cones[cone.vertex as usize] = cone.angle;
        }
        for v in 0..nv {
            if solution.measure[v].abs() > 1e-12 {
                prop_assert_eq!(from_cones[v], solution.measure[v]);
            } else {
                prop_assert_eq!(from_cones[v], 0.0);
            }
        }
        // The report's certificate is reproducible from the returned triple.
        let again = certify(
            &features.problem(&mesh),
            &solution.log_scale_optimal,
            &solution.adjoint,
            &solution.measure,
        ).unwrap();
        prop_assert_eq!(again, solution.report.certificate);
        prop_assert!(again.holds(1e-9));
        prop_assert!(again.bound_violation == 0.0);
        let total = solution.distortion_optimal + solution.cone_cost;
        prop_assert!((again.primal_objective - total).abs() <= 1e-12 * total.abs().max(1.0));
        if mesh.is_closed() {
            let chi = mesh.mesh.euler_characteristic() as f64;
            prop_assert!((solution.measure.iter().sum::<f64>() - 2.0 * PI * chi).abs() < 1e-9);
        }
    }
}

/// A large `λ` on a disk gives no cones at all (paper Fig. 13).
#[test]
fn large_lambda_on_a_disk_gives_no_cones() {
    for seed in 0..20 {
        let mesh = random_mesh(0, seed);
        let solution = solve_ok(&mesh, &Features::plain(100.0)).unwrap();
        assert!(solution.cones.is_empty(), "{}", mesh.name);
        assert!(solution.measure.iter().all(|&m| m == 0.0));
        assert_eq!(solution.cone_cost, 0.0);
        // The adjoint stays inside the threshold.
        assert!(solution.adjoint.iter().all(|p| p.abs() <= 100.0));
    }
}

/// `λ = 0` makes every interior vertex a free site: `µ = Ω`, `u = 0`.
#[test]
fn zero_lambda_reproduces_the_curvature() {
    for kind in 0..N_KINDS {
        let mesh = random_mesh(kind, 7);
        let solution = solve_ok(&mesh, &Features::plain(0.0)).unwrap();
        let omega = angle_defects(&mesh.mesh, &mesh.geometry).unwrap();
        assert!(
            max_abs_diff(&solution.measure, &omega) < 1e-10,
            "{}",
            mesh.name
        );
        assert!(solution.log_scale_optimal.iter().all(|u| u.abs() < 1e-10));
        assert!(solution.objective.abs() < 1e-18);
    }
}

/// A closed convex mesh without a free site gets a cone at every vertex, `µ = Ω` and `u = 0`
/// (paper Sec. 8.1 and Sec. 9); one free site changes that.
#[test]
fn closed_convex_mesh_needs_a_free_site() {
    let base = icosphere(1).unwrap();
    // An ellipsoid: convex, so every angle defect is positive.
    let positions: Vec<[f64; 3]> = base
        .positions
        .iter()
        .map(|p| [1.3 * p[0], 0.8 * p[1], 1.1 * p[2]])
        .collect();
    let mesh = TestMesh::from_positions("ellipsoid", positions, base.faces.clone());
    let omega = angle_defects(&mesh.mesh, &mesh.geometry).unwrap();
    assert!(omega.iter().all(|&o| o > 0.0));
    let nv = mesh.mesh.n_vertices();
    for lambda in [0.01, 0.5, 5.0] {
        let plain = solve_ok(&mesh, &Features::plain(lambda)).unwrap();
        assert_eq!(plain.cones.len(), nv);
        assert!(max_abs_diff(&plain.measure, &omega) < 1e-10);
        assert!(plain.log_scale_optimal.iter().all(|u| u.abs() < 1e-10));
        // All cones are positive, so the cone cost is 4πλ whatever the configuration.
        assert!((plain.cone_cost - 4.0 * PI * lambda).abs() < 1e-10);
        // Every vertex is active: the adjoint is the constant λ.
        assert!(plain.adjoint.iter().all(|p| (p - lambda).abs() < 1e-12));
        assert_eq!(plain.adjoint_shift, (0.0, 0.0));
    }
    // The vertex of greatest curvature becomes a free site.
    let site = (0..nv)
        .max_by(|&a, &b| omega[a].total_cmp(&omega[b]))
        .unwrap();
    let mut costs = vec![1.0; nv];
    costs[site] = 0.0;
    let features = Features {
        cone_costs: Some(costs),
        ..Features::plain(0.5)
    };
    let free = solve_ok(&mesh, &features).unwrap();
    assert!(free.cones.len() < nv / 2, "{} cones", free.cones.len());
    assert_eq!(free.adjoint[site], 0.0);
    // The free site takes most of the curvature, at no cost.
    assert!(free.measure[site] > 2.0 * PI);
    assert!(free.cone_cost < 4.0 * PI * 0.5 - 1.0);
    assert!((free.measure.iter().sum::<f64>() - 4.0 * PI).abs() < 1e-10);
}

/// On a torus with a large `λ` there is no cone, and the adjoint is determined only up to a
/// constant in an interval that the result reports.
#[test]
fn torus_without_cones_reports_the_adjoint_interval() {
    let mesh = random_mesh(4, 3);
    let features = Features::plain(50.0);
    let solution = solve_ok(&mesh, &features).unwrap();
    assert!(solution.cones.is_empty());
    let (lo, hi) = solution.adjoint_shift;
    assert!(lo < 0.0 && hi > 0.0, "{lo} {hi}");
    let max = solution
        .adjoint
        .iter()
        .cloned()
        .fold(f64::NEG_INFINITY, f64::max);
    let min = solution
        .adjoint
        .iter()
        .cloned()
        .fold(f64::INFINITY, f64::min);
    assert!((hi - (50.0 - max)).abs() < 1e-9 && (lo - (-50.0 - min)).abs() < 1e-9);
    // Any constant of the interval gives a triple whose certificate holds.
    let problem = features.problem(&mesh);
    for c in [lo, 0.5 * (lo + hi), hi] {
        let shifted: Vec<f64> = solution.adjoint.iter().map(|p| p + 0.999 * c).collect();
        let certificate = certify(
            &problem,
            &solution.log_scale_optimal,
            &shifted,
            &solution.measure,
        )
        .unwrap();
        assert!(certificate.holds(1e-9), "{certificate:?}");
    }
    // Outside the interval the adjoint is not dual feasible.
    let outside: Vec<f64> = solution.adjoint.iter().map(|p| p + hi + 1e-3).collect();
    let certificate = certify(
        &problem,
        &solution.log_scale_optimal,
        &outside,
        &solution.measure,
    )
    .unwrap();
    assert!(certificate.dual_infeasibility > 5e-4);
    assert!(!certificate.holds(1e-9));
}

/// The scale factor and adjoint that `evaluate_cones` computes from the returned cones are the
/// ones the solver returned, and the adjoint recovers `λ` from the cones.
#[test]
fn evaluate_cones_reproduces_the_solution() {
    for kind in 0..N_KINDS {
        for seed in 0..6 {
            let mesh = random_mesh(kind, 100 + seed);
            let lambda = 0.02;
            let features = Features::plain(lambda);
            let solution = solve_ok(&mesh, &features).unwrap();
            let evaluation = evaluate_cones(
                &features.problem(&mesh),
                &solution.cones,
                Parallelism::Sequential,
            )
            .unwrap();
            let u_scale = solution
                .log_scale_optimal
                .iter()
                .fold(1.0f64, |m, u| m.max(u.abs()));
            assert!(
                max_abs_diff(&evaluation.log_scale, &solution.log_scale_optimal) < 1e-9 * u_scale,
                "{}",
                mesh.name
            );
            assert!((evaluation.distortion - solution.distortion_optimal).abs() < 1e-10);
            // The adjoint agrees up to the constant of a closed mesh.
            let shift = solution.adjoint[0] - evaluation.adjoint[0];
            if !mesh.is_closed() {
                assert!(shift.abs() < 1e-10);
            }
            for v in 0..mesh.mesh.n_vertices() {
                assert!(
                    (evaluation.adjoint[v] + shift - solution.adjoint[v]).abs() < 1e-9,
                    "{} vertex {v}",
                    mesh.name
                );
            }
            // At every cone the adjoint is ±λ, with the sign of the cone.
            for cone in &solution.cones {
                let phi = evaluation.adjoint[cone.vertex as usize] + shift;
                assert!((phi - lambda * cone.angle.signum()).abs() < 1e-9);
            }
        }
    }
}

/// The flat disk has no curvature: nothing to do, whatever the parameters.
#[test]
fn flat_disk_is_already_flat() {
    let base = bumpy_grid_disk(4, 3, 0.0).unwrap();
    let mesh = TestMesh::from_positions("flat", base.positions.clone(), base.faces.clone());
    let solution = solve_ok(&mesh, &Features::plain(0.01)).unwrap();
    assert!(solution.cones.is_empty());
    assert!(solution.log_scale_optimal.iter().all(|u| u.abs() < 1e-12));
    assert!(solution.objective < 1e-20);
}
