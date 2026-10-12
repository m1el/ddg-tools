//! What the numerical certificate may claim when the residuals are not zero
//! (`lean/proofs`: `certificate_of_computed`, `certificate_enclosure`, `dual_lower_bound_of_clip`):
//! the optimal value of the given problem lies in `optimum_enclosure()`, checked against the
//! generic convex solver on triples that are deliberately wrong.

mod common;

use common::{
    FeatureMask, Features, Reference, TestMesh, random_features, random_mesh,
    reference_solve_with_boundary,
};
use ddg_core::generators::SplitMix64;
use ddg_optimal_cones::{BoundaryTreatment, Solution, SolveOptions, Status, certify_with_boundary};

fn solved(mesh: &TestMesh, f: &Features, b: BoundaryTreatment) -> Option<Solution> {
    let s = common::solve_t(&f.problem(mesh), b, &SolveOptions::default()).ok()?;
    (s.report.status == Status::Converged).then_some(s)
}

/// For converged solves the enclosure exists, is tiny and contains the reference value; for
/// perturbed triples it is wider, contains the reference value whenever it exists, and does not
/// exist when a forbidden or prescribed vertex has a residual.
#[test]
fn the_enclosure_contains_the_optimum() {
    let (mut some, mut none, mut forbidden_none) = (0, 0, 0);
    let (mut exact_some, mut exact_none) = (0usize, 0usize);
    let mut widest = 0.0f64;
    for boundary in [
        BoundaryTreatment::Fixed,
        BoundaryTreatment::Free,
        BoundaryTreatment::Penalized,
    ] {
        for kind in 0..5 {
            for seed in 0..25u64 {
                let mesh = random_mesh(kind, seed);
                let f = random_features(&mesh, 31_000 + seed, FeatureMask::UNIQUE);
                let Reference::Solved { objective, .. } =
                    reference_solve_with_boundary(&mesh, &f, boundary)
                else {
                    continue;
                };
                let Some(s) = solved(&mesh, &f, boundary) else {
                    continue;
                };
                let slack = 1e-8 * objective.abs().max(1.0);
                // 1. the exact solution
                let c = s.report.certificate;
                // With a forbidden or prescribed vertex, or an angle bound, the residual of the
                // exact solution is rounding noise, not zero, and the enclosure is withheld
                // there (the repair would have to change a fixed angle); without them it exists.
                let restricted = f.lower_bounds.is_some()
                    || f.upper_bounds.is_some()
                    || !f.prescribed_cones.is_empty()
                    || f.cone_costs
                        .as_ref()
                        .is_some_and(|c| c.iter().any(|x| x.is_infinite()));
                match c.optimum_enclosure() {
                    Some((lo, hi)) => {
                        assert!(
                            lo <= objective + slack && objective <= hi + slack,
                            "{lo} {objective} {hi}"
                        );
                        assert!(hi - lo < 1e-7 * objective.abs().max(1.0));
                        exact_some += 1;
                    }
                    None => {
                        assert!(restricted, "an unrestricted converged solve: {c:?}");
                        exact_none += 1;
                    }
                }
                // 2. perturbed triples
                let mut rng = SplitMix64::new(seed + 17);
                for amplitude in [1e-4, 1e-2, 1e-1] {
                    let u: Vec<f64> = s
                        .log_scale_optimal
                        .iter()
                        .map(|u| u + amplitude * rng.uniform(-1.0, 1.0))
                        .collect();
                    let mu: Vec<f64> = s
                        .measure
                        .iter()
                        .map(|m| m + amplitude * rng.uniform(-1.0, 1.0))
                        .collect();
                    let phi: Vec<f64> = s
                        .adjoint
                        .iter()
                        .map(|p| p + amplitude * rng.uniform(-1.0, 1.0))
                        .collect();
                    // A perturbed boundary value is no value of the fixed-boundary problem.
                    let (u, mu, phi) = if boundary == BoundaryTreatment::Fixed {
                        let zero = |x: &[f64]| -> Vec<f64> {
                            x.iter()
                                .enumerate()
                                .map(|(v, &x)| {
                                    if mesh.mesh.is_boundary_vertex(v) {
                                        0.0
                                    } else {
                                        x
                                    }
                                })
                                .collect()
                        };
                        (zero(&u), zero(&mu), zero(&phi))
                    } else {
                        (u, mu, phi)
                    };
                    let c =
                        certify_with_boundary(&f.problem(&mesh), boundary, &u, &phi, &mu).unwrap();
                    match c.optimum_enclosure() {
                        Some((lo, hi)) => {
                            some += 1;
                            assert!(
                                lo <= objective + slack && objective <= hi + slack,
                                "{} {boundary:?} amplitude {amplitude}: [{lo}, {hi}] vs {objective}\n{c:?}",
                                mesh.name
                            );
                            widest = widest.max(hi - lo);
                            assert_eq!(
                                c.suboptimality_bound(),
                                Some(c.duality_gap + c.repair_cost)
                            );
                        }
                        None => {
                            none += 1;
                            assert!(c.repair_bound_violation > 0.0 || !c.weights_positive);
                            if c.repair_bound_violation > 0.0 {
                                forbidden_none += 1;
                            }
                            assert_eq!(c.suboptimality_bound(), None);
                        }
                    }
                }
            }
        }
    }
    println!(
        "enclosure present in {some} perturbed triples (widest {widest:.2e}), absent in {none} ({forbidden_none} by an inadmissible repair)"
    );
    println!(
        "exact solutions: enclosure in {exact_some}, withheld in {exact_none} (restricted problems)"
    );
    assert!(some > 100 && forbidden_none > 10 && exact_some > 50);
}

/// A residual at an unbounded vertex can be absorbed into its cone angle, at the price
/// `λ wR |r|`; a residual at a forbidden vertex cannot: no enclosure.
#[test]
fn the_repair_cost_and_when_there_is_no_enclosure() {
    let mesh = random_mesh(3, 2); // a sphere: all vertices unknown
    let nv = mesh.mesh.n_vertices();
    let plain = Features::plain(0.05);
    let s = solved(&mesh, &plain, BoundaryTreatment::Fixed).unwrap();
    let problem = plain.problem(&mesh);
    let mut u = s.log_scale_optimal.clone();
    u[1] += 0.01;
    let c = certify_with_boundary(
        &problem,
        BoundaryTreatment::Fixed,
        &u,
        &s.adjoint,
        &s.measure,
    )
    .unwrap();
    assert!(c.state_residual > 1e-4);
    assert_eq!(c.repair_bound_violation, 0.0);
    assert!(c.optimum_enclosure().is_some());
    // The cost is exactly λ wR Σ |r| over the vertices.
    let l = ddg_core::operators::cotan_laplacian(&mesh.mesh, &mesh.geometry).unwrap();
    let omega = ddg_core::operators::angle_defects(&mesh.mesh, &mesh.geometry).unwrap();
    let mut lu = vec![0.0; nv];
    for (r, cc, v) in l.triplets() {
        lu[r as usize] += v * u[cc as usize];
    }
    let cost: f64 = (0..nv)
        .map(|v| 0.05 * (lu[v] - s.measure[v] + omega[v]).abs())
        .sum();
    assert!(
        (c.repair_cost - cost).abs() <= 1e-12 * cost.max(1.0),
        "{} vs {cost}",
        c.repair_cost
    );

    // The same with vertex 0 forbidden: a bump next to it leaves a residual there.
    let mut f = Features::plain(0.05);
    let mut costs = vec![1.0; nv];
    costs[0] = f64::INFINITY;
    f.cone_costs = Some(costs);
    let s = solved(&mesh, &f, BoundaryTreatment::Fixed).unwrap();
    let next_to_zero = mesh.mesh.vertex_neighbors(0).next().unwrap();
    let mut u = s.log_scale_optimal.clone();
    u[next_to_zero] += 0.01;
    let c = certify_with_boundary(
        &f.problem(&mesh),
        BoundaryTreatment::Fixed,
        &u,
        &s.adjoint,
        &s.measure,
    )
    .unwrap();
    assert!(c.repair_bound_violation > 1e-4, "{c:?}");
    assert_eq!(c.optimum_enclosure(), None);
    assert_eq!(c.suboptimality_bound(), None);
}

/// A vertex without distortion weight: the lower bound is not proved.
#[test]
fn zero_weights_void_the_enclosure() {
    let mesh = random_mesh(3, 4);
    let nv = mesh.mesh.n_vertices();
    let mut f = Features::plain(0.05);
    let mut w = vec![1.0; nv];
    w[3] = 0.0;
    f.distortion_weights = Some(w);
    let s = solved(&mesh, &f, BoundaryTreatment::Fixed).unwrap();
    assert!(!s.report.certificate.weights_positive);
    assert_eq!(s.report.certificate.optimum_enclosure(), None);
}
