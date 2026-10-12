//! The factor-once method against the reference method, and its fallback paths.
//!
//! * `agreement_*`: the same problem solved with `SolveMethod::Reference` and
//!   `SolveMethod::FactorOnce` (forced), on generated meshes of every topology, with random
//!   steering (weights, free and forbidden sites, bounds, prescribed cones) and every boundary
//!   treatment. The tolerances are not tuned: both results carry a certificate whose
//!   `suboptimality_bound` `g` is a proved bound on `J − J*`, and since `J` is `½‖u − u*‖²_W`-
//!   strongly convex in `u` (`W = diag(wE m)`), `‖u_a − u_b‖_W ≤ √(2 g_a) + √(2 g_b)`. The
//!   cone angles follow from `µ = L u + Ω`, so `|Δµ_i| ≤ Σ_j |L_ij| |Δu_j|` with
//!   `|Δu_j| ≤ (√(2 g_a) + √(2 g_b)) / √w_j`. The objectives differ by at most the sum of the two
//!   enclosure widths.
//! * `fallback_*`: every way the factor-once method can give up, forced with tiny limits; the
//!   automatic mode still returns the certified optimum and says what happened.

#![allow(clippy::needless_range_loop)]

mod common;

use common::{FeatureMask, Features, N_KINDS, TestMesh, random_features, random_mesh};
use ddg_core::operators::{cotan_laplacian, lumped_mass_normalized};
use ddg_optimal_cones::{
    BoundaryTreatment, ConeError, FactorOnceOptions, FallbackReason, MethodUsed, Solution,
    SolveMethod, SolveOptions, Status, Unsupported, factor_once_support, solve_with_boundary,
};

const TREATMENTS: [BoundaryTreatment; 4] = [
    BoundaryTreatment::Fixed,
    BoundaryTreatment::Free,
    BoundaryTreatment::Penalized,
    BoundaryTreatment::Convex,
];

fn with_method(method: SolveMethod) -> SolveOptions<'static> {
    SolveOptions {
        method,
        ..SolveOptions::default()
    }
}

/// How the two solutions differ, and what their certificates allow.
#[derive(Debug, Default, Clone, Copy)]
struct Difference {
    objective: f64,
    objective_allowed: f64,
    /// Largest `|Δu_j| √w_j / (√(2 g_a) + √(2 g_b))`; at most 1 if the bound holds.
    u_ratio: f64,
    /// Largest `|Δµ_i| / allowed_i`.
    mu_ratio: f64,
    mu: f64,
    u: f64,
}

fn compare(
    mesh: &TestMesh,
    f: &Features,
    boundary: BoundaryTreatment,
    a: &Solution,
    b: &Solution,
) -> Difference {
    let nv = mesh.mesh.n_vertices();
    let mass = lumped_mass_normalized(&mesh.mesh, &mesh.geometry).unwrap();
    let laplacian = cotan_laplacian(&mesh.mesh, &mesh.geometry).unwrap();
    // A proved bound on J - J* when the certificate has an enclosure; otherwise (a prescribed or
    // forbidden vertex with a rounding-level state residual) the gap plus the terms that
    // account for the residual.
    let bound = |s: &Solution| {
        let c = &s.report.certificate;
        c.suboptimality_bound()
            .unwrap_or(c.duality_gap.abs() + c.gap_uncertainty + c.repair_cost)
            .max(0.0)
            // rounding in evaluating the gap
            + 64.0 * f64::EPSILON * c.primal_objective.abs().max(1.0)
    };
    let slack = (2.0 * bound(a)).sqrt() + (2.0 * bound(b)).sqrt();
    let mut d = Difference {
        objective: (a.objective - b.objective).abs(),
        objective_allowed: bound(a) + bound(b),
        ..Difference::default()
    };
    let mut allowed_u = vec![0.0; nv];
    for v in 0..nv {
        let unknown = boundary != BoundaryTreatment::Fixed || !mesh.mesh.is_boundary_vertex(v);
        if !unknown {
            continue;
        }
        let we = f.distortion_weights.as_ref().map_or(1.0, |w| w[v]);
        let w = we * mass[v];
        allowed_u[v] = slack / w.sqrt();
        let du = (a.log_scale_optimal[v] - b.log_scale_optimal[v]).abs();
        d.u = d.u.max(du);
        d.u_ratio = d.u_ratio.max(du * w.sqrt() / slack.max(1e-300));
    }
    let mut allowed_mu = vec![0.0; nv];
    for (r, c, value) in laplacian.triplets() {
        allowed_mu[r as usize] += value.abs() * allowed_u[c as usize];
    }
    for v in 0..nv {
        let unknown = boundary != BoundaryTreatment::Fixed || !mesh.mesh.is_boundary_vertex(v);
        if !unknown {
            continue;
        }
        let dm = (a.measure[v] - b.measure[v]).abs();
        d.mu = d.mu.max(dm);
        d.mu_ratio = d.mu_ratio.max(dm / allowed_mu[v].max(1e-300));
    }
    d
}

struct Tally {
    compared: usize,
    infeasible: usize,
    unsupported: usize,
    reference_failed: usize,
    worst: Difference,
}

fn run_agreement(
    kinds: std::ops::Range<usize>,
    seeds: std::ops::Range<u64>,
    mask: FeatureMask,
) -> Tally {
    let mut t = Tally {
        compared: 0,
        infeasible: 0,
        unsupported: 0,
        reference_failed: 0,
        worst: Difference::default(),
    };
    for kind in kinds {
        for seed in seeds.clone() {
            let mesh = random_mesh(kind, seed);
            let f = random_features(&mesh, seed * 7 + kind as u64, mask);
            let problem = f.problem(&mesh);
            for boundary in TREATMENTS {
                if boundary != BoundaryTreatment::Fixed && mesh.is_closed() {
                    continue;
                }
                let describe = || format!("{} {boundary:?} {f:?}", mesh.name);
                let reference =
                    solve_with_boundary(&problem, boundary, &with_method(SolveMethod::Reference));
                let reference = match reference {
                    Ok(s) => s,
                    Err(ConeError::Infeasible { .. }) => {
                        t.infeasible += 1;
                        continue;
                    }
                    Err(e) => panic!("{}: {e}", describe()),
                };
                if reference.report.status != Status::Converged {
                    t.reference_failed += 1;
                    continue;
                }
                if factor_once_support(&problem, boundary).unwrap().is_some() {
                    {
                        t.unsupported += 1;
                        let err = solve_with_boundary(
                            &problem,
                            boundary,
                            &with_method(SolveMethod::FactorOnce),
                        );
                        assert!(
                            matches!(err, Err(ConeError::InvalidOptions { .. })),
                            "{}: forced factor-once on an unsupported problem must be an error",
                            describe()
                        );
                        continue;
                    }
                }
                let fast =
                    solve_with_boundary(&problem, boundary, &with_method(SolveMethod::FactorOnce))
                        .unwrap_or_else(|e| panic!("{}: {e}", describe()));
                assert_eq!(fast.report.method.used, MethodUsed::FactorOnce);
                assert_eq!(
                    fast.report.status,
                    Status::Converged,
                    "{}: {:?} {:?}",
                    describe(),
                    fast.report.method.fallback,
                    fast.report.certificate
                );
                let d = compare(&mesh, &f, boundary, &reference, &fast);
                assert!(
                    d.objective <= d.objective_allowed + 1e-12,
                    "{}: objectives differ by {:e}, certificates allow {:e}",
                    describe(),
                    d.objective,
                    d.objective_allowed
                );
                assert!(d.u_ratio <= 1.0 + 1e-6, "{}: u ratio {:?}", describe(), d);
                assert!(d.mu_ratio <= 1.0 + 1e-6, "{}: mu ratio {:?}", describe(), d);
                t.compared += 1;
                t.worst.objective = t.worst.objective.max(d.objective);
                t.worst.u_ratio = t.worst.u_ratio.max(d.u_ratio);
                t.worst.mu_ratio = t.worst.mu_ratio.max(d.mu_ratio);
                t.worst.mu = t.worst.mu.max(d.mu);
                t.worst.u = t.worst.u.max(d.u);
            }
        }
    }
    t
}

#[test]
fn agreement_plain_and_weighted() {
    // Unit and positive weights, free and forbidden sites; no bounds, no prescribed cones.
    let mask = FeatureMask {
        bounds: false,
        prescribed: false,
        ..FeatureMask::UNIQUE
    };
    let t = run_agreement(0..N_KINDS, 0..30, mask);
    eprintln!(
        "plain: compared {} infeasible {} unsupported {} reference-failed {} worst {:?}",
        t.compared, t.infeasible, t.unsupported, t.reference_failed, t.worst
    );
    assert!(t.compared > 100);
}

#[test]
fn agreement_bounds_and_prescribed() {
    let t = run_agreement(0..N_KINDS, 100..140, FeatureMask::UNIQUE);
    eprintln!(
        "steered: compared {} infeasible {} unsupported {} reference-failed {} worst {:?}",
        t.compared, t.infeasible, t.unsupported, t.reference_failed, t.worst
    );
    assert!(t.compared > 100);
}

#[test]
fn zero_distortion_weights_are_left_to_the_reference_method() {
    let mesh = random_mesh(0, 3);
    let mut f = Features::plain(0.3);
    let mut w = vec![1.0; mesh.mesh.n_vertices()];
    let v = mesh.interior()[0];
    w[v] = 0.0;
    f.distortion_weights = Some(w);
    let problem = f.problem(&mesh);
    assert_eq!(
        factor_once_support(&problem, BoundaryTreatment::Fixed).unwrap(),
        Some(Unsupported::ZeroDistortionWeight { vertex: v as u32 })
    );
    let auto =
        solve_with_boundary(&problem, BoundaryTreatment::Fixed, &SolveOptions::default()).unwrap();
    assert_eq!(auto.report.method.used, MethodUsed::Reference);
    assert_eq!(
        auto.report.method.fallback,
        Some(FallbackReason::Unsupported(
            Unsupported::ZeroDistortionWeight { vertex: v as u32 }
        ))
    );
    assert!(auto.report.method.factor_once.is_some());
    assert_eq!(auto.report.status, Status::Converged);
    let forced = solve_with_boundary(
        &problem,
        BoundaryTreatment::Fixed,
        &with_method(SolveMethod::FactorOnce),
    );
    assert!(matches!(forced, Err(ConeError::InvalidOptions { .. })));
}

#[test]
fn lambda_zero_is_left_to_the_reference_method() {
    let mesh = random_mesh(0, 4);
    let problem = Features::plain(0.0).problem(&mesh);
    assert_eq!(
        factor_once_support(&problem, BoundaryTreatment::Fixed).unwrap(),
        Some(Unsupported::ZeroLambda)
    );
}

/// A sphere-like closed mesh with enough vertices that tiny limits bind.
fn big_sphere() -> TestMesh {
    let mut m = random_mesh(3, 5);
    for seed in 5.. {
        m = random_mesh(3, seed);
        if m.mesh.n_vertices() >= 40 {
            break;
        }
    }
    m
}

fn assert_fallback(
    options: SolveOptions<'_>,
    expect: impl Fn(&FallbackReason) -> bool,
    mesh: &TestMesh,
    f: &Features,
) {
    let problem = f.problem(mesh);
    let reference = solve_with_boundary(
        &problem,
        BoundaryTreatment::Fixed,
        &with_method(SolveMethod::Reference),
    )
    .unwrap();
    let auto = solve_with_boundary(&problem, BoundaryTreatment::Fixed, &options).unwrap();
    assert_eq!(auto.report.method.used, MethodUsed::Reference);
    let why = auto.report.method.fallback.clone().expect("a reason");
    assert!(expect(&why), "unexpected reason: {why:?}");
    assert!(auto.report.method.factor_once.is_some());
    assert_eq!(auto.report.status, Status::Converged);
    assert!(auto.report.certificate.holds(1e-9));
    assert_eq!(auto.report.method.requested, SolveMethod::Auto);
    let d = compare(mesh, f, BoundaryTreatment::Fixed, &reference, &auto);
    assert!(d.u_ratio <= 1.0 + 1e-6 && d.mu_ratio <= 1.0 + 1e-6, "{d:?}");
}

#[test]
fn fallback_working_set_limit() {
    let mesh = big_sphere();
    let f = Features::plain(0.01);
    let options = SolveOptions {
        factor_once: FactorOnceOptions {
            max_working_set_fraction: 0.0,
            min_working_set: 3,
            ..FactorOnceOptions::default()
        },
        ..SolveOptions::default()
    };
    assert_fallback(
        options,
        |r| matches!(r, FallbackReason::WorkingSetLimit { limit: 3, .. }),
        &mesh,
        &f,
    );
}

#[test]
fn fallback_memory_limit() {
    let mesh = big_sphere();
    let f = Features::plain(0.01);
    let options = SolveOptions {
        factor_once: FactorOnceOptions {
            memory_limit_bytes: 1000,
            ..FactorOnceOptions::default()
        },
        ..SolveOptions::default()
    };
    assert_fallback(
        options,
        |r| {
            matches!(
                r,
                FallbackReason::MemoryLimit {
                    limit_bytes: 1000,
                    ..
                }
            )
        },
        &mesh,
        &f,
    );
}

#[test]
fn fallback_iteration_limit() {
    let mesh = big_sphere();
    let f = Features::plain(0.01);
    let options = SolveOptions {
        factor_once: FactorOnceOptions {
            max_rounds: 1,
            ..FactorOnceOptions::default()
        },
        ..SolveOptions::default()
    };
    assert_fallback(
        options,
        |r| matches!(r, FallbackReason::IterationLimit { rounds: 1 }),
        &mesh,
        &f,
    );
}

#[test]
fn fallback_certificate_failure() {
    // A loose factor-once tolerance leaves cone residuals of about 1e-6 on a noisy sphere with
    // 640 vertices, which the certificate rejects; the reference method meets it.
    let g = ddg_core::generators::icosphere(3).unwrap();
    let mut p = g.positions.clone();
    ddg_core::generators::perturb(&mut p, [0.03; 3], 11);
    let mesh = TestMesh::from_positions("noisy-sphere", p, g.faces);
    let f = Features::plain(0.03);
    let problem = f.problem(&mesh);
    let loose = SolveOptions {
        factor_once: FactorOnceOptions {
            tolerance: 1e-4,
            ..FactorOnceOptions::default()
        },
        ..SolveOptions::default()
    };
    let auto = solve_with_boundary(&problem, BoundaryTreatment::Fixed, &loose).unwrap();
    assert_eq!(auto.report.method.used, MethodUsed::Reference);
    assert!(
        matches!(
            auto.report.method.fallback,
            Some(FallbackReason::CertificateFailed { .. })
        ),
        "{:?}",
        auto.report.method.fallback
    );
    assert_eq!(auto.report.status, Status::Converged);
    // Forced, the same result is returned as not converged, with the reason.
    let forced = SolveOptions {
        method: SolveMethod::FactorOnce,
        ..loose
    };
    let s = solve_with_boundary(&problem, BoundaryTreatment::Fixed, &forced).unwrap();
    assert_eq!(s.report.method.used, MethodUsed::FactorOnce);
    assert_eq!(s.report.status, Status::NotConverged);
    assert!(matches!(
        s.report.method.fallback,
        Some(FallbackReason::CertificateFailed { .. })
    ));
    assert!(s.report.message.is_some());
}

#[test]
fn forced_give_up_returns_a_non_converged_result() {
    let mesh = big_sphere();
    let problem = Features::plain(0.01).problem(&mesh);
    let options = SolveOptions {
        method: SolveMethod::FactorOnce,
        factor_once: FactorOnceOptions {
            max_working_set_fraction: 0.0,
            min_working_set: 3,
            ..FactorOnceOptions::default()
        },
        ..SolveOptions::default()
    };
    let s = solve_with_boundary(&problem, BoundaryTreatment::Fixed, &options).unwrap();
    assert_eq!(s.report.method.used, MethodUsed::FactorOnce);
    assert!(matches!(
        s.report.method.fallback,
        Some(FallbackReason::WorkingSetLimit { .. })
    ));
    assert_eq!(s.report.status, Status::NotConverged);
    assert!(
        s.report
            .message
            .as_deref()
            .unwrap()
            .contains("did not finish")
    );
}

#[test]
fn reference_forced_does_not_try_the_fast_method() {
    let mesh = big_sphere();
    let problem = Features::plain(0.01).problem(&mesh);
    let s = solve_with_boundary(
        &problem,
        BoundaryTreatment::Fixed,
        &with_method(SolveMethod::Reference),
    )
    .unwrap();
    assert_eq!(s.report.method.used, MethodUsed::Reference);
    assert!(s.report.method.fallback.is_none());
    assert!(s.report.method.factor_once.is_none());
}

#[test]
fn invalid_factor_once_options_are_rejected() {
    let mesh = big_sphere();
    let problem = Features::plain(0.01).problem(&mesh);
    for bad in [
        FactorOnceOptions {
            tolerance: f64::NAN,
            ..FactorOnceOptions::default()
        },
        FactorOnceOptions {
            max_working_set_fraction: -1.0,
            ..FactorOnceOptions::default()
        },
        FactorOnceOptions {
            max_rounds: 0,
            ..FactorOnceOptions::default()
        },
    ] {
        let options = SolveOptions {
            factor_once: bad,
            ..SolveOptions::default()
        };
        assert!(matches!(
            solve_with_boundary(&problem, BoundaryTreatment::Fixed, &options),
            Err(ConeError::InvalidOptions { .. })
        ));
    }
}
