//! The public interface: errors that say what is wrong, options, cone extraction and merging,
//! the helpers, and the corner cases of tiny meshes.

// Per-vertex loops that index several arrays read best with the vertex index.
#![allow(clippy::needless_range_loop)]

mod common;

use std::f64::consts::PI;

use common::{
    FeatureMask, Features, TestMesh, max_abs_diff, random_features, random_mesh, summary,
};
use ddg_core::generators::{bumpy_grid_disk, icosphere, torus};
use ddg_core::operators::{angle_defects, cotan_laplacian, lumped_mass_normalized};
use ddg_core::{Geometry, TriMesh};
use ddg_optimal_cones::{
    Cone, ConeError, ConeProblem, FeatureSizeConvention, LinearSolver, Parallelism,
    PrescribedConeIssue, RegularizationMetric, SolveOptions, StageKind, StageStop, Status,
    ThreadPool, UnmergedReason, certify, evaluate_cones, feature_size_weights, soft_threshold,
    solve,
};

fn disk() -> TestMesh {
    let base = bumpy_grid_disk(5, 4, 0.3).unwrap();
    TestMesh::from_positions("disk", base.positions.clone(), base.faces.clone())
}

fn ellipsoid() -> TestMesh {
    let base = icosphere(1).unwrap();
    let positions = base
        .positions
        .iter()
        .map(|p| [1.4 * p[0], 0.7 * p[1], 1.1 * p[2]])
        .collect();
    TestMesh::from_positions("ellipsoid", positions, base.faces.clone())
}

fn solve_default(mesh: &TestMesh, f: &Features) -> ddg_optimal_cones::Solution {
    let solution = solve(&f.problem(mesh), &SolveOptions::default()).unwrap();
    assert_eq!(
        solution.report.status,
        Status::Converged,
        "{}",
        summary(&solution)
    );
    solution
}

/// The reference method, forced: for the tests of its stages, counters and options.
fn solve_reference(mesh: &TestMesh, f: &Features) -> ddg_optimal_cones::Solution {
    let options = SolveOptions {
        method: ddg_optimal_cones::SolveMethod::Reference,
        ..SolveOptions::default()
    };
    let solution = solve(&f.problem(mesh), &options).unwrap();
    assert_eq!(
        solution.report.status,
        Status::Converged,
        "{}",
        summary(&solution)
    );
    solution
}

// ----- input errors -----

#[test]
fn invalid_input_is_reported_not_panicked() {
    let mesh = disk();
    let nv = mesh.mesh.n_vertices();
    let interior = mesh.interior()[0];
    let boundary = (0..nv).find(|&v| mesh.mesh.is_boundary_vertex(v)).unwrap();
    let options = SolveOptions::default();
    let error = |f: Features| solve(&f.problem(&mesh), &options).unwrap_err();

    for lambda in [-1.0, f64::NAN, f64::INFINITY] {
        let e = error(Features::plain(lambda));
        assert!(matches!(e, ConeError::InvalidLambda { .. }), "{e}");
    }
    let e = error(Features {
        distortion_weights: Some(vec![1.0; nv - 1]),
        ..Features::plain(0.5)
    });
    assert_eq!(
        e,
        ConeError::WrongLength {
            what: "distortion_weights",
            expected: nv,
            got: nv - 1
        }
    );
    assert_eq!(
        e.to_string(),
        format!(
            "distortion_weights has {} entries, but the mesh has {nv} vertices",
            nv - 1
        )
    );
    for bad in [-0.5, f64::NAN, f64::INFINITY] {
        let mut w = vec![1.0; nv];
        w[3] = bad;
        let e = error(Features {
            distortion_weights: Some(w),
            ..Features::plain(0.5)
        });
        assert!(
            matches!(e, ConeError::InvalidDistortionWeight { vertex: 3, .. }),
            "{e}"
        );
    }
    for bad in [-0.5, f64::NAN] {
        let mut w = vec![1.0; nv];
        w[4] = bad;
        let e = error(Features {
            cone_costs: Some(w),
            ..Features::plain(0.5)
        });
        assert!(
            matches!(e, ConeError::InvalidConeCost { vertex: 4, .. }),
            "{e}"
        );
    }
    // Bounds must enclose zero.
    for (lower, upper) in [(0.1, 1.0), (-1.0, -0.1), (f64::NAN, 1.0), (-1.0, f64::NAN)] {
        let mut lo = vec![-1.0; nv];
        let mut hi = vec![1.0; nv];
        lo[5] = lower;
        hi[5] = upper;
        let e = error(Features {
            lower_bounds: Some(lo),
            upper_bounds: Some(hi),
            ..Features::plain(0.5)
        });
        assert!(
            matches!(e, ConeError::InvalidBounds { vertex: 5, .. }),
            "{e}"
        );
    }
    // Prescribed cones.
    let prescribed = |cones: Vec<Cone>| {
        error(Features {
            prescribed_cones: cones,
            ..Features::plain(0.5)
        })
    };
    let reason = |e: ConeError| match e {
        ConeError::InvalidPrescribedCone { reason, .. } => reason,
        other => panic!("{other}"),
    };
    let cone = |vertex: usize, angle: f64| Cone {
        vertex: vertex as u32,
        angle,
    };
    assert_eq!(
        reason(prescribed(vec![cone(nv, 1.0)])),
        PrescribedConeIssue::VertexOutOfRange
    );
    assert_eq!(
        reason(prescribed(vec![cone(boundary, 1.0)])),
        PrescribedConeIssue::BoundaryVertex
    );
    assert_eq!(
        reason(prescribed(vec![cone(interior, f64::NAN)])),
        PrescribedConeIssue::NonFiniteAngle
    );
    assert_eq!(
        reason(prescribed(vec![cone(interior, 1.0), cone(interior, 0.5)])),
        PrescribedConeIssue::DuplicateVertex
    );
    let message = prescribed(vec![cone(boundary, 1.0)]).to_string();
    assert!(
        message.contains("boundary") && message.contains(&boundary.to_string()),
        "{message}"
    );
}

#[test]
fn geometry_of_another_mesh_is_rejected() {
    let mesh = disk();
    let other = ellipsoid();
    let problem = ConeProblem::new(&mesh.mesh, &other.geometry, 0.5);
    let e = solve(&problem, &SolveOptions::default()).unwrap_err();
    assert!(matches!(e, ConeError::Geometry(_)), "{e}");
    assert!(certify(&problem, &[], &[], &[]).is_err());
    assert!(evaluate_cones(&problem, &[], Parallelism::Sequential).is_err());
}

#[test]
fn two_single_triangles_have_nothing_to_solve() {
    // Two separate triangles: two components without interior vertices. Defined result: no
    // cones, u = 0, objective 0, a report per component.
    let positions = vec![
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [5.0, 0.0, 0.0],
        [6.0, 0.0, 0.0],
        [5.0, 1.0, 0.0],
    ];
    let mesh = TestMesh::from_positions("two", positions, vec![[0, 1, 2], [3, 4, 5]]);
    for method in [
        ddg_optimal_cones::SolveMethod::Auto,
        ddg_optimal_cones::SolveMethod::Reference,
        ddg_optimal_cones::SolveMethod::FactorOnce,
    ] {
        let options = SolveOptions {
            method,
            ..SolveOptions::default()
        };
        let s = solve(&Features::plain(0.5).problem(&mesh), &options).unwrap();
        assert_eq!(s.report.status, Status::Converged);
        assert!(s.cones.is_empty());
        assert_eq!(s.objective, 0.0);
        assert!(s.log_scale_optimal.iter().all(|&u| u == 0.0));
        assert_eq!(s.report.components.len(), 2);
        assert!(s.report.components.iter().all(|c| c.n_unknowns == 0
            && c.report.stages.is_empty()
            && c.report.status == Status::Converged));
    }
}

#[test]
fn infeasible_closed_problems_say_why() {
    let mesh = ellipsoid();
    let nv = mesh.mesh.n_vertices();
    let options = SolveOptions::default();
    // Cones forbidden everywhere: Σµ = 0, but Gauss–Bonnet needs 4π.
    let forbidden = Features {
        cone_costs: Some(vec![f64::INFINITY; nv]),
        ..Features::plain(0.5)
    };
    match solve(&forbidden.problem(&mesh), &options).unwrap_err() {
        ConeError::Infeasible {
            required,
            min_sum,
            max_sum,
            euler_characteristic,
            n_forbidden,
            ..
        } => {
            assert!((required - 4.0 * PI).abs() < 1e-10);
            assert_eq!((min_sum, max_sum), (0.0, 0.0));
            assert_eq!(euler_characteristic, 2);
            assert_eq!(n_forbidden, nv);
        }
        other => panic!("{other}"),
    }
    // Upper bounds too small.
    let small = Features {
        upper_bounds: Some(vec![0.1; nv]),
        ..Features::plain(0.5)
    };
    let e = solve(&small.problem(&mesh), &options).unwrap_err();
    assert!(matches!(e, ConeError::Infeasible { .. }));
    let text = e.to_string();
    assert!(
        text.contains("2 pi chi") && text.contains("[-inf,"),
        "{text}"
    );
    // Negative cones only.
    let negative = Features {
        upper_bounds: Some(vec![0.0; nv]),
        ..Features::plain(0.5)
    };
    assert!(matches!(
        solve(&negative.problem(&mesh), &options).unwrap_err(),
        ConeError::Infeasible { .. }
    ));
    // Prescribed cones that use up more than 4π while the rest must be positive.
    let too_much = Features {
        lower_bounds: Some(vec![0.0; nv]),
        prescribed_cones: vec![
            Cone {
                vertex: 0,
                angle: 7.0,
            },
            Cone {
                vertex: 1,
                angle: 7.0,
            },
        ],
        ..Features::plain(0.5)
    };
    match solve(&too_much.problem(&mesh), &options).unwrap_err() {
        ConeError::Infeasible {
            n_prescribed,
            prescribed_sum,
            ..
        } => {
            assert_eq!(n_prescribed, 2);
            assert_eq!(prescribed_sum, 14.0);
        }
        other => panic!("{other}"),
    }
    // The same data is feasible on a mesh with boundary, where nothing constrains the sum.
    let disk = disk();
    let n = disk.mesh.n_vertices();
    let forbidden = Features {
        cone_costs: Some(vec![f64::INFINITY; n]),
        ..Features::plain(0.5)
    };
    assert!(solve_default(&disk, &forbidden).cones.is_empty());
    // On a closed mesh the distortion must be measured somewhere.
    let unweighted = Features {
        distortion_weights: Some(vec![0.0; nv]),
        ..Features::plain(0.5)
    };
    assert_eq!(
        solve(&unweighted.problem(&mesh), &options).unwrap_err(),
        ConeError::NoDistortionWeight { component: None }
    );
}

#[test]
fn exactly_feasible_bounds_are_accepted() {
    // A cube: 8 vertices of defect π/2. With 0 ≤ µ ≤ π/2 the only admissible measure is π/2
    // everywhere; the sum of the bounds equals 4π only up to rounding.
    let positions: Vec<[f64; 3]> = (0..8)
        .map(|i| [(i & 1) as f64, ((i >> 1) & 1) as f64, ((i >> 2) & 1) as f64])
        .collect();
    let faces = vec![
        [0, 2, 3],
        [0, 3, 1],
        [4, 5, 7],
        [4, 7, 6],
        [0, 1, 5],
        [0, 5, 4],
        [2, 6, 7],
        [2, 7, 3],
        [0, 4, 6],
        [0, 6, 2],
        [1, 3, 7],
        [1, 7, 5],
    ];
    let mesh = TestMesh::from_positions("cube", positions, faces);
    let features = Features {
        lower_bounds: Some(vec![0.0; 8]),
        upper_bounds: Some(vec![PI / 2.0; 8]),
        ..Features::plain(0.5)
    };
    let solution = solve_default(&mesh, &features);
    assert_eq!(solution.cones.len(), 8);
    assert!(solution.measure.iter().all(|&m| m == PI / 2.0));
    // No vertex is strictly inside its bounds: the adjoint constant is an interval.
    let (lo, hi) = solution.adjoint_shift;
    assert!(lo <= 0.0 && hi == f64::INFINITY, "{lo} {hi}");
    // Slightly smaller bounds are infeasible.
    let features = Features {
        upper_bounds: Some(vec![PI / 2.0 - 1e-6; 8]),
        ..features
    };
    assert!(matches!(
        solve(&features.problem(&mesh), &SolveOptions::default()).unwrap_err(),
        ConeError::Infeasible { .. }
    ));
}

#[test]
fn invalid_options_are_reported() {
    let mesh = disk();
    let problem = Features::plain(0.5).problem(&mesh);
    let bad = |options: SolveOptions<'_>| match solve(&problem, &options).unwrap_err() {
        ConeError::InvalidOptions { reason } => reason,
        other => panic!("{other}"),
    };
    let reason = bad(SolveOptions {
        regularization: vec![1.0, 0.1, 0.1],
        ..SolveOptions::default()
    });
    assert!(
        reason.contains("strictly decreasing") && reason.contains("entry 2"),
        "{reason}"
    );
    assert!(
        bad(SolveOptions {
            regularization: vec![1.0, 0.0],
            ..SolveOptions::default()
        })
        .contains("regularization[1]")
    );
    assert!(
        bad(SolveOptions {
            max_exact_iterations: 0,
            ..SolveOptions::default()
        })
        .contains("max_exact_iterations")
    );
    for field in 0..3 {
        let mut options = SolveOptions::default();
        match field {
            0 => options.newton_tolerance = f64::NAN,
            1 => options.certificate_tolerance = -1.0,
            _ => options.cone_threshold = f64::INFINITY,
        }
        bad(options);
    }
}

// ----- tiny meshes -----

#[test]
fn mesh_without_interior_vertex_has_nothing_to_solve() {
    let mesh = TestMesh::from_positions(
        "triangle",
        vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.5]],
        vec![[0, 1, 2]],
    );
    let solution = solve_default(&mesh, &Features::plain(0.5));
    assert!(solution.cones.is_empty());
    assert_eq!(solution.measure, vec![0.0; 3]);
    assert_eq!(solution.objective, 0.0);
    assert!(solution.report.stages.is_empty());
    assert_eq!(solution.report.linear_solver.ldlt_factorizations, 0);
}

/// One interior vertex: `µ` is the soft threshold of `Ω` in closed form.
#[test]
fn fan_has_the_closed_form_solution() {
    // A cone over a hexagon.
    let mut positions = vec![[0.0, 0.0, 0.6]];
    for k in 0..6 {
        let a = PI / 3.0 * k as f64;
        positions.push([a.cos(), a.sin(), 0.0]);
    }
    let faces = (0..6).map(|k| [0, 1 + k, 1 + (k + 1) % 6]).collect();
    let mesh = TestMesh::from_positions("fan", positions, faces);
    let omega = angle_defects(&mesh.mesh, &mesh.geometry).unwrap()[0];
    let a = cotan_laplacian(&mesh.mesh, &mesh.geometry)
        .unwrap()
        .get(0, 0);
    let m = lumped_mass_normalized(&mesh.mesh, &mesh.geometry).unwrap()[0];
    assert!(omega > 0.5);
    // minimize ½ m ((µ − Ω)/a)² + λ|µ|   ⟹   µ = P_{λ a²/m}(Ω),   φ = −m u / a.
    for lambda in [0.0, 0.01, 0.05, 0.2, 5.0] {
        let solution = solve_default(&mesh, &Features::plain(lambda));
        let mu = soft_threshold(lambda * a * a / m, omega);
        let u = (mu - omega) / a;
        assert!((solution.measure[0] - mu).abs() < 1e-13, "lambda {lambda}");
        assert!((solution.log_scale_optimal[0] - u).abs() < 1e-13);
        assert!((solution.adjoint[0] + m * u / a).abs() < 1e-13);
        assert!((solution.objective - (0.5 * m * u * u + lambda * mu.abs())).abs() < 1e-14);
        if mu != 0.0 {
            assert_eq!(solution.adjoint[0], lambda);
        }
    }
}

// ----- options -----

#[test]
fn every_configuration_finds_the_same_optimum() {
    let pool = ThreadPool::new(2).unwrap();
    for (kind, seed) in [(0, 1), (1, 2), (3, 3), (4, 4), (0, 5), (3, 6)] {
        let mesh = random_mesh(kind, seed);
        let features = random_features(&mesh, 40 + seed, FeatureMask::UNIQUE);
        let problem = features.problem(&mesh);
        let Ok(reference) = solve(&problem, &SolveOptions::default()) else {
            continue;
        };
        assert_eq!(reference.report.status, Status::Converged);
        let variants: Vec<(&str, SolveOptions<'_>)> = vec![
            (
                "mass metric",
                SolveOptions {
                    regularization_metric: RegularizationMetric::Mass,
                    ..SolveOptions::default()
                },
            ),
            (
                "LU",
                SolveOptions {
                    linear_solver: LinearSolver::SaddleLu,
                    ..SolveOptions::default()
                },
            ),
            (
                "LDLt only",
                SolveOptions {
                    linear_solver: LinearSolver::SaddleLdlt,
                    ..SolveOptions::default()
                },
            ),
            (
                "short schedule",
                SolveOptions {
                    regularization: vec![1e-2, 1e-6],
                    ..SolveOptions::default()
                },
            ),
            (
                "no regularization",
                SolveOptions {
                    regularization: vec![],
                    ..SolveOptions::default()
                },
            ),
            (
                "thread pool",
                SolveOptions {
                    parallelism: Parallelism::Pool(&pool),
                    ..SolveOptions::default()
                },
            ),
        ];
        for (name, mut options) in variants {
            // these variants are options of the reference method
            options.method = ddg_optimal_cones::SolveMethod::Reference;
            let solution = solve(&problem, &options).unwrap();
            assert_eq!(
                solution.report.status,
                Status::Converged,
                "{} / {name}: {}",
                mesh.name,
                summary(&solution)
            );
            assert!(
                (solution.objective - reference.objective).abs()
                    < 1e-9 * reference.objective.abs().max(1.0),
                "{} / {name}",
                mesh.name
            );
            assert!(
                max_abs_diff(&solution.measure, &reference.measure) < 1e-8,
                "{} / {name}: {}",
                mesh.name,
                max_abs_diff(&solution.measure, &reference.measure)
            );
            match name {
                "LU" => assert_eq!(solution.report.linear_solver.ldlt_factorizations, 0),
                "LDLt only" => assert_eq!(solution.report.linear_solver.lu_factorizations, 0),
                "no regularization" => {
                    // No Stage I. If the plain Newton iteration from zero did not finish, the
                    // safeguard's regularized problems are reported as such, not as Stage I.
                    let r = &solution.report;
                    assert_eq!(r.regularized_iterations, 0);
                    assert_eq!(r.stages[0].kind, StageKind::Exact);
                    assert!(
                        r.stages
                            .iter()
                            .all(|s| s.safeguard || s.kind == StageKind::Exact)
                    );
                    assert_eq!(r.safeguard_rounds == 0, r.safeguard_iterations == 0);
                    assert_eq!(
                        r.linear_solver.ldlt_factorizations + r.linear_solver.lu_factorizations,
                        r.exact_iterations + r.safeguard_iterations
                    );
                }
                _ => {}
            }
        }
    }
}

#[test]
fn report_counts_what_was_done() {
    let mesh = random_mesh(0, 12);
    let solution = solve_reference(&mesh, &Features::plain(0.01));
    let r = &solution.report;
    let regularized: usize = r
        .stages
        .iter()
        .filter(|s| matches!(s.kind, StageKind::Regularized { .. }))
        .map(|s| s.newton_iterations)
        .sum();
    assert_eq!(
        r.stages.len(),
        SolveOptions::default().regularization.len() + 1,
        "one stage per regularization parameter and the exact one"
    );
    assert_eq!(regularized, r.regularized_iterations);
    assert_eq!(r.stages.last().unwrap().kind, StageKind::Exact);
    assert_eq!(
        r.stages.last().unwrap().stop,
        StageStop::ActiveSetsUnchanged
    );
    assert_eq!(r.safeguard_rounds, 0);
    // One factorization per Newton step, at least one solve each.
    let steps = r.regularized_iterations + r.exact_iterations;
    assert_eq!(
        r.linear_solver.ldlt_factorizations + r.linear_solver.lu_factorizations,
        steps
    );
    assert!(r.linear_solver.solves >= steps);
    assert_eq!(
        r.linear_solver.solves,
        steps + r.linear_solver.refinement_steps
    );
    assert!(r.linear_solver.worst_backward_error < 1e-11);
    assert!(r.timings.total >= r.timings.regularized + r.timings.exact);
    assert!(r.message.is_none());
    assert_eq!((r.bound_correction, r.zeroed_angles), (0.0, 0));
    assert!(r.certificate.holds(1e-9));
    assert!(r.certificate.duality_gap.abs() < 1e-12);
    // Zero distortion weights switch to LU.
    let mut w = vec![1.0; mesh.mesh.n_vertices()];
    w[mesh.interior()[1]] = 0.0;
    let weighted = solve_default(
        &mesh,
        &Features {
            distortion_weights: Some(w),
            ..Features::plain(0.01)
        },
    );
    assert_eq!(weighted.report.linear_solver.ldlt_factorizations, 0);
    assert!(weighted.report.linear_solver.lu_factorizations > 0);
}

/// Without the safeguard and with a tiny iteration limit the result is reported as not
/// converged, with a certificate that says how far off it is, and never as an optimum.
///
/// The problem must be one that one Newton step per stage cannot solve: a disk whose optimum
/// has several cones. (On a convex closed surface, for instance, the optimum `µ = Ω` is reached
/// by one shift of the adjoint constant and one step, and the starved solver rightly reports
/// `Converged`.)
#[test]
fn starved_solver_reports_not_converged() {
    let mesh = random_mesh(0, 21);
    let features = Features::plain(0.001);
    let reference = solve_reference(&mesh, &features);
    assert!(reference.cones.len() >= 2, "{}", summary(&reference));
    assert!(
        reference
            .report
            .stages
            .iter()
            .any(|s| s.newton_iterations > 1),
        "the premise of this test: the default solve needs more than one step in some stage"
    );
    let options = SolveOptions {
        regularization: vec![1e-8],
        max_regularized_iterations: 1,
        max_exact_iterations: 1,
        max_safeguard_rounds: 0,
        method: ddg_optimal_cones::SolveMethod::Reference,
        ..SolveOptions::default()
    };
    let solution = solve(&features.problem(&mesh), &options).unwrap();
    assert_eq!(
        solution.report.status,
        Status::NotConverged,
        "{}",
        summary(&solution)
    );
    assert!(!solution.report.certificate.holds(1e-9));
    assert!(solution.report.message.is_some());
    assert!(
        solution
            .report
            .stages
            .iter()
            .any(|s| s.stop == StageStop::IterationLimit)
    );
    assert_eq!(solution.report.safeguard_rounds, 0);
    assert_eq!(solution.report.safeguard_iterations, 0);
    // With the safeguard the same options converge, to the same optimum.
    let options = SolveOptions {
        max_safeguard_rounds: 200,
        ..options
    };
    let solution = solve(&features.problem(&mesh), &options).unwrap();
    assert_eq!(
        solution.report.status,
        Status::Converged,
        "{}",
        summary(&solution)
    );
    assert!(solution.report.safeguard_rounds > 0);
    assert!(solution.report.safeguard_iterations > 0);
    assert!(solution.report.stages.iter().any(|s| s.safeguard));
    assert!(max_abs_diff(&solution.measure, &reference.measure) < 1e-10);
    assert!((solution.objective - reference.objective).abs() < 1e-12);
}

// ----- cone extraction and merging -----

#[test]
fn threshold_drops_small_cones_and_recomputes_the_scale_factor() {
    let mesh = random_mesh(0, 21);
    let features = Features::plain(0.002);
    let problem = features.problem(&mesh);
    let all = solve_default(&mesh, &features);
    assert!(all.cones.len() >= 2, "{} cones", all.cones.len());
    assert!(!all.extraction.cones_differ_from_measure);
    assert_eq!(all.log_scale, all.log_scale_optimal);
    assert_eq!(all.extraction.distortion_change, 0.0);
    // A threshold between the smallest and the largest cone.
    let mut sizes: Vec<f64> = all.cones.iter().map(|c| c.angle.abs()).collect();
    sizes.sort_by(f64::total_cmp);
    let threshold = 0.5 * (sizes[0] + sizes[1]);
    let options = SolveOptions {
        cone_threshold: threshold,
        ..SolveOptions::default()
    };
    let cut = solve(&problem, &options).unwrap();
    assert_eq!(
        cut.measure, all.measure,
        "the optimal measure does not depend on extraction"
    );
    assert_eq!(cut.cones.len(), all.cones.len() - 1);
    assert_eq!(cut.extraction.n_dropped, 1);
    assert!((cut.extraction.dropped_angle.abs() - sizes[0]).abs() < 1e-15);
    assert!(cut.extraction.cones_differ_from_measure);
    // The returned u belongs to the returned cones, and is no longer the optimal one.
    let evaluation = evaluate_cones(&problem, &cut.cones, Parallelism::Sequential).unwrap();
    assert!(max_abs_diff(&cut.log_scale, &evaluation.log_scale) < 1e-12);
    assert!((cut.distortion - evaluation.distortion).abs() < 1e-14);
    assert!(cut.log_scale != cut.log_scale_optimal);
    assert!(
        (cut.extraction.distortion_change - (cut.distortion - cut.distortion_optimal)).abs()
            < 1e-16
    );
    assert!((cut.extraction.total_angle_change + sizes[0]).abs() < 1e-12);
    assert!(
        (cut.total_cone_angle - cut.cones.iter().map(|c| c.angle.abs()).sum::<f64>()).abs() < 1e-15
    );
}

#[test]
fn clusters_are_merged_at_their_frechet_mean() {
    // A closed convex mesh: every vertex is a cone, all of them form one cluster.
    let mesh = ellipsoid();
    let nv = mesh.mesh.n_vertices();
    let features = Features::plain(0.5);
    let problem = features.problem(&mesh);
    let plain = solve_default(&mesh, &features);
    assert_eq!(plain.cones.len(), nv);
    assert_eq!(plain.extraction.n_clusters, 1);
    assert!(plain.extraction.merged.is_empty() && plain.extraction.unmerged.is_empty());

    // Too large for the default limit: reported, not merged.
    let limited = solve(
        &problem,
        &SolveOptions {
            merge_cone_clusters: true,
            ..SolveOptions::default()
        },
    )
    .unwrap();
    assert_eq!(limited.cones.len(), nv);
    assert_eq!(limited.extraction.unmerged.len(), 1);
    assert_eq!(
        limited.extraction.unmerged[0].reason,
        UnmergedReason::TooLarge
    );
    assert_eq!(limited.extraction.unmerged[0].vertices.len(), nv);
    assert!(!limited.extraction.cones_differ_from_measure);

    // Merged: one cone of 4π.
    let merge_all = SolveOptions {
        merge_cone_clusters: true,
        max_merge_cluster_size: nv,
        ..SolveOptions::default()
    };
    let merged = solve(&problem, &merge_all).unwrap();
    assert_eq!(merged.report.status, Status::Converged);
    assert_eq!(merged.measure, plain.measure);
    assert_eq!(merged.cones.len(), 1);
    assert!(
        (merged.cones[0].angle - 4.0 * PI).abs() < 1e-12,
        "the total angle is kept"
    );
    assert_eq!(merged.extraction.merged.len(), 1);
    assert_eq!(merged.extraction.merged[0].members.len(), nv);
    assert_eq!(merged.extraction.merged[0].cone, merged.cones[0]);
    // The scale factor belongs to the single cone, and the distortion went up from zero.
    let evaluation = evaluate_cones(&problem, &merged.cones, Parallelism::Sequential).unwrap();
    assert!(max_abs_diff(&merged.log_scale, &evaluation.log_scale) < 1e-10);
    assert!(merged.distortion_optimal < 1e-20);
    assert!(merged.distortion > 0.1);
    assert!((merged.extraction.distortion_change - merged.distortion).abs() < 1e-15);
    assert!(merged.extraction.total_angle_change.abs() < 1e-12);

    // The Fréchet mean minimizes Σ |µ_k| dist(v, c_k)² over the members, with distances
    // along the edges of the mesh (every vertex is a member here).
    let target = merged.cones[0].vertex as usize;
    let lengths = mesh.geometry.edge_lengths();
    let mut dist = vec![vec![f64::INFINITY; nv]; nv];
    for v in 0..nv {
        dist[v][v] = 0.0;
    }
    for (e, &[a, b]) in mesh.mesh.edges().iter().enumerate() {
        dist[a as usize][b as usize] = lengths[e];
        dist[b as usize][a as usize] = lengths[e];
    }
    for k in 0..nv {
        for i in 0..nv {
            for j in 0..nv {
                let through = dist[i][k] + dist[k][j];
                if through < dist[i][j] {
                    dist[i][j] = through;
                }
            }
        }
    }
    let variance = |v: usize| -> f64 {
        (0..nv)
            .map(|c| plain.measure[c].abs() * dist[v][c] * dist[v][c])
            .sum()
    };
    let best = (0..nv)
        .min_by(|&a, &b| variance(a).total_cmp(&variance(b)).then(a.cmp(&b)))
        .unwrap();
    assert_eq!(target, best);

    // A merged cone must respect the bounds of its vertex. With µ ≤ 3 everywhere no vertex
    // admits 4π: the cluster stays.
    let bounded = Features {
        upper_bounds: Some(vec![3.0; nv]),
        ..Features::plain(0.5)
    };
    let kept = solve(&bounded.problem(&mesh), &merge_all).unwrap();
    assert_eq!(kept.cones.len(), nv);
    assert_eq!(
        kept.extraction.unmerged[0].reason,
        UnmergedReason::BoundsViolated
    );
    assert!((kept.extraction.unmerged[0].angle - 4.0 * PI).abs() < 1e-12);
    // If only the mean vertex forbids it, the next best member takes the cone.
    let mut upper = vec![f64::INFINITY; nv];
    upper[best] = 3.0;
    let moved = solve(
        &Features {
            upper_bounds: Some(upper),
            ..Features::plain(0.5)
        }
        .problem(&mesh),
        &merge_all,
    )
    .unwrap();
    assert_eq!(moved.cones.len(), 1);
    let second = (0..nv)
        .filter(|&v| v != best)
        .min_by(|&a, &b| variance(a).total_cmp(&variance(b)).then(a.cmp(&b)))
        .unwrap();
    assert_eq!(moved.cones[0].vertex as usize, second);

    // A prescribed cone is never merged or moved; it splits nothing here, the rest merges.
    let fixed = Cone {
        vertex: best as u32,
        angle: 0.25,
    };
    let prescribed = solve(
        &Features {
            prescribed_cones: vec![fixed],
            ..Features::plain(0.5)
        }
        .problem(&mesh),
        &merge_all,
    )
    .unwrap();
    assert_eq!(prescribed.cones.len(), 2);
    assert!(prescribed.cones.contains(&fixed));
    let other = prescribed.cones.iter().find(|c| **c != fixed).unwrap();
    assert!((other.angle - (4.0 * PI - 0.25)).abs() < 1e-10);
}

/// A positive and a negative cone next to each other are a feature of the surface, not one
/// smeared-out cone: they are never merged with each other (summing them would cancel cone
/// angle that the optimum needs).
#[test]
fn cones_of_opposite_sign_are_never_merged() {
    // On a torus with a small λ, positive and negative cones sit next to each other.
    let base = torus(8, 6, 1.0, 0.45).unwrap();
    let mesh = TestMesh::from_positions("torus", base.positions.clone(), base.faces.clone());
    let features = Features::plain(1e-4);
    let nv = mesh.mesh.n_vertices();
    let options = SolveOptions {
        merge_cone_clusters: true,
        max_merge_cluster_size: nv,
        ..SolveOptions::default()
    };
    let solution = solve(&features.problem(&mesh), &options).unwrap();
    assert_eq!(solution.report.status, Status::Converged);
    let raw: Vec<usize> = (0..nv).filter(|&v| solution.measure[v] != 0.0).collect();
    // The premise: some positive cone has a negative cone as a neighbour.
    let mixed_neighbours = raw.iter().any(|&v| {
        mesh.mesh
            .vertex_neighbors(v)
            .any(|w| solution.measure[v] * solution.measure[w] < 0.0)
    });
    assert!(
        mixed_neighbours,
        "the test needs adjacent cones of opposite sign"
    );
    assert!(!solution.extraction.merged.is_empty());
    for cluster in &solution.extraction.merged {
        // One sign per cluster; the merged cone is the sum, at one of the members.
        let positive = cluster.members[0].angle > 0.0;
        assert!(cluster.members.len() >= 2);
        assert!(cluster.members.iter().all(|m| (m.angle > 0.0) == positive));
        let total: f64 = cluster.members.iter().map(|m| m.angle).sum();
        assert!((cluster.cone.angle - total).abs() < 1e-14);
        assert!(
            cluster
                .members
                .iter()
                .any(|m| m.vertex == cluster.cone.vertex)
        );
        assert!(cluster.diameter > 0.0 && cluster.diameter.is_finite());
    }
    // Neither the total angle (zero on a torus) nor the total |angle| changes.
    let sum: f64 = solution.cones.iter().map(|c| c.angle).sum();
    assert!(sum.abs() < 1e-9);
    assert!(solution.extraction.total_angle_change.abs() < 1e-12);
    let positive_before: f64 = solution.measure.iter().filter(|m| **m > 0.0).sum();
    let positive_after: f64 = solution
        .cones
        .iter()
        .map(|c| c.angle)
        .filter(|a| *a > 0.0)
        .sum();
    assert!((positive_before - positive_after).abs() < 1e-12);
    // What merging cost is in the report.
    assert!(solution.distortion >= solution.distortion_optimal);
    assert!(
        (solution.extraction.distortion_change
            - (solution.distortion - solution.distortion_optimal))
            .abs()
            < 1e-16
    );
}

/// The diameter limit, and the count of angles that no flattening can realize.
#[test]
fn wide_clusters_can_be_left_alone_and_impossible_angles_are_counted() {
    // A closed convex mesh: every vertex is a positive cone, all of them one cluster.
    let mesh = ellipsoid();
    let nv = mesh.mesh.n_vertices();
    let problem = Features::plain(0.5).problem(&mesh);
    let merge_all = SolveOptions {
        merge_cone_clusters: true,
        max_merge_cluster_size: nv,
        ..SolveOptions::default()
    };
    let merged = solve(&problem, &merge_all).unwrap();
    assert_eq!(merged.cones.len(), 1);
    let diameter = merged.extraction.merged[0].diameter;
    assert!(
        diameter > 1.0,
        "the cluster spans the ellipsoid: {diameter}"
    );
    // One cone of 4π: not a cone that exists.
    assert_eq!(merged.extraction.n_beyond_two_pi, 1);
    // With a diameter limit below that the cluster stays, and says why.
    let limited = solve(
        &problem,
        &SolveOptions {
            max_merge_cluster_diameter: 0.5 * diameter,
            ..merge_all.clone()
        },
    )
    .unwrap();
    assert_eq!(limited.cones.len(), nv);
    assert_eq!(limited.extraction.unmerged.len(), 1);
    assert_eq!(
        limited.extraction.unmerged[0].reason,
        UnmergedReason::TooWide
    );
    assert_eq!(limited.extraction.unmerged[0].diameter, diameter);
    assert_eq!(limited.extraction.n_beyond_two_pi, 0);
    assert!(!limited.extraction.cones_differ_from_measure);
    // Bounds inside (−2π, 2π) exclude impossible angles, also for merged cones.
    let bounded = Features {
        upper_bounds: Some(vec![6.0; nv]),
        ..Features::plain(0.5)
    };
    let kept = solve(&bounded.problem(&mesh), &merge_all).unwrap();
    assert_eq!(kept.extraction.n_beyond_two_pi, 0);
    assert!(kept.cones.iter().all(|c| c.angle <= 6.0));
    // The option is validated.
    for bad in [-1.0, f64::NAN] {
        let e = solve(
            &problem,
            &SolveOptions {
                max_merge_cluster_diameter: bad,
                ..SolveOptions::default()
            },
        )
        .unwrap_err();
        assert!(e.to_string().contains("max_merge_cluster_diameter"), "{e}");
    }
}

// ----- helpers -----

#[test]
fn feature_size_weights_follow_the_chosen_convention() {
    let mesh = ellipsoid();
    let omega = angle_defects(&mesh.mesh, &mesh.geometry).unwrap();
    let mass = lumped_mass_normalized(&mesh.mesh, &mesh.geometry).unwrap();
    let eps = 0.5;
    let described = feature_size_weights(
        &mesh.mesh,
        &mesh.geometry,
        eps,
        FeatureSizeConvention::AsDescribed,
    )
    .unwrap();
    let printed = feature_size_weights(
        &mesh.mesh,
        &mesh.geometry,
        eps,
        FeatureSizeConvention::AsPrinted,
    )
    .unwrap();
    for v in 0..mesh.mesh.n_vertices() {
        let density = omega[v].abs() / mass[v] + eps;
        assert_eq!(described.feature_size[v], density);
        assert_eq!(described.distortion_weights[v], density);
        assert_eq!(described.cone_costs[v], 1.0 / density);
        assert_eq!(printed.feature_size[v], 1.0 / density);
        assert_eq!(printed.distortion_weights[v], 1.0 / density);
        assert_eq!(printed.cone_costs[v], 1.0 / (1.0 / density));
    }
    // "Large in highly curved regions" is the described convention.
    let curved = (0..omega.len())
        .max_by(|&a, &b| (omega[a] / mass[a]).total_cmp(&(omega[b] / mass[b])))
        .unwrap();
    let flat = (0..omega.len())
        .min_by(|&a, &b| (omega[a] / mass[a]).total_cmp(&(omega[b] / mass[b])))
        .unwrap();
    assert!(described.feature_size[curved] > described.feature_size[flat]);
    assert!(printed.feature_size[curved] < printed.feature_size[flat]);
    // Independent of the length unit.
    let scaled = TestMesh::from_positions(
        "scaled",
        mesh.positions.iter().map(|p| p.map(|x| 8.0 * x)).collect(),
        mesh.mesh.faces().to_vec(),
    );
    let again = feature_size_weights(
        &scaled.mesh,
        &scaled.geometry,
        eps,
        FeatureSizeConvention::AsDescribed,
    )
    .unwrap();
    assert!(max_abs_diff(&again.feature_size, &described.feature_size) < 1e-9);
    // The weights are usable as they are.
    let features = Features {
        distortion_weights: Some(described.distortion_weights.clone()),
        cone_costs: Some(described.cone_costs.clone()),
        ..Features::plain(0.5)
    };
    solve_default(&mesh, &features);
    // Bad arguments.
    for eps in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        let e = feature_size_weights(
            &mesh.mesh,
            &mesh.geometry,
            eps,
            FeatureSizeConvention::AsPrinted,
        )
        .unwrap_err();
        assert!(matches!(e, ConeError::InvalidArgument { .. }), "{e}");
    }
    let other = disk();
    assert!(matches!(
        feature_size_weights(
            &mesh.mesh,
            &other.geometry,
            1.0,
            FeatureSizeConvention::AsPrinted
        ),
        Err(ConeError::Geometry(_))
    ));
}

#[test]
fn evaluate_cones_validates_its_cones() {
    let mesh = disk();
    let problem = Features::plain(0.5).problem(&mesh);
    let nv = mesh.mesh.n_vertices();
    let interior = mesh.interior()[0] as u32;
    let boundary = (0..nv).find(|&v| mesh.mesh.is_boundary_vertex(v)).unwrap() as u32;
    let issue = |cones: &[Cone]| match evaluate_cones(&problem, cones, Parallelism::Sequential) {
        Err(ConeError::InvalidCone { reason, .. }) => reason,
        other => panic!("{other:?}"),
    };
    let cone = |vertex: u32, angle: f64| Cone { vertex, angle };
    assert_eq!(
        issue(&[cone(nv as u32, 1.0)]),
        PrescribedConeIssue::VertexOutOfRange
    );
    assert_eq!(
        issue(&[cone(boundary, 1.0)]),
        PrescribedConeIssue::BoundaryVertex
    );
    assert_eq!(
        issue(&[cone(interior, f64::INFINITY)]),
        PrescribedConeIssue::NonFiniteAngle
    );
    assert_eq!(
        issue(&[cone(interior, 1.0), cone(interior, 1.0)]),
        PrescribedConeIssue::DuplicateVertex
    );
    // No cones at all: the cone-free flattening with u = 0 on the boundary.
    let none = evaluate_cones(&problem, &[], Parallelism::Sequential).unwrap();
    let no_cones = solve_default(&mesh, &Features::plain(1e6));
    assert!(max_abs_diff(&none.log_scale, &no_cones.log_scale_optimal) < 1e-12);
    assert_eq!(none.total_cone_angle, 0.0);
    assert_eq!(none.gauss_bonnet_defect, 0.0);
    // On a closed mesh a cone sum other than 2πχ is reported, not hidden.
    let closed = ellipsoid();
    let problem = Features::plain(0.5).problem(&closed);
    let one = evaluate_cones(
        &problem,
        &[cone(0, 4.0 * PI + 0.5)],
        Parallelism::Sequential,
    )
    .unwrap();
    assert!((one.gauss_bonnet_defect - 0.5).abs() < 1e-10);
    let exact = evaluate_cones(&problem, &[cone(0, 4.0 * PI)], Parallelism::Sequential).unwrap();
    assert!(exact.gauss_bonnet_defect.abs() < 1e-10);
    // The optimal constant: the weighted mean of u vanishes.
    let mass = lumped_mass_normalized(&closed.mesh, &closed.geometry).unwrap();
    let mean: f64 = exact.log_scale.iter().zip(&mass).map(|(u, m)| u * m).sum();
    assert!(mean.abs() < 1e-12);
}

#[test]
fn certify_measures_what_is_wrong_with_a_triple() {
    let mesh = random_mesh(0, 21);
    // The test needs an optimum with a cone. The cone-free configuration is optimal exactly
    // when its adjoint stays within ±λ (the optimality system is necessary and sufficient on a
    // mesh with boundary), so half of its largest adjoint value is a λ that forces a cone.
    let cone_free = evaluate_cones(
        &Features::plain(1.0).problem(&mesh),
        &[],
        Parallelism::Sequential,
    )
    .unwrap();
    let largest = cone_free.adjoint.iter().fold(0.0f64, |m, p| m.max(p.abs()));
    assert!(
        largest > 1e-4,
        "the mesh is too flat for this test: {largest}"
    );
    let features = Features::plain(0.5 * largest);
    let problem = features.problem(&mesh);
    let solution = solve_default(&mesh, &features);
    assert!(!solution.cones.is_empty());
    let nv = mesh.mesh.n_vertices();
    let good = certify(
        &problem,
        &solution.log_scale_optimal,
        &solution.adjoint,
        &solution.measure,
    )
    .unwrap();
    assert!(good.holds(1e-10));
    assert!(good.dual_objective <= good.primal_objective + 1e-12);
    // A perturbed u breaks the state equation and the gap, nothing else.
    let v = mesh.interior()[0];
    let mut u = solution.log_scale_optimal.clone();
    u[v] += 1e-3;
    let c = certify(&problem, &u, &solution.adjoint, &solution.measure).unwrap();
    assert!(c.state_residual > 1e-4 && c.cone_residual == good.cone_residual);
    assert!(!c.holds(1e-9));
    // An adjoint beyond the threshold is not dual feasible.
    let mut phi = solution.adjoint.clone();
    let cone = solution.cones[0].vertex as usize;
    phi[cone] *= 1.5;
    let c = certify(
        &problem,
        &solution.log_scale_optimal,
        &phi,
        &solution.measure,
    )
    .unwrap();
    assert!((c.dual_infeasibility - 0.5 * features.lambda).abs() < 1e-12);
    assert!(c.adjoint_residual > 1e-6);
    // A suboptimal but feasible pair has a positive gap that bounds its suboptimality.
    let none = evaluate_cones(&problem, &[], Parallelism::Sequential).unwrap();
    let c = certify(&problem, &none.log_scale, &none.adjoint, &vec![0.0; nv]).unwrap();
    assert!(c.state_residual < 1e-12 && c.adjoint_residual < 1e-12);
    let suboptimality = c.primal_objective - solution.objective;
    assert!(suboptimality > 0.0);
    assert!(
        c.duality_gap >= suboptimality - 1e-12,
        "{} < {suboptimality}",
        c.duality_gap
    );
    assert!((c.duality_gap - c.complementarity_gap).abs() < 1e-10 * c.primal_objective.max(1.0));
    // Values on the boundary and wrong lengths.
    let boundary = (0..nv).find(|&b| mesh.mesh.is_boundary_vertex(b)).unwrap();
    let mut mu = solution.measure.clone();
    mu[boundary] = 0.25;
    let c = certify(
        &problem,
        &solution.log_scale_optimal,
        &solution.adjoint,
        &mu,
    )
    .unwrap();
    assert_eq!(c.bound_violation, 0.25);
    assert!(matches!(
        certify(&problem, &u[1..], &solution.adjoint, &solution.measure),
        Err(ConeError::WrongLength {
            what: "log_scale",
            ..
        })
    ));
}

/// The identity behind the certificate (`duality-gap-inexact`), for triples that are **not**
/// solutions: with `r = L u − µ + Ω`, `φ̂` the clipped adjoint and `w = wE m`,
///
/// ```text
/// duality_gap − complementarity_gap = Σ_{w_i=0} (L φ̂)_i u_i − Σ φ̂_i r_i
/// ```
///
/// and `gap_uncertainty` is the sum of the magnitudes of these terms. Everything on the right
/// is recomputed here from the assembled Laplacian of `ddg-core`.
#[test]
fn certificate_gap_identity_holds_for_arbitrary_triples() {
    use ddg_core::generators::SplitMix64;
    let mut rng = SplitMix64::new(2018);
    let mut with_zero_weight = 0;
    for case in 0..60u64 {
        let mesh = random_mesh((case % 5) as usize, 100 + case);
        let features = random_features(&mesh, 300 + case, FeatureMask::ALL);
        let problem = features.problem(&mesh);
        let nv = mesh.mesh.n_vertices();
        let interior = mesh.interior();
        if interior.is_empty() {
            continue;
        }
        let mut u = vec![0.0; nv];
        let mut phi = vec![0.0; nv];
        let mut mu = vec![0.0; nv];
        for &v in &interior {
            u[v] = rng.uniform(-1.0, 1.0);
            phi[v] = rng.uniform(-1.5, 1.5) * features.lambda.max(0.05);
            mu[v] = if rng.next_f64() < 0.5 {
                rng.uniform(-2.0, 2.0)
            } else {
                0.0
            };
        }
        let Ok(c) = certify(&problem, &u, &phi, &mu) else {
            // An infeasible closed problem: nothing to certify.
            continue;
        };
        // The right-hand side, from the problem data.
        let laplacian = cotan_laplacian(&mesh.mesh, &mesh.geometry).unwrap();
        let omega = angle_defects(&mesh.mesh, &mesh.geometry).unwrap();
        let mass = lumped_mass_normalized(&mesh.mesh, &mesh.geometry).unwrap();
        let mut fixed = vec![false; nv];
        for cone in &features.prescribed_cones {
            fixed[cone.vertex as usize] = true;
        }
        let phi_hat: Vec<f64> = (0..nv)
            .map(|v| {
                let cost = features.cone_costs.as_ref().map_or(1.0, |c| c[v]);
                if fixed[v] || cost == f64::INFINITY || mesh.mesh.is_boundary_vertex(v) {
                    return phi[v];
                }
                let t = features.lambda * cost;
                let lo = features
                    .lower_bounds
                    .as_ref()
                    .map_or(f64::NEG_INFINITY, |b| b[v]);
                let hi = features
                    .upper_bounds
                    .as_ref()
                    .map_or(f64::INFINITY, |b| b[v]);
                let mut p = phi[v];
                if hi == f64::INFINITY {
                    p = p.min(t);
                }
                if lo == f64::NEG_INFINITY {
                    p = p.max(-t);
                }
                p
            })
            .collect();
        let mut l_u = vec![0.0; nv];
        let mut l_phi_hat = vec![0.0; nv];
        laplacian.mul_vec(&u, &mut l_u).unwrap();
        laplacian.mul_vec(&phi_hat, &mut l_phi_hat).unwrap();
        let (mut extra, mut magnitude) = (0.0, 0.0);
        let mut any_zero = false;
        for &v in &interior {
            let r = l_u[v] - mu[v] + omega[v];
            extra -= phi_hat[v] * r;
            magnitude += (phi_hat[v] * r).abs();
            let w = features.distortion_weights.as_ref().map_or(1.0, |w| w[v]) * mass[v];
            if w == 0.0 {
                any_zero = true;
                extra += l_phi_hat[v] * u[v];
                magnitude += (l_phi_hat[v] * u[v]).abs();
            }
        }
        with_zero_weight += usize::from(any_zero);
        let scale = 1.0 + c.primal_objective.abs() + c.dual_objective.abs() + magnitude;
        assert!(
            (c.duality_gap - c.complementarity_gap - extra).abs() <= 1e-11 * scale,
            "case {case}: gap {} = complementarity {} + {extra}?",
            c.duality_gap,
            c.complementarity_gap
        );
        assert!(
            (c.gap_uncertainty - magnitude).abs() <= 1e-11 * scale,
            "case {case}"
        );
        assert!((c.duality_gap - c.complementarity_gap).abs() <= c.gap_uncertainty + 1e-11 * scale);
        // A random triple is far from feasible: the certificate must not hold, and the gap
        // alone would be misleading (it can have either sign).
        assert!(c.state_residual > 1e-3 && !c.holds(1e-9), "case {case}");
    }
    assert!(
        with_zero_weight >= 3,
        "zero distortion weights were exercised {with_zero_weight} times"
    );
}

/// The certificate in the report belongs to the returned optimum `(log_scale_optimal, adjoint,
/// measure)` exactly as returned (after the rounding-level clean-up of `µ`), and never to the
/// cones that extraction made out of it.
#[test]
fn certificate_is_of_the_returned_optimum_not_of_the_extracted_cones() {
    // A closed mesh, where an active vertex absorbs the rounding of Σ Ω and may be zeroed.
    for (mesh, lambda) in [
        (ellipsoid(), 0.5),
        (random_mesh(4, 9), 0.02),
        (random_mesh(0, 21), 0.001),
    ] {
        let features = Features::plain(lambda);
        let problem = features.problem(&mesh);
        let nv = mesh.mesh.n_vertices();
        let plain = solve_default(&mesh, &features);
        let merged = solve(
            &problem,
            &SolveOptions {
                merge_cone_clusters: true,
                max_merge_cluster_size: nv,
                cone_threshold: 0.05,
                ..SolveOptions::default()
            },
        )
        .unwrap();
        // Extraction changes the cones and their scale factor, and nothing of the optimum.
        assert_eq!(merged.measure, plain.measure);
        assert_eq!(merged.log_scale_optimal, plain.log_scale_optimal);
        assert_eq!(merged.adjoint, plain.adjoint);
        assert_eq!(merged.objective, plain.objective);
        assert_eq!(merged.report.certificate, plain.report.certificate);
        for solution in [&plain, &merged] {
            let again = certify(
                &problem,
                &solution.log_scale_optimal,
                &solution.adjoint,
                &solution.measure,
            )
            .unwrap();
            assert_eq!(again, solution.report.certificate, "{}", mesh.name);
            assert!(again.holds(1e-9));
            assert_eq!(again.primal_objective, {
                let c = &solution.report.certificate;
                c.primal_objective
            });
        }
        if merged.extraction.cones_differ_from_measure {
            // The extracted cones are deliberately not the optimum: certified as a triple of
            // their own (with the optimum's adjoint) they do not pass.
            let mut from_cones = vec![0.0; nv];
            for cone in &merged.cones {
                from_cones[cone.vertex as usize] = cone.angle;
            }
            let c = certify(&problem, &merged.log_scale, &merged.adjoint, &from_cones).unwrap();
            assert!(!c.holds(1e-9), "{}", mesh.name);
            assert!(merged.distortion >= merged.distortion_optimal);
        }
    }
}

#[test]
fn intrinsic_geometry_is_enough() {
    // The solver never sees positions: a geometry built from edge lengths gives the same result.
    let mesh = random_mesh(0, 8);
    let from_lengths =
        Geometry::from_edge_lengths(&mesh.mesh, mesh.geometry.edge_lengths().to_vec()).unwrap();
    let a = solve_default(&mesh, &Features::plain(0.01));
    let problem = ConeProblem::new(&mesh.mesh, &from_lengths, 0.01);
    let b = solve(&problem, &SolveOptions::default()).unwrap();
    assert_eq!(a.measure, b.measure);
    let _ = TriMesh::new(3, vec![[0, 1, 2]]).unwrap();
}
