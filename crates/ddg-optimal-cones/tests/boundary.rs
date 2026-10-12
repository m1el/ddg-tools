//! The boundary treatments (`solve_with_boundary`): optimized boundary values, penalized and
//! convex boundary curvature.
//!
//! * Agreement with the **double** of the mesh: the mesh glued to a mirror copy along its
//!   boundary is a closed mesh, on which the existing closed-mesh path solves the problem with
//!   the seam vertices as free (or penalized, or sign-constrained) cone sites. This route does
//!   not use any of the boundary code. (The comparison with the generic convex solver is in
//!   `convex_reference.rs`.)
//! * Orderings that hold for the problems themselves: the optimized boundary is never worse than
//!   the fixed one, never worse than a penalized or a convex one.
//! * Reporting: `u` and `ν` on the boundary, the certificate, extraction.

#![allow(clippy::needless_range_loop)]

mod common;

use common::{
    FeatureMask, Features, N_KINDS, TestMesh, double_mesh, doubled_features, max_abs_diff,
    random_features, random_mesh, summary,
};
use ddg_core::operators::{boundary_geodesic_curvatures, cotan_laplacian};
use ddg_optimal_cones::{
    BoundaryTreatment, ConeError, SolveOptions, Status, certify, certify_with_boundary, solve,
    solve_with_boundary,
};

const TREATMENTS: [BoundaryTreatment; 3] = [
    BoundaryTreatment::Free,
    BoundaryTreatment::Penalized,
    BoundaryTreatment::Convex,
];

fn l_times(mesh: &TestMesh, x: &[f64]) -> Vec<f64> {
    let l = cotan_laplacian(&mesh.mesh, &mesh.geometry).unwrap();
    let mut out = vec![0.0; x.len()];
    for (r, c, value) in l.triplets() {
        out[r as usize] += value * x[c as usize];
    }
    out
}

/// What the double gives for the same problem, against what the boundary code gives.
#[derive(Debug, Default, Clone, Copy)]
struct Agreement {
    objective: f64,
    measure: f64,
    log_scale: f64,
    boundary_curvature: f64,
}

/// Solves `f` on `mesh` with the boundary treatment and on its double with the corresponding
/// closed-mesh problem; `None` if the double is not a mesh, or the problem is infeasible for both (agreement of verdicts is
/// asserted), an error text on disagreement.
fn compare_with_double(
    mesh: &TestMesh,
    f: &Features,
    boundary: BoundaryTreatment,
) -> Result<Option<Agreement>, String> {
    let describe = || format!("{} with {f:.3?} and {boundary:?}", mesh.name);
    let Some((doubled, geometry, copy)) = double_mesh(&mesh.mesh, &mesh.geometry) else {
        return Ok(None);
    };
    let double = TestMesh {
        name: format!("double of {}", mesh.name),
        positions: Vec::new(),
        mesh: doubled,
        geometry,
    };
    let g = doubled_features(&mesh.mesh, &copy, f, boundary);
    let options = SolveOptions::default();
    let ours = solve_with_boundary(&f.problem(mesh), boundary, &options);
    let theirs = solve(&g.problem(&double), &options);
    let (ours, theirs) = match (ours, theirs) {
        (Err(ConeError::Infeasible { .. }), Err(ConeError::Infeasible { .. })) => return Ok(None),
        (Ok(a), Ok(b)) => (a, b),
        (a, b) => {
            return Err(format!(
                "verdicts differ: boundary code {:?}, double {:?}: {}",
                a.map(|s| s.report.status),
                b.map(|s| s.report.status),
                describe()
            ));
        }
    };
    if ours.report.status != Status::Converged || theirs.report.status != Status::Converged {
        return Err(format!(
            "not converged: ours {}\n  double {}\n  {}",
            summary(&ours),
            summary(&theirs),
            describe()
        ));
    }
    let nv = mesh.mesh.n_vertices();
    let scale = ours.objective.abs().max(1.0);
    let u_scale = ours
        .log_scale_optimal
        .iter()
        .fold(1.0f64, |m, u| m.max(u.abs()));
    let mut a = Agreement {
        objective: (ours.objective - theirs.objective).abs() / scale,
        ..Agreement::default()
    };
    for v in 0..nv {
        a.log_scale = a
            .log_scale
            .max((ours.log_scale_optimal[v] - theirs.log_scale_optimal[v]).abs() / u_scale);
        if mesh.mesh.is_boundary_vertex(v) {
            // A seam cone angle is twice the curvature of the boundary.
            a.boundary_curvature = a
                .boundary_curvature
                .max((ours.boundary_curvature[v] - 0.5 * theirs.measure[v]).abs());
        } else {
            a.measure = a.measure.max((ours.measure[v] - theirs.measure[v]).abs());
            a.measure = a
                .measure
                .max((ours.measure[v] - theirs.measure[copy[v] as usize]).abs());
        }
    }
    if a.objective > 1e-8 || a.measure > 1e-6 || a.log_scale > 1e-6 || a.boundary_curvature > 1e-6 {
        return Err(format!("{a:?}: {}", describe()));
    }
    Ok(Some(a))
}

#[test]
fn boundary_treatments_agree_with_the_double_of_the_mesh() {
    println!(
        "{:<11} {:>5} {:>5} {:>10} {:>10} {:>10} {:>10}",
        "treatment", "cases", "feas", "objective", "measure", "u (rel)", "nu"
    );
    let mut failures = Vec::new();
    for boundary in TREATMENTS {
        let (mut cases, mut feasible) = (0, 0);
        let mut worst = Agreement::default();
        for kind in 0..3 {
            for seed in 0..40u64 {
                let mesh = random_mesh(kind, seed);
                let features =
                    random_features(&mesh, 7000 + seed + 100 * kind as u64, FeatureMask::UNIQUE);
                cases += 1;
                match compare_with_double(&mesh, &features, boundary) {
                    Ok(Some(a)) => {
                        feasible += 1;
                        worst.objective = worst.objective.max(a.objective);
                        worst.measure = worst.measure.max(a.measure);
                        worst.log_scale = worst.log_scale.max(a.log_scale);
                        worst.boundary_curvature =
                            worst.boundary_curvature.max(a.boundary_curvature);
                    }
                    Ok(None) => {}
                    Err(e) => failures.push(e),
                }
            }
        }
        println!(
            "{:<11} {:>5} {:>5} {:>10.2e} {:>10.2e} {:>10.2e} {:>10.2e}",
            format!("{boundary:?}"),
            cases,
            feasible,
            worst.objective,
            worst.measure,
            worst.log_scale,
            worst.boundary_curvature
        );
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

/// `Free` relaxes `Fixed` (a fixed-boundary solution is feasible with `ν = (L u)_B + k_B` at no
/// cost) and `Penalized` and `Convex` (which add a cost, or a constraint, to the same
/// variables), so its optimal objective is the smallest of the four.
#[test]
fn optimized_boundary_is_never_worse() {
    let mut strict = 0;
    let mut cases = 0;
    for kind in 0..3 {
        for seed in 0..40u64 {
            let mesh = random_mesh(kind, seed);
            let mut f = random_features(&mesh, 9000 + seed, FeatureMask::UNIQUE);
            // The bounds of boundary vertices play no role; bounds of interior vertices would
            // make the comparison of the problems unfair only through infeasibility.
            f.prescribed_cones.clear();
            let problem = f.problem(&mesh);
            let options = SolveOptions::default();
            let Ok(fixed) = solve(&problem, &options) else {
                continue;
            };
            let Ok(free) = solve_with_boundary(&problem, BoundaryTreatment::Free, &options) else {
                continue;
            };
            assert_eq!(fixed.report.status, Status::Converged);
            assert_eq!(free.report.status, Status::Converged, "{}", summary(&free));
            cases += 1;
            let scale = fixed.objective.abs().max(1.0);
            assert!(
                free.objective <= fixed.objective + 1e-9 * scale,
                "{}: free {} fixed {}",
                mesh.name,
                free.objective,
                fixed.objective
            );
            if free.objective < fixed.objective - 1e-6 * scale {
                strict += 1;
            }
            for other in [BoundaryTreatment::Penalized, BoundaryTreatment::Convex] {
                if let Ok(s) = solve_with_boundary(&problem, other, &options) {
                    assert_eq!(s.report.status, Status::Converged);
                    assert!(
                        free.objective <= s.objective + 1e-9 * scale,
                        "{} {other:?}: free {} other {}",
                        mesh.name,
                        free.objective,
                        s.objective
                    );
                }
            }
        }
    }
    println!("{cases} cases, free strictly better than fixed in {strict}");
    assert!(strict > cases / 2);
}

fn solved(
    mesh: &TestMesh,
    f: &Features,
    boundary: BoundaryTreatment,
) -> ddg_optimal_cones::Solution {
    let s = solve_with_boundary(&f.problem(mesh), boundary, &SolveOptions::default())
        .unwrap_or_else(|e| panic!("{e}"));
    assert_eq!(s.report.status, Status::Converged, "{}", summary(&s));
    s
}

/// What is reported: `ν = (L u)_B + k_B` at the boundary, `Σ µ + Σ ν = 2πχ`, the adjoint
/// vanishes at free boundary sites, `measure` and `boundary_curvature` agree, `cones` are
/// interior, and the certificate of the returned triple is reproduced bit for bit.
#[test]
fn reported_boundary_quantities() {
    for kind in 0..3 {
        for seed in 0..12u64 {
            let mesh = random_mesh(kind, seed);
            let f = Features::plain(0.02);
            let k = boundary_geodesic_curvatures(&mesh.mesh, &mesh.geometry).unwrap();
            let chi = mesh.mesh.euler_characteristic() as f64;
            for boundary in TREATMENTS.into_iter().chain([BoundaryTreatment::Fixed]) {
                let s = solved(&mesh, &f, boundary);
                let l_u = l_times(&mesh, &s.log_scale_optimal);
                let mut total = 0.0;
                for v in 0..mesh.mesh.n_vertices() {
                    if mesh.mesh.is_boundary_vertex(v) {
                        let nu = l_u[v] + k[v];
                        assert!(
                            (s.boundary_curvature[v] - nu).abs() < 1e-8,
                            "{} {boundary:?} vertex {v}: {} vs {nu}",
                            mesh.name,
                            s.boundary_curvature[v]
                        );
                        if boundary != BoundaryTreatment::Fixed {
                            assert_eq!(s.boundary_curvature[v], s.measure[v]);
                        } else {
                            assert_eq!(s.measure[v], 0.0);
                            assert_eq!(s.log_scale_optimal[v], 0.0);
                        }
                        if boundary == BoundaryTreatment::Free {
                            assert!(s.adjoint[v].abs() < 1e-9, "adjoint on a free boundary");
                        }
                        if boundary == BoundaryTreatment::Convex {
                            assert!(s.boundary_curvature[v] >= 0.0);
                        }
                    } else {
                        assert_eq!(s.boundary_curvature[v], 0.0);
                    }
                    total += s.measure[v]
                        + s.boundary_curvature[v]
                            * f64::from(u8::from(boundary == BoundaryTreatment::Fixed));
                }
                // Gauss–Bonnet: interior angles plus boundary curvatures. (With `Fixed` the
                // boundary curvature is not in `measure`; it is added above.)
                assert!(
                    (total - 2.0 * std::f64::consts::PI * chi).abs() < 1e-7,
                    "{} {boundary:?}: {total} vs {}",
                    mesh.name,
                    2.0 * std::f64::consts::PI * chi
                );
                if boundary != BoundaryTreatment::Fixed {
                    assert!(
                        s.cones
                            .iter()
                            .all(|c| !mesh.mesh.is_boundary_vertex(c.vertex as usize))
                    );
                    let again = certify_with_boundary(
                        &f.problem(&mesh),
                        boundary,
                        &s.log_scale_optimal,
                        &s.adjoint,
                        &s.measure,
                    )
                    .unwrap();
                    assert_eq!(again, s.report.certificate);
                    // The same triple is not an optimum of the fixed-boundary problem.
                    let fixed = certify(
                        &f.problem(&mesh),
                        &s.log_scale_optimal,
                        &s.adjoint,
                        &s.measure,
                    )
                    .unwrap();
                    assert!(!fixed.holds(1e-9) || s.log_scale_optimal.iter().all(|&u| u == 0.0));
                }
            }
        }
    }
}

/// A perturbed triple fails the certificate of the new mode, in each of its residuals.
#[test]
fn certificate_of_the_boundary_modes_detects_errors() {
    let mesh = random_mesh(1, 5);
    let f = Features::plain(0.01);
    for boundary in TREATMENTS {
        let s = solved(&mesh, &f, boundary);
        let problem = f.problem(&mesh);
        let b = (0..mesh.mesh.n_vertices())
            .find(|&v| mesh.mesh.is_boundary_vertex(v))
            .unwrap();
        // A wrong boundary curvature breaks the state equation (and Gauss–Bonnet).
        let mut measure = s.measure.clone();
        measure[b] += 0.1;
        let c = certify_with_boundary(
            &problem,
            boundary,
            &s.log_scale_optimal,
            &s.adjoint,
            &measure,
        )
        .unwrap();
        assert!(
            c.state_residual > 0.05 && !c.holds(1e-9),
            "{boundary:?}: {c:?}"
        );
        assert!(c.gauss_bonnet_residual > 0.05);
        // A wrong boundary value breaks the adjoint equation and the state equation.
        let mut u = s.log_scale_optimal.clone();
        u[b] += 0.1;
        let c = certify_with_boundary(&problem, boundary, &u, &s.adjoint, &s.measure).unwrap();
        assert!(!c.holds(1e-6), "{boundary:?}: {c:?}");
        // A wrong adjoint on the boundary breaks the site condition of `Free`.
        let mut phi = s.adjoint.clone();
        phi[b] += 0.1;
        let c = certify_with_boundary(&problem, boundary, &s.log_scale_optimal, &phi, &s.measure)
            .unwrap();
        assert!(!c.holds(1e-6), "{boundary:?}: {c:?}");
    }
}

/// When extraction changes the cones, the returned scale factor is the one with the boundary
/// optimized again for them, i.e. the solution of "these interior cones, free boundary".
#[test]
fn extracted_cones_get_their_own_optimal_boundary() {
    let mut checked = 0;
    for kind in 0..3 {
        for seed in 0..60u64 {
            let mesh = random_mesh(kind, seed);
            let nv = mesh.mesh.n_vertices();
            let f = Features::plain(0.0005);
            let options = SolveOptions {
                merge_cone_clusters: true,
                max_merge_cluster_size: usize::MAX,
                ..SolveOptions::default()
            };
            let Ok(s) = solve_with_boundary(&f.problem(&mesh), BoundaryTreatment::Free, &options)
            else {
                continue;
            };
            if !s.extraction.cones_differ_from_measure {
                continue;
            }
            assert_eq!(s.report.status, Status::Converged);
            // The independent route: the extracted cones prescribed, everything else forbidden,
            // the boundary free: the minimum-distortion scale factor for these interior cones.
            let mut fixed_cones = Features::plain(0.0);
            fixed_cones.cone_costs = Some(vec![f64::INFINITY; nv]);
            fixed_cones.prescribed_cones = s.cones.clone();
            let again = solved(&mesh, &fixed_cones, BoundaryTreatment::Free);
            let err = max_abs_diff(&again.log_scale_optimal, &s.log_scale);
            assert!(err < 1e-7, "{}: {err}", mesh.name);
            assert!((again.distortion_optimal - s.distortion).abs() < 1e-9);
            // Not optimal for the penalized problem, so never better than the optimum.
            let cost: f64 = s.cones.iter().map(|c| 0.0005 * c.angle.abs()).sum();
            assert!(s.distortion + cost >= s.objective - 1e-9);
            checked += 1;
        }
    }
    assert!(checked >= 3, "only {checked} cases had extraction changes");
    println!("{checked} cases with merged clusters");
}

#[test]
fn closed_meshes_ignore_the_treatment() {
    for seed in 0..5u64 {
        for kind in 3..N_KINDS {
            let mesh = random_mesh(kind, seed);
            let f = Features::plain(0.05);
            let base = solve(&f.problem(&mesh), &SolveOptions::default()).unwrap();
            for boundary in TREATMENTS {
                let s = solve_with_boundary(&f.problem(&mesh), boundary, &SolveOptions::default())
                    .unwrap();
                assert_eq!(s.measure, base.measure);
                assert_eq!(s.log_scale_optimal, base.log_scale_optimal);
                assert!(s.boundary_curvature.iter().all(|&x| x == 0.0));
            }
        }
    }
}
