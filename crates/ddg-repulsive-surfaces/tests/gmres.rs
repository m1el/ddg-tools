//! Linear solves (`SPEC.md` §6, §13.3): the dense nested saddle and the Schur layer against a
//! full saddle LU, the Newton correction and its safeguard, GMRES (both sides) against direct
//! solves, the dense preconditioner on `A′`, and the oracle's `saddle` and `precond` fixtures.

#![allow(clippy::needless_range_loop, clippy::type_complexity)]

#[path = "wp2_fixture/mod.rs"]
mod fixture;

use ddg_core::SparseMatrix;
use ddg_core::generators::{SplitMix64, bumpy_grid_disk, icosphere, perturb, torus};
use ddg_repulsive_surfaces::constraints::{BarycenterMode, Constraint, ConstraintSet};
use ddg_repulsive_surfaces::energy::exact::{energy_and_differential, make_pool};
use ddg_repulsive_surfaces::geometry::{FaceGeometry, mean_edge_length};
use ddg_repulsive_surfaces::metric::dense::dense_metric;
use ddg_repulsive_surfaces::metric::{DenseMatrix, FractionalPairs, SobolevWeights};
use ddg_repulsive_surfaces::solve::gmres::FnOperator;
use ddg_repulsive_surfaces::solve::{
    DensePreconditioner, DenseSaddle, DenseSaddleOperator, GmresOptions, GmresSaddle, GmresSide,
    IdentityPreconditioner, LinearOperator, NewtonSteps, Preconditioner, ProjectionOptions,
    ProjectionStatus, SaddleSolver, SchurComplement, gmres, project,
};
use ddg_repulsive_surfaces::{Exponents, RsError, Surface};
use faer::Mat;

fn perturbed(m: ddg_core::generators::GeneratedMesh, amp: f64, seed: u64) -> (Surface, Vec<[f64; 3]>) {
    let mut p = m.positions;
    let s = Surface::new(p.len(), m.faces).unwrap();
    let l = mean_edge_length(&s, &p);
    perturb(&mut p, [amp * l; 3], seed);
    (s, p)
}

fn flat(v: &[[f64; 3]]) -> Vec<f64> {
    v.iter().flatten().copied().collect()
}

fn norm(x: &[f64]) -> f64 {
    x.iter().map(|v| v * v).sum::<f64>().sqrt()
}

fn dense_rows(c: &SparseMatrix) -> Vec<Vec<f64>> {
    let mut d = vec![vec![0.0; c.ncols()]; c.nrows()];
    for (r, col, x) in c.triplets() {
        d[r as usize][col as usize] += x;
    }
    d
}

/// Solves the full saddle `[[A₃, Cᵀ], [C, 0]] [x; μ] = [a; b]` by a plain LU (all rows inside).
fn full_saddle(a: &DenseMatrix, c: &SparseMatrix, rhs_a: &[f64], rhs_b: &[f64]) -> Vec<f64> {
    use faer::linalg::solvers::Solve;
    let nv = a.n();
    let nd = 3 * nv;
    let k = c.nrows();
    let n = nd + k;
    let cd = dense_rows(c);
    let m = Mat::<f64>::from_fn(n, n, |i, j| {
        if i < nd && j < nd {
            if i % 3 == j % 3 { a.get(i / 3, j / 3) } else { 0.0 }
        } else if i >= nd && j < nd {
            cd[i - nd][j]
        } else if i < nd && j >= nd {
            cd[j - nd][i]
        } else {
            0.0
        }
    });
    let rhs = Mat::<f64>::from_fn(n, 1, |i, _| if i < nd { rhs_a[i] } else { rhs_b[i - nd] });
    let lu = m.partial_piv_lu();
    let x = lu.solve(&rhs);
    (0..n).map(|i| x[(i, 0)]).collect()
}

struct Problem {
    s: Surface,
    p: Vec<[f64; 3]>,
    set: ConstraintSet,
    a: DenseMatrix,
    g: Vec<f64>,
}

fn problem(
    s: Surface,
    p: Vec<[f64; 3]>,
    cons: &[Constraint],
    e: Exponents,
    w: SobolevWeights,
) -> Problem {
    let pool = make_pool(Some(4)).unwrap();
    let set = ConstraintSet::new(&s, &p, cons, BarycenterMode::PerComponent).unwrap();
    let dm = dense_metric(&pool, &s, &p, &e, &w, None).unwrap();
    let geom = FaceGeometry::new(&s, &p).unwrap();
    let (_, g) = energy_and_differential(&pool, &s, &p, &geom, &e).unwrap();
    Problem {
        s,
        p,
        set,
        a: dm.a,
        g: flat(&g),
    }
}

#[test]
fn schur_direction_and_correction_equal_the_full_saddle() {
    let pool = make_pool(Some(4)).unwrap();
    let two = {
        let m = icosphere(1).unwrap();
        let mut p = m.positions.clone();
        let mut f = m.faces.clone();
        let off = p.len() as u32;
        p.extend(m.positions.iter().map(|q| [0.7 * q[0] + 2.2, 0.7 * q[1] + 0.3, 0.7 * q[2]]));
        f.extend(m.faces.iter().map(|t| [t[0] + off, t[1] + off, t[2] + off]));
        perturbed(ddg_core::generators::GeneratedMesh { positions: p, faces: f }, 0.1, 7)
    };
    let cases: Vec<(&str, (Surface, Vec<[f64; 3]>), Vec<Constraint>)> = vec![
        (
            "icosphere: barycenter, area, volume",
            perturbed(icosphere(2).unwrap(), 0.2, 1),
            vec![
                Constraint::Barycenter { component: None, target: None },
                Constraint::Area { component: None, target: Some(13.0) },
                Constraint::Volume { component: None, target: Some(4.0) },
            ],
        ),
        (
            "torus: pins, volume",
            perturbed(torus(10, 6, 1.5, 0.6).unwrap(), 0.1, 2),
            vec![
                Constraint::Pins { vertices: vec![0, 17, 33], targets: None },
                Constraint::Volume { component: None, target: None },
            ],
        ),
        (
            "disk: pins and barycenter, area",
            perturbed(bumpy_grid_disk(5, 5, 0.3).unwrap(), 0.1, 3),
            vec![
                Constraint::Pins { vertices: vec![0, 5], targets: None },
                Constraint::Barycenter { component: None, target: None },
                Constraint::Area { component: None, target: Some(1.1) },
            ],
        ),
        (
            "two components: per-component barycenters, area of 1, volume of 0",
            two,
            vec![
                Constraint::Barycenter { component: None, target: None },
                Constraint::Area { component: Some(1), target: None },
                Constraint::Volume { component: Some(0), target: Some(3.0) },
            ],
        ),
    ];
    for (name, (s, p), cons) in cases {
        let pr = problem(s, p, &cons, Exponents::default(), SobolevWeights::default());
        let q = pr.set.positional_rows();
        let nl = pr.set.nonlinear_rows();
        let cq = pr.set.jacobian(&pr.p, &q).unwrap();
        let cn = pr.set.jacobian(&pr.p, &nl).unwrap();
        let mut saddle = DenseSaddle::new(&pool, &pr.a, &cq).unwrap();
        let schur = SchurComplement::new(&mut saddle, nl.clone(), cn.clone()).unwrap();
        assert!(schur.dropped_rows().is_empty(), "{name}");
        let dir = schur.direction(&mut saddle, &pr.g).unwrap();
        // The full saddle with every row inside.
        let all: Vec<usize> = q.iter().chain(&nl).copied().collect();
        let c_all = pr.set.jacobian(&pr.p, &all).unwrap();
        let nd = 3 * pr.s.n_vertices();
        let full = full_saddle(&pr.a, &c_all, &pr.g, &vec![0.0; all.len()]);
        let err = dir
            .x
            .iter()
            .zip(&full[..nd])
            .fold(0.0f64, |m, (a, b)| m.max((a - b).abs()))
            / norm(&full[..nd]);
        // C x = 0 on all rows.
        let mut cx = vec![0.0; all.len()];
        c_all.mul_vec(&dir.x, &mut cx).unwrap();
        let cxn = cx.iter().fold(0.0f64, |m, v| m.max(v.abs())) / norm(&dir.x);
        // Correction: C_N h = Φ_N, C_Q h = 0, and the full saddle with right side (0; 0, Φ_N).
        let phi = pr.set.values(&pr.p).unwrap();
        let phin: Vec<f64> = nl.iter().map(|&r| phi[r]).collect();
        let h = schur.correction(&phin);
        let mut b = vec![0.0; q.len()];
        b.extend(&phin);
        let fullh = full_saddle(&pr.a, &c_all, &vec![0.0; nd], &b);
        let errh = h
            .iter()
            .zip(&fullh[..nd])
            .fold(0.0f64, |m, (a, b)| m.max((a - b).abs()))
            / norm(&fullh[..nd]).max(1e-300);
        let mut ch = vec![0.0; all.len()];
        c_all.mul_vec(&h, &mut ch).unwrap();
        let mut errc: f64 = 0.0;
        for (i, v) in ch.iter().enumerate() {
            let want = if i < q.len() { 0.0 } else { phin[i - q.len()] };
            errc = errc.max((v - want).abs());
        }
        let scale = phin.iter().fold(1e-300f64, |m, v| m.max(v.abs()));
        eprintln!(
            "{name}: |x - x_full|/|x| {err:.1e}, |Cx|/|x| {cxn:.1e}, |h - h_full|/|h| {errh:.1e}, \
             |Ch - Φ|/|Φ| {:.1e}, pivot {:.2e}, Schur pivot {:.2e}",
            errc / scale,
            saddle.pivot_ratio(),
            schur.min_pivot()
        );
        assert!(err < 1e-10, "{name}: direction {err:e}");
        assert!(cxn < 1e-10, "{name}: C x {cxn:e}");
        if !nl.is_empty() {
            assert!(errh < 1e-10, "{name}: correction {errh:e}");
            assert!(errc < 1e-10 * scale, "{name}: C h");
        }
        // ⟨g, x⟩ = xᵀ A₃ x ≥ 0 (descent).
        let gx: f64 = pr.g.iter().zip(&dir.x).map(|(a, b)| a * b).sum();
        assert!(gx > 0.0, "{name}: <g, x> = {gx:e}");
    }
}

#[test]
fn general_positional_jacobian_uses_the_full_lu() {
    // A positional Jacobian that is not C₁ ⊗ I₃ (one coordinate only) takes the full LU and
    // agrees with the full saddle.
    let pool = make_pool(Some(2)).unwrap();
    let (s, p) = perturbed(icosphere(1).unwrap(), 0.1, 4);
    let pr = problem(
        s,
        p,
        &[Constraint::Barycenter { component: None, target: None }],
        Exponents::default(),
        SobolevWeights::default(),
    );
    let nv = pr.s.n_vertices();
    let mut trip: Vec<(u32, u32, f64)> = (0..nv).map(|v| (0, 3 * v as u32, 1.0)).collect();
    trip.extend((0..nv).map(|v| (1, 3 * v as u32 + 1, 1.0)));
    trip.extend((0..nv).map(|v| (2, 3 * v as u32 + 2, 1.0 + v as f64 * 0.01)));
    let c = SparseMatrix::from_triplets(3, 3 * nv, &trip).unwrap();
    let mut saddle = DenseSaddle::new(&pool, &pr.a, &c).unwrap();
    assert!(!saddle.is_per_axis());
    let sol = saddle.solve(&pr.g, &[0.0; 3]).unwrap();
    let full = full_saddle(&pr.a, &c, &pr.g, &[0.0; 3]);
    let err = sol
        .x
        .iter()
        .zip(&full)
        .fold(0.0f64, |m, (a, b)| m.max((a - b).abs()))
        / norm(&full[..3 * nv]);
    assert!(err < 1e-10, "{err:e}");
    // The per-axis path on the barycenter rows.
    let cq = pr.set.jacobian(&pr.p, &pr.set.positional_rows()).unwrap();
    assert!(DenseSaddle::new(&pool, &pr.a, &cq).unwrap().is_per_axis());
    // No constraint: singular.
    let empty = SparseMatrix::from_triplets(0, 3 * nv, &[]).unwrap();
    assert!(matches!(
        DenseSaddle::new(&pool, &pr.a, &empty),
        Err(RsError::SingularSaddle { .. })
    ));
}

#[test]
fn newton_projection_converges_quadratically_and_the_safeguard_holds() {
    let pool = make_pool(Some(8)).unwrap();
    // A near-sphere with area and volume: ∇A ≈ (2/r) ∇V, S nearly singular.
    let (s, p) = perturbed(icosphere(3).unwrap(), 0.02, 5);
    let geom = FaceGeometry::new(&s, &p).unwrap();
    let a0: f64 = geom.areas.iter().sum();
    let cons = vec![
        Constraint::Barycenter { component: None, target: None },
        Constraint::Area { component: None, target: Some(1.03 * a0) },
        Constraint::Volume { component: None, target: None },
    ];
    let pr = problem(s, p, &cons, Exponents::default(), SobolevWeights::default());
    let q = pr.set.positional_rows();
    let nl = pr.set.nonlinear_rows();
    let cq = pr.set.jacobian(&pr.p, &q).unwrap();
    let cn = pr.set.jacobian(&pr.p, &nl).unwrap();
    let mut saddle = DenseSaddle::new(&pool, &pr.a, &cq).unwrap();
    let schur = SchurComplement::new(&mut saddle, nl.clone(), cn).unwrap();
    let dir = schur.direction(&mut saddle, &pr.g).unwrap();
    eprintln!(
        "near-sphere: Schur pivot {:.3e}, dropped {:?}",
        schur.min_pivot(),
        schur.dropped_rows()
    );
    // A trial that moves the farthest vertex by one mean edge length.
    let l = mean_edge_length(&pr.s, &pr.p);
    let xmax = dir.x.chunks(3).map(norm).fold(0.0f64, f64::max);
    let tau = l / xmax;
    let trial: Vec<[f64; 3]> = pr
        .p
        .iter()
        .enumerate()
        .map(|(v, q)| {
            [
                q[0] - tau * dir.x[3 * v],
                q[1] - tau * dir.x[3 * v + 1],
                q[2] - tau * dir.x[3 * v + 2],
            ]
        })
        .collect();
    let proj = project(&pr.set, Some(&schur), &trial, &ProjectionOptions::default()).unwrap();
    eprintln!(
        "near-sphere projection: status {:?}, rounds {}, history {:?}, max violation {:.2e}",
        proj.status, proj.rounds, proj.history, proj.max_violation
    );
    assert!(proj.history.iter().all(|h| h.is_finite()));
    for w in proj.history.windows(2).take(proj.rounds as usize) {
        assert!(w[1] < w[0], "kept rounds lower the violation");
    }
    assert!(proj.history[proj.rounds as usize] <= proj.history[0]);
    // Without the safeguard, the same frozen-Jacobian rounds (reported for comparison).
    let mut pos = trial.clone();
    pr.set.project_positional(&mut pos).unwrap();
    let mut raw = Vec::new();
    for _ in 0..6 {
        let phi = pr.set.values(&pos).unwrap();
        let phin: Vec<f64> = nl.iter().map(|&r| phi[r]).collect();
        raw.push(
            nl.iter()
                .map(|&r| phi[r].abs() / pr.set.tolerance(r))
                .fold(0.0f64, f64::max),
        );
        let h = schur.correction(&phin);
        for (v, q) in pos.iter_mut().enumerate() {
            for k in 0..3 {
                q[k] -= h[3 * v + k];
            }
        }
        pr.set.project_positional(&mut pos).unwrap();
    }
    eprintln!("near-sphere, unguarded rounds: {raw:?}");

    // A well-conditioned case: one area row, quadratic convergence to tolerance.
    let (s, p) = perturbed(torus(16, 8, 1.5, 0.6).unwrap(), 0.1, 6);
    let geom = FaceGeometry::new(&s, &p).unwrap();
    let a0: f64 = geom.areas.iter().sum();
    let cons = vec![
        Constraint::Barycenter { component: None, target: None },
        Constraint::Area { component: None, target: Some(1.02 * a0) },
    ];
    let pr = problem(s, p, &cons, Exponents::default(), SobolevWeights::default());
    let q = pr.set.positional_rows();
    let nl = pr.set.nonlinear_rows();
    let cq = pr.set.jacobian(&pr.p, &q).unwrap();
    let cn = pr.set.jacobian(&pr.p, &nl).unwrap();
    let mut saddle = DenseSaddle::new(&pool, &pr.a, &cq).unwrap();
    let schur = SchurComplement::new(&mut saddle, nl, cn).unwrap();
    let proj = project(&pr.set, Some(&schur), &pr.p, &ProjectionOptions::default()).unwrap();
    eprintln!(
        "torus area +2 %: status {:?}, rounds {}, history {:?}",
        proj.status, proj.rounds, proj.history
    );
    assert_eq!(proj.status, ProjectionStatus::Converged);
    // Reference mode: exactly one round.
    let one = ProjectionOptions { newton: NewtonSteps::Fixed(1) };
    let proj1 = project(&pr.set, Some(&schur), &pr.p, &one).unwrap();
    assert_eq!(proj1.rounds, 1);
}

#[test]
fn schur_drops_a_dependent_row() {
    // Area and the same area again (as a second constraint on component 0 of a connected
    // mesh): the rank check drops it already; force both into the Schur complement directly.
    let pool = make_pool(Some(2)).unwrap();
    let (s, p) = perturbed(icosphere(1).unwrap(), 0.1, 8);
    let pr = problem(
        s,
        p,
        &[
            Constraint::Barycenter { component: None, target: None },
            Constraint::Area { component: None, target: None },
            Constraint::Area { component: Some(0), target: None },
        ],
        Exponents::default(),
        SobolevWeights::default(),
    );
    assert_eq!(pr.set.dropped_rows(), &[4]);
    let q = pr.set.positional_rows();
    let cq = pr.set.jacobian(&pr.p, &q).unwrap();
    let mut saddle = DenseSaddle::new(&pool, &pr.a, &cq).unwrap();
    let rows = vec![3, 4];
    let cn = pr.set.jacobian(&pr.p, &rows).unwrap();
    let schur = SchurComplement::new(&mut saddle, rows, cn).unwrap();
    assert_eq!(schur.dropped_rows(), vec![4]);
    let dir = schur.direction(&mut saddle, &pr.g).unwrap();
    assert!(dir.max_nonlinear_residual < 1e-10 * norm(&dir.x));
}

/// A random SPD matrix `QᵀQ + I` and a symmetric indefinite one.
fn random_dense(n: usize, seed: u64, indefinite: bool) -> Vec<f64> {
    let mut rng = SplitMix64::new(seed);
    let q: Vec<f64> = (0..n * n).map(|_| rng.uniform(-1.0, 1.0)).collect();
    let mut a = vec![0.0; n * n];
    for i in 0..n {
        for j in 0..n {
            let mut s = 0.0;
            for k in 0..n {
                s += q[k * n + i] * q[k * n + j];
            }
            a[i * n + j] = s / n as f64;
        }
        a[i * n + i] += if indefinite && i % 3 == 0 { -2.0 } else { 1.0 };
    }
    a
}

fn direct(n: usize, a: &[f64], b: &[f64]) -> Vec<f64> {
    use faer::linalg::solvers::Solve;
    let m = Mat::<f64>::from_fn(n, n, |i, j| a[i * n + j]);
    let r = Mat::<f64>::from_fn(n, 1, |i, _| b[i]);
    let x = m.partial_piv_lu().solve(&r);
    (0..n).map(|i| x[(i, 0)]).collect()
}

#[test]
fn gmres_matches_direct_solves_on_both_sides() {
    let n = 120;
    for (label, indefinite) in [("SPD", false), ("indefinite", true)] {
        let a = random_dense(n, 3, indefinite);
        let op = FnOperator::new(n, |x: &[f64], y: &mut [f64]| {
            for i in 0..n {
                y[i] = (0..n).map(|j| a[i * n + j] * x[j]).sum();
            }
            Ok(())
        });
        let jac = FnOperator::new(n, |r: &[f64], z: &mut [f64]| {
            for i in 0..n {
                z[i] = r[i] / a[i * n + i].abs();
            }
            Ok(())
        });
        let b: Vec<f64> = (0..n).map(|i| (i as f64 * 0.37).sin()).collect();
        let xs = direct(n, &a, &b);
        for side in [GmresSide::Right, GmresSide::Left] {
            for (pname, prec) in [
                ("identity", &IdentityPreconditioner { n } as &dyn Preconditioner),
                ("Jacobi", &jac as &dyn Preconditioner),
            ] {
                let opts = GmresOptions {
                    restart: 30,
                    tolerance: 1e-11,
                    max_iterations: 2000,
                    side,
                };
                let (x, rep) = gmres(&op, prec, &b, &opts).unwrap();
                let err = x
                    .iter()
                    .zip(&xs)
                    .fold(0.0f64, |m, (p, q)| m.max((p - q).abs()))
                    / norm(&xs);
                eprintln!(
                    "{label} {side:?} {pname}: {} iterations, {} cycles, residual {:.1e}, error {err:.1e}",
                    rep.iterations, rep.cycles, rep.relative_residual
                );
                assert!(rep.converged);
                assert!(err < 1e-8, "{label} {side:?} {pname}: {err:e}");
                if side == GmresSide::Right {
                    assert!(rep.relative_residual <= 1e-11);
                } else {
                    assert!(rep.preconditioned_relative_residual.unwrap() <= 1e-11);
                }
            }
        }
    }
    // Not converging within the cap is reported (residual below 1e-1) or an error (above).
    let a = random_dense(n, 4, true);
    let op = FnOperator::new(n, |x: &[f64], y: &mut [f64]| {
        for i in 0..n {
            y[i] = (0..n).map(|j| a[i * n + j] * x[j]).sum();
        }
        Ok(())
    });
    let b: Vec<f64> = (0..n).map(|i| (i as f64).cos()).collect();
    let opts = GmresOptions {
        restart: 5,
        tolerance: 1e-12,
        max_iterations: 3,
        side: GmresSide::Right,
    };
    match gmres(&op, &IdentityPreconditioner { n }, &b, &opts) {
        Ok((_, rep)) => assert!(!rep.converged && rep.relative_residual <= 1e-1),
        Err(RsError::SolverDiverged { residual }) => assert!(residual > 1e-1),
        Err(e) => panic!("{e}"),
    }
    // Zero right-hand side.
    let (x, rep) = gmres(&op, &IdentityPreconditioner { n }, &vec![0.0; n], &Default::default())
        .unwrap();
    assert!(rep.converged && x.iter().all(|v| *v == 0.0));
}

#[test]
fn gmres_on_the_saddle_with_the_dense_preconditioner() {
    // The paper expects few iterations with P⁻¹ = L′⁻¹ L_{2−s} L′⁻¹ (Sec. 5.2; SPEC §13.3:
    // fewer than 50 at 1e-4). Meshes up to 1500 faces, barycenter constraint, right side the
    // energy differential; L_{2−s} over face pairs (the SPEC) and over vertex pairs.
    let mut meshes: Vec<(&str, (Surface, Vec<[f64; 3]>))> = vec![
        ("icosphere 320", perturbed(icosphere(2).unwrap(), 0.2, 1)),
        ("bumpy disk 200", perturbed(bumpy_grid_disk(10, 10, 0.3).unwrap(), 0.1, 4)),
    ];
    // The large meshes in optimized builds only (seconds there, minutes unoptimized).
    if !cfg!(debug_assertions) {
        meshes.extend([
            ("icosphere 1280", perturbed(icosphere(3).unwrap(), 0.2, 2)),
            ("torus 1500", perturbed(torus(30, 25, 1.5, 0.6).unwrap(), 0.1, 3)),
            ("bumpy disk 1250", perturbed(bumpy_grid_disk(25, 25, 0.3).unwrap(), 0.1, 4)),
        ]);
    }
    let results: Vec<String> = std::thread::scope(|scope| {
        let handles: Vec<_> = meshes
            .into_iter()
            .map(|(name, (s, p))| {
                scope.spawn(move || {
                    let pool = make_pool(Some(8)).unwrap();
                    let e = Exponents::default();
                    let pr = problem(
                        s,
                        p,
                        &[Constraint::Barycenter { component: None, target: None }],
                        e,
                        SobolevWeights::default(),
                    );
                    let q = pr.set.positional_rows();
                    let cq = pr.set.jacobian(&pr.p, &q).unwrap();
                    let mut saddle = DenseSaddle::new(&pool, &pr.a, &cq).unwrap();
                    let xs = saddle.solve(&pr.g, &vec![0.0; q.len()]).unwrap().x;
                    let xmax = xs.iter().fold(0.0f64, |m, v| m.max(v.abs()));
                    let op = DenseSaddleOperator::new(&pr.a, &cq).unwrap();
                    let mut rhs = pr.g.clone();
                    rhs.extend(vec![0.0; q.len()]);
                    let nd = 3 * pr.s.n_vertices();
                    let mut out = format!("{name} ({} faces):", pr.s.n_faces());
                    for pairs in [FractionalPairs::Faces, FractionalPairs::Vertices] {
                        let prec =
                            DensePreconditioner::new(&pool, &pr.s, &pr.p, &e, &cq, 1e-10, pairs)
                                .unwrap();
                        out += &format!("\n    {pairs:?}:");
                        for side in [GmresSide::Right, GmresSide::Left] {
                            for tol in [1e-4, 1e-8] {
                                let opts = GmresOptions {
                                    tolerance: tol,
                                    side,
                                    ..Default::default()
                                };
                                let (x, rep) = match gmres(&op, &prec, &rhs, &opts) {
                                    Ok(r) => r,
                                    Err(err) => {
                                        out += &format!(" {side:?} {tol:.0e}: {err};");
                                        continue;
                                    }
                                };
                                let err = x[..nd]
                                    .iter()
                                    .zip(&xs)
                                    .fold(0.0f64, |m, (a, b)| m.max((a - b).abs()))
                                    / xmax;
                                out += &format!(
                                    " {side:?} {tol:.0e}: {}{} it (true res {:.0e}, err {err:.0e});",
                                    rep.iterations,
                                    if rep.converged { "" } else { " (not converged)" },
                                    rep.relative_residual
                                );
                                if pairs == FractionalPairs::Vertices && side == GmresSide::Right {
                                    assert!(rep.converged, "{name} {tol}");
                                    assert!(err < 30.0 * tol, "{name} {tol}: {err:e}");
                                    if tol == 1e-4 {
                                        assert!(rep.iterations < 50, "{name}: {}", rep.iterations);
                                    }
                                }
                            }
                        }
                    }
                    let opts = GmresOptions {
                        max_iterations: 300,
                        ..Default::default()
                    };
                    let id = IdentityPreconditioner { n: op.dim() };
                    out += &match gmres(&op, &id, &rhs, &opts) {
                        Ok((_, rep)) => format!(
                            "\n    no preconditioner, right 1e-4: {} it, converged {}",
                            rep.iterations, rep.converged
                        ),
                        Err(err) => format!("\n    no preconditioner: {err}"),
                    };
                    out
                })
            })
            .collect();
        handles.into_iter().map(|h| h.join().unwrap()).collect()
    });
    for l in results {
        eprintln!("{l}");
    }
}

#[test]
fn schur_on_gmres_matches_dense() {
    let pool = make_pool(Some(4)).unwrap();
    let (s, p) = perturbed(icosphere(2).unwrap(), 0.2, 9);
    let e = Exponents::default();
    let pr = problem(
        s,
        p,
        &[
            Constraint::Barycenter { component: None, target: None },
            Constraint::Area { component: None, target: Some(12.0) },
            Constraint::Volume { component: None, target: Some(4.0) },
        ],
        e,
        SobolevWeights::default(),
    );
    let q = pr.set.positional_rows();
    let nl = pr.set.nonlinear_rows();
    let cq = pr.set.jacobian(&pr.p, &q).unwrap();
    let cn = pr.set.jacobian(&pr.p, &nl).unwrap();
    let mut dense = DenseSaddle::new(&pool, &pr.a, &cq).unwrap();
    let sd = SchurComplement::new(&mut dense, nl.clone(), cn.clone()).unwrap();
    let xd = sd.direction(&mut dense, &pr.g).unwrap();
    let op = DenseSaddleOperator::new(&pr.a, &cq).unwrap();
    let prec = DensePreconditioner::new(&pool, &pr.s, &pr.p, &e, &cq, 1e-10, FractionalPairs::Vertices)
        .unwrap();
    for side in [GmresSide::Right, GmresSide::Left] {
        let mut gs = GmresSaddle::new(
            &op,
            &prec,
            3 * pr.s.n_vertices(),
            GmresOptions {
                tolerance: 1e-10,
                side,
                ..Default::default()
            },
        )
        .unwrap();
        let sg = SchurComplement::new(&mut gs, nl.clone(), cn.clone()).unwrap();
        let xg = sg.direction(&mut gs, &pr.g).unwrap();
        let err = xg
            .x
            .iter()
            .zip(&xd.x)
            .fold(0.0f64, |m, (a, b)| m.max((a - b).abs()))
            / norm(&xd.x);
        let its: Vec<usize> = sg.records().iter().map(|r| r.iterations).collect();
        eprintln!(
            "Schur on GMRES ({side:?}): error {err:.1e}, iterations {its:?} + {}, asymmetry {:.1e}",
            xg.record.iterations,
            sg.asymmetry()
        );
        assert!(err < 1e-7, "{side:?}: {err:e}");
    }
}

// ---- the oracle's golden fixtures (relative 1e-9) ----

#[test]
fn oracle_saddle_fixtures() {
    let names = fixture::cases("saddle");
    assert!(names.len() >= 6);
    let pool = make_pool(Some(4)).unwrap();
    for name in names {
        let c = fixture::read_case(&name);
        let o = fixture::read_output(&name);
        let s = Surface::new(c.positions.len(), c.faces.clone()).unwrap();
        let (al, be) = c.exponents.unwrap_or((6.0, 12.0));
        let e = Exponents::new(al, be).unwrap();
        let (wtp, clo, ww) = c.weights.unwrap_or((1.0, 1.0, 0.0));
        let w = SobolevWeights {
            tp_weight: wtp,
            low_order: clo != 0.0,
            willmore_weight: ww,
        };
        let pr = problem(s, c.positions.clone(), &c.constraints, e, w);
        assert!(pr.set.dropped_rows().is_empty(), "{name}");
        let q = pr.set.positional_rows();
        let nl = pr.set.nonlinear_rows();
        let rows = &o.get("rows")[0];
        assert_eq!(q.len(), rows[0].parse::<usize>().unwrap(), "{name}");
        assert_eq!(nl.len(), rows[1].parse::<usize>().unwrap(), "{name}");
        let cq = pr.set.jacobian(&pr.p, &q).unwrap();
        let cn = pr.set.jacobian(&pr.p, &nl).unwrap();
        let mut saddle = DenseSaddle::new(&pool, &pr.a, &cq).unwrap();
        let schur = SchurComplement::new(&mut saddle, nl.clone(), cn).unwrap();
        let dir = schur.direction(&mut saddle, &pr.g).unwrap();
        let phi = pr.set.values(&pr.p).unwrap();
        let phin: Vec<f64> = nl.iter().map(|&r| phi[r]).collect();
        let h = schur.correction(&phin);
        let gx: f64 = pr.g.iter().zip(&dir.x).map(|(a, b)| a * b).sum();
        let mut checks: Vec<(&str, Vec<f64>, Vec<f64>)> = vec![
            ("g", pr.g.clone(), o.vectors("g")),
            ("x", dir.x.clone(), o.vectors("x")),
            ("gx", vec![gx], vec![o.scalar("gx")]),
        ];
        if !nl.is_empty() {
            let phin_o: Vec<f64> = o.get("phin")[0].iter().map(|x| fixture::hex(x)).collect();
            checks.push(("phin", phin, phin_o));
            checks.push(("S", schur.matrix().to_vec(), o.matrix("S").concat()));
            checks.push(("h", h, o.vectors("h")));
        }
        for (what, got, want) in checks {
            let err = fixture::rel_err(&got, &want, 1e-300);
            eprintln!("{name}: {what} relative error {err:.2e}");
            assert!(err < 1e-9, "{name} {what}: {err:e}");
        }
    }
}

#[test]
fn oracle_precond_fixtures() {
    let names = fixture::cases("precond");
    assert!(names.len() >= 3);
    let pool = make_pool(Some(4)).unwrap();
    for name in names {
        let c = fixture::read_case(&name);
        let o = fixture::read_output(&name);
        let s = Surface::new(c.positions.len(), c.faces.clone()).unwrap();
        let (al, be) = c.exponents.unwrap_or((6.0, 12.0));
        let e = Exponents::new(al, be).unwrap();
        let set = ConstraintSet::new(&s, &c.positions, &c.constraints, Default::default()).unwrap();
        let q = set.positional_rows();
        assert_eq!(q.len(), o.get("rows")[0][0].parse::<usize>().unwrap());
        let cq = set.jacobian(&c.positions, &q).unwrap();
        let prec =
            DensePreconditioner::new(
            &pool,
            &s,
            &c.positions,
            &e,
            &cq,
            c.eps_rel.unwrap(),
            FractionalPairs::Faces,
        )
        .unwrap();
        let mut r = flat(&c.r);
        r.extend(vec![0.0; q.len()]);
        let z = prec.apply_full(&r).unwrap();
        let nd = 3 * s.n_vertices();
        let checks: Vec<(&str, Vec<f64>, Vec<f64>)> = vec![
            ("eps", vec![prec.epsilon()], vec![o.scalar("eps")]),
            ("Lfr", prec.fractional().as_slice().to_vec(), o.matrix("Lfr").concat()),
            ("z", z[..nd].to_vec(), o.vectors("z")),
        ];
        for (what, got, want) in checks {
            let err = fixture::rel_err(&got, &want, 1e-300);
            eprintln!("{name}: {what} relative error {err:.2e}");
            assert!(err < 1e-9, "{name} {what}: {err:e}");
        }
        // L′ itself: the operator of the factored saddle against the printed matrix.
        let lp = o.matrix("Lp");
        let n = lp.len();
        let op = prec.laplacian_saddle().operator();
        let mut worst: f64 = 0.0;
        let mut scale: f64 = 0.0;
        let mut e_j = vec![0.0; n];
        let mut col = vec![0.0; n];
        for j in 0..n {
            e_j.iter_mut().for_each(|x| *x = 0.0);
            e_j[j] = 1.0;
            op.apply(&e_j, &mut col).unwrap();
            for i in 0..n {
                worst = worst.max((col[i] - lp[i][j]).abs());
                scale = scale.max(lp[i][j].abs());
            }
        }
        eprintln!("{name}: Lp relative error {:.2e}", worst / scale);
        assert!(worst < 1e-9 * scale, "{name} Lp");
    }
}

#[test]
fn oracle_singular_saddle_is_resolved_by_the_rank_check() {
    // rs-error-singular pins the same vertex twice: the oracle reports SingularSaddle; the
    // rank check drops the duplicate pin and the saddle is regular (SPEC D10).
    let c = fixture::read_case("rs-error-singular");
    let s = Surface::new(c.positions.len(), c.faces.clone()).unwrap();
    let set = ConstraintSet::new(&s, &c.positions, &c.constraints, Default::default()).unwrap();
    eprintln!("rs-error-singular: dropped rows {:?}", set.dropped_rows());
    assert!(!set.dropped_rows().is_empty());
}

/// Symmetric eigen-decomposition (values increasing, vectors as columns).
fn sym_evd(m: &Mat<f64>) -> (Vec<f64>, Mat<f64>) {
    use faer::dyn_stack::{MemBuffer, MemStack};
    use faer::linalg::evd::{ComputeEigenvectors, self_adjoint_evd, self_adjoint_evd_scratch};
    let n = m.nrows();
    let mut s = faer::diag::Diag::<f64>::zeros(n);
    let mut u = Mat::<f64>::zeros(n, n);
    let par = faer::Par::Seq;
    let mut mem = MemBuffer::new(self_adjoint_evd_scratch::<f64>(
        n,
        ComputeEigenvectors::Yes,
        par,
        Default::default(),
    ));
    self_adjoint_evd(
        m.as_ref(),
        s.as_mut(),
        Some(u.as_mut()),
        par,
        MemStack::new(&mut mem),
        Default::default(),
    )
    .unwrap();
    ((0..n).map(|i| s[i]).collect(), u)
}

/// `uᵀ L u = Σ_{v≠w} (u_v − u_w)² a_v a_w |f_v − f_w|^{−(2τ+2)}` (vertex pairs, dual areas).
fn vertex_fractional(p: &[[f64; 3]], mass: &[f64], tau: f64) -> DenseMatrix {
    let n = p.len();
    let mut m = vec![0.0; n * n];
    for v in 0..n {
        for w in 0..n {
            if v == w {
                continue;
            }
            let r2: f64 = (0..3).map(|k| (p[v][k] - p[w][k]).powi(2)).sum();
            let h = 2.0 * mass[v] * mass[w] * r2.powf(-(tau + 1.0));
            m[v * n + w] -= h;
            m[v * n + v] += h;
        }
    }
    DenseMatrix::from_row_major(n, m).unwrap()
}

fn symmetrized(m: &Mat<f64>) -> Mat<f64> {
    Mat::<f64>::from_fn(m.nrows(), m.ncols(), |i, j| 0.5 * (m[(i, j)] + m[(j, i)]))
}

/// The spectrum of the preconditioned scalar operator `G A` on `ker c` (`G = L_C⁻¹ L_{2−s}
/// L_C⁻¹`, one barycenter row), for meshes of growing size: the condition number of `A` alone
/// and with `G`, for `A = B + B⁰` and `A = B`. Diagnostic (`--ignored`).
#[test]
#[ignore]
fn preconditioner_spectrum() {
    use ddg_repulsive_surfaces::geometry::{cotan_stiffness, vertex_areas};
    use ddg_repulsive_surfaces::metric::FaceOperators;
    use ddg_repulsive_surfaces::metric::dense::{fractional_operator, to_dense};
    let pool = make_pool(Some(16)).unwrap();
    for (name, (s, p)) in [
        ("icosphere 80", perturbed(icosphere(1).unwrap(), 0.0, 1)),
        ("icosphere 320", perturbed(icosphere(2).unwrap(), 0.0, 1)),
        ("icosphere 1280", perturbed(icosphere(3).unwrap(), 0.0, 1)),
        ("icosphere 320 perturbed", perturbed(icosphere(2).unwrap(), 0.2, 1)),
        ("torus 600", perturbed(torus(30, 10, 1.5, 0.6).unwrap(), 0.0, 1)),
    ] {
        let n = s.n_vertices();
        let e = Exponents::default();
        let geom = FaceGeometry::new(&s, &p).unwrap();
        let mass = vertex_areas(&s, &geom.areas);
        let tot: f64 = mass.iter().sum();
        let ops = FaceOperators::new(&s, &p).unwrap();
        let w = to_dense(&cotan_stiffness(&s, &p).unwrap());
        let lfr = fractional_operator(&pool, &ops, 2.0 - e.s()).unwrap();
        // L_C⁻¹: top-left block of [[W, cᵀ], [c, 0]]⁻¹.
        let k = Mat::<f64>::from_fn(n + 1, n + 1, |i, j| {
            if i < n && j < n {
                w.get(i, j)
            } else if i == n && j < n {
                mass[j] / tot
            } else if j == n && i < n {
                mass[i] / tot
            } else {
                0.0
            }
        });
        let kinv = {
            use faer::linalg::solvers::Solve;
            k.partial_piv_lu().solve(Mat::<f64>::identity(n + 1, n + 1))
        };
        let lci = Mat::<f64>::from_fn(n, n, |i, j| kinv[(i, j)]);
        let lf = Mat::<f64>::from_fn(n, n, |i, j| lfr.get(i, j));
        let g = symmetrized(&(&lci * &lf * &lci));
        let (gv, gu) = sym_evd(&g);
        let gmax = gv[n - 1];
        let d = Mat::<f64>::from_fn(n, n, |i, j| {
            if i == j && gv[i] > 1e-12 * gmax {
                gv[i].sqrt()
            } else {
                0.0
            }
        });
        let gh = &gu * &d * gu.transpose();
        // Π = I − c cᵀ/|c|² (Euclidean projector onto ker c).
        let c: Vec<f64> = mass.iter().map(|m| m / tot).collect();
        let c2: f64 = c.iter().map(|x| x * x).sum();
        let pi = Mat::<f64>::from_fn(n, n, |i, j| f64::from(u8::from(i == j)) - c[i] * c[j] / c2);
        // Spectral M Δ^t = M^{1/2} (M^{-1/2} W M^{-1/2})^t M^{1/2}: how far are B and L_{2−s}
        // from their continuous orders (generalized condition numbers on the complement of
        // the constants)?
        let ms: Vec<f64> = mass.iter().map(|m| m.sqrt()).collect();
        let (lv, lu) = sym_evd(&symmetrized(&Mat::<f64>::from_fn(n, n, |i, j| {
            w.get(i, j) / (ms[i] * ms[j])
        })));
        let gen_cond = |op: &DenseMatrix, t: f64| -> (f64, f64) {
            // (M Δ^t)^{-1/2} op (M Δ^t)^{-1/2} on the nonzero modes of Δ.
            let d = Mat::<f64>::from_fn(n, n, |i, j| {
                if i == j && lv[i] > 1e-10 * lv[n - 1] { lv[i].powf(-0.5 * t) } else { 0.0 }
            });
            let h = &lu * &d * lu.transpose();
            let o = Mat::<f64>::from_fn(n, n, |i, j| op.get(i, j) / (ms[i] * ms[j]));
            let (v, _) = sym_evd(&symmetrized(&(&h * &o * &h)));
            let top = v[n - 1];
            let lo = v.iter().copied().find(|x| *x > 1e-10 * top).unwrap();
            (lo, top)
        };
        {
            let dm = dense_metric(&pool, &s, &p, &e, &SobolevWeights::default(), None).unwrap();
            let (b0, b1) = gen_cond(&dm.b, e.s());
            let (l0, l1) = gen_cond(&lfr, 2.0 - e.s());
            let lw = Mat::<f64>::from_fn(n, n, |i, j| w.get(i, j));
            let wd = DenseMatrix::from_row_major(n, (0..n * n).map(|k| lw[(k / n, k % n)]).collect()).unwrap();
            let (w0, w1) = gen_cond(&wd, 1.0);
            let lvert = vertex_fractional(&p, &mass, 2.0 - e.s());
            let (v0, v1) = gen_cond(&lvert, 2.0 - e.s());
            eprintln!(
                "{name:>24}: B vs M Δ^s: [{b0:.2e}, {b1:.2e}] cond {:.1}; L_(2-s) vs M Δ^(2-s): \
                 [{l0:.2e}, {l1:.2e}] cond {:.1}; vertex L_(2-s): [{v0:.2e}, {v1:.2e}] cond {:.1}; \
                 W vs M Δ: cond {:.1}",
                b1 / b0,
                l1 / l0,
                v1 / v0,
                w1 / w0
            );
            // cond(G A) with the vertex-based L_{2−s}.
            let lvm = Mat::<f64>::from_fn(n, n, |i, j| lvert.get(i, j));
            let gv2 = symmetrized(&(&lci * &lvm * &lci));
            let (gvv, gvu) = sym_evd(&gv2);
            let gvmax = gvv[n - 1];
            let dv = Mat::<f64>::from_fn(n, n, |i, j| {
                if i == j && gvv[i] > 1e-12 * gvmax {
                    gvv[i].sqrt()
                } else {
                    0.0
                }
            });
            let ghv = &gvu * &dv * gvu.transpose();
            let a = Mat::<f64>::from_fn(n, n, |i, j| dm.a.get(i, j));
            let (tv, _) = sym_evd(&symmetrized(&(&ghv * &a * &ghv)));
            let tmax = tv[n - 1];
            let tmin = tv.iter().copied().find(|x| *x > 1e-10 * tmax).unwrap();
            eprintln!("{name:>24}: cond(G_vertex A) {:.2e}", tmax / tmin);
        }
        for (label, low) in [("A = B + B0", true), ("A = B", false)] {
            let w = SobolevWeights {
                low_order: low,
                ..Default::default()
            };
            let dm = dense_metric(&pool, &s, &p, &e, &w, None).unwrap();
            let a = Mat::<f64>::from_fn(n, n, |i, j| dm.a.get(i, j));
            let (av, _) = sym_evd(&symmetrized(&(&pi * &a * &pi)));
            let amax = av[n - 1];
            let amin = av.iter().copied().find(|x| *x > 1e-10 * amax).unwrap();
            let (tv, _) = sym_evd(&symmetrized(&(&gh * &a * &gh)));
            let tmax = tv[n - 1];
            let nz: Vec<f64> = tv.iter().copied().filter(|x| *x > 1e-10 * tmax).collect();
            let tmin = nz[0];
            let q = |f: f64| nz[((nz.len() - 1) as f64 * f) as usize] / tmin;
            eprintln!(
                "{name:>24} {label:<10}: cond(Π A Π) {:.2e}; cond(G A) {:.2e} (spectrum / min at \
                 10%, 50%, 90%, 99%: {:.1}, {:.1}, {:.1}, {:.1}; zero eigenvalues {})",
                amax / amin,
                tmax / tmin,
                q(0.1),
                q(0.5),
                q(0.9),
                q(0.99),
                n - nz.len()
            );
        }
    }
}
