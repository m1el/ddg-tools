# WP2 notes (dense metric, constraints, saddle + Schur, GMRES, dense preconditioner)

Build: `CARGO_TARGET_DIR=$PWD/target-wp2 cargo test --release -p ddg-repulsive-surfaces`
(`--test metric|constraints|gmres`; `-- --ignored preconditioner_spectrum` for the spectral diagnosis).

## Done

- `metric/mod.rs`: `Metric`, `L2Kind`, `SobolevWeights` (`tp_weight`, `low_order`, `willmore_weight`),
  `df_coefficients`/`df_face` (`D_f`, plain edges), `FaceOperators` (per-face `a`, `X`, `N`, `D_f`
  coefficients; `average`, `gradient`, `row(s, Channel)`), `Channel::{Average (U), Gradient (V)}`,
  `FaceDiagonal` (obstacle hook), `FractionalPairs::{Faces, Vertices}`. `DenseMatrix` is the curves crate's.
- `metric/dense.rs`: kernel matrices (`high_order_kernel`, `low_order_kernel` symmetrized,
  `fractional_kernel`), `kernel_matrix_form` (App. A, any symmetric zero-diagonal `H`; used for the
  production assembly), `kernel_matrix_apply` (3 columns; WP6's reference product), `pair_sum_matrix`
  (Alg. 1 literally, any kernel; test reference), `fractional_operator`, `high_order_operator`,
  `low_order_operator`, `vertex_fractional_operator`, `dense_metric` (B, B⁰, W, M, W M⁻¹ W, A),
  `add_face_diagonal`, `willmore_matrix`, `comparison_metric` (L² mass/identity, H¹, H²; sparse), `to_dense`.
- `constraints/{mod,rank}.rs`: `Constraint`, `BarycenterMode`, `RowKind`, `Row`, `ConstraintSet`
  (`new`, `analyze`, `values`, `jacobian(positions, rows)`, `area_gradient`, `volume_gradient`,
  `project_positional`, `positional_rows`, `nonlinear_rows`, `dropped_rows`, `tolerance`,
  `set_target`), greedy rank check (pins > barycenters > area/volume, `RANK_TOLERANCE` 1e-10),
  translation check.
- `solve/mod.rs`: `SaddleSolver` trait, `SaddleSolution`, `SolveRecord`, `LinearReport`.
- `solve/dense.rs`: `DenseSaddle` (equilibrated LU + one refinement step; per-axis LU of size
  `|V| + k_Q/3` when `C_Q = C₁ ⊗ I₃`, full LU otherwise; `SingularSaddle` below pivot ratio 1e-13),
  `DenseSaddleOperator` (`A′` product), `DensePreconditioner` (§6.3, `FractionalPairs` option),
  `preconditioner_epsilon`.
- `solve/schur.rs`: `SchurComplement` (Z, S = −C_N Z, equilibrated pivoted Cholesky of −S with row
  dropping at 1e-10, `direction`, `correction` = −Z S⁻¹ Φ), `NewtonSteps`, `ProjectionOptions`,
  `project` (positional, Newton rounds with the safeguard), `ProjectionStatus::{Converged, MaxRounds,
  NotImproving}`.
- `solve/gmres.rs`: `LinearOperator`, `Preconditioner`, `IdentityPreconditioner`, `FnOperator`,
  `GmresSide`, `GmresOptions` (30, 1e-4, 200, Right), `GmresReport`, `gmres`, `GmresSaddle`
  (SaddleSolver over any operator/preconditioner of `A′`).
- Tests: `tests/metric.rs` (10), `tests/constraints.rs` (7), `tests/gmres.rs` (10 + 1 ignored);
  fixture reader `tests/wp2_fixture/mod.rs` (a directory module, not a test target).

## Measurements

- Oracle fixtures (df, metric, constraints, saddle, precond): max relative error 2e-11 (`⟨g, x⟩` of
  `rs-icosa-saddle`), typically 1e-15. `rs-error-singular` (pin twice): the rank check drops the
  duplicate (rows 3–5), so no `SingularSaddle`.
- Schur direction/correction vs full saddle LU: 1e-13 or better; `C x`, `C h − Φ` at round-off.
- Kernel of the metrics: B per-component constants; **A with B⁰ only the global constants** on
  multi-component meshes (B⁰, L_τ couple all face pairs); smallest nonzero/largest eigenvalue 9e-4..6e-2.
- Preconditioner (GMRES iterations to 1e-4, right side; `tests/gmres.rs`):

  | mesh | faces (SPEC) right / left | vertices right / left | none |
  |---|---|---|---|
  | icosphere 320 | 38 / 17 (true res 4e-2) | 5 / 6 | 58 |
  | icosphere 1280 | 74 / 13 (true res 1e-1) | 5 / 4 | 269 |
  | torus 1500 | >200 / 19 (true res 2e-1) | 7 / 7 | >300 |
  | bumpy disk 1250 | stalls at 4e-2 / 18 (true res 2e-1) | 27 / 16 | >300 |

  Spectral diagnosis (`preconditioner_spectrum`): B vs `M Δ^s` cond 1.5–2.3 at every size; face-pair
  `L_{2−s}` vs `M Δ^{2−s}` cond 12.6 / 68 / 224 (icosphere 80 / 320 / 1280) — the face averages nearly
  annihilate mesh-scale oscillations; vertex-pair `L_{2−s}` cond 1.3–2.1, preconditioned metric ≈ 2.
  The reference's left-preconditioned stop hides this (true residual 0.1–0.2 at its tolerance).
- Safeguard on a near-sphere (icosphere 1280, area +3 %, volume kept, a step of one mean edge):
  equilibrated Schur pivot 1.1e-5 (not dropped), the first Newton round raises the violation ×12 →
  undone, `NotImproving`; unguarded rounds reach 7e71. One area row on a torus: converged in 3 rounds
  (2e6 → 9e3 → 82 → 0.75 in units of the tolerance).

## Decisions

- Production assembly of B, B⁰, L_τ by the kernel-matrix form (exactly App. A; O(|F|²)), Alg. 1 kept
  as `pair_sum_matrix` for tests (agree to 1e-12).
- Rank check: true linear dependencies only. A barycenter plus some pins on a component is independent
  and kept; `project_positional` then repeats translation + pin reset until the barycenter holds
  (`POSITIONAL_TOLERANCE` 1e-14 × coordinate scale, ≤ 100 rounds).
- `BarycenterMode::Global`: one global group; the translation check accepts it for every component
  (reference convention); dense solves on several components are then `SingularSaddle` only if B⁰ is off
  (with B⁰ the global constants are the only kernel — see above).
- Safeguard applies in both `NewtonSteps` modes; `Fixed(n)` takes n rounds unless a round is not improving.
- GMRES: convergence is decided on the explicit residual at the start of each cycle; `SolverDiverged`
  only when not converged and the true relative residual > 1e-1. Default `max_iterations` 200.
- Volume rows require `Consistent` components (also `Orientable` ones with the caller's order are refused).

## For WP5 / WP6

- One step: `set.jacobian(positions, &set.positional_rows())` → `DenseSaddle::new(pool, &A, &C_Q)`
  (or `GmresSaddle::new(op, prec, 3|V|, opts)`); `SchurComplement::new(&mut solver, set.nonlinear_rows(),
  C_N)`; `schur.direction(&mut solver, &g)` → `x` (step `f − τx`); after the line search
  `project(&set, Some(&schur), &trial, &opts)`. `schur.fill_report(&mut report)` + `direction.record`.
  Vectors are interleaved `3 v + axis`; `A′` vectors append the `k_Q` multipliers.
- WP6: implement `LinearOperator` for the hierarchical `A′` and `Preconditioner` for the sparse
  `L′`-based `P⁻¹` (input/output of length `3|V| + k_Q`); check products with `kernel_matrix_apply`
  (same `H`) and `DenseSaddleOperator`, the preconditioner with `DensePreconditioner`.
  For `FractionalPairs::Vertices` a vertex BVH (points `f_v`, masses `a_v`, kernel `ρ^{2s−6}`) is needed.
- Remeshing changes the vertex set: rebuild the `ConstraintSet` with explicit targets.

## Remaining / open (for the orchestrator)

- SPEC decisions proposed in the WP2 report: vertex-pair `L_{2−s}` (D32), the kernel statement of
  §4.4/§13.3, the D10 wording, GMRES `max_iterations`.
- `D7` (Willmore term in the preconditioner's Laplacian) is not in `DensePreconditioner` (WP6 measures it).
