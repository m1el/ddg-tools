# WP9 notes (comparison optimizers: AQP, H1 L-BFGS, BCQN)

Build/test: `CARGO_TARGET_DIR=$PWD/target-wp9 cargo test --release -p ddg-repulsive-surfaces --test optimize --test session --test flow`
(optimize: 9 tests, 1 s). Sanity run: `target-wp9/release/examples/compare_methods OBJ METHOD [STEPS] [THREADS]`.

## Done
- `src/optimize/`: `soft.rs` (`soft_problem`: Area -> `AreaDeviation` weight 100, Volume -> `VolumeDeviation`
  weight 100, positional constraints kept, pin schedules re-indexed; per-component area/volume and growth
  schedules -> `Unsupported`), `state.rs` (`MethodState`, `MethodDiagnostics`, `aqp_theta`, `two_loop_reference`,
  `spectral_norm`, `apply_k3`, `blend_beta`), `mod.rs` (docs).
- Session hook (small): `Session` holds `method_state: Mutex<Option<MethodState>>` (the objective borrows
  `self` during the step); `step()` calls `begin_step` (AQP: positions become the Nesterov point `y`),
  `end_step`/`abort_step`; in `step_inner` the direction is `MethodState::direction` instead of
  `SchurComplement::direction` when a method is set; the line search starts at `tau0 = 1` for quasi-Newton
  methods that hold pairs. `FlowOptions::effective()` sets `metric = H1` for every method but projected gradient
  (so the saddle is the sparse `[[W (x) I3, C_Q^T],[C_Q, 0]]`; `StepSolver::as_sparse`, `SparseSaddle::metric_matrix`
  added). `Session::method_diagnostics()`.
- Options: `Method::Bcqn { memory, variant: BcqnVariant }` (new enum), `RsError::MethodNeedsSoftProblem` gains
  `constraint: String` ("Area (constraint 1)"); `validate` rejects area/volume constraints for the three methods and
  `memory == 0`, `kappa < 1`. Scene defaults: memory 20 (was 10 for BQN).
- Tests: `tests/optimize.rs` (9) + 1 unit test (momentum resets), `tests/session.rs` adapted.

## Decisions
- AQP: gradient and H1 step at `y = (1+theta) x_k - theta x_{k-1}`, Armijo against `F(y)`; reset (prev := none)
  when `F(x_{k+1}) > F(x_k)` (energies of two steps' Barnes-Hut lists) or on remeshing. The iterate is kept after
  a reset (no rollback).
- L-BFGS standard: the curves crate's `Lbfgs` (`H0 = gamma W'^-1`, `gamma = s^T y / y^T W'^-1 y`, N&W 7.20);
  a rejected pair (`s^T y <= 1e-12 |s||y|`) clears the memory; a non-descent direction clears it and takes the H1
  gradient. Reference: own recursion with all inner products in `W` (no pair rejection, memory reset on non-descent).
- BCQN: `y_hat = (1-beta) y + beta W s`, `beta = clamp(norm(W) y^T W s / B, 0, 1)`, `B = A^{4/3}` in both variants
  (formula as coded by the reference; Zhu et al. not re-read, see below); standard uses the power-iteration spectral
  norm of `W` (<= 200 products, rel 1e-6), reference the Frobenius norm of the constrained Laplacian saddle
  (about 70x larger on the test mesh: beta saturates at 1 more often).
- Memory reset on splits/collapses only (the dimension changes); flips and smoothing keep the pairs (the
  reset-on-every-flip variant wiped the memory on 37-44 of 50 steps of the sanity run).
- Positional pins/barycenters stay exact: `s` and every direction lie in the nullspace of `C_Q`.

## Sanity run (EntangledDoubleTorus/doubletorus10k.obj, 10 506 faces, 50 steps, 4 threads each, 9 runs in parallel)
Remeshing and CCD on (defaults). `E` = unweighted tangent-point energy (alpha, beta = 6, 12) after the step. Hard =
barycenter + area constraints; soft = barycenter + area penalty (weight 100). The soft runs let the area grow by
26-29 % (E scales with the size, so E is not comparable across the hard/soft groups, only within them).

| method | problem | E at step 10 | step 25 | step 50 | seconds (50 steps) | faces | area / A0 |
|---|---|---|---|---|---|---|---|
| Hs (ours) | hard | 138.2 | 116.3 | 95.1 | 56.4 | 10634 | 1.000 |
| L2 | hard | 335.8 | 323.9 | 313.7 | 12.2 | 10506 | 1.000 |
| H1 | hard | 198.9 | 165.8 | 149.0 | 13.3 | 10508 | 1.000 |
| H2 | hard | 369.8 | 365.9 | 360.1 | 17.4 | 10506 | 1.000 |
| AQP | soft | 258.9 | 145.3 | 118.7 | 11.9 | 10820 | 1.269 |
| H1 L-BFGS | soft | 159.6 | 126.8 | 111.5 | 12.5 | 10970 | 1.294 |
| H1 L-BFGS (ref) | soft | 220.8 | 140.7 | 115.8 | 12.2 | 10858 | 1.279 |
| BCQN | soft | 231.5 | 145.9 | 116.9 | 13.5 | 10886 | 1.278 |
| BCQN (ref) | soft | 259.7 | 148.9 | 122.4 | 14.2 | 10748 | 1.256 |

Per-step logs: `runs_tmp/wp9/*.log` (not committed). Timings are with 9 runs sharing the machine, so only the
Hs row (hierarchical solves) is inflated by it. This is a 50-step sanity check, not the evaluation (WP16).

## Not verified
- Zhu et al. 2018 Eq. 13 was not re-read: the formula is taken from the reference's comment and code; `B = A^{4/3}`
  is the reference's value for d = 3. If the paper's `B(T)` or norm differ, change `blend_beta`/`direction` only.
- Kovalsky et al. 2016 uses a fixed rest-pose proxy, adaptive kappa and a flip-avoiding bound; D28 deliberately does not.

## Proposed SPEC edits
- §9.1: `Bcqn { memory: usize, variant: BcqnVariant }` with `enum BcqnVariant { Standard, Reference }`;
  `MethodNeedsSoftProblem { method, constraint }`.
- D26: "memory reset on splits and collapses (not on flips/smoothing), on a rejected pair, on a non-descent
  direction"; quasi-Newton methods start the line search at tau = 1 once pairs exist.
- D28: "reset on an energy increase between consecutive iterates (no rollback)".
- §14.2: soft problem helper `optimize::soft_problem`.
