# WP5 notes (session, step, line search, projection, schedules, remeshing, stop rules)

Build/test: `CARGO_TARGET_DIR=$PWD/target-wp5 cargo test --release -p ddg-repulsive-surfaces --test flow --test session --test oracle`.

## Done
- `src/step/`: `options.rs` (`FlowOptions`, `LineSearch`, `ArmijoRule`, `StartRule` (+ `Absolute`),
  `StopRule` (+ `Off`), `Schedule`, `PinSchedule`, `ScheduleTiming`, `Method`, `LbfgsVariant`,
  `effective()` = `reference_compat`), `objective.rs` (frozen BH lists + per-obstacle
  `CrossEvaluator`, Coulomb/Willmore energies, obstacle diagonal, Willmore weight),
  `direction.rs` (`StepSolver`: dense / hierarchical (reuses the BH face tree) / `SparseSaddle`
  for L2/H1/H2), `session.rs` (`Session`, `flow`, reports), `scene_run.rs` (`scene_options`,
  `run_scene`, CSV log).
- `tests/oracle.rs`: the `step` task is replayed through `Session` (fixtures and fuzz).
- `RsError::Internal { what, detail }`.
- §16 unification: `FractionalPairs::{Faces, Vertices}` → `FractionalForm::{FacePairs, VertexPairs}`
  (default `VertexPairs`) in every file that used it; `DEFAULT_DENSE_MAX_FACES` 2000 → 1000 (§6.3).

## Decisions
- `barycenter_mode: Option<BarycenterMode>`, `None` = Auto (SPEC §9.1 lists the enum without
  `Auto`; the constraints module was not edited).
- `F = w_TP ℰ + Σ penalties` and the metric × `w_TP` (the reference's reading of D6; the oracle
  weights only the metric — the replay rescales `τ`).
- Constraints are rebuilt every step from explicit targets (per-component barycenters expanded),
  so remeshing only remaps components (through face parents) and pins.
- CCD (WP8): `LineSearch::ccd: Option<CcdOptions>`; `first_contact` with `velocity = −x`, `τ_end = τ₀`;
  `τ₀ ← min(τ₀, ⅔ τ_max)` only when a contact was found; collapse guard during remeshing;
  obstacles from `FlowOptions::ccd_obstacles`.
- Stuck steps still project, advance schedules and remesh (SPEC §7.2).

## Remaining
- see the bottom of this file (updated as work goes on).
