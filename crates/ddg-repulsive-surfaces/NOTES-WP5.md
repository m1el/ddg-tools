# WP5 notes (session, step, line search, projection, schedules, remeshing, stop rules)

Build/test: `CARGO_TARGET_DIR=$PWD/target-wp5 cargo test --release -p ddg-repulsive-surfaces --test flow --test session --test oracle`
(release ~10 s; debug ~1 min). Smoke runs: `scripts/rs-fetch-scenes.sh`, then
`target-wp5/release/examples/run_scene data/repulsive-surfaces/scenes/<Scene>/scene.txt --steps N --threads 4 [--ccd] [--no-stop] --csv X.csv --obj X.obj`
(one summary line at the end; logs/CSVs of the runs below in `runs_tmp/wp5/{smoke,timing}/`, not committed).

## Done
- `src/step/`: `options.rs` (`FlowOptions`, `LineSearch`, `ArmijoRule`, `StartRule` (+ `Absolute`),
  `StopRule` (+ `Off`), `Schedule`, `PinSchedule`, `ScheduleTiming`, `Method`, `LbfgsVariant`,
  `effective()` = `reference_compat`), `objective.rs` (frozen BH lists + per-obstacle `CrossEvaluator`,
  Coulomb/Willmore energies, obstacle diagonal → `FaceDiagonal` (`d = c·a/2`), Willmore weight),
  `direction.rs` (`StepSolver`: dense / hierarchical (reuses the BH face tree) / `SparseSaddle` for
  L2/H1/H2: sparse LDLT of `K + εM` (ε = 1e-10 tr K/tr M for H1/H2) + refinement against the exact
  saddle), `session.rs` (`Session`, `flow`, `StepReport`, `FlowReport`, `Status`, `Control`, `Trial`,
  `Timings`), `scene_run.rs` (`scene_options`, `run_scene`, `CSV_HEADER`, `csv_line`).
- `examples/run_scene.rs`; `tests/flow.rs` (12), `tests/session.rs` (7); `tests/oracle.rs`: the `step`
  task replayed through `Session` (7 fixtures, 0 disagreements, worst ratio 0.2 of tolerance; fuzz
  cases too, rescaled for `w_TP`, `w_W > 0` cases skipped).
- `RsError::Internal { what, detail }`; lib.rs re-exports the §9.1 API.
- §16 unification: `FractionalPairs::{Faces, Vertices}` → `FractionalForm::{FacePairs, VertexPairs}`
  (default `VertexPairs`) in every file that used it; `DEFAULT_DENSE_MAX_FACES` 2000 → 1000 (§6.3).

## Decisions
- `barycenter_mode: Option<BarycenterMode>`, `None` = Auto (Global with B⁰, PerComponent otherwise and
  for the comparison metrics); the constraints module was not edited.
- `F = w_TP ℰ + Σ penalties`, metric × `w_TP` (the reference's D6). The oracle weights only the metric:
  same trajectory with τ scaled by `1/w_TP`.
- Constraints rebuilt every step from explicit targets (per-component barycenters expanded); remeshing
  remaps components through face parents (error `Internal` if inconsistent) and pins through
  `Remeshed::pinned`; squared-error targets `transport_vec3`, attractor sets `transport_set(½)`;
  `Session::face_origin()` composes face parents back to the initial faces.
- Schedules: target(k) = T₀ + min(k, steps)/steps (final − T₀); `BeforeProjection` uses k,
  `AfterProjection` (reference) k − 1. Convergence is not declared while a schedule moves.
- Barycenter shift: global mean (Global mode) or per-component means added to x; the barycenter targets
  follow the line-search result (not projected).
- CCD (WP8): `first_contact(velocity = −x, τ_end = τ₀)`; `τ₀ ← min(τ₀, ⅔ τ_max)` only when a contact was
  found (else `τ_max = τ₀` would always shrink the start); collapse guard in remeshing.
  **Default ON** (D20 measurement below); `scene_options` turns it off for subcritical exponents.
- Stuck: `τ < τ_min = 1e-12 ℓ̄₀/‖x‖_∞` or `max_halvings`, or `⟨g,x⟩ ≤ 0`; projection/schedules/remesh
  still run; `stuck_limit` consecutive → `Status::Stuck`.
- `StopRule::SobolevNorm`: `√(⟨g,x⟩_k/⟨g,x⟩_1) < eps`.
- `reference_compat` sets D8, D11, D12, D13 (+ face pairs D32, printed form D33, ε_multi 1e-2), D15,
  D16, D21, D25; not reproducible: D9 3× volume row, D14 absolute ε, absolute τ_min 1e-20, D24 values.

## Measurements (2026-11-02, release)
Smoke (10 runs in parallel, 4 threads each, `--no-stop`); timing: 30 steps each run alone, 4 threads.

| scene | faces → | steps | F start → end | ms/step (alone) | GMRES its/step | split/coll/flip (total) |
|---|---|---|---|---|---|---|
| HopfLink | 9216 → 9078 | 200 | 847 → 308 | 563 | 32 | 0/69/4026 |
| EntangledDoubleTorus | 10506 → 10552 | 200 | 377 → 8.25 | 609 | 21 | 947/924/10423 |
| KleinBottle (3.5, 7) | 1800 → 1546 | 200 | 793 → 156 | 147 | 46 | 79/206/983 |
| HandcuffsPole | 10688 → 12288 | 200 | 6.5e4 → 6.9e3 | 703 | 31 | 2224/1424/8608 |
| Bunny | 28578 → 34298 | 100 | 1.6e5 → 56 | 2357 | 34 | 4659/1798/5545 |

No step stuck; final meshes intersection-free except the Klein bottle (an immersion: 10 pairs; CCD run 45).
Step breakdown (Bunny, ms): factor 475 (BCT + L′ + fractional tree), solve 1326 (GMRES, 3 solves), energy
154, remesh 145, projection 53, line search 48, BVH 38. CCD share of a step (median, alone): 0.5–2.0 %
(Hopf 1.0, EDT 1.8, Handcuffs 2.0, Bunny 1.9, Klein 0.5). With CCD the EDT flow is slower (F 42 vs 8.2
after 200 steps, τ capped ≈ 0.5 vs 3.7–7.4; no run intersected at any step end) and the Klein bottle
stalls (F 791 after 30 steps, its immersion has touching pairs).

## Remaining
- WP9 plugs the comparison optimizers into `Method` (Session::new returns `Unsupported` now).
- Normal constraints of the scene format (`constrain normals`) are `Unsupported`.
- A metric-only Willmore term (oracle `w_W` without a Willmore energy) has no option.
