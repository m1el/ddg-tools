# SOLVERFIX notes: robustness failures of NOTES-WP15 ("Failures and reproducers")

Branch: the worktree branch of this task, on top of `repulsive-surfaces` 51b1963. Build:
`CARGO_TARGET_DIR=$PWD/target-solverfix cargo build --release -p ddg-repulsive-surfaces --examples`.
Diagnostics: `examples/solverfix.rs` (`spectrum` = eigenvalues of the preconditioned saddle operator, dense
`A′` times `P⁻¹` column by column; `gmres` = residual histories of the hierarchical solve; `product` =
hierarchical vs dense `A′` at χ = 0; `cobflow` = the cobblestones flow with one line per step). Meshes
`h:L:pins|bary`, `torus:N:M:pins|bary`, `spindle:NM:NN:R`, `cobble:N:L[:pc]`. Logs and per-step CSVs in
`runs_tmp/solverfix/` (gitignored). The step CSV got two columns: `fallbacks` (escalations of the solver
chain) and `projection_blocked`. Regression tests: `tests/solverfix.rs` (6 tests, ≈ 12 s release).

## 1. Hard pins: GMRES diverged at step 0 (H of Fig. 22, 3840 faces, 24 pins)

**Not a product bug.** Hierarchical `A′` at χ = 0 against the dense one on this mesh: relative difference
2.0e-15. **Root cause: the scale of the augmented preconditioner (D33).** `P⁻¹ = L′⁻¹ [[L_{2−s}, C_Qᵀ],
[C_Q, 0]] L′⁻¹` reproduces the constraint rows exactly, so `P⁻¹A′` has the eigenvalue 1 on the `2 k_Q`
constraint directions, while on the rest it is `≈ W⁻¹ L_{2−s} W⁻¹ B`, a constant of about 200 for
`s = 5/3` (spectra, dense, χ = 0: torus 900 faces bulk 120–300, H 960 faces 32–3342; the unit cluster:
6 eigenvalues with a barycenter, 144 with the 24 pins). With a barycenter only six outliers do not matter
(WP6 measured that case); with 72 pin rows restarted GMRES(30) must resolve two clusters two orders of
magnitude apart and stagnated (true residual 0.965 after 200 iterations in the Schur column solve; the
fallback chain then hid it with restart 60). The constant is the planar symbol ratio
`4 I(σ) I(2−s)`, `I(t) = π Γ(1−t)/(t 4^t Γ(1+t))` (199.86 at `s = 5/3`, checked against quadrature;
`rs-preconditioner-scale`, cited in `lean/proofs/CITED.md`).

**Fix** (`solve/precond.rs`): `FractionalScale::Symbol` (default) divides `L_{2−s}` by
`w_TP · 4 I(σ) I(2−s)` in the augmented form (`SparsePreconditioner::with_metric_weight`; `Unit` keeps
WP6's form; the congruence form is unaffected, a global factor). Spectra after: torus with 3 pins
0.59–1.94, H 960 faces 0.16–16.7 (sharp creases; the pin-free H is the same).

| H, 3840 faces, pins + volume | Schur column + direction iterations | fallbacks / step |
|---|---|---|
| before (`Unit`) | diverged (0.965); retry restart 60: 210–248 its per step | 1 |
| after | 81 + 54; in the flow 135–161 its per step | 0 |

Test: `hard_pins_on_the_h_converge` (scaled converges in ≤ 120 per solve, `Unit` returns
`SolverDiverged`). `h_dense_volume` (all sharp-edge vertices pinned; WP15: 295 of 300 steps unconverged
at 200 iterations): see the table at the end.

## 3. Many components: ≈ 2800 GMRES iterations per step (cobblestones, 14 spheres)

**Root cause:** with one global barycenter (`BarycenterMode::Auto` = Global with `B⁰`) the components'
translations have no rows of their own, `L′ = [[W + ε_multi M, C_Qᵀ], [C_Q, 0]]` sees them only through
`ε_multi`, and `P⁻¹A′` has exactly `3 (n_c − 1)` outliers (39 for 14 spheres) between 1e2 and 1.4e5,
scaling as `ε_multi⁻²` (1.4e9 at `ε_multi = 1e-6`, while the rest lies in 0.36–1). Plus the scale of §1.

**Fix** (`solve/coarse.rs`, attached in `hierarchical_saddle`): a balancing (deflation) coarse correction
on `Z` = translations of the free components (+ the multiplier coordinates of the rows touching them),
`P_c⁻¹ = Z E⁺ Zᵀ + (I − Z E⁺ ZᵀA′) P⁻¹ (I − A′Z E⁺ Zᵀ)`, `E = ZᵀA′Z` (pseudo-inverse); `A′Z` from
`⌈n_c/3⌉` products (three components in the three axes of one product), kept dense (at most 2²⁵
entries, else skipped); `PreconditionerOptions::coarse` (default on, augmented form only). After: the 39
outliers are gone (spectrum 0.36–1.40, independent of `ε_multi`).

| 14 spheres, global barycenter + 14 areas | per solve | per step (15 solves) | solve time / step |
|---|---|---|---|
| before | ≈ 200 (cap; residual up to 1.6e-3, unconverged) | 2768–2963 | 13.8–17.6 s |
| scale only | 52 | 774 | |
| scale + coarse | 7–9 | ≈ 128 (first 200 steps) | |
| one sphere / per-component barycenters | 4 / 9 | | |

Test: `many_components_with_a_global_barycenter` (≤ 15 per solve with the coarse space, > 3× more without).

## 2. Divergence on self-intersecting inputs (subcritical step 0, Klein step 74, envelopes level 3)

Not one cause; three different ones, none a bug in how the operators are assembled:

* **Klein bottle** (one component, barycenter, no CCD): the failure at step 74 does not recur after the
  scale fix (120 steps, no fallback, ≈ 109 iterations per step). The full 700-step run with the default
  back end still breaks down later (step ≈ 184: a remeshing pass in the intersection band, where the
  collapse guard exempts pairs, raises `F` from 3.3e3 to 5.3e8; afterwards dense fallbacks every step,
  area lost by step 220, `SolverDiverged` at 218–225). With `Backend::Dense` the same 700 steps run clean
  (area kept to 6e-7): the trajectories separate through the GMRES/χ error, and the remeshing of an
  immersed surface is the fragile part. Open.
* **Subcritical spindle torus** (R = 1, r = 1.02): the "self-intersection" is a pinch: the ring `v = π`
  has radius 0.02 with 64 vertices, i.e. needle triangles of aspect ratio ≈ 100 around the axis. There the
  discrete `W` and the vertex-pair `L_{2−s}` model the face-based `B` poorly: `P⁻¹A′` reaches 4.3e4
  (1024 faces) against ≤ 16 for r = 0.9 (embedded, aspect ≈ 18). Softening the vertex-pair kernel
  (`|f_v − f_w|² + δ²`) lowers the top to 800 but GMRES on 4096 faces got worse; not adopted. Even the
  dense LU has residuals up to 2e2 late in WP15's 2000-step run: the metric itself is near-singular on
  this degenerating mesh. Remedy as before: the fallback chain (dense LU at ≤ 6000 vertices). After the
  fixes, 20 steps at 4096 faces with the default back end: (3.5, 7) 508 iterations per step, 4 fallbacks,
  18 steps above the 1e-4 tolerance; (3, 6) 1046 per step, 10 fallbacks in 11 steps. Open.
* **Envelopes / walnut / obstacles in general**: the mesh-obstacle term `O_diag` is in `A` but not in the
  preconditioner (§16 open item). An envelope starts 0.3 ℓ̄ from its obstacle, so `O_diag` dominates `B`
  for every wavelength above that distance; `P⁻¹A′` spreads by `(size/d)^{2s−2}`. WP15's level 1 ran all
  150 steps at the 400-iteration cap with residuals 4–9e-2 (below the 0.1 divergence threshold, so it was
  never reported as a failure); same after these fixes. A symbol analysis shows that adding `O` to `L′`'s
  Laplacian block, to the fractional middle factor, or both, does not narrow the spread (the low
  frequencies are off by `w ℓ/k` either way); it needs `(K L̃⁻¹ K + O)⁻¹`, i.e. an inner solve. Open.
  **Side effect of the scale fix here:** the primary solve now ends above 0.1 from step 3 on, so the
  retry stage (restart 60) runs: 1000 iterations per step at residual 4e-2 instead of 400 at 4–9e-2,
  ≈ 3× the solve time on envelope level 1 (16.8k faces, ≈ 30 s per step). Proposed: an adaptive scale
  (Rayleigh quotient of `P⁻¹A′` on a few vectors) or the obstacle term above; and escalation to the
  dense fallback when GMRES ends above 1e-2.

## 4. Remeshing explosion (`RS_REMESH=1 apps subcritical`)

**Root cause** (per-step CSV, `runs_tmp/solverfix/out/o_remesh`): not a needle loop of the remesher.
(1) Without CCD (subcritical runs turn it off), on the self-intersecting spindle, every remeshing pass
raises `F` (×1.2 typical, ×79 at step 5): `F` climbs from 74 to 1e4–1e5 by step 80 (a crumpled mesh).
(2) At steps 88–91 the Newton projection fails (`NotImproving`, or the area row dropped from the Schur
complement for being near-singular) and the line search accepts steps that leave the area constraint
(violation 7e-2 → 28 → 320 → 1e4: the area grows 250×). (3) The remesher keeps `L₀` fixed, so the face
count follows the area: 5.3k → 6.7k → 10.9k → 23.8k → 62.3k faces; (4) with `Backend::Dense` the
assembly at 31k vertices took 30 GB RSS and the next step asked for 118 GB.

**Fixes:** (a) the projection guard (§5) also requires that the projection leaves every nonlinear row
(kept or dropped by the Schur complement) at most `max(1, m_base)` tolerances from its target, `m_base`
the value at the step's base positions; otherwise the trial is halved (the curves crate's
`ProjectionFailed`). (b) `RemeshOptions::max_growth` (default 2): the splits of one pass may at most
double the face count; beyond, only the longest edges within the budget are split
(`RemeshReport::splits_deferred`). (c) `dense_metric` refuses more than `DENSE_MAX_VERTICES = 16384`
vertices with `RsError::InvalidOption` ("… GB per |V| x |V| matrix …; use Backend::Auto or
Backend::Hierarchical") before allocating. Tests: `remesh_growth_is_bounded`,
`dense_metric_refuses_huge_meshes`. (d) `RemeshOptions::guard_collapses` (default on): without CCD the
collapses are still checked by the CCD collapse guard (default gap, exempt band of the faces that
intersect at the start of the pass); it addresses (1). After (a)–(d), the same run
(`RS_REMESH=1 apps subcritical`, (3.5, 7)): `F` stays 48–66 for 100 steps (before: 1e4 by step 40), then
climbs to ≈ 1e4 by step 160 with 4–5k faces, the area constraint held to ≤ 1e-7 throughout, 198 steps
without explosion (before: 62k faces and 30 GB at step 93).

## 5. Obstacles breached by the projection (cobblestones)

**Mechanism confirmed.** CCD guarded only the line search; the Newton projection onto the growing area
targets moves the surface after it, unchecked. Baseline (51b1963 code, `RS_NO_REMESH=1`): Stuck at step
320, 152 vertices outside the slab, max |z| 1.387 (half 0.9); WP15: 218 with remeshing. Minimal case
(`projection_respects_obstacles`): a sphere 0.02 inside a weak slab, area target ×2 in one step: max |z|
0.672 unguarded against the slab at 0.52.

**Rule** (`LineSearch::guard_projection`, default on; `session.rs` `ProjectionGuard`): a trial that passes
Armijo is kept only if its projection is admissible: `F` finite after it (implicit obstacles, degenerate
faces), the constraint test of §4(a), and, with CCD on, the linear motion from the trial to the projected
positions passes the line search's CCD (self pairs with the band rule, and the CCD obstacles). Otherwise
`TrialOutcome::ProjectionBlocked` and τ is halved — unless the projection of the *unmoved* positions
meets a contact too (halving cannot help): then the trial is kept and only the admissible fraction of its
projection is applied (`⅔ τ_max` of the CCD along the projection's motion, halved while `F` is not
finite), `ProjectionStatus::Blocked`, `StepReport::{projection_blocked, projection_fraction}`; a stuck step
does the same with its own projection. The constraint then lags behind its target, reported in
`max_violation`, never hidden. Remeshing: the pass's motion (new vertices from their vertex-map positions
on the old surface to their new positions) is checked against the CCD obstacles (`IntersectionBand::all`
exempts self pairs); a pass that meets one is not applied this step (`StepReport::remesh_rejected`).

**Consequence, open: freezing at the obstacle.** With the guard the surface can no longer tunnel, and
the flow slows down where the descent direction itself pushes a vertex against the slab: the line
search's CCD cap is a single τ for the whole surface, so one vertex approaching the gap shrinks every
step geometrically (`cobflow`: τ_ccd 4.8e-2, 3.5e-2, 3.0e-3, …, 8.7e-9 at steps 41–61, vertex 673 against
slab triangle 172), while the area targets keep growing (violation 0.2 → 2.7 by step 220). The baseline
froze too (`Stuck` at step 320), but only after pushing 152 vertices through the slab. Making the CCD target
`min(gap, close_factor · d₀)` for every pair (tried) does not help (the approach is still geometric) and
was reverted. A real fix needs a contact-aware step: per-vertex limiting or the contacts as inequality
rows in the direction and the projection.

## Before / after (full runs)

`RS_AUTO_BACKEND=1` = default back end (hierarchical above 1000 faces); 4 threads, 9–13 jobs in parallel on
the 64-core machine (times not clean). "Before" = WP15's runs (showcase CSVs) or this task's baseline runs
on 51b1963.

| run | before | after |
|---|---|---|
| `h_sparse_volume` (auto), 30 steps | primary diverges every step, fallback retry; 210–248 its/step | 0 fallbacks, 202 its/step (2 solves), pins exact, volume ×1.1 on schedule |
| `h_dense_volume` (auto), 276 pins | WP15: 295/300 steps unconverged at 200 its | 30 steps: 0 unconverged, 0 fallbacks, 84.5 its/step |
| `cobblestones` (auto), 400 steps | 2768–2963 its/step, 30–60 s/step (WP15); breach: 218 vertices outside (dense, remeshing) | 134 its/step, 511 s for 400 steps (1.3 s/step); 0 vertices outside, max \|z\| 0.893 (half 0.9); area ×1.23 instead of ×2.5 (violation 3.9, guard blocked 359 projections) |
| `cobblestones` no remeshing (dense) | Stuck at 320, 152 outside, max \|z\| 1.387 | Stuck at 54 (CCD freeze, §5), 0 outside, max \|z\| 0.893, violation 0.22 |
| `klein` (auto), 700 steps | `SolverDiverged` at step 74 (0.29) | clean to step ≈ 183, then remeshing-induced blow-up and `SolverDiverged` at 218 (§2); dense: 700 clean steps |
| `subcritical` (auto), 20 steps | fallback at step 1, 1000–1400 its/step | same (pinch, §2) |
| `RS_REMESH=1 subcritical` | 62k faces + 30 GB at step 93, 118 GB request next | 198 steps, ≤ 5k faces, area ≤ 1e-7, F 48–66 to step 100, ≈ 1e4 by 160 |
| envelope level 1 (auto) | 400 its/step at 4–9e-2 (WP15, never converged) | 1000 its/step at 4e-2 (retry stage), ≈ 3× slower (§2) |

## Proposed SPEC edits

* **§6.3** (D33, after "The augmented form"): "The fractional factor is divided by `w_TP · 4 I(σ) I(2−s)`,
  `I(t) = π Γ(1−t)/(t 4^t Γ(1+t))` (`rs-preconditioner-scale`, ≈ 200 at `s = 5/3`): without it the
  augmented form's constraint part (eigenvalue 1) and fractional part (≈ 200) are two clusters, and
  restarted GMRES stagnates with many positional rows (24 pins: residual 0.965). `FractionalScale::Unit`
  keeps the unscaled form." And: "With a global barycenter on several components the preconditioner adds a
  balancing coarse correction on the translations of the components without rows of their own
  (`PreconditionerOptions::coarse`); `ε_multi` then no longer affects the iteration count."
* **§6.2** (new paragraph "Admissibility"): "A trial that passes Armijo is kept only if its projection is
  admissible: `F` finite after it; every nonlinear row at most `max(1, m_base)` tolerances from its target
  (`m_base` at the step's base positions; also rows dropped from the Schur complement); with CCD on, the
  linear motion from the trial to its projection passes the CCD of §7.3. Otherwise τ is halved; if the
  projection of the unmoved positions meets a contact too, the trial is kept and only `⅔ τ_max` of its
  projection is applied (`ProjectionStatus::Blocked`, `projection_fraction`). Constraints may then lag their
  targets; this is reported in `max_violation`. `LineSearch::guard_projection` (default on)."
* **§7.3**: replace "CCD guards the line search" by "CCD guards the line search, the projection (§6.2) and
  the remeshing (§10)"; note the open freeze: one vertex pushed against an obstacle caps τ for the whole
  surface.
* **§10**: "A pass splits at most `(max_growth − 1)·|F|` faces' worth of edges (default `max_growth = 2`;
  longest first beyond it, `splits_deferred`). Collapses are guarded also without CCD
  (`guard_collapses`, default on: default gap, exempt band). With CCD obstacles, a pass whose motion
  (vertex-map positions on the old surface → new positions) meets an obstacle is not applied
  (`remesh_rejected`)."
* **§8.4 / §6.3 Dense**: "`dense_metric` refuses more than 16384 vertices (`InvalidOption`), before
  allocating."
* **§16 open items**: the obstacle term of the preconditioner (envelopes: never converged, now 3× slower
  through the retry stage); contact-aware stepping (freeze at obstacles, lagging area in cobblestones);
  remeshing of immersed surfaces inside the exempt band (Klein bottle, default back end); needle/pinch
  meshes (subcritical spindle) for the vertex-pair preconditioner.
