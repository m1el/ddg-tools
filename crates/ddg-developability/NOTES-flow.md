# Flow driver, cut and flatten: progress notes

Done (SPEC §5.1, §6, §7.1, §10.7, §10.8): `src/flow.rs`, `src/cut.rs`, `examples/develop.rs`
(CLI: OBJ in, flow, OBJ + per-vertex energy out, optional cut and flatten with UVs),
`examples/study.sh` (the parallel study below), `tests/flow.rs`, `tests/cut.rs`, `benches/flow.rs`.

## How the pieces behave

* `flow`: per round subdivide (pinned ids by `Subdivided::pinned_after`, new vertices get
  `vertex_source = NONE`), remesh once, then `optimize` in segments. The progress callback checks
  every `check_interval` accepted iterations for triangles with a corner angle below the small
  angle; if there are more than the count that was left after the last remesh (those are
  unfixable by the legality rules), it breaks, remeshes, and restarts the optimizer with the
  remaining budget of the round (fresh L-BFGS memory). A segment that ends with `LineSearchFailed`
  or `NonDescent` also triggers one remesh and a restart if the remesh changed something.
  At the end of the budget a pending check is still remeshed, so the returned mesh has no fixable
  bad triangle. `RoundReport.optimize` is the merge of `segments` (histories concatenated, a
  restart's first entry is the energy right after the remesh). `RemeshReport.energy_before/after`
  are filled in.
* `FlowOptions::threads = Some(n)` runs the flow in a local pool; a user progress callback (not
  `Send`) stays on the calling thread and is served over channels (positions are copied per call).
  Results are bitwise independent of the thread count (tested with 1 and 3).
* `Progress.iteration` is cumulative within the round over restarts; `optimize.pinned` must be
  empty (error), use `FlowOptions::pinned`.
* `cut_and_flatten` follows §6. Closed genus-0 components with fewer than 2 terminals get their
  highest-energy vertices added up front; a retry loop on `BffError::SphereNeedsTwoCones` covers
  closed pieces created by `must_cut`. Faces with zero plane area get `face_log_scale = NaN`
  (excluded from the vertex mean and from `mean_abs_log_scale`); flipped faces use |A_uv| and are
  counted in `charts[i].distortion.flipped_triangles`. Extra output beyond the SPEC struct:
  `mean_abs_log_scale` (area-weighted mean of |u_f|), helpers `n_charts`, `flipped_triangles`,
  `is_injective`. `CutOptions` has `length_weight` (c_len), `alignment_weight` (c_align, default 0),
  `extension`, `overlap_samples_per_triangle`, `threads` in addition to the SPEC fields.
* No new `DevError` variant was needed.

## Tests

`cargo test -p ddg-developability --test flow --test cut` (fast, < 2 s in release).
Long ones: `cargo test --release -p ddg-developability --test flow --test cut -- --ignored --nocapture`
(about 15 s in release). Exact/analytic checks: a flat disk, a tube and an annulus with seams unroll
isometrically (mean |u| < 1e-6, 0 flips); `c_align` makes the tube cut follow a ruling; torus cut has
2 loop edges; closed spheres get 2 added terminals; noisy cylinder and cone sheets: E_λ falls from
about 2 to 5e-7 in 600 iterations; icosphere: top 20 % of vertices hold 78 % of the energy after
400 iterations (22 % before); subdivision preserves the round-to-round energy to 1e-9 for all six
energies; remeshing in a flow removes an injected needle and keeps `vertex_source` consistent.
Models (release, `--ignored`): crumpled E_λ 61.9 -> 9.1e-4 in 500 iterations (-99.99 %), dihedral
angles above 30 degrees 354 -> 35 of 800 edges, mean |angle defect| 0 -> 0.025 (the flow does not
preserve isometry; it only makes normals coplanar per star).

## Findings

* **(Resolved by D30, see "Cut study" below.)** **Spec expectation not met on spot.** SPEC §10.7 expects at most 0.1 % flipped UV triangles after
  flow + cut on spot. Measured (control mesh, 2 subdivisions, E_λ, ε = 1e-3): 300 iterations per
  round: 99 of 5112 (1.9 %), mean |u| 0.90; 1000 per round: 86 of 4568 (1.9 %), mean |u| 3.4.
  Lower ε helps (300 iterations: ε = 1e-4: 25 flips, 1e-5: 6 flips = 0.12 %, mean |u| 0.085, with
  1086 terminals). Cause (measured with Σ|Ω|): the flow leaves a lot of canceling positive and
  negative angle defect at vertices whose energy is below ε (spot, 300 it: Σ|Ω| 80 at terminals,
  138 at the other interior vertices). E_λ is quadratic in the normal spread while the defect is
  not, so a small energy does not mean a small defect. BFF with zero boundary scale must absorb
  that as area distortion. The tests assert <= 5 % for spot/bob/blub (blub 15 %) and print the
  numbers.
* **(Resolved by D30.)** **Over-running the flow gives spikes.** Spot, 3 rounds, 3000 iterations each: E_λ = 1e-8, no
  terminal above ε, and the cut/flatten has 195 flipped triangles and mean |u| 6. The minimizers
  include non-embedded stars (spikes, double covers; SPEC D28, paper Fig. 9). The same shows on
  blub (E < 1e-3 everywhere after 2 x 500 iterations but Σ|Ω| = 99). Stopping early matters; an
  embeddedness check would be a feature request (out of scope).
* Remeshing is active all the time: 10 to 60 collapses per round on spot/blub/bob control meshes
  (needles keep forming), flips are rare. Face counts drop 10 to 15 % after the first round on
  control meshes. No stalls and no line-search failures were seen in the study with the 1 degree
  threshold and the 0.5 step cap; a tuning was not needed.
* Crumpled developable with E_λ: `cut_and_flatten` finds no terminals after the flow (E_λ ~ 6e-4),
  one chart, 0 flips, mean |u| 0.01.

## Cut study (D30): flipped UV triangles after flow + cut

Question: why does the paper's terminal rule (`E_λ > ε`) leave 1.9 % flipped UV triangles on spot
(SPEC §10.7 expected 0.1 %), 5 % when the flow is over-run, and what fixes it at what cut length.

**Hypothesis (non-embedded stars carry the defect): confirmed.** `star_embeddedness` (smallest
cosine between a face normal and the vertex normal; <= 0 means the star is not a graph over its
tangent plane) on spot control + 2 subdivisions, `E_λ`, ε = 1e-3, paper rule:

| run | terminals | Σ\|Ω\| at terminals | Σ\|Ω\| at non-terminals | of which on non-embedded stars | non-terminals with \|Ω\| > 0.1 (non-embedded) | flipped | mean \|u\| |
|---|---|---|---|---|---|---|---|
| 3 x 300 it | 78 | 80.4 | 138.2 | 134.9 | 87 (86) | 99 (1.9 %) | 0.90 |
| 3 x 3000 it | 0 (+2 added) | 9.4 | 222.2 | 222.2 | 77 (77) | 195 (5.0 %) | 6.08 |

Non-embedded stars are 14 % of the interior vertices (360 of 2558), yet hold 98 % of the defect of
the non-terminals; only 68 % of the flipped triangles touch one (the distortion spreads through
the chart), so the flips are not located exactly at the spikes, but the defect that causes them is.

**Study.** `examples/cut_study.sh` (flows once per mesh x flow energy x budget, written to OBJ, then
`develop --cut-only` for 3 terminal energies x ε in {1e-2, 1e-3, 1e-4} x `defect_tolerance` in
{off, 0.1, 0.02}: 1188 cut jobs, `xargs -P 56`, 3 min wall), `examples/cut_study.py` aggregates
(`runs_tmp/dev-cut-study/{summary.md,all.csv}`). Flows: spot, bob, blub control meshes and
CrumpledDevelopable (control: 2 subdivisions, 300 and 3000 it/round; crumpled 500 and 3000), seven
dense meshes (1 subdivision, 300 it/round), flow energies `E_λ`, intrinsic, `E_P`. 4 of 45 flows are
unusable for the cut or the flow failed (see "Problems found").

Control meshes, flow `E_λ`, ε = 1e-3. Cells: terminals / flipped UV triangles / mean |u| / cut length
as % of the total edge length.

| mesh, it/round | paper: E_λ, off | intrinsic, off (a) | E_λ + tol 0.1 (b, default) | intrinsic + tol 0.1 (c) | E_λ + tol 0.02 |
|---|---|---|---|---|---|
| spot 300 | 78 / 99 (1.94 %) / 0.90 / 2.06 | 309 / 6 (0.12 %) / 4.2e-2 / 4.60 | 165 / 0 / 4.0e-3 / 4.06 | 329 / 0 / 2.4e-3 / 4.76 | 301 / 0 / 1.0e-3 / 4.93 |
| spot 3000 | 2 / 195 (5.0 %) / 6.08 / 0.13 | n/a (energy undefined) | 79 / 0 / 1.3e-3 / 4.05 | n/a | 89 / 0 / 6.2e-5 / 4.13 |
| bob 300 | 255 / 0 / 4.9e-2 / 2.18 | 339 / 0 / 1.1e-2 / 2.77 | 276 / 0 / 1.7e-2 / 2.37 | 345 / 0 / 1.0e-2 / 2.80 | 352 / 0 / 6.1e-3 / 2.70 |
| bob 3000 | 2 / 145 (1.75 %) / 0.59 / 0.96 | 115 / 1 / 0.14 / 2.18 | 90 / 0 / 4.6e-3 / 2.10 | 175 / 0 / 3.3e-3 / 2.49 | 168 / 0 / 9.3e-4 / 2.56 |
| blub 300 | 2 / 140 (4.88 %) / 2.49 / 0.76 | 165 / 8 (0.28 %) / 0.10 / 5.10 | 72 / 0 / 6.2e-3 / 3.73 | 194 / 0 / 5.1e-3 / 5.54 | 129 / 0 / 8.7e-4 / 4.79 |
| blub 3000 | 2 / 158 (5.67 %) / 7.93 / 0.13 | n/a | 59 / 0 / 1.8e-3 / 4.05 | n/a | 64 / 0 / 6.3e-6 / 4.28 |

Summary over all usable flows (worst case over meshes; "cut" is the mean over meshes):

* control meshes, flow `E_λ`, 300 it: paper rule ε = 1e-2 / 1e-3 / 1e-4: worst 4.9 / 4.9 / 4.9 % flipped,
  cut 0.8 / 1.7 / 4.0 %. Tolerance 0.1: 0 flipped for every ε, worst mean |u| 3.3e-2 / 1.7e-2 / 5.8e-3,
  cut 3.1 / 3.4 / 5.3 %. Tolerance 0.02: 0 flipped, worst mean |u| 1.1e-2 / 6e-3 / 8.6e-4,
  cut 4.0 / 4.1 / 5.8 %.
* control meshes, flow `E_λ`, 3000 it: paper rule 5.7 % flipped, mean |u| up to 7.9, cut 0.4 %;
  tolerance 0.1: 0 flipped, mean |u| <= 4.6e-3, cut 3.4 %.
* control meshes, flow intrinsic or `E_P`: the paper rule is worse for coarser ε: worst flipped
  0.1 to 0.6 % at ε = 1e-4, 0.4 to 3.4 % at 1e-3, 0.5 to 5.5 % at 1e-2 (worst with `E_λ` terminals
  on the intrinsic flow); tolerance 0.1: 0 flipped in every case.
* dense meshes (seven, 23k to 57k faces, three flow energies): 0 flipped in every variant;
  tolerance 0.1 changes neither mean |u| nor the cut length (3.88 % -> 3.89 % at ε = 1e-3); mean |u|
  is set by ε alone (medians 9e-5 at ε = 1e-4, 8e-3 at 1e-3, 0.2 at 1e-2).
* CrumpledDevelopable: 0 flipped everywhere; the tolerance adds one terminal (1 cut edge, 0.17 %
  of the length) and lowers mean |u| from 1e-2 to 1.4e-4 (500 it) or 1e-6 (3000 it).
* The terminal energy matters little once the tolerance is set (`E_λ`, intrinsic and `E_P` give the
  same 0 flips and mean |u| within a factor 2); intrinsic costs 20 to 30 % more cut length.

**Decision (D30):** `CutOptions::energy` stays `Covariance`, new `defect_tolerance` defaults to
`Some(0.1)` rad (`None` = the paper's rule). Reason: it removes all flips on all control meshes and
over-runs, it is free on dense meshes, it cannot error (the intrinsic energy is undefined when a
face normal is antipodal to the vertex normal, which happens exactly on the spiky stars), and the
cut stays at 2.1 to 4.1 % of the edge length on control meshes (the paper rule is 0.1 to 2 % there
only because it leaves the defect uncut). Tolerance 0.02 is the choice for lowest distortion
(mean |u| below 1e-3 on spot and blub at +20 % cut length).

Problems found in other modules: fixed, see "Robustness study" below (D31, D32).

## Robustness study (D31, D32, D33)

Two failures of the cut study, both found in the flow itself. Logs: `runs_tmp/dev-robust/`
(`repro/`, `sweep/`, `sweep2/`, `over/`, `study/`); scripts `examples/robust_study.sh`,
`robust_summary.py`, `robust_compare.py`.

**1. Degenerate faces with `E_P` (and the intrinsic energies).** Mechanism, measured:

* The energies are scale free: squeezing a face to zero area costs (almost) nothing. Bare
  `optimize` (no remeshing) on small meshes (perturbed icospheres level 1 and 2, tori, bumpy grids,
  6 seeds, `E_P` and intrinsic, 800 iterations): 40 of 48 runs end with a face of angle
  1e-9 to 1e-12; the weak Wolfe search accepts every such step because the energy goes down a bit.
* In the flow the remeshing does not stop it: the blub control mesh with `E_P`, round 0, iteration
  90: one needle appears, its collapse is refused ("normal turn of the merged star"), the face is
  counted as `known_unfixable`, and from then on the flow only re-checks for *more* bad faces than
  that. The unfixable face shrinks from 1.2e-2 rad through 1e-4 (iteration 120), 1e-7 (iteration
  157) to 1e-9 (iteration 175) while the energy changes in the 7th digit. Refusal reasons over the
  three failing flows (29 + 1, 19 + 16 + 40 + 5, 4 + 3 + 8 + 18 + 8 + 8 refused operations):
  collapse "normal turn of the merged star" and "link condition", flip "normal turn", "degenerate
  new face", "smallest angle does not increase", "opposite vertices adjacent": all legal refusals,
  not bugs, but they leave faces of angles (0, 0, 180) and (0, 90, 90) in place.
* The final faces have angles like (2e-8, 180, 1.8e-8) degrees-scale; `Geometry::from_positions`
  rejects them in the cut ("violates the triangle inequality", Heron with rounding).

*Fix (D31):* admissibility floor `min_face_angle = 1e-4` rad in the optimizer (an inadmissible
trial is `+∞` like an undefined energy; the starting point is an error), the same floor as a collapse
rule of the remeshing, and the vertices of faces within a factor 2 of the floor are blocked
(gradient and direction zeroed, like pinned) so that one squeezed face does not stop the whole
optimization. The first version of the floor (no blocking) made the line search fail at the wall:
spot control mesh, intrinsic, 300 it/round, ended after 164 iterations at E = 16.7 (budget 900;
E = 6.2 before the floor), and the next rounds made 0 iterations; with the blocking the same flow
reaches E = 5.05 in 900 iterations. The default flow (E_λ, 2 rounds x 500 it) on spot, bob and blub
(triangulated and control meshes) runs the whole budget (`runs_tmp/dev-robust/default/`). Cost of
the floor on the energy: on 24 small runs the final energy with floors 1e-5, 1e-4, 1e-3 is within
1 % of the unconstrained run in most cases (see SPEC D31).

**2. Intrinsic energy undefined (blub, 3000 it/round).** The error is raised at the start of
round 1: `evaluate_energy` of the *subdivided* mesh. Mechanism: the end state of round 0 (93
vertices, E = 0.79) is a flat pancake, many edges with two big faces of nearly opposite normals
(normal cosines -1.000000000, -0.999999486, -0.999998921 ..., areas 0.45 / 0.41); the midpoint of
such an edge has two (planar-arc) half-stars with nearly opposite normals, vector area nearly
zero, so D12 applies there although it did not at the old vertices. `E_λ` flows reach the same
pancakes (E = 1e-11, smallest fold gap 4e-10, spot / blub / bob covariance 3000). So the line search
treats undefined trials as `+∞` correctly; the undefined state appears only in the subdivision
(an operation outside the line search). Candidates measured:

* A second wall on the fold angle (`min_fold_angle` = 1e-3, 1e-4, 1e-5, 1e-6 rad). It does not
  work: a hard wall in L-BFGS ends the run when one edge reaches it (`LineSearchFailed`; blub `E_λ`
  with 1e-3: E = 1.4 after 0.05 s, unconstrained 2.7e-11; spot `E_λ` 1e-5: 0.021, unconstrained
  5.5e-8), and even at 1e-6 spot intrinsic-max still hit the undefined state (sweeps in
  `runs_tmp/dev-robust/sweep/`, `sweep2/`). Removed again.
* Decision (D32): evaluate the energy of the next round on the candidate state before committing. If
  it is undefined, the flow does not run that round, returns the state of the previous round and
  reports `FlowReport::stopped = Some(FlowStop { round, after_subdivision, vertex, reason })` (an
  error for round 0). A remeshing result with an undefined energy or a face below the floor is not
  applied (`RemeshReport::rejected`). Nothing is repaired.
* In the over-run set below this happens in 2 of 18 flows (spot intrinsic: before round 2; spot
  intrinsic-max: before round 1); the other 16 complete all rounds. Blub intrinsic 3000 now
  completes (E = 1.07, was an error).

**3. Found by the oracle fuzz (D33):** a boundary star with at most two faces contributed noise to
the `E_λ` gradient (the eigenvector of a rank-deficient matrix is arbitrary); it now contributes
exactly nothing.

### Over-runs: 3 rounds of 3000 iterations, spot / bob / blub control meshes, six energies

All 18 flows complete without an error; none has a face below 1e-4 rad (`min angle`), the result of
each builds a `Geometry`. "faces < 1 deg" are the thin faces left by the remeshing (reported as
`FlowReport::bad_faces`). Before (cut study binary): blub `E_P`, spot `E_P`, bob `E_P`/intrinsic left
degenerate faces (angles 3e-10 to 2e-11 rad), blub intrinsic and spot intrinsic-max failed with
errors. `min fold gap` (pi minus the largest angle between adjacent normals) is a measurement, not a
constraint. Seconds are single-thread under a machine load of 60 to 150.

| run | faces | status | E_final | s | inadmissible | min angle | faces < 1 deg | min fold gap |
|---|---|---|---|---|---|---|---|---|
| blub_ctl-combinatorial-max | 3324 | max/max/max | 2.270763e0 | 20.905 | 4 | 2.13e-02 | 0 | 1.12e-03 |
| blub_ctl-combinatorial | 2144 | stalled/stalled/stalled | 1.129796e0 | 0.923 | 54 | 3.30e-02 | 0 | 3.41e-08 |
| blub_ctl-covariance-max | 3012 | max/max/max | 2.084945e-2 | 43.724 | 14 | 1.10e-04 | 9 | 4.17e-08 |
| blub_ctl-covariance | 2832 | max/lsf/max | 9.052823e-9 | 11.960 | 0 | 1.44e-04 | 9 | 9.51e-10 |
| blub_ctl-intrinsic-max | 3192 | max/max/max | 1.532855e0 | 31.878 | 16 | 1.15e-04 | 3 | 1.98e-04 |
| blub_ctl-intrinsic | 2860 | stalled/lsf/lsf | 1.070741e0 | 0.571 | 2 | 3.03e-02 | 0 | 3.13e-03 |
| bob_ctl-combinatorial-max | 10192 | max/max/max | 2.529696e0 | 86.508 | 0 | 2.00e-02 | 0 | 5.43e-01 |
| bob_ctl-combinatorial | 8596 | max/max/max | 4.286885e-1 | 86.455 | 2 | 1.90e-04 | 5 | 1.79e-01 |
| bob_ctl-covariance-max | 10462 | max/max/max | 4.346150e-1 | 123.544 | 0 | 1.93e-02 | 0 | 1.85e-05 |
| bob_ctl-covariance | 8366 | max/max/max | 2.417610e-2 | 71.141 | 3 | 1.11e-04 | 1 | 1.14e-06 |
| bob_ctl-intrinsic-max | 10152 | max/max/max | 2.002256e0 | 83.406 | 0 | 2.57e-02 | 0 | 1.11e-02 |
| bob_ctl-intrinsic | 9040 | max/lsf/stalled | 2.444759e0 | 49.355 | 13 | 1.18e-02 | 2 | 3.26e-05 |
| spot_ctl-combinatorial-max | 5442 | max/max/max | 1.008967e1 | 38.432 | 8 | 1.20e-04 | 2 | 2.18e-04 |
| spot_ctl-combinatorial | 3952 | max/max/max | 1.037054e0 | 35.798 | 47 | 1.03e-04 | 19 | 2.52e-05 |
| spot_ctl-covariance-max | 4790 | max/max/max | 2.640672e-1 | 61.758 | 0 | 1.78e-02 | 0 | 2.66e-08 |
| spot_ctl-covariance | 3808 | max/max/max | 7.550502e-4 | 21.302 | 20 | 1.05e-04 | 11 | 7.10e-09 |
| spot_ctl-intrinsic-max | 324 | lsf | 1.315386e1 | 1.728 | 0 | 4.80e-02 | 0 | 1.64e-03 |
| spot_ctl-intrinsic | 1260 | max/lsf | 3.322428e0 | 7.378 | 19 | 1.95e-04 | 2 | 4.56e-07 |

### Before / after on the flows of the cut study (45 flows)

Same flows as `examples/cut_study.sh` (control meshes 300 and 3000 it/round, dense meshes 300
it/round), then the cut with covariance terminals, eps 1e-3, defect tolerance 0.1 (D30). Old =
`runs_tmp/dev-cut-study`, new = `runs_tmp/dev-robust/study`. "flipped uv / status": flipped UV
triangles after the cut, then the optimizer status per round (lsf = line search failed). All 45
flows are listed; the dense-mesh flows (300 it/round) are identical in energy, iterations and faces,
the crumpled ones differ in the 4th digit or by 2 faces (the floor changes a few accepted steps).
Seconds differ because of the machine load (about 60 to 150 during the new runs, less during the
old), not because of the code (the check is O(|F|) per evaluation, measured below 1 % of an
evaluation). The cut succeeds on all 45 (3 failed before; the 4th failure was the flow itself) and
flips 0 UV triangles everywhere.

| flow | E_final old | E_final new | iterations old | new | s old | s new | faces old | new | flipped uv / status old | flipped uv / status new | new min angle | new bad faces |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| CrumpledDevelopable-combinatorial-3000 | 4.404057e-6 | 4.404057e-6 | 3000 | 3000 | 1.787 | 1.835 | 507 | 507 | 0 / max | 0 / max | 5.884e-2 | 0 |
| CrumpledDevelopable-combinatorial-500 | 1.478583e-3 | 1.478583e-3 | 500 | 500 | 0.322 | 0.322 | 510 | 510 | 0 / max | 0 / max | 4.100e-2 | 0 |
| CrumpledDevelopable-covariance-3000 | 9.034545e-8 | 1.189913e-7 | 3000 | 3000 | 1.262 | 1.368 | 497 | 497 | 0 / max | 0 / max | 4.890e-2 | 0 |
| CrumpledDevelopable-covariance-500 | 5.724026e-4 | 8.912045e-4 | 500 | 500 | 0.253 | 0.221 | 499 | 501 | 0 / max | 0 / max | 1.954e-2 | 0 |
| CrumpledDevelopable-intrinsic-3000 | 2.723743e-7 | 2.326685e-7 | 3000 | 3000 | 1.342 | 1.351 | 506 | 508 | 0 / max | 0 / max | 5.655e-2 | 0 |
| CrumpledDevelopable-intrinsic-500 | 6.472441e-5 | 7.021902e-5 | 500 | 500 | 0.239 | 0.194 | 508 | 508 | 0 / max | 0 / max | 5.745e-2 | 0 |
| blub_ctl-combinatorial-300 | 1.733067e0 | 1.678362e0 | 807 | 900 | 2.992 | 2.022 | 3330 | 3292 | cut fails / max/max/stalled | 0 / max/max/max | 8.917e-3 | 2 |
| blub_ctl-combinatorial-3000 | 3.733598e0 | 1.129796e0 | 602 | 2795 | 0.389 | 0.792 | 3204 | 2144 | cut fails / lsf/lsf/lsf | 0 / stalled/stalled/stalled | 3.296e-2 | 0 |
| blub_ctl-covariance-300 | 1.319488e-2 | 1.319488e-2 | 900 | 900 | 1.231 | 0.999 | 2868 | 2868 | 0 / max/max/max | 0 / max/max/max | 2.020e-2 | 0 |
| blub_ctl-covariance-3000 | 2.719293e-11 | 9.052823e-9 | 6056 | 7252 | 2.491 | 7.702 | 2788 | 2832 | 0 / max/max/lsf | 0 / max/lsf/max | 1.444e-4 | 9 |
| blub_ctl-intrinsic-300 | 2.142529e0 | 2.142529e0 | 900 | 900 | 1.811 | 1.979 | 3040 | 3040 | 0 / max/max/max | 0 / max/max/max | 1.767e-2 | 0 |
| blub_ctl-intrinsic-3000 | ERR error: flow failed: the energy is undefi | 1.070741e0 |  | 1360 |  | 0.383 |  | 2860 |  | 0 / stalled/lsf/lsf | 3.030e-2 | 0 |
| blub_quadrangulated-combinatorial-300 | 4.146214e0 | 4.146214e0 | 600 | 600 | 29.408 | 30.676 | 56806 | 56806 | 0 / max/max | 0 / max/max | 2.348e-2 | 0 |
| blub_quadrangulated-covariance-300 | 5.186006e0 | 5.186006e0 | 600 | 600 | 27.035 | 27.720 | 56598 | 56598 | 0 / max/max | 0 / max/max | 1.906e-2 | 0 |
| blub_quadrangulated-intrinsic-300 | 6.198339e0 | 6.198339e0 | 600 | 600 | 33.605 | 33.360 | 56604 | 56604 | 0 / max/max | 0 / max/max | 1.734e-2 | 4 |
| blub_triangulated-combinatorial-300 | 3.724887e0 | 3.724887e0 | 600 | 600 | 30.003 | 29.304 | 56822 | 56822 | 0 / max/max | 0 / max/max | 2.206e-2 | 0 |
| blub_triangulated-covariance-300 | 5.632107e0 | 5.632107e0 | 600 | 600 | 33.106 | 33.190 | 56614 | 56614 | 0 / max/max | 0 / max/max | 1.801e-2 | 0 |
| blub_triangulated-intrinsic-300 | 6.474141e0 | 6.474141e0 | 600 | 600 | 35.852 | 36.584 | 56634 | 56634 | 0 / max/max | 0 / max/max | 1.759e-2 | 0 |
| bob_ctl-combinatorial-300 | 1.229317e0 | 1.229317e0 | 900 | 900 | 5.960 | 5.657 | 10422 | 10422 | 0 / max/max/max | 0 / max/max/max | 1.813e-2 | 0 |
| bob_ctl-combinatorial-3000 | 9.014284e-1 | 4.286885e-1 | 2362 | 9000 | 6.799 | 45.902 | 9278 | 8596 | 0 / lsf/stalled/stalled | 0 / max/max/max | 1.899e-4 | 5 |
| bob_ctl-covariance-300 | 9.265637e-1 | 9.265637e-1 | 900 | 900 | 4.505 | 4.268 | 10116 | 10116 | 0 / max/max/max | 0 / max/max/max | 1.833e-2 | 0 |
| bob_ctl-covariance-3000 | 1.194679e-2 | 2.417610e-2 | 9000 | 9000 | 30.704 | 32.823 | 8292 | 8366 | 0 / max/max/max | 0 / max/max/max | 1.108e-4 | 1 |
| bob_ctl-intrinsic-300 | 1.639766e0 | 1.639766e0 | 900 | 900 | 6.247 | 5.448 | 10254 | 10254 | 0 / max/max/max | 0 / max/max/max | 2.454e-2 | 0 |
| bob_ctl-intrinsic-3000 | 4.165154e0 | 2.444759e0 | 2078 | 6576 | 6.214 | 30.871 | 9126 | 9040 | 0 / lsf/lsf/lsf | 0 / max/lsf/stalled | 1.185e-2 | 2 |
| bob_isotropic-combinatorial-300 | 2.863124e0 | 2.863124e0 | 600 | 600 | 14.504 | 14.575 | 24638 | 24638 | 0 / max/max | 0 / max/max | 2.761e-2 | 0 |
| bob_isotropic-covariance-300 | 4.627136e0 | 4.627136e0 | 600 | 600 | 14.970 | 14.220 | 24592 | 24592 | 0 / max/max | 0 / max/max | 1.890e-2 | 0 |
| bob_isotropic-intrinsic-300 | 4.626109e0 | 4.626109e0 | 600 | 600 | 15.563 | 14.990 | 24610 | 24610 | 0 / max/max | 0 / max/max | 1.895e-2 | 0 |
| bob_quad-combinatorial-300 | 3.491903e0 | 3.491903e0 | 600 | 600 | 23.304 | 23.207 | 42734 | 42734 | 0 / max/max | 0 / max/max | 2.393e-2 | 0 |
| bob_quad-covariance-300 | 7.352449e0 | 7.352449e0 | 600 | 600 | 22.278 | 22.826 | 42660 | 42660 | 0 / max/max | 0 / max/max | 2.074e-2 | 0 |
| bob_quad-intrinsic-300 | 8.249370e0 | 8.249370e0 | 600 | 600 | 21.949 | 24.263 | 42630 | 42630 | 0 / max/max | 0 / max/max | 2.951e-2 | 0 |
| bob_tri-combinatorial-300 | 3.076487e0 | 3.076487e0 | 600 | 600 | 23.108 | 23.044 | 42750 | 42750 | 0 / max/max | 0 / max/max | 2.283e-2 | 0 |
| bob_tri-covariance-300 | 7.024482e0 | 7.024482e0 | 600 | 600 | 22.363 | 23.569 | 42662 | 42662 | 0 / max/max | 0 / max/max | 3.613e-2 | 0 |
| bob_tri-intrinsic-300 | 7.126420e0 | 7.126420e0 | 600 | 600 | 23.987 | 25.133 | 42656 | 42656 | 0 / max/max | 0 / max/max | 2.476e-2 | 0 |
| spot_ctl-combinatorial-300 | 1.530048e0 | 1.530048e0 | 900 | 900 | 3.072 | 3.235 | 5512 | 5512 | 0 / max/max/max | 0 / max/max/max | 2.640e-2 | 0 |
| spot_ctl-combinatorial-3000 | 2.485319e0 | 1.037054e0 | 975 | 9000 | 1.213 | 23.911 | 5216 | 3952 | cut fails / lsf/lsf/lsf | 0 / max/max/max | 1.027e-4 | 19 |
| spot_ctl-covariance-300 | 3.154606e-1 | 3.154606e-1 | 900 | 900 | 2.819 | 2.431 | 5112 | 5112 | 0 / max/max/max | 0 / max/max/max | 1.805e-2 | 0 |
| spot_ctl-covariance-3000 | 1.092982e-8 | 7.550502e-4 | 9000 | 9000 | 12.045 | 14.817 | 3888 | 3808 | 0 / max/max/max | 0 / max/max/max | 1.047e-4 | 11 |
| spot_ctl-intrinsic-300 | 6.202401e0 | 5.054394e0 | 900 | 900 | 3.373 | 2.347 | 5584 | 5448 | 0 / max/max/max | 0 / max/max/max | 1.024e-4 | 7 |
| spot_ctl-intrinsic-3000 | 4.892556e0 | 3.322428e0 | 6313 | 5359 | 23.363 | 5.157 | 5484 | 1260 | 0 / stalled/max/max | 0 / max/lsf | 1.951e-4 | 2 |
| spot_quadrangulated-combinatorial-300 | 8.276170e0 | 8.276170e0 | 600 | 600 | 13.451 | 12.570 | 23392 | 23392 | 0 / max/max | 0 / max/max | 2.351e-2 | 0 |
| spot_quadrangulated-covariance-300 | 9.867035e0 | 9.867035e0 | 600 | 600 | 14.241 | 14.118 | 23044 | 23044 | 0 / max/max | 0 / max/max | 1.869e-2 | 0 |
| spot_quadrangulated-intrinsic-300 | 1.351481e1 | 1.351481e1 | 600 | 600 | 16.480 | 15.103 | 23008 | 23008 | 0 / max/max | 0 / max/max | 1.799e-2 | 0 |
| spot_triangulated-combinatorial-300 | 8.062792e0 | 8.062792e0 | 600 | 600 | 13.707 | 12.577 | 23388 | 23388 | 0 / max/max | 0 / max/max | 2.784e-2 | 0 |
| spot_triangulated-covariance-300 | 1.143797e1 | 1.143797e1 | 600 | 600 | 12.757 | 12.467 | 23168 | 23168 | 0 / max/max | 0 / max/max | 2.085e-2 | 0 |
| spot_triangulated-intrinsic-300 | 1.664337e1 | 1.664337e1 | 600 | 600 | 17.093 | 16.104 | 23278 | 23278 | 0 / max/max | 0 / max/max | 2.856e-2 | 0 |

Reading: where the floor is not reached nothing changes (identical energies). Where it is reached
the flows run further (bob `E_P` 3000: 9000 iterations, E 0.43 instead of 0.90 after 2362; spot `E_P`
3000: E 1.04 instead of 2.49 after 975), because the old runs ended in `LineSearchFailed` on
squeezed faces. `E_λ` over-runs of the control meshes end at a higher energy than before (spot
7.6e-4 instead of 1e-8, bob 2.4e-2 instead of 1.2e-2) because they no longer collapse faces to
reach the flat pancake; blub `E_λ` 3000 ends at 9e-9 instead of 2.7e-11. These are the
over-runs; with 300 it/round the energies are identical.

## Study (SPEC task 5)

`sh crates/ddg-developability/examples/study.sh 60 "100 400"` after
`cargo build --release -p ddg-developability --example develop`. Every OBJ of `data/models`
(11 files; polygons fan-triangulated) x 6 energies x iteration budgets 100 and 400 per round, two
rounds (coarse + one subdivision), 1 thread per job, 60 jobs in parallel, ε = 1e-3, default options
(remeshing on). 132 jobs in 3 min 7 s wall (machine load average about 60 from other agents; the
times below are single-thread seconds under that load). Logs: `runs_tmp/dev-flow-study/job-*.log`.
Result: all 132 jobs finished; no error, no non-finite energy or position, no line-search failure,
no `NonDescent`; 1 stall (`EnergyStalled` in round 0 of spot_control_mesh Intrinsic, budget 400);
undefined trials: 1 in total; cut and flatten succeeded in all 132.
Columns: `flips/collapses` are remesh operations over both rounds; `cut` is terminals(+added)/cut
edges/charts; `flipped uv` and `mean |u|` are from the cut at ε = 1e-3 on the final mesh. E0 is
the energy of the input (energy-specific). Energies of different kinds are not comparable.

Budget 400 iterations per round:

| mesh | energy | faces in->out | E0 | E_final | s | status | flips/collapses | cut: term(+add)/edges/charts | flipped uv | mean abs u |
|---|---|---|---|---|---|---|---|---|---|---|
| CrumpledDevelopable | Combinatorial | 512->2040 | 33.9 | 0.00173 | 4.116 | max_iter/max_iter | 0/1 | 0(+0)/0/1 | 0 | 0.000831 |
| CrumpledDevelopable | CombinatorialMax | 512->2044 | 126 | 0.0224 | 1.774 | max_iter/max_iter | 0/1 | 1(+0)/6/1 | 0 | 0.00115 |
| CrumpledDevelopable | Covariance | 512->2000 | 61.9 | 0.000639 | 3.534 | max_iter/max_iter | 1/8 | 0(+0)/0/1 | 0 | 0.011 |
| CrumpledDevelopable | CovarianceMax | 512->2024 | 18.3 | 0.00209 | 3.505 | max_iter/max_iter | 0/3 | 0(+0)/0/1 | 0 | 0.00149 |
| CrumpledDevelopable | Intrinsic | 512->2032 | 79.6 | 6.68e-05 | 1.305 | max_iter/max_iter | 0/2 | 0(+0)/0/1 | 0 | 2e-05 |
| CrumpledDevelopable | IntrinsicMax | 512->2040 | 91.3 | 0.00518 | 1.676 | max_iter/max_iter | 0/1 | 0(+0)/0/1 | 0 | 0.000933 |
| blub_control_mesh | Combinatorial | 220->812 | 23.7 | 2.11 | 1.526 | max_iter/max_iter | 6/25 | 125(+0)/152/1 | 57 | 1.45 |
| blub_control_mesh | CombinatorialMax | 220->852 | 98.2 | 4.62 | 1.637 | max_iter/max_iter | 5/5 | 280(+0)/301/1 | 3 | 0.0873 |
| blub_control_mesh | Covariance | 220->700 | 37.6 | 0.00983 | 0.988 | max_iter/max_iter | 7/33 | 2(+0)/16/1 | 36 | 2.52 |
| blub_control_mesh | CovarianceMax | 220->788 | 10.3 | 0.0641 | 2.324 | max_iter/max_iter | 0/13 | 16(+0)/41/1 | 57 | 2.13 |
| blub_control_mesh | Intrinsic | 220->746 | 89 | 1.48 | 1.002 | max_iter/max_iter | 2/19 | 137(+0)/146/1 | 35 | 1.22 |
| blub_control_mesh | IntrinsicMax | 220->832 | 76 | 3.5 | 2.103 | max_iter/max_iter | 0/6 | 227(+0)/250/1 | 19 | 0.113 |
| blub_quadrangulated | Combinatorial | 14208->56798 | 19 | 3.88 | 167.064 | max_iter/max_iter | 0/5 | 663(+0)/1690/1 | 0 | 0.0378 |
| blub_quadrangulated | CombinatorialMax | 14208->56828 | 61.8 | 15.8 | 140.425 | max_iter/max_iter | 0/2 | 3150(+0)/5354/1 | 0 | 0.00633 |
| blub_quadrangulated | Covariance | 14208->56542 | 45.7 | 4.79 | 169.838 | max_iter/max_iter | 0/49 | 1041(+0)/2104/1 | 0 | 0.0117 |
| blub_quadrangulated | CovarianceMax | 14208->56824 | 10.5 | 1.91 | 165.817 | max_iter/max_iter | 0/1 | 214(+0)/746/1 | 0 | 0.252 |
| blub_quadrangulated | Intrinsic | 14208->56536 | 49.1 | 5.91 | 179.217 | max_iter/max_iter | 2/49 | 1279(+0)/2512/1 | 0 | 0.00919 |
| blub_quadrangulated | IntrinsicMax | 14208->56832 | 44.4 | 8.22 | 105.901 | max_iter/max_iter | 0/0 | 1722(+0)/3110/1 | 0 | 0.00812 |
| blub_triangulated | Combinatorial | 14208->56812 | 19.6 | 3.56 | 152.299 | max_iter/max_iter | 0/4 | 587(+0)/1495/1 | 0 | 0.0374 |
| blub_triangulated | CombinatorialMax | 14208->56832 | 62.2 | 17.1 | 161.897 | max_iter/max_iter | 0/0 | 2958(+0)/4912/1 | 0 | 0.00788 |
| blub_triangulated | Covariance | 14208->56534 | 45.6 | 5.13 | 181.483 | max_iter/max_iter | 0/59 | 1168(+0)/2406/1 | 0 | 0.00964 |
| blub_triangulated | CovarianceMax | 14208->56824 | 10.7 | 2.04 | 167.955 | max_iter/max_iter | 0/1 | 259(+0)/1006/1 | 0 | 0.0919 |
| blub_triangulated | Intrinsic | 14208->56580 | 48.4 | 6.81 | 186.133 | max_iter/max_iter | 1/45 | 1375(+0)/2647/1 | 0 | 0.00867 |
| blub_triangulated | IntrinsicMax | 14208->56824 | 44.4 | 8.48 | 142.241 | max_iter/max_iter | 0/1 | 1865(+0)/3271/1 | 0 | 0.00878 |
| bob_controlmesh | Combinatorial | 668->2572 | 16.3 | 1.51 | 6.612 | max_iter/max_iter | 0/17 | 327(+0)/348/1 | 0 | 0.00112 |
| bob_controlmesh | CombinatorialMax | 668->2656 | 58 | 10.3 | 6.478 | max_iter/max_iter | 0/2 | 753(+0)/785/1 | 0 | 0.000141 |
| bob_controlmesh | Covariance | 668->2490 | 47.2 | 1.39 | 5.863 | max_iter/max_iter | 0/31 | 255(+0)/306/1 | 1 | 0.029 |
| bob_controlmesh | CovarianceMax | 668->2640 | 11.1 | 1.13 | 7.905 | max_iter/max_iter | 0/4 | 386(+0)/440/1 | 1 | 0.0143 |
| bob_controlmesh | Intrinsic | 668->2546 | 67.7 | 2.79 | 6.529 | max_iter/max_iter | 0/21 | 349(+0)/383/1 | 0 | 0.00122 |
| bob_controlmesh | IntrinsicMax | 668->2606 | 58.5 | 8.15 | 8.086 | max_iter/max_iter | 1/9 | 670(+0)/696/1 | 0 | 0.00984 |
| bob_isotropic | Combinatorial | 6174->24628 | 17.4 | 2.79 | 90.926 | max_iter/max_iter | 0/10 | 579(+0)/1136/1 | 0 | 0.0515 |
| bob_isotropic | CombinatorialMax | 6174->24696 | 56.8 | 11.5 | 74.044 | max_iter/max_iter | 0/0 | 2374(+0)/3028/1 | 0 | 0.00185 |
| bob_isotropic | Covariance | 6174->24554 | 43.8 | 4.3 | 105.770 | max_iter/max_iter | 2/26 | 1096(+0)/1558/1 | 0 | 0.00532 |
| bob_isotropic | CovarianceMax | 6174->24696 | 11.9 | 2.23 | 112.171 | max_iter/max_iter | 0/0 | 448(+0)/960/1 | 0 | 0.0488 |
| bob_isotropic | Intrinsic | 6174->24590 | 44.4 | 4.48 | 106.043 | max_iter/max_iter | 4/23 | 1163(+0)/1588/1 | 0 | 0.00482 |
| bob_isotropic | IntrinsicMax | 6174->24672 | 48.5 | 9.07 | 80.680 | max_iter/max_iter | 0/3 | 2103(+0)/2729/1 | 0 | 0.00202 |
| bob_quad | Combinatorial | 10688->42736 | 13.3 | 3.35 | 135.340 | max_iter/max_iter | 0/2 | 420(+0)/856/1 | 0 | 0.162 |
| bob_quad | CombinatorialMax | 10688->42752 | 38.2 | 13.3 | 91.447 | max_iter/max_iter | 0/0 | 3017(+0)/3895/1 | 0 | 0.0172 |
| bob_quad | Covariance | 10688->42640 | 42.1 | 6.96 | 142.175 | max_iter/max_iter | 3/23 | 2784(+0)/3217/1 | 0 | 0.00943 |
| bob_quad | CovarianceMax | 10688->42752 | 8.58 | 2.53 | 132.234 | max_iter/max_iter | 0/0 | 342(+0)/842/1 | 0 | 0.163 |
| bob_quad | Intrinsic | 10688->42584 | 42.9 | 7.58 | 144.470 | max_iter/max_iter | 2/36 | 3034(+0)/3434/1 | 0 | 0.00913 |
| bob_quad | IntrinsicMax | 10688->42752 | 34.5 | 10.3 | 131.256 | max_iter/max_iter | 0/0 | 2504(+0)/3180/1 | 0 | 0.0185 |
| bob_tri | Combinatorial | 10688->42738 | 13.6 | 3.01 | 136.723 | max_iter/max_iter | 0/4 | 400(+0)/838/1 | 0 | 0.159 |
| bob_tri | CombinatorialMax | 10688->42752 | 38 | 12.7 | 125.328 | max_iter/max_iter | 0/0 | 2784(+0)/3606/1 | 0 | 0.0159 |
| bob_tri | Covariance | 10688->42636 | 42.9 | 6.47 | 134.003 | max_iter/max_iter | 1/22 | 1765(+0)/2510/1 | 0 | 0.0103 |
| bob_tri | CovarianceMax | 10688->42752 | 8.81 | 2.32 | 131.083 | max_iter/max_iter | 0/0 | 291(+0)/775/1 | 0 | 0.158 |
| bob_tri | Intrinsic | 10688->42638 | 43.8 | 6.59 | 110.589 | max_iter/max_iter | 1/24 | 1701(+0)/2528/1 | 0 | 0.0102 |
| bob_tri | IntrinsicMax | 10688->42752 | 35.5 | 9.52 | 130.692 | max_iter/max_iter | 0/0 | 2195(+0)/2987/1 | 0 | 0.0163 |
| spot_control_mesh | Combinatorial | 372->1352 | 40.2 | 2.09 | 1.211 | max_iter/max_iter | 0/23 | 196(+0)/234/1 | 15 | 0.206 |
| spot_control_mesh | CombinatorialMax | 372->1424 | 170 | 17.5 | 1.069 | max_iter/max_iter | 0/8 | 535(+0)/544/1 | 0 | 0.00581 |
| spot_control_mesh | Covariance | 372->1236 | 99 | 0.294 | 2.772 | max_iter/max_iter | 3/42 | 70(+0)/123/1 | 55 | 0.871 |
| spot_control_mesh | CovarianceMax | 372->1382 | 24.3 | 1.02 | 5.375 | max_iter/max_iter | 1/14 | 146(+0)/177/1 | 12 | 0.336 |
| spot_control_mesh | Intrinsic | 372->1394 | 172 | 9.4 | 2.406 | stalled/max_iter | 4/17 | 437(+0)/451/1 | 0 | 0.000628 |
| spot_control_mesh | IntrinsicMax | 372->1438 | 168 | 21.8 | 3.795 | max_iter/max_iter | 0/10 | 533(+0)/549/1 | 0 | 0.000436 |
| spot_quadrangulated | Combinatorial | 5856->23362 | 42 | 8.22 | 71.721 | max_iter/max_iter | 0/13 | 1249(+0)/1866/1 | 0 | 0.00374 |
| spot_quadrangulated | CombinatorialMax | 5856->23416 | 126 | 31.6 | 77.321 | max_iter/max_iter | 0/4 | 3065(+0)/3839/1 | 0 | 0.000839 |
| spot_quadrangulated | Covariance | 5856->23000 | 111 | 8.78 | 80.551 | max_iter/max_iter | 7/83 | 1499(+0)/2009/1 | 2 | 0.00263 |
| spot_quadrangulated | CovarianceMax | 5856->23388 | 25 | 4.93 | 83.495 | max_iter/max_iter | 0/6 | 914(+0)/1398/1 | 0 | 0.0136 |
| spot_quadrangulated | Intrinsic | 5856->22902 | 115 | 13.1 | 76.713 | max_iter/max_iter | 7/99 | 1862(+0)/2457/1 | 0 | 0.00128 |
| spot_quadrangulated | IntrinsicMax | 5856->23334 | 104 | 23.2 | 67.407 | max_iter/max_iter | 0/15 | 2627(+0)/3200/1 | 0 | 0.00087 |
| spot_triangulated | Combinatorial | 5856->23360 | 39.1 | 7.64 | 60.386 | max_iter/max_iter | 1/11 | 1383(+0)/1860/1 | 0 | 0.00794 |
| spot_triangulated | CombinatorialMax | 5856->23406 | 123 | 30 | 55.407 | max_iter/max_iter | 0/3 | 3449(+0)/3961/1 | 0 | 0.000835 |
| spot_triangulated | Covariance | 5856->23086 | 112 | 10.6 | 74.507 | max_iter/max_iter | 5/64 | 1609(+0)/2078/1 | 0 | 0.00275 |
| spot_triangulated | CovarianceMax | 5856->23378 | 25.5 | 4.8 | 74.532 | max_iter/max_iter | 0/8 | 879(+0)/1363/1 | 0 | 0.0158 |
| spot_triangulated | Intrinsic | 5856->23256 | 118 | 15.6 | 79.991 | max_iter/max_iter | 1/27 | 2014(+0)/2596/1 | 0 | 0.00176 |
| spot_triangulated | IntrinsicMax | 5856->23314 | 105 | 24.8 | 57.766 | max_iter/max_iter | 0/16 | 2926(+0)/3500/1 | 0 | 0.000963 |

Reading the table: the dense meshes (5.8k to 14k faces, 23k to 57k after the subdivision) run in
60 to 190 s single-threaded for 800 iterations; the control meshes in 1 to 8 s. The cut
distortion is small (mean |u| ~ 1e-3 to 1e-2, 0 flips) on dense meshes, large on the coarse control
meshes (spot_control_mesh Covariance 0.87, blub_control_mesh 1 to 2.5 with 20 to 60 flips of ~800
faces): coarse meshes with few vertices cannot hide the residual curvature. Budget 100 per round
gives the same picture (`runs_tmp/dev-flow-study/table100.md`, not reproduced; e.g. blub_control_mesh
Covariance 85 flips, spot_control_mesh Covariance 97 flips, all dense meshes 0 flips).

## Benchmarks (`cargo bench -p ddg-developability --bench flow`)

Measured while the machine had a load average of 55 to 60 from other agents, so the "64 threads"
column shows contention rather than scaling, and the single-thread numbers are probably slightly
pessimistic. Remeshing on, default options.

| job | faces | iterations | evaluations | 1 thread | 64 threads |
|---|---|---|---|---|---|
| spot_triangulated, E_λ, 50 it | 5852 | 50 | 62 | 0.25 s (4.96 ms/it, 4.0 ms/eval) | 0.32 s |
| spot_triangulated, intrinsic | 5850 | 50 | 67 | 0.29 s (5.77 ms/it) | 0.33 s |
| bob_tri, E_λ | 10688 | 50 | 58 | 0.41 s (8.2 ms/it, 7.1 ms/eval) | 0.40 s |
| bob_tri, intrinsic | 10688 | 50 | 57 | 0.44 s (8.9 ms/it) | 0.43 s |
| blub_triangulated, E_λ | 14208 | 50 | 61 | 0.59 s (11.8 ms/it, 9.7 ms/eval) | 0.52 s |
| blub_triangulated, intrinsic | 14208 | 50 | 53 | 0.56 s (11.2 ms/it) | 0.47 s |
| spot control + 2 subdivisions, E_λ, 3 x 300 it | 5112 | 900 | 1427 | 1.68 s (1.9 ms/it) | 3.53 s |
| crumpled, E_λ, 300 it | 503 | 300 | 452 | 0.12 s | 0.57 s |
| bob control + 2 subdivisions, E_λ, 3 x 300 it | 10116 | 900 | 1397 | 3.62 s (4.0 ms/it) | 4.87 s |

Per evaluation (value + gradient + line search overhead) about 0.7 us per face single-threaded
(7 ms for 10.7k faces). The paper: "seconds to a few minutes" single-threaded on a 3.4 GHz i7
laptop for its examples; our spot run (the SPEC target: <= 60 s single-threaded) takes 1.7 s for 900
iterations. Multithreading does not help at these sizes under load (thread pool overhead, ~5k to 14k
vertices); a quiet-machine measurement of the speed-up remains to be done.

## Open / for others

* Bugs found since: see "Robustness study" (D31 to D33).
* `Evaluation` for the intrinsic energies with `omit_valence_three = false` errors on undefined
  vertices (cut stage propagates the error); not seen on any data model.
* Not done: gluing cracks and overlap avoidance after the cut (D23, out of scope), Python/Blender
  bindings.
