# Notes: the two defects from the WP18b black-box comparison

Logs and scratch under `runs_tmp/defects/` (gitignored). Release builds, `CARGO_TARGET_DIR=target-defects`.

## Defect 1: `SolverDiverged { residual: 1.466e3 }` (compat, splits/collapses off, blocky Klein bottle)

**Cause.** Not the CCD band (compat has no CCD), not a singular `L'`. In step 299 the flow shrinks the edge
(211, 419) from 0.048 to 7e-4 (splits/collapses are off, nothing repairs it; the immersion is subcritical, so the
energy is no barrier), and the Delaunay flip at the end of the step produces two needle faces
(area 1e-4 against a median 5.6e-2, min angle 0.1 degrees). `F` jumps from 24.7 (the value `energy_after` reports,
before the flip) to 4.4e9 at the start of step 300. On that mesh the printed preconditioner of compat (left
GMRES, face-pair `L_{2-s}`, congruence form) is useless: 200 iterations end at relative residual 1.47e3. Dump
of step 299 and the reproducer: `tests/fixtures/solver_fallback/`.

**Fix** (`src/step/direction.rs`, `RobustHierarchical`; `SPEC.md` 6.3). A hierarchical solve that ends in
`SolverDiverged` is repeated, for that solve and the rest of the step, by (1) right GMRES with the augmented
preconditioner on vertex pairs (the library defaults), restart 60, at least 800 iterations, on a rebuilt
preconditioner; if that fails too, (2) a dense LU of the metric (at most 6000 vertices). Each escalation is
recorded in `LinearReport::fallbacks` (a new `Vec<String>`, empty normally). The error stays only when neither
works. Justification for repairing instead of marking the step `Stuck`: a `Stuck` step leaves the positions
unchanged, and the next step sees the same needle mesh and the same failure; the dense solve gives a real
descent direction and the flow repairs itself (steps 2-5 from the dump: F 4.4e9, 1.1e9, 30.7, 27.8, 25.3, 26.8,
no further fallback). The dense path is limited by memory, so larger meshes depend on stage (1) only.

**Test.** `tests/flow.rs::solver_fallback_recovers_the_klein_bottle_needle_state` (4 s, release): fails
with `SolverDiverged` before the change, now step 1 reports a fallback and `F` drops below 1e3 in four steps.

## Defect 2: Reconstruction self-energy (point-cloud bunny)

**Checked and equal to the reference**: the point-cloud energy (one-sided, mesh face barycenter and normal
against the points, unit point weights, `|v.n|^alpha / |v|^beta`; `tp_pointcloud_obstacle_barnes_hut_0.cpp`),
the weight 1, the objective at step 0 (3.0064e5 against 3.0023e5, the Barnes-Hut trees differ), the volume
schedule (`0.11x` over 200 steps; volume at step 150 3.246 against 3.332 — the reference lags by D9), the
barycenter shift.

**Found: the reference adds the obstacle's metric diagonal for point clouds too.** `hs.h` builds the obstacle
BCT for any obstacle energy (`TPPointCloudObstacleBarnesHut0` has zero projectors and a zero first-order
operator), so `hi_diag`, `lo_diag` of the cloud are added to the metric. We refused point clouds
(`obstacle_diagonal` returned an error; the test asserted it). With the term (`K_2(T,S) = 0`,
`a_T = 1`) our metric is the reference's: **the first line-search step is `tau` 3.40e-4 against the
reference's 3.37e-4** (it was 1.9e-6: the unscaled `A` was 170 times too small), `F` after step 1 288990
against 288808 (0.06 %), and the objective follows the reference to step 125 (1.182e6 against 1.168e6;
before: 6 % off from step 10 on). Implemented in `CrossEvaluator::metric_diagonal`; it is **on in
`reference_compat` only** (`objective.rs::obstacle_diagonal(.., point_clouds)`): the preconditioner has no
obstacle term, and with the term in default mode GMRES takes 1600 iterations (300 s) per step on this scene.
Test: `tests/penalties.rs::obstacle_diagonal_matches_the_brute_force_formula` (point-cloud part).

**Not explained: the final self-energy.** Exact-kernel energies (`compare_reference stats`, `bh0.25`) of the
compat run with the term: step 50 / 100 / 150 = 11.1 / 1408 / 4267 against the reference's 2.8 / 5.5 / 14.9
(before the term: 449 at 150). Shape statistics (sd, extents, area, volume, mesh quality: min angle 28.6 against
32.0 degrees, edge range the same) agree. The energy is concentrated: at step 100, 83 % of ours sits in
100 faces (adjacent faces at 12-45 degrees at a distance of 0.015, i.e. creases) near the cloud, whose
nearest-cloud distance is 0.13-0.15 (the reference's faces with the highest energy are 0.31 from the cloud
and sit on its outer, unshrunk part; its minimum distance to the cloud is 0.24 against ours 0.13).
So our surface wraps a sharp feature of the cloud more tightly and bends around it, while the reference stays
away from it. The objective, dominated by the obstacle term (1e5-1e6), hardly sees this. Candidates not tested:
the reference's left GMRES with a large true residual smoothing the direction, `fr` part of the obstacle in its
preconditioner, remesher differences (D16), the CCD-free line search accepting different steps late in the run.
Open.
