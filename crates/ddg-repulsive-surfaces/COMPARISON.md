# Comparison with the authors' program (WP18, `SPEC.md` §14.7)

Black-box comparison of `ddg-repulsive-surfaces` with the authors' published program (Yu,
Brakensiek, Schumacher, Crane, *Repulsive Surfaces*, 2021; reference `c51b464`, MIT). The reference
is built outside the repository (`/var/data/discrete-geometry/RepulsiveSurfaces-build/`, GCC 15,
MKL `mkl_gnu_thread`, headless driver, one thread; its README lists the one numerically neutral
patch). Everything below was produced on 2026-11-02 with the tools of this commit:

* `tests/reference_fixtures.rs` + `tests/fixtures/reference/` (parity fixtures, replayed by `cargo test`);
* `examples/compare_reference.rs` (`run`: our flow on a scene file with mode `compat` / `default` /
  `noedge` / `compat-noedge`; `stats`: area, volume, extents, exact energy, self-intersection of OBJs;
  `swept`: exact swept test between two states of equal connectivity);
* `scripts/rs-compare-reference.sh` (one scene, both programs, in parallel), `scripts/rs_compare_tables.py`
  and `scripts/rs_compare_timing.py` (the tables below). Logs, CSVs, OBJs and PNG plots of this study are in
  `runs_tmp/wp18b/` (not committed).

## 1. Parity fixtures (exact agreement of the formulas)

Three small closed meshes, dumped by the reference's `--dump-initial --dump-metric`: `tube_regular`
(the reference's own `Lowres/tube_regular.obj`, 369 vertices / 734 faces) and two generated lumpy
tori (`torus192`: 96 / 192, `torus504`: 252 / 504; `tests/fixtures/reference/*.obj` is the input both
programs read, 17 digits). `(α, β) = (6, 12)`, no `reference_compat` switch is needed for these
quantities (they do not involve a step). Measured by `cargo test --release --test reference_fixtures -- --nocapture`:

| quantity | tube_regular | torus192 | torus504 |
|---|---|---|---|
| all-pairs energy, relative difference | 7.1e-16 | 4.2e-16 | 1.8e-16 |
| all-pairs differential, max component error / max component | 1.1e-14 | 1.4e-15 | 5.7e-15 |
| dense `A = B + B⁰` (the reference's `M_small`, rows checked), max abs error / largest entry | 1.1e-14 (6.5e-11 / 5.9e3) | 4.1e-15 (4.3e-12 / 1.1e3) | 5.1e-15 (5.1e-12 / 1.0e3) |
| barycenter, area rows (`a_v/A`, `∇A`), max abs error | 9e-19, 2e-17 | 3e-18, 2e-16 | 2e-18, 1e-16 |
| volume row, `3×` ours vs reference (D9) | 1.6e-17 | 2.2e-16 | 2.6e-16 |
| Barnes–Hut θ = 0.5 energy: ours (leaf 8) vs reference | 2.4e-5 | 3.5e-6 | 3.0e-5 |
| same, ours with leaf size 1 | 5.0e-4 | 6.7e-5 | 7.1e-5 |
| BH energy: reference vs exact / ours (leaf 8) vs exact | 8.6e-5 / 1.1e-4 | 6.9e-7 / 4.2e-6 | 1.2e-5 / 4.2e-5 |
| BH differential, 2-norm relative: ours vs ref / ref vs exact | 1.6e-4 / 1.0e-4 | 5.4e-6 / 2.9e-7 | 5.5e-6 / 4.2e-6 |

* The energy and its exact differential, the dense Sobolev matrix (all 1107² entries of the tube
  were compared at generation time against the block structure `A ⊗ I₃`: exactly block-triplicated,
  symmetric to 5e-13) and the constraint rows agree to round-off (the metric to 1e-13 entry-wise
  relative to the entry or to `10⁻⁶` of the largest entry: 7e-14, 1.7e-13, 2.3e-13). The reference's metric block is
  `FillMatrixHigh + FillMatrixLow` with the same `w_TP = 1` weight and no separate scaling: no
  convention fix is needed. Only the **volume row** differs, by the exact factor 3 (D9: the
  reference's `Σ aN` against our exact `⅓ Σ aN`). The committed metric fixtures store the
  `x`-channel block `A[i][j] = M[3i][3j]` (all rows for `torus192`, the first 12 for the other two)
  and the five constraint rows in full, < 1 MB in all.
* **Barnes–Hut** agrees to the approximation error, not to round-off: the two trees split
  differently (ours: leaves of ≤ 8 faces, per-leaf admissibility `max(r_ℓ, r_C) < θ d_box`; the reference:
  its own cluster tree, splitting policy and leaf size), so the sets of lumped far-field pairs differ.
  Both are within ≈ 1e-4 of the exact value on the tube, ours within 4.2e-5 on the other two,
  and our differential is within 1.6e-4 (relative 2-norm) of the reference's, of the same size as the
  reference's own distance to the exact differential (1.0e-4). With leaf size 1 ours is *less*
  accurate at the same θ on the tube (5e-4), i.e. the difference is the tree shape and not a bug;
  the BH energy of the reference is not "more correct", only differently approximate (D4/D5: both
  differentiate `Ẽ` exactly).

## 2. Trajectories

### Method

* **Runs.** Per scene, the reference (`rsurfaces_headless`, `--threads 1`, remeshing as in the
  scene) and ours twice: **compat** = `reference_compat: true` with CCD off (the reference has none),
  **default** = our defaults (CCD on; off for the subcritical Klein exponents, as `scene_options`
  does). All eleven reference runs and 24 of ours ran at the same time on the 64-core machine
  (ours with 2 threads), for part of the time together with a ~15-core job that was not ours:
  **the ms/step of these runs are contended and not comparable** (§2.3 has the timed pass).
  Step counts: the scene's `iteration_limit` where affordable (HopfLink 200, Bunny 200, Walnut 220),
  else 300 (EntangledDoubleTorus, HandcuffsPole, KleinBottle), 350 (BorromeanRings, the scene's own
  limit), 150 (Reconstruction, limit 220) and 220 of 600 (LungGrowing; the reference's steps 232 and 233
  took 1335 s and 559 s, see §3, so the comparison stops at the last common dump). The energy-stall
  stop rule is off in ours.
* **Measures.** `F` = the run's own objective (energy + penalties, Barnes–Hut, θ = 0.5, each program's own
  tree) at every step; and, for the final mesh of every run (the reference's OBJs have 6 digits), our **exact
  all-pairs tangent-point energy** (no penalties, so it is not `F` on scenes with penalties), |V|, |F|, area, volume
  (signed, `⅙Σ f·(g×h)`; meaningless for the Klein bottle), bounding-box extents and the principal
  standard deviations `sd1 ≥ sd2 ≥ sd3` of the area-weighted face centroids (a rotation-invariant shape
  measure), the self-intersection of the dumped states (`ccd::intersects`; the reference every 10 or 20 steps, ours at every step
  end), the drift of area and volume against the initial values (the scenes' area/volume targets are the
  initial values unless the scene grows them, Walnut ×8 area, Lung ×3 area and ×2.8 volume,
  Reconstruction volume; there the drift column is the growth).
* **Not step for step.** Remeshing amplifies round-off; the interesting comparison is the energy curve and
  the shape. Reference dumps are every 10 steps (20 for Bunny, Lung and the blocky Klein bottle).

### 2.1 Tunnelling

The reference writes only the states after remeshing, whose connectivity differs from step to step, so
a swept test of the real motion is not possible. Two substitutes (`compare_reference swept`:
exact-time CCD, gap 1e-6 ℓ̄, 10⁶ advances, between consecutive dumps of equal connectivity):

* **Runs without remeshing**, every step dumped (`--no-remesh`; EntangledDoubleTorus 14 steps, HopfLink
  and BorromeanRings 30): the reference, ours compat and ours default have **no contact** in 14 + 30
  steps on EntangledDoubleTorus and HopfLink (ours compat takes the reference's small first steps:
  0.32 ℓ̄ per step, so it does not repeat the step-6 jump of §16 that our *default* start rule makes
  without CCD). **BorromeanRings is different**: after the **first step the reference's mesh
  self-intersects and stays so for steps 1–9** (nine consecutive dumps; min. distance 1e-4 ℓ̄ at step 1),
  then separates (30 steps: clean again from step 10; its first step is `τ = 2.1e-5` with a direction of
  0.6 ℓ̄ per step). Ours compat does the same at step 1 only (one intersecting state, clean from step 2,
  and the swept test finds the contact in steps 0→1 and 1→2); **ours default (CCD) is clean in all 31
  states**. The swept test of the reference's dumps finds a contact in nine of its first ten steps (all
  but 5→6), between states that are intersecting or nearly touching; its line-search steps there are
  small (`τ` 2e-5 to 3e-6, with directions of up to 0.6 ℓ̄).
* **Full runs**, state after each step (ours) or each 10th/20th step (reference) tested for
  self-intersection: none intersects, **except** the input of Walnut (`inside.obj` itself intersects, the reference's step-0 dump) and ours compat
  on BorromeanRings step 1 (1 of 351 step ends; the reference's dumps at steps 10, 20, … are clean, step 1
  was not dumped). Ours default (CCD on) is clean at every step end of every scene. This cannot detect a
  crossing that heals within a step, which is what the EntangledDoubleTorus step-6 jump of
  `SPEC.md` §16 looks like in a no-CCD run; the reference does not take such a jump in the
  14 steps without remeshing checked above.

### 2.2 Results per scene

Tables: `F` is each program's own objective at its last step; "exact TP energy" is our exact all-pairs
evaluation of the final mesh (comparable between programs); "drift" = final/initial − 1 for area and volume
(Klein bottle: the volume is not defined). ms/step in these tables is contended
(see above); the timed pass is §2.3. A step of the reference includes its remeshing.


#### EDT

| run | steps | ms/step | F (own, BH) | exact TP energy of the final mesh | V / F | area | volume | extents | sd1,sd2,sd3 | area drift, volume drift | dumped states intersecting |
|---|---|---|---|---|---|---|---|---|---|---|---|
| reference | 300 | 3935 | 20.061 | 19.939 | 5261 / 10526 | 158.33 | 88.477 | 7.49x7.16x7.21 | 2.132,2.017,1.276 | +2.3e-06, +1.3e+00 | 0/31 |
| ours compat | 300 | 4900 | 19.953 | 19.807 | 5289 / 10582 | 158.33 | 88.674 | 7.45x7.23x7.15 | 2.140,2.007,1.283 | +7.4e-07, +1.3e+00 | 0/301 step ends |
| ours default | 300 | 5039 | 8.426 | 8.3338 | 5434 / 10872 | 158.33 | 106.69 | 7.31x8.58x4.94 | 2.344,1.947,1.178 | +1.6e-06, +1.7e+00 | 0/301 step ends |

Objective F against step (each run's own Barnes-Hut value, penalties included):

| run | 0 | 1 | 5 | 10 | 25 | 50 | 100 | 150 | 200 | 300 |
|---|---|---|---|---|---|---|---|---|---|---|
| reference | 377.05 | 358.83 | 307.4 | 260.19 | 169.42 | 114.85 | 79.992 | 65.771 | 59.219 | 20.061 |
| ours compat | 376.52 | 359.74 | 306.48 | 259.25 | 168.97 | 114.49 | 79.529 | 65.788 | 58.626 | 19.953 |
| ours default | 376.52 | 328.87 | 156.38 | 138.22 | 116.28 | 95.05 | 69.589 | 61.896 | 41.8 | 8.426 |

#### Hopf

| run | steps | ms/step | F (own, BH) | exact TP energy of the final mesh | V / F | area | volume | extents | sd1,sd2,sd3 | area drift, volume drift | dumped states intersecting |
|---|---|---|---|---|---|---|---|---|---|---|---|
| reference | 200 | 4769 | 306.26 | 304.93 | 4572 / 9144 | 29.511 | 5.1407 | 2.97x2.97x3.83 | 0.961,0.635,0.635 | +2.6e-08, +4.1e-01 | 0/21 |
| ours compat | 200 | 8206 | 305.69 | 304.02 | 4568 / 9136 | 29.511 | 5.1678 | 2.96x2.95x3.83 | 0.960,0.635,0.633 | +2.3e-07, +4.1e-01 | 0/201 step ends |
| ours default | 200 | 5575 | 307.75 | 306.13 | 4539 / 9078 | 29.511 | 5.1665 | 2.96x2.95x3.84 | 0.961,0.635,0.631 | +1.6e-07, +4.1e-01 | 0/201 step ends |

Objective F against step (each run's own Barnes-Hut value, penalties included):

| run | 0 | 1 | 5 | 10 | 25 | 50 | 100 | 150 | 200 |
|---|---|---|---|---|---|---|---|---|---|
| reference | 849.34 | 776.16 | 586.61 | 414.96 | 312.3 | 306.34 | 306.28 | 306.24 | 306.26 |
| ours compat | 846.96 | 787.52 | 581.87 | 412.54 | 307.95 | 306.11 | 305.78 | 305.89 | 305.69 |
| ours default | 846.96 | 735.36 | 318.16 | 311.44 | 308.46 | 308.23 | 307.83 | 307.66 | 307.75 |

#### Borromean

| run | steps | ms/step | F (own, BH) | exact TP energy of the final mesh | V / F | area | volume | extents | sd1,sd2,sd3 | area drift, volume drift | dumped states intersecting |
|---|---|---|---|---|---|---|---|---|---|---|---|
| reference | 350 | 4067 | 7302.6 | 7324.7 | 5796 / 11592 | 12.129 | 0.72888 | 2.44x1.23x2.27 | 0.603,0.576,0.222 | +2.8e-07, +2.2e-01 | 0/36 |
| ours compat | 350 | 6932 | 4534.6 | 4549 | 5658 / 11316 | 12.129 | 0.85869 | 1.67x1.67x1.78 | 0.422,0.419,0.419 | +1.6e-06, +4.3e-01 | 1/351 step ends |
| ours default | 350 | 6213 | 4543 | 4557.2 | 5643 / 11286 | 12.129 | 0.85865 | 1.66x1.71x1.82 | 0.421,0.421,0.419 | +6.7e-07, +4.3e-01 | 0/351 step ends |

Objective F against step (each run's own Barnes-Hut value, penalties included):

| run | 0 | 1 | 5 | 10 | 25 | 50 | 100 | 150 | 200 | 300 | 350 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| reference | 4.6938e+05 | 2.6643e+05 | 46650 | 22569 | 14496 | 7343.1 | 7303.5 | 7302.8 | 7302.1 | 7301.9 | 7302.6 |
| ours compat | 4.6938e+05 | 2.8928e+05 | 13908 | 9726.3 | 6372.8 | 6227.6 | 5638.4 | 4714.1 | 4538 | 4534.4 | 4534.6 |
| ours default | 4.6938e+05 | 1.9186e+05 | 32737 | 14951 | 6403.8 | 6272 | 6083.5 | 5210.3 | 4590.5 | 4545.4 | 4543 |

#### Bunny

| run | steps | ms/step | F (own, BH) | exact TP energy of the final mesh | V / F | area | volume | extents | sd1,sd2,sd3 | area drift, volume drift | dumped states intersecting |
|---|---|---|---|---|---|---|---|---|---|---|---|
| reference | 200 | 10920 | 56.134 | 53.252 | 16990 / 33976 | 9.4692 | 1.5976 | 2.18x2.45x1.77 | 0.645,0.642,0.265 | +5.7e-08, -4.5e-06 | 0/11 |
| ours compat | 200 | 12431 | 56.424 | 53.306 | 16847 / 33690 | 9.4692 | 1.5976 | 2.12x2.45x1.87 | 0.646,0.641,0.265 | +8.9e-08, -1.1e-07 | 0/201 step ends |
| ours default | 200 | 12531 | 56.473 | 53.282 | 17122 / 34240 | 9.4692 | 1.5976 | 2.14x2.45x1.83 | 0.646,0.641,0.265 | +4.4e-08, -3.8e-08 | 0/201 step ends |

Objective F against step (each run's own Barnes-Hut value, penalties included):

| run | 0 | 1 | 5 | 10 | 25 | 50 | 100 | 150 | 200 |
|---|---|---|---|---|---|---|---|---|---|
| reference | 1.6245e+05 | 53032 | 4766.8 | 557.12 | 71.613 | 56.43 | 56.117 | 56.1 | 56.134 |
| ours compat | 1.6245e+05 | 36748 | 1949.8 | 375.07 | 68.458 | 56.797 | 56.418 | 56.468 | 56.424 |
| ours default | 1.6245e+05 | 70796 | 13416 | 1128.5 | 73.815 | 57.005 | 56.464 | 56.417 | 56.473 |

#### Handcuffs

| run | steps | ms/step | F (own, BH) | exact TP energy of the final mesh | V / F | area | volume | extents | sd1,sd2,sd3 | area drift, volume drift | dumped states intersecting |
|---|---|---|---|---|---|---|---|---|---|---|---|
| reference | 300 | 3834 | 6387.8 | 2783.6 | 6086 / 12176 | 2.1903 | 0.10829 | 1.48x1.13x0.34 | 0.426,0.332,0.070 | +3.6e-07, +8.6e-02 | 0/31 |
| ours compat | 300 | 4938 | 6395.4 | 2816.6 | 6081 / 12166 | 2.1903 | 0.10802 | 1.48x1.13x0.35 | 0.428,0.332,0.069 | +8.1e-07, +8.3e-02 | 0/301 step ends |
| ours default | 300 | 5393 | 6404 | 2816.8 | 6127 / 12258 | 2.1903 | 0.10796 | 1.48x1.13x0.32 | 0.428,0.333,0.069 | +2.2e-07, +8.3e-02 | 0/301 step ends |

Objective F against step (each run's own Barnes-Hut value, penalties included):

| run | 0 | 1 | 5 | 10 | 25 | 50 | 100 | 150 | 200 | 300 |
|---|---|---|---|---|---|---|---|---|---|---|
| reference | 64686 | 52266 | 33351 | 21476 | 16433 | 16101 | 15101 | 12237 | 8122.7 | 6387.8 |
| ours compat | 64702 | 52233 | 33352 | 21461 | 16491 | 16119 | 14907 | 11869 | 7752.9 | 6395.4 |
| ours default | 64702 | 56999 | 27149 | 16869 | 16180 | 15851 | 13260 | 9661.2 | 6625.6 | 6404 |

#### Walnut

| run | steps | ms/step | F (own, BH) | exact TP energy of the final mesh | V / F | area | volume | extents | sd1,sd2,sd3 | area drift, volume drift | dumped states intersecting |
|---|---|---|---|---|---|---|---|---|---|---|---|
| reference | 220 | 4452 | 4136.1 | 2338.9 | 3682 / 7360 | 122.03 | 27.17 | 4.51x4.46x4.47 | 1.076,1.071,1.063 | +7.0e+00, +4.0e+00 | 1/23 |
| ours compat | 220 | 4908 | 4082.7 | 2360.8 | 3625 / 7246 | 122.03 | 27.814 | 4.48x4.44x4.43 | 1.073,1.067,1.062 | +7.0e+00, +4.1e+00 | 0/221 step ends |
| ours default | 220 | 9056 | 4059.6 | 2188.1 | 3573 / 7142 | 122.03 | 28.301 | 4.49x4.50x4.51 | 1.082,1.074,1.072 | +7.0e+00, +4.2e+00 | 0/221 step ends |

Objective F against step (each run's own Barnes-Hut value, penalties included):

| run | 0 | 1 | 5 | 10 | 25 | 50 | 100 | 150 | 200 |
|---|---|---|---|---|---|---|---|---|---|
| reference | 24.132 | 9.5456 | 9.2685 | 13.8 | 47.79 | 188.48 | 609.05 | 1714.9 | 4273.7 |
| ours compat | 24.168 | 9.2048 | 9.3143 | 13.839 | 44.542 | 140.99 | 654.3 | 1741.9 | 4123.2 |
| ours default | 24.168 | 8.8078 | 9.983 | 15.07 | 47.191 | 155.88 | 540.34 | 2167 | 4238 |

#### Lung

| run | steps | ms/step | F (own, BH) | exact TP energy of the final mesh | V / F | area | volume | extents | sd1,sd2,sd3 | area drift, volume drift | dumped states intersecting |
|---|---|---|---|---|---|---|---|---|---|---|---|
| reference | 220 | 2492 | 70.893 | 48.559 | 3477 / 6950 | 5.8858 | 0.93016 | 1.18x2.02x0.78 | 0.625,0.332,0.237 | +8.8e-01, +7.8e-01 | 0/12 |
| ours compat | 220 | 3259 | 68.606 | 45.531 | 3567 / 7130 | 5.8858 | 0.93409 | 1.20x2.09x0.78 | 0.628,0.329,0.239 | +8.8e-01, +7.9e-01 | 0/221 step ends |
| ours default | 220 | 3219 | 68.77 | 45.404 | 3577 / 7150 | 5.8983 | 0.93597 | 1.21x2.09x0.78 | 0.627,0.331,0.240 | +8.8e-01, +7.9e-01 | 0/221 step ends |

Objective F against step (each run's own Barnes-Hut value, penalties included):

| run | 0 | 1 | 5 | 10 | 25 | 50 | 100 | 150 | 200 |
|---|---|---|---|---|---|---|---|---|---|
| reference | 19.97 | 19.93 | 20.833 | 22.001 | 25.009 | 29.309 | 36.987 | 46.61 | 64.597 |
| ours compat | 20.174 | 19.556 | 20.146 | 21.382 | 24.627 | 28.923 | 36.63 | 46.358 | 63.923 |
| ours default | 20.174 | 21.062 | 20.386 | 21.566 | 24.777 | 29.091 | 36.757 | 46.682 | 64.11 |

#### Recon

| run | steps | ms/step | F (own, BH) | exact TP energy of the final mesh | V / F | area | volume | extents | sd1,sd2,sd3 | area drift, volume drift | dumped states intersecting |
|---|---|---|---|---|---|---|---|---|---|---|---|
| reference | 150 | 18362 | 3.2113e+06 | 14.666 | 40925 / 81846 | 11.508 | 3.332 | 2.06x2.11x1.64 | 0.670,0.524,0.447 | -4.7e-01, -6.5e-01 | 0/16 |
| ours compat | 150 | 17634 | 3.1761e+06 | 449.09 | 40374 / 80744 | 11.794 | 3.246 | 2.06x2.07x1.73 | 0.672,0.523,0.445 | -4.6e-01, -6.6e-01 | 0/151 step ends |
| ours default | 150 | 19020 | 3.4378e+06 | 161.06 | 40552 / 81100 | 11.648 | 3.2029 | 2.03x2.07x1.70 | 0.669,0.521,0.442 | -4.7e-01, -6.7e-01 | 0/151 step ends |

Objective F against step (each run's own Barnes-Hut value, penalties included):

| run | 0 | 1 | 5 | 10 | 25 | 50 | 100 | 150 |
|---|---|---|---|---|---|---|---|---|
| reference | 3.0023e+05 | 2.8881e+05 | 2.6141e+05 | 2.5064e+05 | 2.2346e+05 | 2.6305e+05 | 6.0933e+05 | 3.2113e+06 |
| ours compat | 3.0064e+05 | 2.8616e+05 | 2.5697e+05 | 2.3621e+05 | 2.0848e+05 | 2.6285e+05 | 6.0394e+05 | 3.1761e+06 |
| ours default | 3.0064e+05 | 2.8349e+05 | 1.9833e+05 | 1.6838e+05 | 1.9559e+05 | 2.6352e+05 | 6.1741e+05 | 3.4378e+06 |

#### kbb

| run | steps | ms/step | F (own, BH) | exact TP energy of the final mesh | V / F | area | volume | extents | sd1,sd2,sd3 | area drift, volume drift | dumped states intersecting |
|---|---|---|---|---|---|---|---|---|---|---|---|
| reference | 300 | 1893 | 24.695 | 23.697 | 900 / 1800 | 100.01 | 84.789 | 5.64x5.63x5.63 | 1.629,1.625,1.623 | +1.5e-05, +1.3e+00 | n/a |
| ours compat | 300 | 1301 | 24.759 | 23.63 | 817 / 1634 | 100 | 88.507 | 5.64x5.66x5.63 | 1.629,1.625,1.622 | +7.7e-06, +1.4e+00 | n/a (immersion) |
| ours default | 300 | 1338 | 25.104 | 23.974 | 733 / 1466 | 100 | 82.735 | 5.62x5.74x5.56 | 1.659,1.611,1.604 | -1.7e-05, +1.2e+00 | n/a (immersion) |
| ours compat-noedge | 299 | 1477 | 24.739 | - | 900 / 1800 | 100.01 | 84.282 | - | - | +5.3e-05, +1.2e+00 | n/a (immersion) |
| ours noedge | 300 | 1234 | 24.735 | 23.563 | 900 / 1800 | 100.01 | 82.475 | 5.64x5.67x5.64 | 1.630,1.627,1.623 | +1.9e-05, +1.2e+00 | n/a (immersion) |

Objective F against step (each run's own Barnes-Hut value, penalties included):

| run | 0 | 1 | 5 | 10 | 25 | 50 | 100 | 150 | 200 | 300 |
|---|---|---|---|---|---|---|---|---|---|---|
| reference | 793.91 | 661.25 | 318.31 | 277.09 | 194.69 | 186.98 | 154.89 | 144 | 55.212 | 24.695 |
| ours compat | 793.23 | 598.07 | 288.6 | 269.43 | 179.19 | 170.37 | 60.084 | 29.959 | 27.304 | 24.759 |
| ours default | 793.23 | 413.25 | 245.14 | 234.81 | 181.33 | 167.52 | 152.43 | 137.87 | 155.73 | 25.104 |
| ours compat-noedge | 793.23 | 598.07 | 288.6 | 271.69 | 189.44 | 161.29 | 150.32 | 127.5 | 76.209 | - |
| ours noedge | 793.23 | 413.25 | 248.15 | 215.25 | 144.9 | 126.87 | 78.701 | 42.964 | 24.702 | 24.735 |

#### kb2

| run | steps | ms/step | F (own, BH) | exact TP energy of the final mesh | V / F | area | volume | extents | sd1,sd2,sd3 | area drift, volume drift | dumped states intersecting |
|---|---|---|---|---|---|---|---|---|---|---|---|
| reference | 40 | 40721 | 161.67 | 160.97 | 8192 / 16384 | 27.14 | 10.432 | 2.97x3.22x2.61 | 0.907,0.750,0.672 | +6.3e-05, +2.1e-01 | n/a |
| ours default | 40 | 23018 | 117.91 | 117.62 | 8267 / 16534 | 27.139 | 11.089 | 3.26x2.68x2.90 | 0.862,0.751,0.742 | -6.4e-06, +2.8e-01 | n/a (immersion) |
| ours noedge | 40 | 28213 | 149.65 | 149.04 | 8192 / 16384 | 27.139 | 10.69 | 3.29x2.85x2.76 | 0.882,0.756,0.696 | +6.8e-06, +2.4e-01 | n/a (immersion) |

Objective F against step (each run's own Barnes-Hut value, penalties included):

| run | 0 | 1 | 5 | 10 | 25 |
|---|---|---|---|---|---|
| reference | 387.2 | 257.46 | 208.99 | 181.19 | 165.48 |
| ours default | 387.37 | 259.81 | 236.82 | 198.45 | 142.36 |
| ours noedge | 387.37 | 259.81 | 249.37 | 212.98 | 151.83 |

#### kb1

| run | steps | ms/step | F (own, BH) | exact TP energy of the final mesh | V / F | area | volume | extents | sd1,sd2,sd3 | area drift, volume drift | dumped states intersecting |
|---|---|---|---|---|---|---|---|---|---|---|---|
| reference | 1 | 9451 | 323.29 | - | 8192 / 16384 | - | - | - | - | - | n/a |
| ours default | 40 | 14214 | 112.53 | 112.5 | 6662 / 13324 | 24.19 | 7.7632 | 3.01x3.29x1.88 | 0.955,0.806,0.508 | +1.8e-06, +3.5e-01 | n/a (immersion) |
| ours noedge | 40 | 25516 | 123.22 | 120.95 | 8192 / 16384 | 24.19 | 7.5122 | 3.09x3.38x1.78 | 1.005,0.799,0.485 | +7.0e-06, +3.1e-01 | n/a (immersion) |

Objective F against step (each run's own Barnes-Hut value, penalties included):

| run | 0 | 1 | 5 | 10 | 25 |
|---|---|---|---|---|---|
| reference | 309.77 | 323.29 | - | - | - |
| ours default | 309.93 | 248.94 | 186.01 | 152.34 | 122.64 |
| ours noedge | 309.93 | 248.94 | 247.49 | 177.85 | 147.45 |

#### Reading the tables

* **Compat tracks the reference.** With `reference_compat` and no CCD our energy curve follows the
  reference's to about 1 % on EntangledDoubleTorus (F 377.05/376.52 at step 0, 79.99/79.53 at step
  100, 59.22/58.63 at 200, 20.06/19.95 at 300; exact final energies 19.94/19.81; extents 7.49×7.16×7.21 /
  7.45×7.23×7.15; sd within 0.5 %), HopfLink (306.26/305.69 at 200; 304.93/304.02 exact; sd equal to
  2–3 digits), HandcuffsPole (exact 2783.6/2816.6; extents 1.48×1.13×0.34 / 1.48×1.13×0.35), Bunny
  (53.25/53.31; sd 0.645,0.642,0.265 / 0.646,0.641,0.265), Walnut (2338.9/2360.8; extents 4.51×4.46×4.47 /
  4.48×4.44×4.43), the blocky Klein bottle (23.70/23.63; sd 1.629,1.625,1.623 / 1.629,1.625,1.622) and
  LungGrowing at step 220 (48.6/45.5). The residual differences are trajectory noise (different BH trees,
  the D14/D9/D24 items below, remeshing round-off). Area is held to 1e-7–2e-6 relative by both on the orientable scenes (Klein bottle: reference 1.5e-5 and
  6.3e-5, ours 2e-6–7.7e-5 over the variants and meshes).
* **Default vs reference.** Our defaults reach the same states with a different path: the
  start rule (D15) takes larger first steps (EntangledDoubleTorus `F` 156 at step 5 against 307; HopfLink
  318 against 587; Walnut/Bunny differ in the first 25 steps, Bunny's `F` is higher at step 1–10 and
  equal from step 50) and ends at the same plateau on HopfLink (307.75 vs 306.26, within 0.5 %), Bunny
  (56.47 vs 56.13), HandcuffsPole (6404 vs 6388). On EntangledDoubleTorus the default is *ahead*: at 300
  steps `F = 8.4` (reference 20.1), volume 106.7 against 88.5, i.e. it has already finished the
  untangling that the reference completes later (the 200-step values 41.8 against 59.2). On Walnut the default's
  final exact energy is 6.5 % lower (2188 vs 2339), on LungGrowing 6.5 % lower than the reference (45.4 vs 48.6).
* **BorromeanRings is the one scene where the final states differ clearly:** the reference
  stagnates at `F ≈ 7302` from step 50 (extents 2.44×1.23×2.27, sd 0.603,0.576,0.222), ours reaches
  4535 (compat) and 4543 (default) with a round symmetric shape (sd 0.422,0.419,0.419; extents
  ~1.7×1.7×1.8); both our runs end in the same state although they
  take different paths. The reference's line search takes `τ ≈ 1e-5` for the first
  steps and the reference mesh passes through itself in steps 1–9 (§2.1); whether its plateau is a local
  minimum reached after that or an effect of the stalled line search was not investigated.
* **Constraints.** Area: both hold it to 1e-7–1e-6 relative on every scene. **Volume (Bunny): the reference's
  volume lags its target, ours does not** (final drift −4.5e-6 for the reference, −1.1e-7 compat, −3.8e-8
  default; D9 and one Newton step), 40 to 120 times tighter. With growth schedules
  (Walnut area ×8: both reach 122.03 = 8 × 15.25 to the printed digits).
* **Reconstruction (point-cloud obstacle, 82k → 41k vertices)** is the one case where ours is not
  close in the exact tangent-point energy of the final mesh: 14.7 (reference) against 449 (compat) and 161
  (default), while the *objective* is within 1.1 % (compat 3.176e6 against 3.211e6) or 7 % (default 3.438e6)
  and the shape statistics agree (sd 0.670,0.524,0.447 / 0.672,0.523,0.445 / 0.669,0.521,0.442; extents within 5 %).
  The objective is dominated by the obstacle term at these steps (it rises from 2.1–2.6e5 at step 50 to
  6.1e5 at step 100 and 3.2–3.4e6 at step 150 in all three runs as the growing volume approaches the cloud), and
  its model differs (D23: the reference's one-sided obstacle with metric diagonal against our two-sided
  kernel; not reproduced by `reference_compat`), so the balance with the self-repulsion differs. Not
  investigated further.

### 2.3 Time per step (each run alone, one thread, 30 steps; Bunny 15, Reconstruction 6)

The timed pass ran after all other runs of this study had finished, one process at a time (no process of
ours was running next to it; the other user's job was not running at the start). Mean ms/step over the steps
(median in brackets); the reference includes its remeshing and `compress()`, ours the whole `Session::step`.

| scene | faces | reference | ours compat | ours default | ours noedge | steps timed |
|---|---|---|---|---|---|---|
| EDT | 10506 | 1034 (med 1010) | 1901 (med 1869) | 2110 (med 2118) | - | 30 |
| Hopf | 9216 | 604 (med 595) | 2641 (med 2381) | 1450 (med 1461) | - | 30 |
| Handcuffs | 10688 | 830 (med 831) | 1467 (med 1544) | 1879 (med 1929) | - | 30 |
| Walnut | 1280 | 234 (med 139) | 160 (med 154) | 160 (med 139) | - | 30 |
| Lung | 5120 | 441 (med 453) | 653 (med 650) | 613 (med 611) | - | 30 |
| Bunny | 28576 | 3597 (med 3329) | 5934 (med 5684) | 6411 (med 6374) | - | 15 |
| Recon | 81920 | 8668 (med 8625) | 9828 (med 9669) | 10769 (med 10735) | - | 6 |
| kbb | 1800 | 378 (med 387) | 298 (med 294) | 395 (med 329) | 387 (med 337) | 30 |

Ours is **slower per step at one thread** than the reference on the large orientable scenes: 1.1× to 1.8×
(compat) and 1.2× to 2.3× (default), and on HopfLink 4.4× (compat) / 2.4× (default); it is faster on Walnut
(1280 faces, 0.7×) and equal or faster on the blocky Klein bottle (compat 0.8×, default 1.0×). Hints from the
logs, not profiled further: (a) GMRES iterations per step of ours (compat / default): EntangledDoubleTorus
20 / 35, HopfLink 59 / 29, Bunny 27 / 43; the reference's log shows about 23 per step on EntangledDoubleTorus
(equal to compat) and about 11 on HopfLink (340 iteration markers in 30 steps, counted from its stdout), so on
HopfLink our `reference_compat` solver needs five times the iterations of the reference although the
preconditioner has the reference's form; the difference to the reference is its absolute `ε` (D14), which is
not reproduced, and the left-preconditioned residual; (b) the default's right-preconditioned GMRES takes more
iterations than the left one on EntangledDoubleTorus and Bunny; (c) the paper's design target is many threads
(WP6/WP17), so this table only says that single-thread work per step has not been tuned against the reference.

## 3. Klein bottle

* **The published `master` cannot load any of the three meshes** of its own `scenes/KleinBottle/`:
  `rsurfaces_headless --scene KleinBottle/scene.txt` aborts in geometry-central's `ManifoldSurfaceMesh`
  (`GC_SAFETY_ASSERT … duplicate edge in list 249 -- 364`; `runs_tmp/wp18b/kb_master.log`). Fig. 26 (bottom) of the
  paper cannot have been made with it.
* **The `nonorientable` branch's commit `ce89c2c`** (general `SurfaceMesh`; edge splits and
  collapses commented out, flips and smoothing stay) loads them:
  * `klein_bottle_blocky.obj` (900 / 1800, exponents 3.5, 7, scene `KleinBottle/scene.txt`): 300 steps,
    `F` 793.91 → 24.70, |V|, |F| constant. All of ours end near the same state: `F` 24.7–25.1,
    exact energy 23.56–23.97 (reference 23.70), area 100.0, extents 5.62–5.74 on all axes, sd 1.60–1.66
    (nearly round and symmetric; sd within 3 % of each other). Ours with splits and collapses off
    ("noedge", the variant's remeshing) and the variant itself agree to 0.2 % in `F` (24.735 / 24.695) and 0.6 %
    in the exact energy (23.56 / 23.70) at 300 steps; with full remeshing ("default") the mesh
    coarsens to 733 / 1466 and the exact energy is 1.2–1.7 % higher (23.97) — a remeshing effect.
    The descent is faster in ours (F at step 100: 154.9 reference, 60.1
    compat, 152.4 default, 78.7 noedge; at step 150: 144.0 / 30.0 / 137.9 / 43.0), with plateaus at different steps in each run — trajectory noise
    of a flow that crosses its own sheets.
  * `klein_bottle_2.obj` (8192 / 16384): the variant takes 41 s per step in the contended run (23 s for ours
    default); 40 steps: `F` 387.2 → 161.7 reference, 117.9 default, 149.7 noedge (exact energy 161.0 / 117.6 / 149.0).
  * `klein_bottle.obj` (8192 / 16384): the variant's energy **rises** in step 1 (309.77 → 323.29) and
    step 2 does not finish (> 8 min earlier, killed after 30 min here by `timeout`). Ours: 309.93 → 112.5 in 40
    steps with full remeshing (6662 / 13324 faces at the end), 123.2 with splits and collapses off (exact
    energy 112.5 / 121.0). The stall is the reference's (no CCD, its line search at `τ` it cannot
    reduce, not investigated); ours has CCD off here because the exponents are subcritical.
* Ours accepts all three (D1: orientation optional), including the full-remeshing flow that the published
  code cannot run, and `compat-noedge` reproduces the variant's remeshing. Self-intersection and volume
  are meaningless for these immersions and not reported.

## 4. What explains the gaps; what `reference_compat` cannot reproduce

From `SPEC.md` §11 and §16, as seen in the tables:

* **D9** the reference's volume gradient is 3× ours (volume row of the parity fixtures): its volume lags the target
  (Bunny −4.5e-6, ours −1e-7 compat, −4e-8 default); `reference_compat` keeps our exact row.
* **D14** its absolute preconditioner `ε` (1e-8 / 1e-2) against our relative `1e-10 tr W / tr M`;
  **τ_min = 1e-20** absolute; **D24** its penalty value/gradient conventions; **D23** its obstacle model (Walnut, Lung,
  Reconstruction: both programs descend the same objective only up to those conventions — F at step 0 differs
  by 0.15 % on Walnut (24.13 against 24.17) and the same at the end by a few %).
* **D4/D5** Barnes–Hut trees differ (Parity table): every `F` in the tables is a different approximation of
  the same energy, within ≈ 1e-4 relative; the exact energies of the final meshes are the fair numbers.
* **D15** the line-search start rule, **D12** Newton steps, **D16** remeshing (`reference_compat`:
  smooth→flip, but our remesher is f64 with collapse/flip/smoothing checks where the reference uses
  float angles and a flip cap), **D26/D27** only matter for the comparison optimizers (not run here).
* **No CCD in the reference**: it self-intersects in steps 1–9 of BorromeanRings (§2.1); ours default never
  does at a step end.
* **Machine**: the reference was compiled with GCC 15 (it only builds with clang unpatched, README of the build) and MKL's
  `mkl_gnu_thread`, `-O3 -march=native`; one-thread timings may differ from a clang build.

## 5. Where ours is better, worse, and open items

**Better:** all 11 scenes run, including the Klein bottle that `master` cannot load and the one that stalls the
variant (§3); no self-intersection at any step end with CCD (BorromeanRings: the reference intersects for nine
steps and ends at `F = 7302` against 4535–4543); a tighter volume constraint; lower final energy than the
reference at equal step counts on EntangledDoubleTorus (8.3 against 19.9 at 300 steps), Walnut (−6.5 %),
LungGrowing (−6.5 %), BorromeanRings (−38 %); exact-tolerance parity of the formulas (§1).

**Worse:**
* Single-thread time per step is 1.1–2.3× the reference's on the large orientable scenes and 4.4× on HopfLink in
  compat mode (§2.3).
* **Possible bug (ours): GMRES divergence in compat mode on the Klein bottle at step 300.**
  `compare_reference run data/repulsive-surfaces/scenes/KleinBottle/scene.txt --mode compat-noedge --steps 300 --threads 1`
  ends with `error: the linear solver diverged: relative residual 1.466e3 (above 1e-1)` in step 300
  (deterministic: the same number with 2 threads; the run was at `F = 24.74`, `τ = 3.59`, no remeshing
  edge operations). `reference_compat` solves with left-preconditioned GMRES and the face-pair form;
  the other three Klein variants complete their 300 steps. A hard error ends the run: a step that fails
  the residual check could instead be retried with the right-preconditioned solver.
* Reconstruction: our final self-energy is 11–30× the reference's (449/161 against 14.7) at a similar objective
  (§2.2); the obstacle terms are different (D23), not investigated further.
* The `compat` mode intersects at step 1 of BorromeanRings like the reference (no CCD by construction).
* Not done: reference-side swept tests of runs with remeshing (the dumps do not allow it), BorromeanRings
  plateau analysis, LungGrowing beyond step 232 (the reference's steps 232 and 233 took 1335 s and 559 s
  of its own gradient step at one thread while the machine was loaded; the same scene had
  taken ~2 s per step before, not reproduced alone), a Fig. 12 timing curve (WP17).
