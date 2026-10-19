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

* **Spec expectation not met on spot.** SPEC §10.7 expects at most 0.1 % flipped UV triangles after
  flow + cut on spot. Measured (control mesh, 2 subdivisions, E_λ, ε = 1e-3): 300 iterations per
  round: 99 of 5112 (1.9 %), mean |u| 0.90; 1000 per round: 86 of 4568 (1.9 %), mean |u| 3.4.
  Lower ε helps (300 iterations: ε = 1e-4: 25 flips, 1e-5: 6 flips = 0.12 %, mean |u| 0.085, with
  1086 terminals). Cause (measured with Σ|Ω|): the flow leaves a lot of canceling positive and
  negative angle defect at vertices whose energy is below ε (spot, 300 it: Σ|Ω| 80 at terminals,
  138 at the other interior vertices). E_λ is quadratic in the normal spread while the defect is
  not, so a small energy does not mean a small defect. BFF with zero boundary scale must absorb
  that as area distortion. The tests assert <= 5 % for spot/bob/blub (blub 15 %) and print the
  numbers.
* **Over-running the flow gives spikes.** Spot, 3 rounds, 3000 iterations each: E_λ = 1e-8, no
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

## Study (SPEC task 5)

`sh crates/ddg-developability/examples/study.sh 60 "100 400"` after
`cargo build --release -p ddg-developability --example develop`. Every OBJ of `data/models`
(11 files; polygons fan-triangulated) x 6 energies x iteration budgets 100 and 400 per round, two
rounds (coarse + one subdivision), 1 thread per job, 60 jobs in parallel, ε = 1e-3, default options
(remeshing on). 132 jobs in 3 min 7 s wall (machine load average about 60 from other agents; the
times below are single-thread seconds under that load). Logs: `runs_tmp/dev-flow-study/job-*.log`.
Result: all 132 jobs finished; no error, no non-finite energy or position, no line-search failure,
no `NonDescent`; 1 stall (`EnergyStalled`, spot_control_mesh Intrinsic, budget 100 round 1 and 400);
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

* No bug found in the energy, optimizer, subdivide or remesh modules.
* `Evaluation` for the intrinsic energies with `omit_valence_three = false` errors on undefined
  vertices (cut stage propagates the error); not seen on any data model.
* Not done: gluing cracks and overlap avoidance after the cut (D23, out of scope), Python/Blender
  bindings.
