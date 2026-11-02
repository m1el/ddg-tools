# WP17 notes (performance: benches and the optimizations they showed were needed)

Build: `CARGO_TARGET_DIR=$PWD/target-wp17 cargo build --release -p ddg-repulsive-surfaces --bench step --bench product`
(or `cargo bench -p ddg-repulsive-surfaces --bench step -- --header --suite fig13`). Scripts and raw logs of
the runs below: `runs_tmp/wp17/` (not committed): `study_tol.sh` (tolerance study), `exp3.sh` (obstacles),
`timed.sh` (the timed pass: `ab45k.sh`, Fig. 13 suite, `scaling.sh`, `ref1t.sh` + `ref1t_table.py`).

## Benches

* `benches/step.rs`: whole steps of `Session` with the breakdown of `StepReport::seconds` (BVH, energy +
  gradient, factor = metric tree + fractional tree + `L′`, solve = GMRES with products and preconditioner,
  line search and its trials, CCD, projection, remeshing), GMRES iterations, max residual, final energy
  (17 digits, for bitwise comparisons). `--suite fig13` (the 14 Fig. 13 models or analogues at the paper's
  sizes, barycenter + area as the paper's runs), `--suite wp18b` (the scenes of `COMPARISON.md` §2.3),
  `--mesh SPEC` (generators, OBJ, scene), `--threads`, `--warmup`, `--steps`, `--gmres-tol`, `--restart`,
  `--chi`, `--no-lo`, `--obstacle-prec`, `--compat`, `--log` (per-step CSV), `--cache DIR` (generated meshes).
* `benches/product.rs`: one hierarchical product with 1–4 vectors, the block tree build, the vertex tree,
  one preconditioner application with 1–2 vectors, per thread count.
* `StepReport::seconds` (`Timings`) gained `metric_tree`, `fractional_tree`, `lprime` (parts of `factor`)
  and `preconditioner` (part of `solve`).

## What was slow (45k-face figure-8 tunnel, 4 threads, before; 2 steps, `perf`)

Before WP17 a step on the 44,848-face figure-8 tunnel (remeshed to 53.5k faces in the first step) took
7.0–8.8 s at 4 threads (load average 4–11 from other sessions): GMRES 5.2 s of it (94 iterations per
step, 2 solves: products 32 ms, preconditioner 19 ms each), the block cluster tree build 0.68 s, energy
0.34 s, remeshing 0.2 s. The near-field loop of the products was 48 % of all samples; it ran in SSE2
(the crate is built for generic x86-64) over faces in their original order (random gathers of 96-byte
arguments, 4-byte column index per entry); the curves crate's BCT build flattened 18M near entries
through per-row `Vec`s on one thread and evaluated every kernel twice. The preconditioner solved each
axis column separately with a refinement step (4 back substitutions per column per application).

## Changes (all bitwise equal to before unless stated)

1. **`spatial::block_tree::BlockTree`** (replaces the curves crate's `BlockClusterTree` in
   `FaceMetricTree`, `FaceFractionalTree`, `VertexFractionalTree`; that crate is not edited, D30): the
   same blocks, far coefficients and near entries, stored in BVH order as dense `|A|×|B|` blocks (no
   column indices, contiguous arguments); each unordered pair's kernel evaluated once (the mirror block
   is the transpose), the build parallel over blocks; products parallel over leaves, any channel count,
   kernel split known at compile time (vectorized), AVX (no FMA, so the same roundings) with run-time
   detection. **Bitwise equal products** to the curves crate's (test `tests/block_tree.rs`, 2 meshes × 3
   leaf sizes × 4 χ × 1/3/8 threads), and the 45k step's energies are bitwise those of the base build.
2. **Batched solves**: GMRES is now a resumable computation (`GmresRun`) and `gmres_batch` runs several
   right-hand sides in lockstep, sending the products they need at the same time through one
   `LinearOperator::apply_batch` / `Preconditioner::apply_batch`. `SaddleSolver::solve_many`;
   `SchurComplement::with_direction` solves the `k` Schur columns and the direction's `y = A′⁻¹ g` as one
   batch (the session uses it; the quasi-Newton path is unchanged). One pass over the 18M near entries
   serves 2–3 right-hand sides (App. B's "3 columns at once" extended to the Schur columns). Every
   solution is bitwise the one of a single solve. The fallback chain (`RobustHierarchical`, merged from
   `fd010c7`) batches through its primary solver and falls back to one-by-one solves on a divergence.
3. **Preconditioner**: the three axis columns of a vector in one pass over the `LDLᵀ` factor
   (`solve_many_in_place`), vectors in parallel; **no refinement step inside the preconditioner**
   (`PreconditionerOptions::refine`, default `false`; `LaplacianSaddle::solve` still refines). Measured
   (5 scenes, 2–5 steps): iterations 98.8 → 97.6 (Borromean), 36.4 → 36.2 (EDT), 108.5 → 107.0
   (figure-8 tunnel), 23.2 → 22.6 (Hopf), 105.6 → 106.0 (Klein bottle); energies after the steps equal
   to 5 digits. Two tests that check the exact preconditioner algebraically (sparse = dense to 1e-15,
   the D33 behaviour of the printed form) set `refine: true`.
4. **Barnes–Hut**: the face data (barycenter, normal, area) packed in BVH order, so a leaf's faces and
   an exact block's partners are contiguous (bitwise equal).
5. **GMRES vector operations** (modified Gram–Schmidt, norms) in parallel on the step's pool in fixed
   4096-entry chunks, each update fused with the next dot product; results independent of the pool
   (differ from the sequential sums by round-off).
6. **Defaults changed (flow only, `solve::iterative::FLOW_GMRES`): GMRES tolerance `1e-4 → 1e-3`,
   restart `30 → 60`.** `GmresOptions::default()` (the solver primitive) and `reference_compat`
   (`1e-4`, 30, left) are unchanged. Justification below.

## The GMRES tolerance and restart (default change)

Study: 11 scenes/meshes × tolerance `1e-4` / `1e-3` / `1e-2` (and restart 100 at `1e-4`; restart
60/100 at `1e-3` on Walnut and the punctured torus), 200 steps (Bunny and the trefoil tunnel 100), 2
threads each, all in parallel (only counts and energies are compared). Every step of every run was
accepted; max constraint violations equal across tolerances (projection); no run failed except the
punctured torus, which fails in **every** variant and in the base build (below).

| scene | its/step 1e-4 / 1e-3 / 1e-2 | F at the end, 1e-4 / 1e-3 / 1e-2 |
|---|---|---|
| EntangledDoubleTorus | 33.4 / 27.4 / 21.3 | 41.80 / 40.26 / 45.38 |
| HopfLink | 31.7 / 25.7 / 20.7 | 307.75 / 307.87 / 307.76 |
| HandcuffsPole | 30.6 / 26.1 / 17.0 | 6626 / 6842 / 6842 |
| Walnut (restart 30) | 80.9 / 61.5 / 25.9 | 4238 / 4352 / 4180 |
| Walnut, restart 60 | – / 49.7 / – | – / 4369 / – |
| LungGrowing | 26.3 / 19.6 / 15.7 | 64.11 / 64.25 / 64.09 |
| BorromeanRings | 55.9 / 46.1 / 36.7 | 4591 / 4579 / 4550 |
| KleinBottle (blocky) | 34.0 / 23.9 / 20.1 | 25.29 / 24.71 / 25.22 |
| Bunny (100 steps) | 34.4 / 29.8 / 23.2 | 56.48 / 56.38 / 56.36 |
| trefoil tunnel 23k (100) | 50.8 / 41.7 / 32.6 | 159.90 / 159.91 / 159.85 |
| genus-3 slab 6k | 20.3 / 16.0 / 11.1 | 56.30 / 56.34 / 56.63 |

The energies after 200 steps differ between tolerances only as much as trajectories do (remeshing
amplifies round-off; the same spread as WP5's CCD on/off runs); no scene is systematically worse at
`1e-3`. Why the direction needs no more: (a) the compressed operator itself is only accurate to
`1e-3` (random) – `2e-2` (smooth) at `χ = 0.5`, and hierarchical directions differ from the dense ones
by 4e-3–2e-2 (`NOTES-WP6.md`), so `1e-4` solves an approximate system ten times more accurately than it
is approximated; (b) `C_N x = 0` holds to round-off whatever the accuracy of the Schur columns (the
direction uses the `N = C_N Z` of the same `Z`), and the positional rows are restored by the
projection; (c) the reference stops on the preconditioned residual with true residuals of 0.1–0.2.
`1e-3` keeps a factor 100 below the divergence threshold. `1e-2` was as good on these runs and 20 %
cheaper still, but is not taken (robustness first). Restart 30 → 60 halves the iterations where GMRES
needs several cycles (Walnut 61.5 → 49.7 at `1e-3`; 100 gives the same as 60) and changes nothing on
the others (one cycle). Restart 100 at `1e-4`: identical on the 9 one-cycle scenes, Walnut 80.9 → 57.8.

**Punctured torus (open surface, barycenter + area, no boundary terms)**: the energy falls to 0.01, the
step size doubles every step and GMRES diverges at step 38–80 in every variant, also in the build before
WP17 (step 61, residual 0.73): a property of that problem (free boundary), not of these changes; the
Fig. 13 bench times only its first steps.

## Obstacles in the preconditioner (orchestrator request)

Default mode, iterations per step (60 steps; Hamburger 40), with `O_diag` in the operator only, and
with the obstacle term added to the preconditioner's Laplacian block as the sparse vertex matrix
`O = 2 Σ_S a_S (hi_S g_c·g_c' + lo_S/9)` (`solve::iterative::obstacle_vertex_matrix`,
`HierarchicalOptions::obstacle_in_preconditioner`):

| scene | without | with |
|---|---|---|
| Walnut (mesh obstacle) | 22.4 (200 steps: 80.9) | 22.6 (200 steps: 81.3) |
| LungGrowing (mesh obstacle) | 23.9 | 29.6 |
| Hamburger (mesh obstacle) | 23.8 | 30.3 |
| Reconstruction (point cloud, `O_diag` off by default since `fd010c7`) | 19.5 (4 steps) | – |

The obstacle term in `L′` does not help (worse on two of three) and stays off. The reference's way —
adding the obstacle's `fr` diagonal to the middle factor `L_{2−s}` — makes `P⁻¹ = L′⁻¹ (L + D) L′⁻¹`
larger where `(A + O_diag)⁻¹` is smaller (the wrong direction for a positive `O_diag`) and was not
ported. The Walnut slowdown over 200 steps follows the growth of the mesh (1.3k → 6.6k faces, area ×8)
with or without the obstacle term; restart 60 halves it. Reconstruction in default mode: 19.5 GMRES
iterations per step, no slowdown (the 1600 iterations came from the point-cloud `O_diag`, now
`reference_compat` only).

## Results (timed pass, 2026-11-02, `runs_tmp/wp17/timed/`)

Machine: AMD Ryzen Threadripper 3970X (32 cores, Zen 2 — the paper's 3990X is the 64-core part of the
same family). One run at a time; outside load (another session's application runs, `load.txt`): 1-minute
load average 16–17 during the A/B, 6.5–12 during Fig. 13 and scaling, 8–13 during the 1-thread
comparison, 6–11 during the base 1-thread runs. "Before" = the branch at `31f9057` (with only the
bench and the timing fields added), "after" = this branch.

**45k figure-8 tunnel, 4 threads, barycenter + area (44,848 faces, 53.5k after the first remeshing),
interleaved A B A B (`ab45k.csv`):**

| | step | GMRES its | solve | products | preconditioner | metric tree | energy | F after 4 steps |
|---|---|---|---|---|---|---|---|---|
| before | 9.30 / 9.48 s | 94.3 | 7.32 s | 4.5 s (98) | 2.6 s | 0.78 s | 0.37 s | 886.175015813004 |
| after | 4.16 / 3.82 s | 62.3 | 2.47 s | 1.46 s (64) | 0.90 s | 0.31 s | 0.40 s | 886.843721243789 |

(The energies differ by 0.08 % through the tolerance change; with changes 1–2 alone they were bitwise
equal.)

**Fig. 13, 4 threads, warm-up 2, mean of 5 steps (`fig13.csv`)** — paper values read off Fig. 13:

| # | model (analogue) | faces (input → after remeshing) | paper s/step | ours s/step | GMRES its/step | solve | metric tree | energy | remesh |
|---|---|---|---|---|---|---|---|---|---|
| 1 | triple torus (genus-3 slab) | 6,200 → 6,758 | ≈0.25 | 0.23 | 17 | 0.08 | 0.04 | 0.03 | 0.03 |
| 2 | moving bars (HandcuffsPole scene) | 10,688 | ≈0.45 | 0.40 | 23 | 0.14 | 0.05 | 0.04 | 0.05 |
| 3 | broken figure-8 (reference mesh) | 11,740 → 11,416 | ≈0.37 | 0.50 | 37 | 0.23 | 0.07 | 0.04 | 0.04 |
| 4 | punctured torus | 9,886 → 10,046 | ≈0.65 | 0.43 | 37 | 0.22 | 0.04 | 0.03 | 0.03 |
| 5 | handcuffs (k = 1) | 10,960 → 12,340 | ≈0.45 | 0.41 | 22 | 0.14 | 0.06 | 0.04 | 0.05 |
| 6 | 2×-wound handcuffs | 16,020 → 18,270 | ≈0.75 (14k) | 0.70 | 31 | 0.29 | 0.10 | 0.09 | 0.07 |
| 7 | point cloud (sphere 20k + bunny cloud, volume) | 20,480 | ≈1.6 | 0.73 | 17 | 0.21 | 0.10 | 0.09 | 0.05 |
| 8 | nested cages (sphere 20k + bunny mesh obstacle) | 20,480 | ≈1.75 | 0.83 | 19 | 0.22 | 0.10 | 0.12 | 0.06 |
| 9 | walnut (genus-0 slab) | 23,328 → 25,634 | ≈1.4 | 0.76 | 13 | 0.21 | 0.12 | 0.11 | 0.09 |
| 10 | trefoil tunnel | 23,344 → 26,586 | ≈1.4 | 1.16 | 39 | 0.59 | 0.13 | 0.13 | 0.09 |
| 11 | 6×-wound handcuffs | 21,640 → 24,372 | ≈1.35 | 1.41 | 63 | 0.79 | 0.14 | 0.13 | 0.10 |
| 12 | double torus (`scene-bench.txt`, 35k) | 35,454 → 35,520 | ≈1.8 | 1.70 | 31 | 0.62 | 0.22 | 0.23 | 0.17 |
| 13 | keys (genus-4 slab) | 37,312 → 39,772 | ≈2.05 | 1.47 | 19 | 0.48 | 0.22 | 0.25 | 0.13 |
| 14 | figure-8 tunnel | 44,848 → 53,542 | ≈2.7 | 3.15 | 59 | 1.85 | 0.27 | 0.31 | 0.22 |

At or below the paper on 12 of 14; above on #3 (0.50 vs ≈0.37 s, 37 GMRES iterations/step) and #14
(3.15 s vs ≈2.7 s). #14 is our analogue of the paper's mesh: it is remeshed from 44.8k to 53.5k faces
in the first step (the generator's edge lengths are uneven) and needs 59 iterations/step for 2
solves; per face it is 0.059 ms (paper ≈0.060 ms/face at 45.2k), i.e. on a 45k mesh of the same
difficulty about 2.65 s. The reference program on the same mesh, same machine, 4 threads, before
WP17's runs (`ref_f8/`, load ≈5): 4.0–5.2 s/step (it remeshes to 57–58k faces; 13 + 26 GMRES
iterations per step, left preconditioned).

**Thread scaling, 40k torus (40,050 faces, barycenter + area, 14 iterations/step), warm-up 1, mean of
3 (`scaling.csv`):** 1 / 2 / 4 / 8 / 16 / 32 / 64 threads: 3.32 / 2.08 / 1.44 / 1.14 / 0.97 / 0.93 /
1.03 s (speed-up 1.6 / 2.3 / 2.9 / 3.4 / 3.5 / 3.2). Parts that stop scaling at 8–16 threads: the
preconditioner (≈150–160 ms: the `L′` back substitutions are sequential per vector), remeshing
(≈150 ms, sequential), the energy (≈140 ms: the Barnes–Hut gradient's sequential passes), the metric
tree build (≈105–115 ms). Products scale to 16 threads (621 → 149 ms for 16).

**1 thread against the reference (WP18b scenes; 30 steps, Bunny 15, Reconstruction 6; mean ms/step,
median in brackets; `ref1t/`, `base1t/`; the blocky Klein bottle with the reference's non-orientable
build):**

| scene | faces | reference | ours before | ours after | after / reference | GMRES its/step before → after |
|---|---|---|---|---|---|---|
| EntangledDoubleTorus | 10,506 | 1149 (1139) | 2367 (2381) | 1257 (1238) | 1.09 | 36.0 → 29.8 |
| HopfLink | 9,216 | 667 (666) | 1617 (1614) | 967 (954) | 1.45 | 30.0 → 24.0 |
| HandcuffsPole | 10,688 | 955 (973) | 2126 (2163) | 1127 (1138) | 1.18 | 32.0 → 27.4 |
| Walnut | 1,280 | 259 (149) | 186 (166) | 140 (132) | 0.54 | 14.1 → 9.6 |
| LungGrowing | 5,120 | 506 (508) | 694 (688) | 482 (483) | 0.95 | 23.8 → 18.7 |
| Bunny | 28,576 | 4205 (4026) | 7657 (7701) | 3522 (3463) | 0.84 | 45.8 → 38.1 |
| Reconstruction | 81,920 | 9764 (9672) | 11496 (11546) | 7638 (7706) | 0.78 | 20.0 → 18.0 |
| KleinBottle (blocky) | 1,800 | 403 (412) | 464 (420) | 237 (225) | 0.59 | 69.3 → 42.8 |

Faster than the reference on 5 of 8 scenes; slower on EntangledDoubleTorus (1.09×), HandcuffsPole
(1.18×) and HopfLink (1.45×). Why: those three spend 29.8 / 27.4 / 24 GMRES iterations per step on a
true residual of `1e-3` where the reference stops after 20.7 / 13.0 / 9.3 iterations per step (counted
in its stdout of these runs) on its preconditioned residual (true residuals 0.1–0.2, `NOTES-WP6.md`); and at 9–11k faces a 1-thread step of ours is
≈ 1/4 tree builds (metric tree 0.15–0.2 s, vertex tree), ≈ 1/3 solve, the rest energy, line-search
trials, CCD (Hopf: factor 212, solve 345, energy 122, line search 120, projection 66, CCD 50 ms).

## Tried and not kept

* Rows outer / blocks inner in the near-field loop (accumulator kept for a whole row): 1.3–2× slower
  (the argument blocks are no longer reused across the rows of a leaf).
* `cbrt`-based kernel powers (`ρ^{−5/3}`, `ρ^{−4/3}`) in the tree builds: slower than `powf` (Rust's
  software `cbrt`), reverted.
* The obstacle term in the preconditioner (above).

## Proposed SPEC §14.5 "Results of record" (WP17)

> Benches `benches/step.rs` (whole steps, breakdown; `--suite fig13`, `--suite wp18b`) and
> `benches/product.rs`. Threadripper 3970X, 4 threads (Fig. 13's setting), barycenter + area: the 14
> Fig. 13 models (analogues where the paper's meshes are not published) at 0.23–1.70 s/step up to
> 40k faces, at or below the paper on 12 of 14; the figure-8 tunnel analogue (44.8k faces, 53.5k after
> its first remeshing, 59 GMRES iterations/step) 3.15 s against the paper's ≈2.7 s (0.059 ms/face, the
> paper's per-face rate), the reference program on the same mesh and machine 4.0–5.2 s. Thread
> scaling (40k torus): 3.3 s at 1 thread, 1.44 s at 4, 0.97 s at 16, no gain past 16 (sequential `L′`
> back substitutions, remeshing, Barnes–Hut gradient passes). One thread against the reference on the
> WP18b scenes: 0.54–0.95× on Walnut, LungGrowing, Bunny, Reconstruction, blocky Klein bottle; 1.09×,
> 1.18×, 1.45× on EntangledDoubleTorus, HandcuffsPole, HopfLink (more GMRES iterations: ours stops on
> the true residual). WP17 changed the defaults of the flow's GMRES to tolerance `1e-3` and restart
> 60 (§6.3; measured: same energies after 200 steps on 10 scenes, 20 % fewer iterations; `1e-4`/30
> stays the solver default and the `reference_compat` setting) and the preconditioner's `L′` solves to
> no refinement step (iterations unchanged); all other changes are bitwise or round-off.

## Not done / possible next steps

* The products are memory-bound (288 MB of near-field values per pass at 45k faces); halving that needs
  `f32` storage of the near field (a result change beyond round-off; not done) or a higher-order far field
  with a larger admissibility (fewer near entries).
* Remeshing (serial) and the Barnes–Hut gradient's sequential parts limit the scaling past 8 threads.
* The three `needless_range_loop` clippy warnings of `examples/compare_reference.rs` (WP18's file,
  not in WP17's paths) remain; the library, tests and benches are clippy-clean.
* The preconditioner's 3-column block solves (one pass over the factor per vector) were chosen without
  a clean A/B (the machine was loaded at that time); with more than ~8 threads one task per column may
  scale better.
* Allocation churn: fresh `|F| × C` buffers per product (≈4 % of a 4-thread step in page faults).
