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

## Results (timed pass)

(filled from `runs_tmp/wp17/timed/`)

## Not done / possible next steps

* The products are memory-bound (288 MB of near-field values per pass at 45k faces); halving that needs
  `f32` storage of the near field (a result change beyond round-off; not done) or a higher-order far field
  with a larger admissibility (fewer near entries).
* Remeshing (serial) and the Barnes–Hut gradient's sequential parts limit the scaling past 8 threads.
* The two `needless_range_loop` clippy warnings of `examples/compare_reference.rs` (WP18's file) remain.
