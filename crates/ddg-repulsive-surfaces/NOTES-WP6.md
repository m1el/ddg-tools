# WP6 notes (hierarchical metric, sparse preconditioner, iterative solves)

Build: `CARGO_TARGET_DIR=$PWD/target-wp6 cargo test --release -p ddg-repulsive-surfaces --test hierarchical`
(12 tests, ~2 s release). Measurements: `cargo build --release -p ddg-repulsive-surfaces --example
wp6_precond`, then `target-wp6/release/examples/wp6_precond <solve|error|breakdown|dense|modes> <mesh> …`
(see the example's docs; one CSV line per run; logs of the runs below in `runs_tmp/wp6/`).

## Done

- `src/metric/hierarchical.rs`: `ClusterTreeOptions { chi 0.5, leaf_size 8, leaf_pairs 64 }`,
  `separated` (SPEC §8.3, `max(r)² ≤ χ² d²`, `d > 0`), `FaceMetricTree` (`B̃`, `B̃⁰` on one curves-crate
  `BlockClusterTree<2>` over WP3's `FaceTree`; far kernels from cluster centroids and averaged projectors,
  near kernels with exactly the dense expressions; `H̃a` once per tree; one 12-channel pass for 3 columns
  × (3 gradient + 1 average); obstacle diagonal hook), `FaceFractionalTree` (face-form `L_{2−s}`),
  `VertexFractionalTree` (D32 vertex pairs: BVH over vertices, mass `a_v`). Product form
  `2 Mᵀ diag(a) (diag(H̃a) M u − H̃ (a ⊙ M u))`, so `L 1 = 0` is bitwise exact for `B⁰`, `L_{2−s}`.
  Tags: `rs-kernel-matrix-form`, `rs-hierarchical-null-space`, `rs-projector-average`, `rs-box-distance`.
- `src/solve/precond.rs`: `LaplacianSaddle` (sparse unpivoted `LDLᵀ` of `L′`, order = AMD of the
  Laplacian block then the multipliers; per-axis factor when `C_Q = C₁ ⊗ I₃`; one refinement step per
  solve; the three axes solved in parallel, bitwise independent of the pool), `PreconditionerOptions
  { eps_rel 1e-10, eps_multi_rel 1e-4, fractional Vertices, cluster, form Augmented }`,
  `translations_fixed_per_component` (decides `ε` vs `ε_multi`), `SparsePreconditioner`
  (`rs-fractional-preconditioner`; `usage()`, `factor_seconds()`, `epsilon()`), `PreconditionerForm`.
- `src/solve/iterative.rs`: `MetricExtras { willmore: Option<SparseMatrix>, obstacle: Option<FaceDiagonal> }`
  (WP7 hooks), `HierarchicalOperator` (`A′`, `LinearOperator`, `apply_metric`, `usage()`),
  `HierarchicalOptions { cluster, gmres, preconditioner, willmore_in_preconditioner: false }`
  (`::reference()` = face pairs + left), `Backend::{Dense, Hierarchical(opts), Auto { dense_max_faces,
  hierarchical }}` (`Default` = Auto 2000; `resolve(n_faces)`), `hierarchical_saddle` (→
  `GmresSaddle<HierarchicalOperator, SparsePreconditioner>`, label "gmres (hierarchical)"),
  `dense_metric_with_extras`, `saddle_solver` (→ `Box<dyn SaddleSolver + Send>`).
- `tests/hierarchical.rs` (12): χ = 0 exactness of `B`, `B⁰`, face and vertex `L_{2−s}` (4 meshes × leaf
  1/3/8 × 2 exponent pairs; 2e-16..1.1e-15); error decay in χ and < 5 % at 0.5 (B, B⁰); `B⁰ 1 = L 1 = 0`
  bitwise, `B 1` ~ round-off; symmetry 1e-17..8e-17; operator with `w_TP`, `c_lo`, Willmore, obstacle
  diagonal vs `dense_metric_with_extras` at χ = 0 (≤ 7e-16); sparse vs dense preconditioner (congruence
  form, ≤ 1.6e-15, `C z = 0`) and the augmented form's properties; `L′` solve residual ≤ 7e-16;
  ε_multi selection; GMRES + Schur (direction and correction, area + volume, pins) = dense direct solve
  (χ = 0: ≤ 4e-12; χ = 0.5: ≤ 8e-4); bitwise equality at 1, 3, 8 threads (products and Schur
  directions, both fractional forms); `Backend::Auto` switching; vertex pairs beat face pairs;
  augmented vs congruence on a right side with net force and with pins.

## Decisions (and proposed SPEC changes)

- **D33 (proposed): the preconditioner as printed is singular on `A′`.** `P⁻¹ = L′⁻¹ blockdiag(L,0) L′⁻¹`
  has no multiplier component in its range (`blockdiag(·, 0)` zeroes it and `L 1 = 0` removes the
  translations). Right GMRES cannot produce the multipliers and stalls; left GMRES stops on its
  preconditioned residual with a large true residual. The multipliers are nonzero whenever the
  positional rows carry a force: pins (reaction forces) and per-component barycenters on interacting
  components (Hopf: right GMRES stalls at 8.6e-4 after 200 iterations; with 3 pins on a torus right
  stalls at 0.12, left "converges" at a true residual of 0.04; random right side: 0.08–0.11). Fix,
  same cost: `P⁻¹ = L′⁻¹ + L′⁻¹ blockdiag((L_{2−s} − K) ⊗ I₃, 0) L′⁻¹` with `K = W + εM (+ w_W WM⁻¹W)`.
  With `L′⁻¹ = [[E, F], [Fᵀ, G]]`, `L′⁻¹ blockdiag(K, 0) L′⁻¹ = [[E, 0], [0, −G]]`, so this is
  `[[E L E, F + E L F], [Fᵀ + Fᵀ L E, 2G + Fᵀ L F]]`: the same position block `E L E` (order
  `−2s`), plus the coupling `F` of `L′⁻¹` — the structure of `A′⁻¹`. Symmetric, nonsingular,
  `C_Q P⁻¹ [r; b] = b`. Proposed `rs-fractional-preconditioner` statement for WP10: the identity
  `L′⁻¹ blockdiag(K,0) L′⁻¹ = blockdiag(E, −G)`, symmetry, and `C_Q (P⁻¹ [r; b])_x = b`.
  `PreconditionerForm::Congruence` keeps the printed form (tests, reference mode).
- `FractionalPairs` is WP2's enum (default `Faces` there); the WP6 options default to `Vertices`
  (D32). SPEC §6.3 names it `FractionalForm::{VertexPairs, FacePairs}`: rename at merge (WP2's file).
- `Backend` is defined here (SPEC §8.4 sketch: `Hierarchical { bh, chi, gmres, eps_rel }`); actual:
  `Hierarchical(HierarchicalOptions)` and `Auto { dense_max_faces, hierarchical }` — θ of Barnes–Hut
  stays in `BhOptions` (the energy's).
- `ε` choice is automatic: `ε_multi` iff some component has no positional row of its own.
- `ClusterTreeOptions::leaf_pairs = 64`; the curves crate's block split (both sides) instead of the
  reference's "split the larger radius".
- Willmore: `MetricExtras::willmore` (WP7's matrix, unweighted) or computed `W M⁻¹ W`; weighted by
  `SobolevWeights::willmore_weight`; added to the preconditioner's Laplacian block when
  `willmore_in_preconditioner` (D7; measurement below).
- The obstacle `fr` diagonal of the reference (obstacle term in the preconditioner's `L_{2−s}`) is
  not implemented (only `hi`/`lo` of `FaceDiagonal`, in `A′`).

## Measurements (release, 2026-11-02)

Sweeps (a)–(c) ran in parallel (14 jobs × 4 threads; their times are under load, iterations are
exact). Timings (d), (e) ran alone, one job at a time; outside load during them: another session's
`rustc` (1–2 cores, load average 6–15 on 64 cores). Raw CSV: `runs_tmp/wp6/{solve,em,sweep,dense,
dense_chi,bd,bd_chi}.csv` (not committed).

**(a) Product error vs χ** (relative 2-norm against the χ = 0 tree, i.e. the exact operator; random
input / smooth input `(sin x + z, xy, cos 2z)`):

| mesh | χ | B | B⁰ | L_{2−s} faces | L^V_{2−s} vertices | near entries (B) |
|---|---|---|---|---|---|---|
| torus 2.1k | 1 / 0.5 / 0.25 | 3.2e-3 / 8.2e-4 / 1.7e-4 | 4.9e-3 / 1.2e-3 / 1.6e-4 | 1.3e-2 / 4.8e-3 / 1.3e-3 | 2.7e-2 / 9.6e-3 / 3.2e-3 | 0.36M / 0.85M / 2.1M |
| torus 10k | 1 / 0.5 / 0.25 | 3.5e-3 / 9.3e-4 / 2.1e-4 | 2.1e-3 / 5.3e-4 / 1.1e-4 | 1.3e-2 / 4.9e-3 / 1.7e-3 | 2.4e-2 / 8.2e-3 / 2.8e-3 | 1.9M / 4.1M / 11.9M |
| icosphere 5k | 1 / 0.5 / 0.25 | 3.4e-3 / 7.8e-4 / 1.0e-4 | 2.3e-3 / 5.4e-4 / 8.4e-5 | 1.2e-2 / 3.5e-3 / 6.6e-4 | 2.2e-2 / 6.7e-3 / 1.4e-3 | |
| hopf 9.2k | 1 / 0.5 / 0.25 | 2.4e-3 / 7.9e-4 / 2.0e-4 | 4.8e-3 / 1.1e-3 / 8.4e-5 | 1.3e-2 / 6.1e-3 / 2.1e-3 | 2.2e-2 / 8.1e-3 / 3.1e-3 | 2.0M / 4.6M / 9.9M |
| smooth, torus 10k | 1 / 0.5 / 0.25 | 2.4e-2 / 9.6e-3 / 3.4e-3 | 6.0e-2 / 2.4e-2 / 7.7e-3 | 1.0e-1 / 5.1e-2 / 1.9e-2 | 2.1e-1 / 1.2e-1 / 5.3e-2 | |
| smooth, hopf | 1 / 0.5 / 0.25 | 4.8e-2 / 1.9e-2 / 5.0e-3 | 5.2e-2 / 1.2e-2 / 1.6e-3 | 1.2e-1 / 6.5e-2 / 2.5e-2 | 2.1e-1 / 1.2e-1 / 5.3e-2 | |

`O(χ²)`-ish; at χ = 0.5 the metric is within 1e-3 (random) / 2.4e-2 (smooth: `L u` is a small
difference there) — §13.3's 5 % holds. Directions (GMRES to 1e-4, vertex pairs) vs the dense direct
solve: χ = 1 / 0.5 / 0.25: torus 3k 1.1e-2 / 4.2e-3 / 2.2e-3, stretched 3k 1.2e-2 / 1.0e-3 / 5.9e-4,
hopf 1.5e-2 / 2.2e-2 / 5.7e-3 (at 1e-4 the GMRES tolerance, not χ, dominates on Hopf).

**(b) GMRES iterations to 1e-4** (one solve `A′ x = [g; 0]`, `g` the Barnes–Hut TPE differential,
per-component barycenters, χ = 0.5, augmented form unless noted; "L" = left, its *true* residual in
brackets; ">200" = not converged, true residual in brackets):

| mesh | faces | vertex R | vertex L | face R | face L | vertex R congruence | face L congruence (reference) |
|---|---|---|---|---|---|---|---|
| icosphere | 1280 | 4 | 4 (1e-3) | 62 | 10 (0.19) | 4 | 10 (0.19) |
| icosphere | 5120 | 4 | 4 (2e-3) | 144 | 8 (0.25) | 4 | 8 (0.25) |
| icosphere | 20480 | 4 | 4 (6e-4) | >200 (3e-3) | 7 (0.29) | 4 | 7 (0.29) |
| icosphere | 81920 | 4 | 3 (2e-4) | >200 (1e-2) | 4 (0.52) | 4 | 4 (0.52) |
| torus | 2112 | 6 | 6 (7e-4) | >200 (2e-4) | 12 (0.24) | 6 | 12 (0.24) |
| torus | 10000 | 6 | 5 (2e-3) | >200 (6e-3) | 11 (0.26) | 6 | 11 (0.26) |
| torus | 40000 | 6 | 6 (2e-3) | >200 (1e-2) | 8 (0.44) | 6 | 8 (0.44) |
| torus | 102400 | 7 | 6 (1e-2) | >200 (1e-2) | 6 (0.63) | 7 | 6 (0.63) |
| stretched bumpy torus | 2112 | 13 | 11 (7e-4) | >200 (2e-3) | 16 (0.20) | 13 | 16 (0.20) |
| stretched | 10000 | 14 | 11 (9e-4) | >200 (7e-3) | 13 (0.32) | 14 | 13 (0.32) |
| stretched | 40000 | 15 | 10 (4e-3) | >200 (8e-3) | 12 (0.35) | 15 | 12 (0.35) |
| stretched | 102400 | 16 | 10 (4e-3) | >200 (2e-2) | 11 (0.48) | 16 | 11 (0.48) |
| doubletorus10k | 10506 | 16 | 13 (3e-3) | 108 | 13 (0.51) | 16 | 13 (0.51) |
| doubletorus35k | 35454 | 16 | 14 (1e-3) | 103 | 14 (0.40) | 16 | 14 (0.40) |
| hopf (2 comp.) | 9216 | 12 | 10 (5e-3) | 57 | 14 (2.2) | **>200 (9e-4)** | 6 (0.67) |
| borromean (3) | 11520 | 37 | 10 (9e-2) | >200 (3e-3) | 10 (0.20) | 40 | 10 (0.20) |
| double torus + ring (2) | 14464 | 16 | 11 (1e-2) | >200 (6e-2) | 15 (0.74) | 16 | 14 (0.76) |

Vertex pairs, right: 4–16 iterations at every size (no growth 2k → 100k), borromean 35–39 (also at
χ = 0 and 0.25: a property of the scene, not of the compression). Face pairs (D32 confirmed at
10k–100k): right GMRES does not reach 1e-4 in 200 iterations on most meshes; left GMRES (the
reference) stops after 4–16 iterations with true residuals 0.2–2.2. The congruence form fails on
Hopf (per-component barycenters carry a force: D33). With area + volume rows (Schur, 3 solves):
torus 10k 9+9+6, hopf 14+14+12, torus 40k 10+10+6, stretched 40k 26+27+15, doubletorus35k 20+20+16.
χ = 0.25 / 1 at 40k: the same iterations (torus 6/6, stretched 15/15).

**(c) Barycenter modes** (one direction; `Global` = `Auto` with `B⁰` on; iterations to 1e-4):

| scene | comps | PerComponent | Global ε_multi 1e-2 / 1e-3 / 1e-4 / 1e-6 (rel.) | ‖x_G − x_PC‖/‖x_G‖ | ‖x_G − (x_PC + shift)‖/‖x_G‖, cos | max net force/‖g‖ |
|---|---|---|---|---|---|---|
| hopf | 2 | 12 | 52 / 13 / 8 / 11 | 0.46 | 0.46, 0.889 | 2e-4 |
| borromean | 3 | 37 | 49 / 40 / 38 / 46 | 0.67 | 1.9e4, 0.040 | 1.23 |
| double torus + ring | 2 | 16 | 24 / 12 / 16 / – | 0.092 | 0.15, 0.989 | 8e-4 |

ε_multi 1e-1: no convergence in 200 on all three; ε 1e-10 (no ε_multi): diverged (true residual 49 on
Hopf). The direction differs materially between the modes (Hopf 46 %, borromean 67 %): Global lets
the components move relative to each other through the metric; PerComponent freezes the relative
barycenters, and the paper's shift (Sec. 7.1.1) adds an L² mean force whose scale is unrelated to the
`Hs` direction (borromean: the shift is 2e4 × ‖x‖, the direction becomes a rigid translation).
With ε_multi = 1e-4 Global costs about the same as PerComponent → `Auto` (Global with `B⁰`) is
affordable; default ε_multi_rel changed to 1e-4.

**D7 (Willmore term in the preconditioner)**, torus 10k / stretched 10k, iterations: w_W = 0.01:
with 154, without 6; w_W = 1: with diverged (0.57), without 5 / 11; w_W = 100: with diverged, without
13 / 20. → `willmore_in_preconditioner: false` (default); the term stays in `A′`.

**(d) Step breakdown, 4 threads, alone** (one direction with barycenter + area + volume: 3 GMRES
solves; best of 3; seconds):

| mesh | faces | BH build | energy+grad | BCT build | L′ factor | L^V tree | GMRES (iters) | products (n, ms each) | precond (ms each) | total |
|---|---|---|---|---|---|---|---|---|---|---|
| torus | 40000 | 0.05 | 0.21 | 0.46 | 0.03 | 0.08 | 1.04 (10+10+6) | 29, 22 | 13 | 1.91 |
| torus | 45000 | 0.06 | 0.26 | 0.57 | 0.04 | 0.10 | 1.15 (10+10+6) | 29, 25 | 14 | 2.21 |
| torus, 1 thread | 45000 | 0.11 | 0.59 | 1.10 | 0.04 | 0.18 | 3.30 | 29, 80 | 34 | 5.33 |
| torus, 8 threads | 45000 | 0.06 | 0.19 | 0.46 | 0.04 | 0.09 | 0.92 | 29, 17 | 14 | 1.78 |
| torus χ = 1 | 45000 | 0.06 | 0.27 | 0.26 | 0.04 | 0.06 | 0.83 (11+11+7) | 32, 12 | 13 | 1.57 |
| torus χ = 0.25 | 45000 | 0.06 | 0.25 | 1.46 | 0.04 | 0.26 | 2.40 | 29, 65 | 18 | 4.56 |
| stretched | 40000 | 0.06 | 0.27 | 0.58 | 0.03 | 0.10 | 2.92 (26+27+15) | 71, 26 | 14 | 3.91 |
| stretched χ = 1 | 40000 | 0.06 | 0.27 | 0.26 | 0.03 | 0.05 | 1.80 | 71, 12 | 12 | 2.49 |
| doubletorus35k | 35454 | 0.05 | 0.25 | 0.63 | 0.04 | 0.11 | 2.56 (20+20+16) | 59, 34 | 9 | 3.67 |

The paper's floor is ≈ 2.5 s per *full* step at 45k faces on 4 threads (Fig. 13; with line search and
remeshing). Our direction alone is 2.2 s at 45k on a torus (1.9 s at 40k) with two Schur rows, so the
full step will be at or slightly above the floor; with one solve (barycenter only) the GMRES part is
~0.4 s. The BCT build (curves crate: near-field rows, single-threaded flattening) and the products
(memory-bound: 18M near entries × 12 channels, faces in original order) dominate and scale poorly
past 4 threads. Levers, not done: renumber faces/vertices in BVH order (locality of the near-field
gathers), χ = 1 (−30 % time, 2.4–3.5e-3 product error, 1–1.5 % direction error), or a parallel
flattening in the curves crate's BCT build.

**(e) Dense vs hierarchical direction** (dense = `dense_metric` + per-axis LU + solve; hierarchical =
build + GMRES; 4 threads, alone; seconds): 512 faces 0.073 / 0.030; 960 0.109 / 0.034; 1500 0.161 /
0.089; 2000 0.45 / 0.17; 3072 0.78 / 0.19; 5120 1.66 / 0.09–0.2. Hierarchical is faster from ~500
faces on; the dense path is exact (no χ, no GMRES tolerance: hierarchical directions differ by
1e-3–2e-2).

**Recommended defaults (proposed SPEC text, §6.3 / §8.3 / §14.3):**
- `dense_max_faces = 1000` (dense costs ≤ 0.11 s there and is exact; above, hierarchical is 3–8×
  faster and the dense memory grows as `|F|²`). Code default is still 2000 (`DEFAULT_DENSE_MAX_FACES`)
  until the SPEC decides.
- `χ = 0.5` (product error ≤ 1e-3, direction error ≤ 2 %; χ = 1 is 30 % faster per direction at
  2–3× the error; χ = 0.25 doubles the time for no fewer iterations).
- GMRES: right side, tolerance 1e-4 on the true residual, restart 30, cap 200 (vertex pairs need 4–16,
  borromean 37); left only as the reference mode (its stop hides true residuals up to 0.1).
- `ε_rel = 1e-10` (single-component and per-component cases: iterations unchanged, `L′` solve
  residual ≤ 1e-15 after one refinement step); `ε_multi_rel = 1e-4` (not the reference's 1e-2).
- `FractionalForm::VertexPairs` (D32 confirmed to 100k faces) and the augmented form (D33).
- D7: no Willmore term in the preconditioner.

## For WP5 (one step, and after remeshing)

Everything is rebuilt per step from the current positions; nothing carries over (remeshing needs no
special handling beyond rebuilding the `ConstraintSet` with explicit targets, WP2):

```rust
let pool = Arc::new(make_pool(threads)?);
let cq = set.jacobian(&positions, &set.positional_rows())?;
let mut solver = saddle_solver(pool.clone(), &surface, &positions, &exponents, &weights,
                               &MetricExtras { willmore, obstacle }, &cq, &backend)?;   // Box<dyn SaddleSolver + Send>
let schur = SchurComplement::new(&mut *solver, set.nonlinear_rows(), cn)?;  // then direction / correction (WP2)
```

To report the hierarchical details (products, preconditioner time, `ε`), call `backend.resolve(nf)`
and `hierarchical_saddle(…, Some(bh.face_tree().clone()))` directly (concrete type; reusing the
Barnes–Hut face tree saves its build when the leaf sizes agree). `weights = SobolevWeights { tp_weight,
low_order, willmore_weight }`; `willmore`: WP7's `W M⁻¹ W` (interior rows) or `None`; `obstacle`: WP7's
`FaceDiagonal` or `None`. `BarycenterMode::Auto` is resolved by WP5 when building the `ConstraintSet`.

## Remaining

- WP10: restate `rs-fractional-preconditioner` for the augmented form (above) if D33 is accepted.
- Obstacle `fr` diagonal in the preconditioner (reference) — not done; only matters with obstacles.
