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
  { eps_rel 1e-10, eps_multi_rel 1e-2, fractional Vertices, cluster, form Augmented }`,
  `translations_fixed_per_component` (decides `ε` vs `ε_multi`), `SparsePreconditioner`
  (`rs-fractional-preconditioner`; `usage()`, `factor_seconds()`, `epsilon()`), `PreconditionerForm`.
- `src/solve/iterative.rs`: `MetricExtras { willmore: Option<SparseMatrix>, obstacle: Option<FaceDiagonal> }`
  (WP7 hooks), `HierarchicalOperator` (`A′`, `LinearOperator`, `apply_metric`, `usage()`),
  `HierarchicalOptions { cluster, gmres, preconditioner, willmore_in_preconditioner: true }`
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

## Measurements

(filled below)

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
