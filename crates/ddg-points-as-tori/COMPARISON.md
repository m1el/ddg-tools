# Points as Tori: ours vs the authors' program (WP14, SPEC.md §16.1, §16.7)

This compares our implementation (`ddg-points-as-tori`) with the authors' program (points-as-tori,
MIT, Nicole Feng), run as a black box on identical inputs. Their build and venv are at
`/var/data/discrete-geometry/PointsAsTori-build/` (recipe in `runs_tmp/pat-ref/README.md`; no
patches). Our driver is `scripts/pat_eval_ref.py`, which calls their public API
`PointsAsTori(points, normals)`, `.signed_distance(G, accelerate=True|False)`, `_tdf.get_tori()`, and
`model(P[idx], N[idx])` for the raw coefficients. Everything ran on the CPU on both sides (JAX CPU,
f32 network; their C++ in f64). The GPU belonged to the training run.

## Inputs

| set | clouds | points | source |
|---|---|---|---|
| WP14 evaluation clouds | 366 | 512 (342 clouds) and 2048 (24) | every cloud of EVALUATION.md: ABC 100, blobs 100, Thingi10K 45, analytic 70, reference meshes 3, saddle sets 2 × 11, App. E saddle shapes 2 × 13 |
| golden set | 12 | 512–12,917 | `runs_tmp/pat-ref/cpu_run1` (the reference agent's run, WP1's parity test) |

Both programs get the same `points.npy`/`normals.npy` (float64) and the same 64³ grid,
`meshgrid(linspace(−1, 1, 64), …, indexing='ij')`. On our side these are the crate defaults
(`SignRule::Reference`, 32-NN `PerQuery`, and `Neighborhood::All` for `accelerate=False`), which is
what reproduces the reference.

## Results (all 378 clouds)

Maximum over the clouds of each set, with the median of the per-cloud medians in brackets:

| set | clouds | coefficients |Δa|/(1+|a|) | torus center (rel.) | radii (rel.) | axis (sin of angle) | SDF 32-NN |Δφ| | SDF all |Δφ| | kNN rows differing | minor-radius sign flips |
|---|---|---|---|---|---|---|---|---|---|
| ABC | 100 | 3.3e-6 (6e-8) | 1.1e-5 (5e-8) | 1.9e-5 (4e-8) | 2.4e-4 (4e-7) | 1.4e-6 (9e-9) | 9.0e-7 (6e-9) | 0 | 0 |
| blobs | 100 | 1.7e-6 (5e-8) | 1.0e-4 (1e-7) | 1.3e-4 (1e-7) | 2.3e-4 (4e-7) | 8.0e-7 (1e-8) | 7.5e-7 (8e-9) | 0 | 0 |
| Thingi10K | 45 | 2.4e-6 (7e-8) | 8.1e-6 (5e-8) | 4.8e-6 (4e-8) | 2.9e-4 (4e-7) | 1.1e-6 (1e-8) | 9.6e-7 (8e-9) | 0 | 0 |
| analytic | 70 | 2.4e-6 (5e-8) | 4.4e-4 (8e-8) | 4.4e-4 (1e-7) | 9.2e-4 (5e-7) | 1.8e-7 (6e-9) | 1.3e-7 (5e-9) | 0 | 0 |
| reference meshes | 3 | 1.1e-6 (6e-8) | 9.5e-7 (5e-8) | 7.3e-7 (5e-8) | 1.7e-5 (3e-7) | 9.7e-8 (1e-8) | 8.3e-8 (7e-9) | 0 | 0 |
| saddle set 512 / 2048 | 22 | 2.8e-6 (6e-8) | 2.0e-5 (1e-7) | 2.3e-5 (1e-7) | 1.2e-3 (6e-7) | 1.0e-6 (5e-9) | 6.7e-7 (3e-9) | 0 | 0 |
| App. E saddle 512 / 2048 | 26 | 3.0e-6 (6e-8) | 8.7e-4 (2e-7) | 8.8e-4 (2e-7) | 9.2e-4 (7e-7) | 1.8e-7 (5e-9) | 1.2e-7 (3e-9) | 0 | 0 |
| golden 12 | 12 | 2.0e-5 (6e-8) | 1.1e-3 (1e-7) | 1.3e-3 (1e-7) | 8.6e-4 (6e-7) | 2.6e-7 (4e-9) | 1.5e-7 (2e-9) | 0 | 0 |

Center and radii are compared relative to 1 + R + |r|. Tori with R or |r| > 10³ (near-flat
neighborhoods: 2 per App. E cloud set, 4 on the golden set) are compared only through the SDF, as the
reference README recommends. The largest SDF difference on any cloud is 1.43e-6 (`abc55t3_134`); p99.9 is
≤ 1e-6 everywhere. The 64-NN sets agree on every row of every cloud. No torus flips the sign of its
minor radius.

**Consequence for EVALUATION.md:** the "authors' program" rows equal the PAT rows to 4 digits in
every table, with the same catastrophic-failure counts. What EVALUATION.md reports about the
method's accuracy, the near-sphere tori, the saddle statistics and the continuity of the blend
therefore describes the authors' program as released, not an artifact of the port.

## Where the differences come from

- **Coefficients** (≤ 3e-6 relative at 512–2048 points, 2e-5 on the golden set). Our CPU network
  (burn 0.22, f32) and their JAX CPU network (XLA, f32) sum the matmuls in different orders. WP1
  measured raw outputs within 6.8e-6 and identical σ to 2 ulp. The golden set's 2e-5 is femur,
  which has the largest activations (WP1, "default feature order").
- **Tori** (center and radii up to ~1e-3 relative). A torus radius is 1/|κ|, so near-flat
  neighborhoods amplify the coefficient round-off. The worst cases are points with |κ|σ ~ 1e-3
  (R ~ 10²–10³). The axis angle reaches 1e-3 rad where κmin ≈ κmax (near-umbilic: the principal
  direction is ill-conditioned). Most tori are near-umbilic (EVALUATION.md "Sign rule"), which is why
  the axis has the largest relative spread. None of this reaches the SDF beyond 1.4e-6: the blend
  weights put such points far away, and an ill-defined axis of a near-sphere barely changes its
  distance.
- **SDF** (≤ 1.4e-6). Coefficient round-off propagated through the fit. Given the same tori, our
  blend equals theirs to ≤ 1.8e-12 on the grid (WP1, `tests/parity.rs`). The f64 arithmetic is the
  same; only the summation order differs for `All`, where we sum in index order as they do.

No discrepancy is unexplained.

## Timings (CPU only, measured under load)

These were measured on the 315 512-point clouds of the four datasets, in the same time window. Ours: 4 threads
(rayon). Theirs: `OMP_NUM_THREADS=4`, `XLA_FLAGS=--xla_cpu_multi_thread_eigen=false
intra_op_parallelism_threads=4`. Their process used about 4–5 cores at times (JIT and Eigen
helpers), so this is not a strict single-core comparison. The machine had an outside load of 40–100
(other agents, the GPU training run's host threads), so absolute numbers are inflated for both. All
values are medians over clouds.

| stage | ours | authors' | ratio |
|---|---|---|---|
| build: kNN + network + fit (theirs: warm second constructor call, no JIT) | 1.84 s | 1.97 s | 1.1× |
| 64³ grid (262,144 queries), 32-NN | 0.32 s | 0.73 s | 2.3× |
| 64³ grid, all points | 1.54 s | 6.71 s | 4.4× |

The network dominates the build on both sides. Our CPU network ran at ~280 neighborhoods/s at
4 threads under this load. WP1 measured 147/s per thread on a quiet host, so the load cost about
half. Quiet-host and GPU timings belong to WP16 (§16.5). The golden set's own timings (theirs, all
64 cores, load ~20–30) are in `runs_tmp/pat-ref/tables.md`. They are not comparable with ours here
because the thread counts and load differ.

## Reproduction

```
JAX_PLATFORMS=cpu OMP_NUM_THREADS=4 XLA_FLAGS="--xla_cpu_multi_thread_eigen=false intra_op_parallelism_threads=4" \
  uv run --no-project --python /var/data/discrete-geometry/PointsAsTori-build/.venv/bin/python -- \
  python -I crates/ddg-points-as-tori/scripts/pat_eval_ref.py <file listing case directories>
target-pat14/release/examples/eval case <case-dir> --threads 4     # writes result.json with "black_box"
```

`scripts/pat_eval.sh ref case report` runs both for every case. The `black_box` object of each
`result.json` holds the per-cloud statistics above, and `report/tables.md` has the aggregate and the
largest differences. The golden 12 are compared from `runs_tmp/pat-ref/cpu_run1` (copied into
`cases/reffull/<name>/ref/` by `eval prepare reffull`).
