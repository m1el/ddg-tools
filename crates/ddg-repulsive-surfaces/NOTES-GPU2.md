# GPU-2 notes (the preconditioner on the GPU)

Status: in progress (2026-11-03).

## Done / remaining

* Done: `src/solve/ldlt.rs`: faer's `LDLᵀ` of `L′` read out as a simplicial factor (order: AMD +
  multipliers as before, then a postorder of the elimination tree), triangular solves in a fixed
  chunked order (`CHUNK = 32` terms summed from 0, chunk sums subtracted in sequence), the GPU
  schedule (chains of the elimination tree, levels) and its CPU emulation (`solve_scheduled`, bitwise
  equal to `solve_rows`, unit tests). `LaplacianSaddle` uses it by default
  (`TriangularSolve::Chunked`; `TriangularSolve::Faer` keeps the old path).
* Remaining: CUDA kernels, the preconditioner pipeline on the device, GMRES vectors, measurements.

## Measurements so far (`examples/gpu2_ldlt.rs`, logs `runs_tmp/gpu2/struct*.log`)

| | 45k (fig12:3) | 135k (tunnel8:135000) |
|---|---|---|
| saddle n (per axis) | 22,425 | 67,491 |
| nnz(L) (incl. explicit zeros of relaxed supernodes) | 849k | 3.07M |
| etree height | 994 | 1,471 |
| chains / chain levels / Σ longest chain per level | 13,123 / 30 / 1,207 | 39,711 / 35 / 2,293 |
| factor: chunked (faer + read-out + postorder) / faer | 37 / 28 ms | 131 / 93 ms |
| solve, 3 columns, 1 thread: chunked / faer | 2.30 / 3.31 ms | 7.9 / 16.8 ms |
| max relative difference chunked vs faer | 1.9e-11 | 2.0e-11 |

The chunked CPU solve is faster than faer's supernodal one (row-major `n × 3`, one pass for the
three axes), so it can be the CPU default and determinism level 2 keeps holding.
