# Speed-ups of the dense saddle solve and the multigrid projections (2026-11-03)

Working notes of the work package; the results of record go to `SPEC.md` §17t and the README.
Logs: `/var/data/discrete-geometry/ddg-tools/runs_tmp/rcspeed/` (`diag/`, `ab*/`, `fig20-*.csv`).
Binaries of the state before the package: `runs_tmp/rcspeed/bin/fig20-before`, `components-before`.

## 1. What dominated before (confirmed, branch tip 4fb4c8e)

`Backend::Auto` and the Fig. 20 "Barnes–Hut + dense" configuration already use the Schur
complement solve (`DenseMethod::Schur`, §17l item 5): Cholesky of `M = A + ρuuᵀ`, its explicit
inverse, `S = C M̄⁻¹ Cᵀ`, Cholesky of `S`, the LU only as fallback. The LU numbers of §17k
(87–96 % of an iteration) are history.

Traces (`--trace 1`, 1 thread, machine otherwise idle, 4 jobs in parallel, diagnostic only):

* **BH + dense, Freedman 4000:** 6.46 s/iteration; metric 4.22 s (65 %: Gram 1.25, Schur
  factorization 2.57 = Cholesky 0.49 + inverse 1.0 + `S` assembly 0.41 + Cholesky of `S` 0.46),
  projections 0.86 s (3 solves, 0.29 s each), direction solve 0.30 s, energy 0.52 s, line-search
  objective 0.33 s. A dense solve was 0.3 s: the products with `M⁻¹` and `A` (strided
  interleaved dot products, latency bound) and the refinement step.
* **Multigrid, Freedman 16000:** 5.33 s/iteration (8 steps); direction solve 1.39 s (26 %),
  energy 1.09 s (21 %), projections 1.39 s (26 %; most steps project with Newton-polish steps
  only and cost 0.25–0.8 s, step 1 6.5 s), line-search objective 0.50 s, CCD 0.23 s. Since
  WP14b the projections are no longer 60–80 % on Freedman.
* **Multigrid, knot 96 at 8192:** 29–41 s/iteration, **projections 59–65 %** (17–27 s: 4–5
  multigrid solves of 4–5 s each per step; the trial at `2 τ_prev` is given up after one chord
  solve, the accepted one at `τ_prev/2` takes 3 chord solves and a Newton step).

## 2. Dense solve

* `(M⁻¹ ⊗ I_d) x` and `(A ⊗ I_d) x`: de-interleave `x`, one pass of contiguous dot products per
  row with 8 fixed lanes (vectorizable, same order at every thread count). One solve with
  refinement at 4000: 0.30 → 0.10 s; projections of a step 0.90 → 0.33 s.
* `S` assembly: `g = M̄⁻¹ cᵢᵀ` per coordinate (contiguous axpys), lower triangle only (the
  Cholesky reads only that; the upper starts as NaN so a reader would fall back to the LU,
  checked by `tests/dense_schur.rs`). Factorization at 4000, 1 thread: 2.57 → 2.25 s; 16
  threads 0.78 → 0.73 s. Result hash identical at 1 and 16 threads.
* Not changed: the explicit inverse (needed: every entry of `M⁻¹` enters `S`), one refinement
  step per solve (robustness).

## 3. Multigrid projections

(in progress)
