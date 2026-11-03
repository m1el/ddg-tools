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

What a projection solve costs now: one full-multigrid start (24 smoothing steps per level, about
one V-cycle) plus 1–2 FCG iterations, each a V-cycle; on knot 96 the level-0 product's near
field (`bct.rs`, `multiply_unchecked` step 4) is 43 % of the whole run (perf, knot 2048).

Paired A/B runs (1 thread each, all jobs at the same time, diagnostic only; `ab1/`–`ab4/`):

* **Recycled warm starts** (`MgOptions::recycle_projections`, kept): the projection solves of
  one step keep `(b, x)` orthonormalized; a new solve starts from `x₀ + C†(b − C x₀)` and
  stops on the cold solve's residual and energy targets (`vcycle::Reference`). Cycles per step
  fall about 2× (knot 2048 step 4: 12 → 5), time less, because every solve still pays the
  full-multigrid start: projections per step knot 8192 18.9 → 15.5 s (−18 %), knot 2048 −14 %,
  Freedman 2000 0.61 → 0.46 s (−25 %); Freedman 16000 unchanged (its projections are mostly
  Newton-polish steps since WP14b). Energies after 5–15 steps agree to 5–7 digits.
* **Near field of the block cluster tree** (kept): the kernel values of an entry spread over the
  12 channels, then a contiguous multiply-add (AVX via runtime detection, no FMA, same order):
  product bitwise identical, 70 → 60 ms (knot 8192, `theta_t 0.5`), solve 5.06 → 4.44 s.
* **Growth back-off** (`FlowOptions::growth_backoff`, kept): constant-τ phases had a rejected
  `2 τ_prev` trial in every step (`[1:1:MaxIterations 0.5:4:ok]` on the knot, `[1:2:armijo
  0.5:2:ok]` on Freedman). s/iteration (paired): knot 2048 2.17 → 1.93, knot 4096 8.32 → 6.89,
  knot 8192 24.8 → 23.2, Freedman 2000 1.90 → 1.68, BH + dense 2000 2.04 → 1.89, Freedman
  16000 5.98 → 5.55; energies at the last step equal or lower except Freedman 16000 (+1 % at
  step 15, one held step; the same energy at equal time).
* **Warm start of the direction solve from the previous direction** (tried, removed): the
  previous direction, projected and Galerkin-scaled, has a larger Euclidean residual than zero on
  Freedman 8000/16000 and knot 4096 (its rough part changes every step), so the cold path ran;
  forcing it (energy-optimal scale) left the cycles at 2 and the time unchanged.
* Not done: predicting the chord rate to skip a trial before its first solve (the back-off
  removes most of those trials with no model); a cheaper correction for linear rows (the
  default rows are all nonlinear: the barycenter is length-weighted, D6).

## 4. Results of record

In `SPEC.md` §17t (Fig. 20 before/after, 24-knot subset, thread scaling, bitwise check) and the
README. Summary: every Fig. 20 configuration 6–21 % faster on one thread (Barnes–Hut + dense 4000
5.76 → 4.91 s, knot 16384 28.5 → 22.9 s), 5–6 % at 4–64 threads; 24 knots × 200 steps: no stuck
step or isotopy violation, trials per step 1.93 → 1.45, lower energy at equal wall time on 23/24.
Defaults: `recycle_projections = true`, `growth_backoff = true`.
