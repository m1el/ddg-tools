# WP3 notes: Lean oracle for Points as Tori

Code: `lean/oracle/DdgOracle/PointsAsTori/{Basic,Eval,Format,Gen,Selftest,Fixtures}.lean`; command `ddg-oracle pat`;
fixtures `lean/oracle/fixtures/pat` (35 cases); format and rules in `lean/oracle/README.md` ("Points as Tori").

## Agreement with the authors' program (`runs_tmp/pat-ref/cpu_run1`, 12 clouds, 34,222 tori)

Inputs: their `points.npy`, `normals.npy`, `net_coeffs.npy` (float32 widened), `queries_random.npy`. Scripts in the scratchpad
(not committed). Fit: oracle tori vs `tori_*`: all tori agree: centers, R, r to <= 6e-15 relative, signs identical, axes to 4e-10 (their
(theta, phi) round trip). No near-flat disagreement showed: only 4 tori (holey-bunny) have R or r above 1e3 and they agree too.
Blend (their per-query shift/lambda, `accelerate` = knn 32, and all points): oracle with their tori vs `query_sdf_*`:
|dSDF| <= 1.3e-14, |dgrad| <= 1.3e-12; with the oracle's own fitted tori the same bounds. So the oracle reproduces the reference to rounding.

## SPEC issues and choices

1. **PaperHeuristic overflow.** With constant `lambda = 1e3/D` and the 32-NN fallback (queries far from the surface), `exp(-lambda(rho - sigma_x))`
   overflows (`lambda sigma_x` ~ 1300 > 709): NaN. SPEC §7.4 says "sigma_x is still the shift". The oracle uses the shift `min rho` for the constant
   rules (nearest weight = 1, no overflow); `sigma_x` is still printed. Rust should do the same (or log-sum-exp).
2. **Radius sets with one point at distance 0.** `sigma_x = 0`: with PerQuery, `lambda = C/0` is a division by zero (error `DegenerateQuery`);
   with constant / paper lambda it is fine. SPEC §7.4 only says "W underflow cannot occur".
3. **SHC sum exactly 0** (e.g. the only set member is the query's own data point, skipped): `phi` undefined; the oracle reports `DegenerateQuery`.
4. **Laplacian at singular points** is not specified: the oracle returns 0 for kind 1 and 2 (axis, core).
5. **Counters.** SPEC lists `clamped_disc, clamped_eigvec, clamped_curvature, umbilic`; the denominators `ev.I.ev` have no counter and `umbilic` is undefined.
   The oracle prints five flags (`disc, eig, denom, curv, umbilic`); `umbilic` is defined as the disc clamp (`arg < eps`).
6. **Exact gradient ties.** `argmax` of rho is the last member of the `(rho^2, index)`-sorted set.
7. **Paper rule D.** `D` = mean over all points of the mean distance to the `min(k_mean, n)` nearest points, self included (the reference's
   `compute_optimal_radius` is racy and not used; SPEC leaves the definition open). `radius auto` = `R_eval = 2C/lambda`; the SPEC enum has `r` and
   `PaperHeuristic` separate, the oracle lets the fixture say which.
8. `pat-blend-laplacian` (App. D divergence) is not written out in the SPEC; the oracle derives it (formula in the README), checked by differences.
9. `sgn(0)`: JAX gives 0 at exact zeros of `ev` components; SPEC says use C++ `sgn+`; the oracle follows the SPEC (only matters at exact umbilics).
