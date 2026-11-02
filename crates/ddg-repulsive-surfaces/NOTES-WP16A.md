# WP16a notes (consistency study, paper Sec. 8.1 / Fig. 10; SPEC 14.1)

Files: `examples/consistency.rs`, `scripts/rs_consistency_plot.py` (tables, fits, PNG; `uv run`). Logs (not committed):
`runs_tmp/wp16a/*.log`, tables `runs_tmp/wp16a/tables.txt`, plot `runs_tmp/wp16a/consistency.png`.
Build: `cargo build -p ddg-repulsive-surfaces --release --example consistency --target-dir target-wp16a`.
Modes: `check`, `smooth SURF A B XRULE DRULE THREADS`, `discrete SURF FACES THREADS [THETAS]`; `SURF` = `pert` | `torus`.

## Test surface
`f = f0 + Phi(f0)`, `f0` the standard torus R = 1, r = 1/3 (the paper's printed `f0` has typos), `Phi` = WP12's
`PerturbedTorus` with degree 2, amplitude 0.2, seed 16 (our choice; the paper gives no `Phi`). The example rebuilds the
terms of `Phi` exactly as the library does (asserted against `PerturbedTorus::point`) to get **analytic** `f_u, f_v`
(`f_u = (I + DPhi) p_u`); checked against central differences (2e-12). Embedded (brute-force sweep: 0 intersections on
a 300x120 sampling), Jacobian |f_u x f_v| in [0.118, 0.606], area 13.163. The plain torus `Phi = 0` is the sanity case.

## Smooth energy (cubature)
`E = int int |N(x).(f(x)-f(y))|^a |f(x)-f(y)|^-b dA dA`, `dA = |f_u x f_v| du dv`, `N = f_u x f_v/|f_u x f_v|`, over
`[0,2pi]^4` (ordered pairs, both orders, same convention as the discrete energy). Substitute `y = x + delta`,
`delta in [-pi,pi]^2`: for fixed `delta != 0` the integrand is smooth and periodic in `x`; the only non-smoothness is at
`delta = 0` (bounded, direction-dependent limit).
* x rules: `trap:M` (M x M trapezoid, spectrally accurate for smooth periodic data), `gl:P:n` (P x P panels of n-point
  Gauss-Legendre; the independent x rule), `traps:M` (torus of revolution only: u-translation symmetry, `u = 0`, weight 2pi).
* delta rules, both parallel (rayon, explicit pool) over x:
  * **Duffy** `duffy:Pr:nr:nw`: each quadrant = 2 triangles `(r, r w)`, `(r w, r)`; Jacobian `r` removes the direction
    dependence, the integrand is smooth in `(r, w)`; GL panels in r.
  * **Cells** `cells:D:n`: the cells touching the diagonal are subdivided geometrically (D levels, L-shape of 3 squares per
    level per quadrant plus the corner cell of size pi/2^D), plain tensor n x n GL on every cell, no change of variables.
    The error is the corner cell (area ~ 4^-D times a bounded integrand) plus GL error on the other cells.
  The two delta rules share nothing but the integrand function.

Convergence (the (6,12) energy of the perturbed torus; x = trap:72 unless stated):

| x rule | delta rule | E |
|---|---|---|
| trap:64 | duffy 4:12:12 / 6:16:16 / 8:20:20 | 460.7158077 / 460.7162350 / 460.7162243 |
| trap:48/64/80/96 | duffy 6:16:16 | 460.71620635 / 460.71623496 / 460.71623506 / 460.71623506 |
| trap:72 | duffy 10:22:22 / 12:24:24 | 460.71622437 / **460.71622436** |
| trap:72 | cells 14:12 / 18:16 / 22:20 | 460.71802 / 460.71573 / 460.7162070 |
| gl:6:12 | cells 22:20 / duffy 12:24:24 | 460.7161985 / 460.7162156 |

Results of record (value; agreement of the two independent delta rules, relative):

| surface | (a, b) | E | Duffy vs cells |
|---|---|---|---|
| perturbed torus | (6, 12) | **460.7162244** | 3.8e-8 (cells 460.71620703; with x = GL: 460.7161985 vs 460.7162156) |
| perturbed torus | (4, 8) | **96.40066290** | 1.3e-9 (x = trap:96 vs trap:72: < 1e-15) |
| torus of revolution | (6, 12) | **103.7704318369** | 4e-10 (Duffy 8:20:20 = 12:24:24 to 1e-12; traps:128 = 256; full 24x24 trapezoid grid = symmetric rule to 1e-15) |
| torus of revolution | (4, 8) | **61.64509301** | 2e-11 |

So the (6,12) value is good to at least 7 significant digits, SPEC 14.1's 6 digits are met with margin. Cost: a few
minutes of CPU per big run (10^8..10^9 integrand evaluations, 8 threads, 1-2 min). The Duffy value is the reference.

## Discrete energies
Sample `PerturbedTorus::sample(3 nv, nv, 0.4, 0.4)` (parameter grid warped by the sine squeeze), `L0` from the target
face count (equilateral triangles of area A/F), then 4 x (WP4 `remesh` with defaults; per-vertex parameters transported
as (cos, sin) pairs through `vertex_map`; projection back by damped Gauss-Newton on `|f(u,v) - p|^2` in parameter
space, gradient < 1e-14). Largest vertex move in projection 4.8e-2 (coarsest) .. 4.6e-4 (finest). Final face counts are
about 1.34x the target. `h` = longest edge. Exact energy `energy_value` (WP1), 6 threads at most: 85k faces 5.8 s.

Perturbed torus, relative error `|E_h - E|/E` (all errors positive; `E_h` converges from above):

| faces | h | E_h (6,12) | err (6,12) | err (4,8) |
|---|---|---|---|---|
| 648 | 0.2925 | 1374.36 | 1.983 | 0.3256 |
| 1332 | 0.2166 | 947.32 | 1.056 | 0.1618 |
| 2646 | 0.1575 | 716.53 | 0.5553 | 0.08389 |
| 5360 | 0.1176 | 610.13 | 0.3243 | 0.04795 |
| 10796 | 0.0856 | 535.80 | 0.1630 | 0.02524 |
| 21498 | 0.0640 | 502.35 | 0.09037 | 0.01400 |
| 42812 | 0.04486 | 482.98 | 0.04832 | 0.007636 |
| 85554 | 0.03257 | 473.05 | 0.02677 | 0.004184 |

Fitted rate (log error vs log h): (6,12) 1.97 over all eight, 1.86 over the last four; (4,8) 1.97 / 1.84. Plain torus:
2.05 / 1.93 (6,12), 2.00 / 1.92 (4,8). So between O(h) and O(h^2), at the upper end, as in the paper. The rate drifts
slightly below 2 for the finest meshes (the remeshed meshes are not a nested family; no sign of O(h)).

## Barnes-Hut error vs theta
Perturbed torus (6,12), leaf size 8, signed error `(E~ - E_exact)/E_exact` of the same mesh (theta = 0 is bit-identical to the exact energy):

| theta | 10796 faces | 42812 faces | 85554 faces |
|---|---|---|---|
| 0.0625 | -3.2e-7 | -1.5e-5 | -3.8e-5 |
| 0.125 | -3.8e-5 | -2.0e-4 | -1.9e-4 |
| 0.25 | -4.3e-4 | +1.2e-5 | +9.5e-4 |
| 0.3536 | -6.6e-4 | +1.5e-3 | +3.5e-3 |
| 0.5 | +1.1e-4 | +4.8e-3 | +8.2e-3 |
| 0.7071 | +2.3e-3 | +1.3e-2 | +1.7e-2 |
| 1 | +1.1e-2 | +2.3e-2 | +2.9e-2 |
| 1.414 | +2.1e-2 | +4.3e-2 | +5.7e-2 |

Fitted order on theta >= 0.35 (only where all errors are positive): (6,12) 2.41 (42812), 1.99 (85554); (4,8) 2.03 (42812), 1.88 (85554); 3.89 for (4,8) on 10796 faces, which has a sign change at 0.35 (not meaningful). In short: **O(theta^2) for theta >= 0.35**
(0.25 on the 42k mesh). At small theta the error is *negative*, grows with the mesh resolution, and the curve is not
monotone: a second (negative) error term, ~ -2e-4 at theta = 0.125 on the 42k and 85k meshes, which does not vanish as fast as
theta^2 (at theta = 0.125: 3.8e-5, 2.0e-4, 1.9e-4 for 11k/43k/86k faces). The absolute error at theta = 0.125..0.25 is therefore
not a clean O(theta^2) and the signed error crosses zero (which makes `|error|` plots, and the fit, depend on theta range).
Cause not investigated (both terms are midpoint-rule effects of lumping near-field clusters).

## Surprises
* The theta error has a sign change; the pure O(theta^2) regime starts only at theta >= ~0.35.
* The Duffy substitution made the "diagonal" singularity a non-issue: bounded cells with geometric refinement need depth 22
  to reach 4e-8 where Duffy has 1e-10 with 12 x 24 nodes.
* (4,8) converges at the same rate as (6,12) with 5x smaller errors.

## Proposed text for SPEC 14.1 "Results of record" (not applied)

> **Results of record (WP16a).** Perturbed torus: R = 1, r = 1/3, `Phi` = `PerturbedTorus(degree 2, amplitude 0.2, seed 16)`
> (embedded, area 13.1631). Smooth energy by our own cubature over `[0,2pi]^4` (offset coordinates, trapezoid in the base
> point, Duffy / geometric-cell rules for the diagonal; two delta rules agree to 4e-8): `E(6,12) = 460.71622`, `E(4,8) = 96.400663`;
> torus of revolution `E(6,12) = 103.770432`, `E(4,8) = 61.645093`. Discrete energy on squeezed grids remeshed 4 passes with
> projection, 648 .. 85 554 faces (h 0.29 .. 0.033): relative error 2.0 .. 2.7e-2, fitted rate 1.97 (1.86 over the last four),
> torus 2.05; (4,8) 1.97. Barnes-Hut (leaf 8): error ~ theta^2 for theta >= 0.35 (fits 1.9 .. 2.4); for theta <= 0.25 the signed error
> changes sign and is not monotone (|error| <= 2e-4 at theta = 0.125, 4e-5..4e-4 at 0.0625..0.25 on 11k..86k faces).
