# Investigating the camelshotnoise reproduction

Research experiment, 2026-10-13. This does not change either production solver.
The authors' actual recipe is still unknown. These experiments distinguish
reproducing their cone file from obtaining comparable UV area distortion.

## Objectives and evaluation conventions

The existing cone optimizer minimizes `E = 1/2 sum_i m_i u_i^2` with
`L u = mu - Omega`. It is a convex discrete problem with a fixed cotan matrix.

The new experiment instead evaluates a BFF map `X(mu)` on the **original** mesh:

```
q_t(mu) = 1/2 log(|signed_area_t(X(mu))| / A_t)
w_t = A_t / sum_t A_t
D*(mu) = sum_t w_t (q_t(mu) - sum_s w_s q_s(mu))^2
minimize D*(mu), subject to sum_i |mu_i| <= Phi
```

This is a local, nonlinear outer optimization with a BFF solve at each objective
evaluation. It does not inherit the cone paper's convexity guarantee. It uses
SLSQP, finite differences of step `1e-4`, fixed cone signs and positions,
angle magnitudes between `0.01` and `6.1`, and the published total angle as `Phi`.
It permits flips; absolute area is a distortion statistic, not an injectivity
guarantee. Every candidate reports the flip count and rejects degenerate faces.

Two scores are reported:

* **All faces:** the objective above, including the absolute area of reversed
  faces. This is the primary optimization objective.
* **Excluding flips:** the existing `published_uvs` experiment's convention:
  discard reversed faces, optimize the mean over the remainder, but retain
  `A_t / total_original_area` as the weights. This is useful for comparison
  with the previous README tables, but is unsafe as a sole quality criterion.

## What happens to the needles in the published UVs

The mesh has 22 needle tips (interior angle defect > 6 radians), touching 132
triangles. These triangles occupy **0.25871%** of the original area. In the
published UV map:

* 68 of those 132 triangles are reversed; there are 70 reversed faces in total.
* The median ratio of principal stretches on the needle triangles is **6.591**;
  the maximum is **354.57**. The median elsewhere is **1.094**.
* Their face log-area scales range from **-6.294 to -3.259**, median **-4.408**.
  This differs radically from the Poisson log scales at their tips, about
  **-12,581 to -17,030** for the published cone angles and fixed boundary.
* Including all faces gives `D* = 0.09531245`, versus `0.06813801` when reversed
  faces are omitted. Needle faces contribute `0.05421298` to the former score;
  all reversed faces contribute `0.02714015` using the all-face mean.

Exponentiating the Poisson state into edge lengths produces 315 invalid
triangles, including **all 132 needle faces**. BFF does not realize that metric
triangle by triangle: it constructs boundary data and harmonically extends
coordinates. A small BFF area score therefore cannot be interpreted as a small
energy of the Poisson state, or as exact affine conformality on these faces.

The NumPy/SciPy operators were assembled independently from position cross
products. Their published-cone energy, `84784.6524257`, agrees with the existing
Rust diagnostic, as do the fixed/free-boundary angle refits below.

## Directly optimizing the measured UV score

First hold the published six positions fixed. This tests the angle objective,
not cone placement. Use `Phi = 20.06489` for camelshotnoise, `20.39531` for camel.

| mesh / candidate | all-face D* | excluding-flips D* | flips |
|---|---:|---:|---:|
| camelshotnoise: published UV | 0.09531245 | 0.06813801 | 70 |
| camelshotnoise: BFF with published angles | 0.09421870 | 0.06602276 | 72 |
| camelshotnoise: optimize all-face D* | **0.09175994** | **0.06344621** | 72 |
| camelshotnoise: optimize excluding-flips D* | 0.09176633 | 0.06343996 | 72 |
| camel: published UV | 0.03739924 | 0.03739924 | 0 |
| camel: BFF with published angles | 0.03731138 | 0.03731138 | 0 |
| camel: optimize all-face D* | **0.03410536** | **0.03410536** | 0 |

SLSQP reports successful termination in each run; this is not a global certificate.
The spiky all-face run used 217 objective calls. Its angles, in the published
file's order (zero-based vertices `246, 101, 854, 1575, 191, 864`), are:

```
4.77191691, 4.76826042, -2.39131868, 2.88397334, 2.37641678, -2.87300386
```

Their total absolute angle is exactly the budget to floating-point accuracy.
The largest deviation from the published angles is 0.22042 radians. Thus
optimizing the UV score obtains good results but does **not** reproduce the
published angle file, even when given its positions. On the clean control it
also finds a lower score and different angles.

## A strategy that does not use the published positions

1. For cone selection only, replace each detected needle tip by the mean of its
   original one-ring positions. Keep connectivity and all other positions.
2. Run the existing cone solver on this proxy at normalized `lambda = 0.0045`,
   with free boundary, same-sign merging and no cluster-size limit.
3. Refit the merged cones on the proxy with free boundary and no cone penalty.
   This produces seven cones. Retain the six largest in absolute angle.
4. On the **untouched original mesh**, refine these six angles with the all-face
   BFF objective above, under the published total-angle budget.

The selected zero-based vertices are `246, 101, 233, 2408, 2585, 1962`.
None is a needle tip. Two equal published positions exactly. Neither the
published positions nor individual angles initialize this experiment; the
comparison supplies only the target count of six and total-angle budget.
The lambda was selected from a sweep to obtain a small candidate set, and is
not a generally validated default.

The initial six-cone map scores `D* = 0.10300466` including all faces and
`0.07437954` excluding flips. SLSQP terminated successfully after 450 evaluations:

| candidate | cones | total absolute angle | all-face D* | excluding-flips D* | flips |
|---|---:|---:|---:|---:|---:|
| published | 6 | 20.06489 | 0.09531245 | 0.06813801 | 70 |
| proxy selection, UV-refined angles | 6 | **17.16845** | **0.09685539** | **0.06854184** | 74 |

This uses **14.4% less total angle**, with all-face distortion **1.62% higher**
and the earlier excluding-flips score **0.59% higher**. It is a comparable
six-cone result without knowing the published sites, not an exact reproduction.
The area objective alone does not eliminate folds: 71 of the 74 reversed faces
touch needle tips. No claim of local or global injectivity is made.

Final angles in the vertex order above:

```
4.1533507934, 4.1927511326, 2.3036643897, 2.9045275070, -2.2111913031, -1.4029651297
```

The result is a local solution; the angle cap is inactive. The original published
sites lead to a better score when directly optimized, as the preceding table
shows. Position refinement remains an opportunity for improvement.

## Alternatives that did not explain the file

For these experiments, hold the six published sites fixed, allow their angles
to vary freely, and solve both fixed and free boundary variants. Results below
are maximum absolute differences from the published angles, in radians.

| model | fixed boundary | free boundary |
|---|---:|---:|
| original quadratic energy | 39.21 | 44.07 |
| zero distortion weight at tips | 1.044 | 0.5514 |
| clamp negative **assembled edge** cotan weights to zero | 0.6610 | 0.3808 |
| flatten tip positions in the proxy, original quadratic energy | 0.6151 | 0.2599 |

Area-weighted Huber, soft-L1, Cauchy and squared-asinh losses were also tested,
each at scale `0.1`, `1`, `10`, using IRLS. The best free-boundary discrepancies
remain about 0.55 radians. IRLS runs report iteration count and iterate change;
the redescending losses are nonconvex, so these are local experiments, not
exclusion proofs for all scales/initializations. One fixed-boundary Cauchy run
at scale 0.1 hit 100 iterations without reaching the iterate tolerance.

Clamping cotans also worsens the clean camel's free-boundary fit, from 0.00343
to 0.30863 radians. It is not a convincing reconstruction of the authors' method.
Simply downweighting extreme `u` values cannot undo all of the effect the
needles have on the **state equation** and on scale factors elsewhere.

## A candidate intrinsic objective, before UV layout

There is a useful analytic counterexample to treating the Poisson state as an
accurate predictor of face-area distortion at extreme aspect ratios.

Take a regular six-triangle cone with base radius `rho`, height `h`, slant
length `l = sqrt(h^2 + rho^2)`. In a discrete vertex-scaled metric
`l'_ij = exp((s_i+s_j)/2) l_ij`, leave the base scale factors at zero and set

```
s_tip = 2 log(rho/l).
```

The new tip edges have length `rho`: the fan is a flat disk, has positive
triangle areas, and its tip curvature is zero. Its face area ratio is
asymptotically proportional to `rho/h`. Consequently its contribution to a
face-area objective behaves as

```
D_area,needle = O(rho h log^2(h/rho)),
```

which tends to zero as `rho -> 0` at fixed height. By contrast the Poisson
quadratic state energy grows like `h^3/rho`.

The `needle` experiment verifies zero flattened apex curvature to `4e-15`.
With height 0.5 and one additional unit of unchanged area:

| h/rho | Poisson vertex D (no 1/2 factor) | exact metric face-area D* |
|---:|---:|---:|
| 10 | 8.4459 | 0.097717 |
| 100 | 106.7891 | 0.041662 |
| 1000 | 1093.7103 | 0.0093095 |
| 10000 | 10963.3107 | 0.0016404 |

This is **discrete metric conformality**, not angle preservation by the affine
map on each original triangle. In fact those affine maps become strongly
anisotropic. It is not a counterexample to the smooth cone-energy calculation.

A possible optimization is therefore

```
minimize_s  D_area*(areas of l'(s)) + lambda sum_interior |K_i(s)|,
```

with valid intrinsic geometry and suitable boundary treatment. Triangle areas
and curvature are computable from lengths before UV layout. The actual
nonlinear `K(s)` replaces `Omega + L s`. There is no established convexity
guarantee for this new objective.

A fixed-connectivity Newton prototype for prescribed published cone curvature
was tested on both camels. Neither converged: it approached degenerate
triangles, with maximum curvature residuals 4.881 radians (spiky) and 3.660
radians (clean). Directional finite-difference checks of its initial curvature
Jacobian agree to relative errors `1.0e-6` and `2.6e-9`, respectively. The intermediate energies are
**not feasible flattening results**. A credible implementation needs a method
that handles triangulation changes or another reliable nonlinear metric solver,
plus area integration across any changed triangulation. This remains a research
direction rather than a working general optimizer in this experiment.

## Reproduction

From `ddg-tools`, build only the new example:

```sh
CARGO_TARGET_DIR=target-camel-investigation cargo build --release -p ddg-uv --example camel_objective
```

The tested Python environment has NumPy 2.2.6 and SciPy 1.15.3. For example:

```sh
OPENBLAS_NUM_THREADS=1 UV_CACHE_DIR=/tmp/uv-camel-investigation uv run --no-project --with numpy==2.2.6 --with scipy==1.15.3 python scripts/investigate_camel.py probe camelshotnoise
```

Replace `probe` with `refit`, `uv-opt`, `proxy-scan`, `proxy-opt`, `nonlinear`,
or `needle`. Use `camel` as the name for the clean control. After saving the
`proxy-opt camelshotnoise` output to `runs_tmp/camel-objective-proxy-opt.log`,
`export camelshotnoise` writes `runs_tmp/camelshotnoise-investigation.obj` and
checks unchanged positions/connectivity and agreement with independent NumPy
area measurements. The research worker
uses the build directory above; neither program edits input dataset files.

Logs from this investigation are under `runs_tmp/camel-objective-*.log`.
`proxy-scan` includes both the coarse and finer grids used during investigation.
The UV refinement uses repeated full BFF solves and takes minutes; its purpose
is to test an objective, not establish a performance baseline.
