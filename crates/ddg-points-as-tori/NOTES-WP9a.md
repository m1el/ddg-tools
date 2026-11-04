# WP9a notes: normal orientation (Fig. 27) and the classical 2D variant (Fig. 10)

Files: `src/apps/orient.rs`, `src/pat2d/{mod,fit,field,grid,output}.rs`, `tests/{orient,pat2d}.rs`, `examples/pat2d.rs`.
Build dir `target-pat9`, logs and pictures in `runs_tmp/pat-wp9a/` (`orient-sweep.tsv`, `sweep.log`, `pca.log`, `pat2d-measure.log`,
`pat2d-{heart,box,moon,pie}.png/.svg`). Everything was run with at most 16 threads while the machine had an outside load of 40-80.

## Orientation (SPEC §8.4)

API: `apps::orient::{orient_normals, OrientOptions, Epsilon, OrientReport, orientation_metric, orientation_votes, pca_normals_random_sign}`.
`orient_normals(positions, normals, model, opts)` returns the oriented normals and the report (flips per iteration, `negative_per_iteration`,
final count of `<n_i, grad phi(p_i)> < 0`, converged, total flipped, undecided). Each iteration builds the field with the network (the
facade), evaluates `grad phi` at all points in one batch, takes the 64 nearest neighbors (fixed, computed once), weights them with
`exp(-M_i(p_i - p_j))` (the smallest `M` is factored out: a positive factor, only the sign of the vote matters), and flips all points whose
vote disagrees (Jacobi). `kappa_min, kappa_max, v_min, v_max` come from `fit_torus` on the stored coefficients (same fit options as the build),
`sigma_i` from the network scale. `[formula: pat-orientation-metric]` on `orientation_metric`.

### Result (Fig. 27 protocol)

Clean clouds with true normals; flip a seeded random fraction; 10 iterations; defaults (k = 64, eps = 1/sigma_i, reference weights, 32-NN
`PerQuery` blend). 11 clouds: five WP8 analytic solids (seeds 3, 12, 26, 40, 52; 2048 points), and the golden clouds synth_torus, tooth,
knot, block, femur, snail. "Residual" = fraction of points whose final normal is opposite to the true normal.

| flipped | eps | mean residual | median | worst cloud |
|---|---|---|---|---|
| 30 % | 1/sigma (default) | 1.18 % | 0.00 % | femur 6.6 % |
| 50 % | 1/sigma | 4.5 % | 3.0 % | 12.9 % |

Paper: 30 % -> 0.3 % after 10 iterations (its clouds are not specified). We match that on 8 of 11 clouds (0.00-0.05 %: knot, block,
box/hexagon/pentagon revolutions, torus, circle extrusion converge in 1-2 iterations: `[633, 0]`) and are worse on the thin or
high-curvature ones: femur 6.6 % (13 points with `<n, grad phi> < 0` at the end), snail 2.0 %, arc extrusion 3.2 %, tooth 1.2 %. These do not
converge in 10 iterations (flips per iteration for femur: 1242, 88, 41, 19, 18, 22, 14, 12, 7, 6), so more iterations would help; it was not tried.
Failure mode: small patches that are consistently wrong. Inside such a patch every vote confirms the wrong sign, so it is a fixed point;
the final count of `<n_i, grad phi(p_i)> < 0` (femur 13, snail 9) is far below the number of wrong normals (257, 80). That count is a poor
measure of success.
Per cloud and per eps, all runs: `runs_tmp/pat-wp9a/orient-sweep.tsv`.

### eps sweep (SPEC §16.3), mean residual over the 11 clouds

| eps | 0 | 0.25/s | 0.5/s | **1/s** | 2/s | 4/s | 16/s |
|---|---|---|---|---|---|---|---|
| 30 % flipped | 1.27 % | 1.21 % | 1.23 % | 1.18 % | 1.22 % | 1.17 % | 1.55 % |
| 50 % flipped | 4.89 % | 4.86 % | 4.74 % | 4.52 % | 4.13 % | 4.16 % | 4.97 % |

The result is flat within the noise from 0.25/sigma to 4/sigma; `eps = 0` (the metric is then only a seminorm, no weight along the normal) is
slightly worse and 16/sigma is worse. 2/sigma-4/sigma is marginally better at 50 %, 1/sigma-4/sigma equal at 30 %. **Recommendation: keep
1/sigma_i** (the SPEC default); no single value is measurably better. One seed (11) per cell; differences below ~0.3 % are not significant.

### Unoriented estimates (PCA normals with random signs)

`pca_normals_random_sign(positions, 16, seed)`: smallest principal axis of the 16 nearest points, random sign (so about 50 % wrong). Then
`orient_normals` with defaults. Residual is counted up to a global sign (the data has no global orientation), 10 iterations:

| cloud | PCA axes more than 25 deg off | wrong before | after (global sign free) |
|---|---|---|---|
| synth_torus | 3 | 261 / 512 | 0 |
| knot | 346 | 981 / 2048 | 0 |
| box revolution | 161 | 1018 | 0 |
| pentagon revolution | 130 | 1051 | 0 |
| block | 577 | 1081 | 1 |
| arc extrusion | 332 | 1012 | 50 (2.4 %) |
| hexagon revolution | 124 | 1039 | 32 (1.6 %) |
| tooth | 132 | 1000 | 89 (4.3 %) |
| snail | 335 | 2052 / 4096 | 230 (5.6 %) |
| femur | 956 | 2008 / 3897 | 586 (15 %) |
| circle extrusion | 245 | 1009 | 354 (17 %) |

So the method copes with lines that are roughly right, but not as well as with clean lines plus flipped signs: where PCA is off at creases and
thin parts (femur, circle-extrusion caps), the tori are fitted to wrong normals, the votes are noisy, and flips continue to the last iteration.
It never re-estimates a line; a normal that is 60 deg off stays 60 deg off (or 120 deg after a flip). A normal refinement (re-estimate the
line from the fitted `n*`) is a natural extension, not tried. The paper only claims the flip case.

## Classical 2D (SPEC §8.7)

API: `pat2d::{Cloud2, ClassicalFit, QuadraticPredictor, Predicted2, Quad2, fit_circles, Fit2dOptions, Pat2d, BuildOptions2, Pat2dField,
Eval2dOptions, evaluate_grid, marching_squares, chain_segments, Raster, write_png, write_svg_contours}`.
`[formula: pat2d-circle-fit]` on `ClassicalFit` and `fit_circles`, `pat2d-circle-sdf` on `circle_sdf`/`circle_sdf_gradient`, `pat2d-blend` on the
blend in `Pat2dField` (`field.rs`). Conventions match `Pat2d.lean`: frame `(t, n)` with `t = (-n_y, n_x)`, `n* = (-a1 t + n)/sqrt(1+a1^2)`,
`kappa = 2 a2 / (1+a1^2)^{3/2}`, centre `o + n*/kappa`, sign `+1` iff `kappa < 0`; tested: `g(o) = 0`, `grad g(o) = n*`.

- **Fit.** Weighted least squares of `h` against `s` over the `k = 16` nearest points (self included), Gaussian weights
  `exp(-d^2/(2 b^2))`, `b = 0.6 d_k`; solved in scaled coordinates by Cholesky; a singular system falls back to a line (`a2 = 0`), counted
  in `fit.degenerate_fits`. Defaults by measurement (`pat2d-measure.log`: k 8/16/32/64, weights, n = 200/800 on 8 shapes): `k = 8` is slightly
  best on clean data, 16 is the safer default (k = 64 is clearly worse: 0.004-0.024 mean error on 200 points, because a quadratic cannot follow
  a curve over a long arc); weights help a little.
- **Blend.** As §7.4 in 2D: same `Neighborhood`, `LambdaRule` (PerQuery, Constant, PaperHeuristic with `Pat2dField::paper_options`) and
  `GradientMode` types as 3D, same shift rules, frozen and exact gradients, `DegenerateQuery` for sigma = 0.
- **Accuracy** (`tests/pat2d.rs`, 120x120 grid over the box padded by 25 %, exact SDFs of WP8, 600 uniform boundary points, default options):
  mean |phi - exact| / bounding width: circle 0.0000, ellipse 0.0000, pill 0.00001, hexagon 0.00057, box 0.00107, heart 0.00036; within 0.1 of the
  curve 0.0000-0.0004; max 0.002-0.017 (at creases). With 200 points: 0.0000-0.009 (box 0.0036, hexagon 0.0031). Asserted: mean < 0.3 % of the
  width, band < 0.1 %. Gradient: exact against central differences (`All` + `PerQuery`/`Constant` with `Exact`, `Constant` with `Frozen`) to
  1e-5 relative at 60 random queries each; with the default 32-NN, median agreement 1e-6, p90 1e-3 (set changes between the two sides of a
  difference). Bitwise equal for 1, 3, 8, 16 threads (kappa, grid values, gradients).
- **Far field inside a concave bite.** The crescent ("moon", params 0.5, 0.6, 0.35) has large errors (mean 0.07, max 1.0) with any k: beyond the
  concave arc the nearest circles are the bite circle, `g = r - |x - c|`, which turns negative beyond distance `r` from its centre, while the true
  distance stays positive. Near the curve it is exact (band error 0.004). This is a property of blending local circles (and tori), not a bug;
  `pat2d-moon.png` shows it. WP8's `Shape2d` closed form for the moon agrees with the traced polyline (mean 4e-5), so the reference is right.
- **Picture** (`examples/pat2d.rs`): input points with normals (red), fitted circle arcs (green), SDF colours with level lines every 4 % of the
  width, zero contour (black), for heart, box, moon, pie with 70 points.

### Interface for the learned 2D predictor (WP9b)

`trait QuadraticPredictor: Sync { fn neighborhood_size(&self) -> usize; fn predict(&self, cloud: &Cloud2, nbhds: &Neighborhoods, threads: Option<usize>) -> Result<Predicted2, PatError> }`
with `Predicted2 { coefficients: Vec<Quad2 {a0, a1, a2}>, sigma: Vec<f64>, fallbacks: usize }`. `Pat2d::build(&cloud, &predictor, &opts)` builds
the fitting neighborhoods of `neighborhood_size()` tokens (self first, `crate::Neighborhoods` rows, computed with the 3D kd-tree on the plane
`z = 0`), calls `predict`, and `fit_circles` turns the coefficients (world units, frame `(t, n)`) into circles. A network should apply its own
scale `sigma_i` and output scaling (as `pat-output-scaling` does) and return the `sigma` it used (the line threshold uses it). `Quad2` is the
2D analogue of the six coefficients `(a00, a01, a10, a11, a02, a20)`; its three outputs are `(a00, a10, a20)` in that notation. The test
`a_custom_predictor_plugs_in` shows a predictor of three lines.

## SPEC proposals

1. **§8.7 `kappa_floor` is dimensionless here.** `Fit2dOptions::line_threshold` (default 1e-6): a point is a line iff `|kappa| sigma_i < threshold`
   (`sigma_i` = the predictor's scale; if NaN, 1). An absolute floor depends on the cloud's units (a cloud of size 1000 and one of size 0.001
   differ by 10^6). At 1e-6 the radius is above a million neighborhood radii and the circle-versus-line difference over the neighborhood is
   below 1e-6 of it.
2. **§8.7:** add the interface above and the defaults `k = 16`, bandwidth `0.6 d_k` (SPEC says "weighted least-squares quadratic over the k nearest
   points", unspecified).
3. **§8.4:** state that the final count `<n_i, grad phi(p_i)> < 0` is not a measure of correctness (consistent wrong patches), and prefer the
   residual against ground truth in tests. Keep `eps = 1/sigma_i`; §16.3 result above (flat).
4. **§8.4/§16.3:** the paper's 0.3 % at 30 % is reached on most, not all clouds; thin and high-curvature clouds need more than 10 iterations. An
   option `max_iterations` exists; whether to default to more is a decision for the orchestrator (not measured).

## Open

- Orientation sweep used one seed; the 10-iteration cap leaves femur/snail unconverged (flips still 6/46 at the end).
- No Python/Blender wiring (WP12/13). `apps/mod.rs` does not re-export `orient`; use `apps::orient::...`.
- `Raster`/`write_png` is an uncompressed PNG writer for 2D slices; WP11 may replace it with its own.
- The 2D field uses the kd-tree of the crate on `z = 0`; a 2D tree would be faster by a constant, not needed.
