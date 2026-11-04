# WP11: the sphere-tracing renderer (SPEC §8.5, Fig. 12)

## What was built

| item | where |
|---|---|
| Tracer, camera, options, report, PNG | `src/render/mod.rs` |
| Shading (Lambert + IBL-GGX + Fresnel rim, colour map, isobands) | `src/render/shade.rs` |
| Cube map (demo's `ForestBlurred0..5.png`) and constant sky | `src/render/cubemap.rs`, `assets/` |
| CUDA backend | `src/render/gpu.rs`, kernels `trace_knn_{32,64,128}`, `trace_all` in `src/gpu/kernels.cu`, `GpuField::trace` in `src/gpu/mod.rs`, `TorusField::gpu_trace` in `src/field.rs` |
| Examples | `examples/render.rs`, `examples/render_study.rs` |
| Tests | `tests/render.rs` (CPU, 5 tests + 4 unit tests), `tests/render_gpu.rs` (CUDA) |
| Python | `pat_render` / `points_as_tori.render` in `crates/ddg-py` (`src/points_as_tori.rs`, `src/python/points_as_tori.rs`, `ddg_tools.pyi`, class `Rendering` in `python-tests/test_points_as_tori.py`) |
| Showcase | `showcase/points-as-tori/render/` (11 images at 1024², `render.sh`, `safety-study.txt`) |
| Formulas | `pat-render-step` (proved: `sphere_trace_step_safe` in `Apps.lean`), `pat-render-ggx`, `pat-render-fresnel` (cited, `CITED.md`); `check-formulas.sh`: 0 missing |

The assets are the demo's six 64x64 RGB PNGs, copied unchanged (15.9 KB in all, so no downsampling
question arises) with `assets/LICENSE-points-as-tori-MIT.txt` and `assets/README.md`.

### Design

`Tracer::render(&dyn ImplicitField, &Camera)` works on any field (`TorusField`, `Csg` with offsets and
Booleans, analytic fields). It runs in phases over a `Backend` (CPU: rayon pool of
`options.threads`; CUDA: `render::gpu::GpuBackend`): (0) per-pixel rays and clip ranges, (A) primary
rays, (B) shadow rays of the ground pixels nothing else covers, (D) `phi` at visible cut-plane
points, (C) shading. Every pixel is independent and computed in f64 in a fixed order.

Traversal (`trace_ray`): clip to the bounding sphere (centre of the field's bbox, radius
`0.8*diag + max(level, 0)`), then `t += safety*|phi - level| / max(|grad|, g_min)`, hit at
`|phi - level| < eps`, at most `max_steps` evaluations. Defaults: safety 0.9, `g_min` 0.1,
128 steps, `eps = 1e-3*diag`; `TraceOptions::demo()` is the demo's (1.0, 0, 40, `1e-2 s`, shadow 80
steps and `1e-4 s`, ray started at the camera), with `s = diag/(2 sqrt 3)` so a `[-1,1]^3` scene is
the demo's unit scene.

Shading is the demo's: diffuse `max(n.l, 0.3)` with light `(0,10,2) s`; `0.3 * normalShade`
(IBL-GGX: metallic 0.07, roughness 0.1, `F0` 0.2, GGX roughness 0.07, specular direction
`(-0.5,1,0)`, `pow(., 2.2)` and `pow(., 1/2.2)`); `0.6 * (1 - |w.n|)^80`; then `min(., 1)` and `sqrt`.
Colour `(68,74,145)/255`, white background, ground plane shadow-only (a shadowed ground pixel is
0.8). Cut plane: the demo's `getCutPlaneBasis`, colour map over `[-2s, 2s]`, isobands
`0.9(1 - 0.1 smoothstep(-0.01, 0.01, sin(-60 phi/s)))`.

## Determinism

`tests/render.rs::determinism_across_threads`: a fitted field with an offset and a Boolean, cut plane
and shadows; the image bytes and nine report counters are identical at 1 (twice), 3 and 8 threads.
Python: `Rendering.test_bitwise_at_any_thread_count` (1, 3, 4 threads). The CUDA kernel repeats
`trace_ray` operation for operation on the f64 blend (no FMA through `--fmad=false`, the crate's
`pat_exp`); shading stays on the host, so only hit positions and gradients can differ from the CPU,
by the blend's measured `exp` tolerance.

## Measurements (CPU, 8 threads, machine at load average ~12 from other jobs)

| image (1024²) | render | steps mean / max | hits |
|---|---|---|---|
| tooth, level 0 | 2.8 s | 4.3 / 77 | 73,740 |
| tooth, offset 0.16 | 3.4 s | 4.9 / 93 | 118,796 |
| snail, level 0 | 3.5 s | 4.6 / 78 | 66,245 |
| bunny, level 0 (12,917 points) | 4.6 s | 4.9 / 96 | 83,109 |
| analytic 5 / subtract / union (CSG of two fields) | 3.7 / 8.7 / 7.9 s | 5.6-6.6 / 123-128 | 140,201 / 102,649 / 150,779 |

Build of the fields (2048 points) 2.6 s; 12,917 points 13 s. Most pixels miss: the mean steps are
over traversed rays only (`clipped` rays are not traversed). Shadow rays cost about as much as
primary rays. The paper's 40-step shader runs in real time on a laptop GPU; ours on the CPU takes
3 to 9 s per megapixel. The CUDA path was **not run** (see below).

**Safety study** (`showcase/points-as-tori/render/safety-study.txt`, 256², vs a reference render with
safety 0.5 and 2048 steps), three clouds: safety 0.9 / 128 steps differs from the reference in 2 to 12
pixels of 3500-4750 hits (silhouette pixels at thin parts), mean steps 4.3-4.6 against 7.9-8.6 for
safety 0.5; safety 1.0 is as good on tooth and bunny (0 and 7 pixels), 7 on the snail. `g_min` in
{0, 0.1, 0.3} changes nothing here (gradient norms stay near 1). 64 steps gives 1-6 step-limit rays;
the demo's 40 steps with safety 1 gives 20-56 step-limit rays per 65k pixels (silhouettes are cut).
The default 0.9 / 128 stays; the evidence does not distinguish 0.9 from 1.0.

## Deviations from the demo

1. Rays start at the bounding-sphere entry (the demo tests the sphere and then restarts at `t = 0`).
   `start_at_sphere = false` is the demo's behaviour (`TraceOptions::demo`).
2. Step `safety` 0.9 and the `g_min` floor (SPEC D17); 128 steps; relative `eps`.
3. The ground plane's normal is `up` (default `+y`) and its height is below the bbox along `up`; the
   demo puts it at `bbox_min.z - 1` with normal `+y` (an axis mix-up, NOTES-REFERENCE §14).
4. The demo's cut-plane map is black outside `[-range, range)` (its `segment` is 0 there, and
   `pow` of a negative base is undefined in GLSL); `CutPlane::clamp_colormap` (default true) saturates
   instead. `false` reproduces the black.
5. `pow` bases are clamped at 0 where `refl -= spec` can go negative (undefined in GLSL).
6. `textureBlurred`'s six directional averages are computed once per environment (they do not depend
   on the argument).
7. Shadow rays are traced only for ground pixels nothing else covers; ground pixels the demo
   colours by the object are the same.
8. Cube-map faces use the OpenGL face table with bilinear, clamped lookups, rows in PNG order.
   Seams are not blended (GL's non-seamless default).

## Decisions for the user

- **Defaults** `safety = 0.9`, `g_min = 0.1`, `max_steps = 128`, `eps_hit = 1e-3 * diag` (SPEC D17 values;
  `g_min` is ours, SPEC leaves it open). The study above does not favour a change.
- **Up axis**: the ground plane and the camera default to `+y` as the demo; clouds in Blender are
  z-up, so the Python/Blender caller passes `ground_up=(0,0,1)` and `up`. We did not guess a global default.
- **Python return value**: `pat_render` returns a flat `array('B')` of RGBA bytes (alpha 0 on the plain
  background) rather than an ndarray, per the repo's "numpy never required" convention (SPEC §12.1 says
  `-> ndarray`); `numpy.asarray(a).reshape(h, w, 4)` gives the array.
- **GPU scope**: the CUDA tracer handles a single `TorusField` on a CUDA device (kNN with `k <= 128`, or
  all points; level and field isovalue offsets). Radius neighborhoods and Booleans of several fields
  render on the CPU. A CSG program kernel (several leaves) is possible if wanted.

## GPU tests written, not run (the GPU is busy with the training run)

- `tests/render_gpu.rs::gpu_trace_matches_cpu_and_repeats` (`#[ignore = "run after training (WP16)"]`):
  CPU vs GPU image (at most 0.5% differing bytes, hits within 2), and GPU run to run bitwise.
- `tests/render_gpu.rs::gpu_tracer_rejects_cpu_field` runs without a device (an error before any device is opened).

The kernels were compiled with `nvcc -ptx -arch=sm_80 --fmad=false` on the host only (no GPU use), and the
Rust side builds and passes clippy with `--features cuda`. NVRTC compilation of the new kernels at run time
has not been exercised.

## Open issues

- Evaluation errors along a ray (`DegenerateQuery`) abort the render with the pixel index (`EVAL_ERROR`);
  never seen in practice.
- Cube-map seams are visible only in the rim/IBL terms, which are blurred; not blended.
- `examples/render.rs` reads PLY (6 floats), OBJ with `vn`, or an OBJ mesh (uniform sampling); no
  other formats.
- Animated sweeps over the offset level (Fig. 12's series) are a shell loop over `level=`.

## SPEC proposals

1. §8.5: state `g_min = 0.1` default and the bounding sphere `0.8*diag + max(level, 0)`; say rays start at the sphere entry.
2. §12.1: `render(...) -> array('B')` RGBA (or ndarray where numpy is present), options as in `ddg_tools.pyi`.
3. §17 WP11: mention `TraceOptions::demo()` and the CUDA scope (single field).
4. §16.3 (tracer safety): record the study above as the measurement; keep 0.9.
