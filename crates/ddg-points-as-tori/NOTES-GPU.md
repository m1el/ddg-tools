# WP4 notes: the CUDA back end of Points as Tori

Branch `points-as-tori`, build dir `target-pat4`, logs in `runs_tmp/pat-wp4/`. Feature `cuda` (off by default).
RTX 4080 16 GB (`sm_89`), driver 595.99, CUDA 12.4 NVRTC; Threadripper 3970X (32 cores / 64 threads).

## What runs where

| stage | CPU | GPU (`Device::Cuda`) | determinism on the GPU | GPU vs CPU |
|---|---|---|---|---|
| kd-tree build | ours (`spatial::KdTree`) | — (CPU always) | — | — |
| fitting neighborhoods (64-NN, self first) | kd-tree, pool (used by `build`: faster) | `gpu::fitting_neighborhoods` (flattened tree, insertion list) | level 2 | bitwise equal |
| kNN batch | `spatial::knn_batch` | `gpu::knn_batch`, `k ≤ 128` | level 2 | bitwise equal |
| features | pool | — (CPU always, overlapped with the GPU forward pass) | — | — |
| network | burn-flex, `CHUNK = 1` workers | burn-cuda + fusion, `GPU_BATCH = 256` | level 1 (tested), batch-shape independent | ≤ 4e-7 raw (`Split3`), ≤ 2e-4 (`Tf32`) |
| torus fit | pool (used by `build`: as fast) | `gpu::fit_tori` (both sign rules, both eigen rules) | level 2 | bitwise equal, same report |
| blend φ, ∇φ (frozen/exact), Δφ, f64 | pool | `TorusField::evaluate` with `EvalOptions::device` | level 1 | see "Tolerances"; 85–99 % of the values bitwise |
| blend, f32 stable form | — | `EvalPrecision::F32Stable` | level 1 | see "Tolerances" |
| single-point queries (`distance`, `distance_gradient`, `laplacian`, `blend_at`) | always the CPU (latency) | — | — | — |

Reports: `PredictReport::device` (the forward pass), `BuildReport::{knn_device, fit_device}`, `EvalReport::device`,
`TorusField::device()`.

## Device selection

`Device::{Cpu, Cuda { ordinal }, Auto}` as in ddg-repulsive-surfaces:
- `Cuda` is `PatError::Gpu` when the driver, the ordinal or NVRTC is missing, and (for evaluation) when the neighborhood
  does not run on the GPU (`Knn { k }` with `128 < k < n`, `Radius` with `fallback_k > 128`).
- `Auto` uses device 0 when it opens and the neighborhood is supported, else the CPU (`F32Stable` with `Auto` and no
  usable GPU is an error: it has no CPU path). `F32Stable` with `Device::Cpu` is `PatError::InvalidOption`.
- Without the feature: `Cuda` and `F32Stable` are `PatError::Gpu`, `Auto` is the CPU.
- `PointsAsTori::build` takes the precompute device from `BuildOptions::predict.device` (network, and kNN / fit when
  they pay, below) and the evaluation device from `BuildOptions::eval.device`.
- `gpu::probe`, `gpu::devices`, `gpu::resolve` open the driver only; `gpu::GpuContext::new` compiles our kernels
  (NVRTC, ~1 s, memoized per compute capability: the PTX is a pure function of the source and the flags).

## The network on the GPU

**Framework.** The same burn 0.22.1 as the CPU path (`=0.22.1`), on burn's dispatch device `Device::cuda(ordinal)`;
CubeCL compiles its kernels with NVRTC at run time (nothing CUDA at build time; cudarc 0.19.10 with `cuda-12040`, the
same as ddg-repulsive-surfaces and CubeCL). The forward pass is written out in `network/cuda.rs` (the same
architecture as `network/burn.rs`; burn's `LayerNorm`, `softmax`, `gelu_approximate`; attention as in burn's
`MultiHeadAttention`, with q, k and v as one `E → 3E` matmul) because the matmuls need the treatment below.

**Finding 1: CubeCL's f32 matmul is TF32.** On `sm_80+` the default strategy (`SimpleCyclicCmma`) promotes f32
operands to TF32 (`adjust_dtypes`, cubek-matmul 0.3.1; PTX `mma.sync…m16n16k8.f32.tf32.tf32.f32`, dumped with
`CUBECL_CUDA_DUMP_PTX`). Raw outputs then differ from the CPU by up to **1.5e-4** relative (3,000 ellipsoid
neighborhoods), 75× the CPU-vs-JAX difference. There is no configuration switch (no precision knob in burn or CubeCL;
`autotune` would pick by timing, i.e. TF32 again, and breaks run-to-run determinism).

**Fix: exact 3×TF32 splitting** (`GpuMatmul::Split3`, the default). `a = a_h + a_l` with `a_h` exact in TF32;
`a·b ≈ a_h·b_h + (a_h·b_l + a_l·b_h)` (the dropped `a_l·b_l` is ~2⁻²² relative). Weights are split on the host by
truncating the mantissa (`lo` is exact), activations on the device by Veltkamp's split written FMA-proof:
`c = 8192·a + a` (the power-of-two product is exact, so contraction cannot change the one rounding), `a_h = c − (c − a)`,
`a_l = a − a_h`. Raw outputs: **3.8e-7** relative of the CPU. (The textbook `c = 8193·a` breaks under fusion: NVRTC
contracts `8193·a − a` into an FMA and the error returns to 1.6e-4; measured.) `GpuMatmul::Tf32` stays as an option.

**Finding 2: fusion.** Without burn's fusion every element-wise op is its own kernel and the pass is memory bound:
Split3 5.9k/s, TF32 15.8k/s (40,000 neighborhoods, batch 512). With fusion (turned on through `burn-cubecl`'s feature
`fusion`, because `burn/fusion` pulls `burn-remote` into the resolution, which fails on `gloo-timers ^0.4`):
Split3 15.4k/s, TF32 42k/s at batch 256 (machine idle, `exp4-fusion-split.log`). Fusion keeps both properties below.

**Determinism (D4).**
- Run to run: bitwise, in the same process and across processes (hash `7241133e162a8ee1` for TF32 in three
  separate runs; Split3 hashes equal in every run of the same build).
- Batch shape: bitwise identical outputs for batches of 1, 7, 64, 256, 1000 and 4096 neighborhoods (both modes, with
  and without fusion). The matmul heuristic picks the kernel from the shape and no kernel uses atomics or split-K.
  `GPU_BATCH = 256` is still a constant (the last batch is padded with copies of its first neighborhood), so the
  property holds by construction, and 256 is the fastest shape (the activations of 256 × 64 tokens × 512 floats =
  33 MB stay in the 64 MB L2): 128 / 256 / 512 / 1024 gave 30k / 42k / 32k / 24k per s (TF32).
- No candle: burn-cuda is repeatable; the §5.2 fallback is not needed.

**Tolerance against the CPU.** Raw outputs ≤ 3.8e-7 (Split3). Coefficients in the pipeline (5,000-point ellipsoid,
our kNN, both on GPU and CPU): `max |Δa|/(1 + |a|)` = **3.9e-6** (SPEC proposal 1e-5: met); SDF at 4,000 queries:
median 2.3e-9, p99.9 8.9e-9, max 1.0e-8. TF32: raw 1.5e-4 (coefficients ~1e-4).

## Our kernels (`src/gpu/kernels.cu`, cudarc 0.19 + NVRTC)

Workspace GPU rule: one thread per output in the CPU's operation order, gathers only, no atomics, no cross-thread
reductions; compiled with `--fmad=false --prec-div=true --prec-sqrt=true --ftz=false` (one module, no level-1/level-2
split: the kernels are not FMA bound and level 2 needs the flag).

- **kNN** (`knn_batch_{32,64,128}`, `fitting_{32,64,128}`): the CPU's kd-tree flattened (`KdTree::flat`, a
  `pub(crate)` accessor added to WP1's `kdtree.rs`); stack traversal with the CPU's strict pruning; a sorted insertion
  list of `k ≤ KMAX` by `(d², index)`; `d²` as the CPU (`dx² + dy² + dz²`, no FMA). The result is unique, so it is
  bitwise the CPU's for any traversal order (tested with grid clouds full of ties and duplicates, k = 1…128).
- **Fit** (`fit_tori`): `fit_torus` line by line, both `SignRule`s and both `EigenRule`s, clamp flags as bits; tori
  and the whole `FitReport` bitwise the CPU's (4,000 points with random and degenerate coefficients: zeros, umbilics,
  planes, cylinders, `v2 = 0`, saddles, tiny curvatures).
- **Blend** (`eval_knn_{32,64,128}`, `eval_all`, `eval_radius` after `radius_counts`): the CPU's `blend` line by
  line: the shift, the three λ rules, the frozen and exact gradients, the App. D Laplacian, area weights, the
  isovalue, the counters (axis, core, at data point, radius fallbacks, terms) and `DegenerateQuery`. Radius sets: a
  count pass, a host prefix sum, a fill pass with an in-place heap sort by `(d², index)` per query.
- **F32Stable** (`eval_*_f32`): below.

**Level 2 and `exp`.** Everything but `exp` is IEEE-exact in both (`+ − × /`, `sqrt`), so kNN, fitting
neighborhoods, the fit and the per-torus SDF are bitwise equal. The blend's weights use CUDA's `exp`, the CPU uses
glibc's: on 2,000,000 arguments uniform in [−64, 64], **6.2 % differ, by 1 ulp at most** (`exp_against_host_libm`).
So the blend is level 1 only, bitwise equal to the CPU on 85–99 % of the values. Measured (3,000 queries in the
bounding box plus 50 data points, 3,000-point ellipsoid with the network's tori), `max |Δ|/(1 + |a|)`:

| case | φ | ∇φ | Δφ | bitwise |
|---|---|---|---|---|
| 32-NN, PerQuery C = 64, frozen | 1.4e-16 | 3.7e-14 | 2.4e-11 | 87.2 % |
| 32-NN, exact gradient | 1.4e-16 | 3.7e-14 | — | 92.8 % |
| 64-NN, constant λ = 40 | 1.5e-16 | 6.4e-15 | 2.8e-13 | 84.9 % |
| 100-NN, area weights, isovalue | 1.9e-16 | 3.4e-14 | 1.2e-11 | 86.3 % |
| kNN with k ≥ n (index order) | 1.5e-16 | 1.3e-14 | 1.6e-12 | 84.6 % |
| All, exact gradient | 1.5e-16 | 1.2e-14 | — | 91.4 % |
| Radius 0.1, fallback 32 | 1.4e-16 | 3.7e-14 | 2.4e-11 | 87.2 % |
| paper heuristic (radius, λ = 10³/D) | 1.4e-16 | 1.2e-12 | 1.5e-8 | 99.1 % |

The gradient and Laplacian amplify the 1-ulp weight differences by the cancellation of App. D (`λ ~ 10³–10⁵`), which is
the CPU's own conditioning. Counters are identical in every case. Asserted: 1e-14 / 1e-10 / 1e-7.

## F32Stable (SPEC §7.3, D8)

Coordinates relative to the bounding-box center in f32; per torus the expansion point is its data point `p` (not
`o = Q*(0,0)`, see proposal 3), with `w = p − c`, `v_p = ‖w × u‖`, `d_p = (v_p − R, w·u)` and
`K₂ = ‖d_p‖² − r² = (‖d_p‖ − |r|)(‖d_p‖ + |r|)` computed on the host in f64 and rounded. With `y = x − p`:
- `v² − v_p² = |y|² + 2 y·w − (y·u)² − 2 (y·u)(w·u)`, `δ₁ = (v² − v_p²)/(v + v_p)`, `δ₂ = y·u`, `d = d_p + δ`
  (no `O(R)` cancellation);
- `‖d‖² − r² = K₂ + 2 d_p·δ + |δ|²`, `g = s (‖d‖² − r²)/(‖d‖ + |r|)` (no `O(r)` cancellation; SPEC §7.3 only treats
  the first);
- the blend in the cancellation-free arrangement: `∇φ = Σ w ∇g/W + λ Σ w r̂ (φ − g)/W`,
  `Δφ = Σ w [Δg − 2λ r̂·(∇g − ∇φ) + (g − φ)(λ² − 2λ/ρ)]/W` (three passes, torus terms recomputed). Equal to App. D in
  exact arithmetic; in f32 the App. D form `A + φB` loses `ε·λ|g| ≈ 1e-3`.
- kNN stays f64 (the CPU's neighbor sets).

The rearrangement is exact: written in f64 (`stable_form_in_f64_equals_cpu`) it matches the CPU's Eq. 23 blend to
3e-16 (φ), 6e-14 (∇φ), 3e-11 (Δφ). F32Stable against the CPU (f64), absolute errors (`blend_f32_stable_tolerance`):

| cloud | neighborhood | φ | ∇φ | Δφ / (1 + |Δφ|) |
|---|---|---|---|---|
| ellipsoid 0.6 × 0.45 × 0.3 | 32-NN | 3.8e-7 | 7.9e-5 | 2.3e-2 |
| same | All | 9.4e-7 | 6.6e-5 | 2.2e-3 |
| plane, curvatures ≤ 2e-4 and 0 (R, r up to 1e10) | 32-NN | 9.5e-7 | 2.2e-4 | 0.33 |
| same | All | 1.4e-6 | 1.0e-4 | 8.8e-3 |
| ellipsoid moved to (1000, −500, 250) | 32-NN | 3.8e-7 | 1.1e-4 | 4.7e-2 |

Without the stable form an f32 torus with `R = 1e10` would carry an error of `ε·R ≈ 600`. The Laplacian in f32 is
poor (the `λ²(g − φ)` terms, `λ² ~ 10⁷`, amplify the f32 position quantization); use f64 for Laplacians.
At queries within f32 resolution of a data point the f32 path sees `ρ = 0` (the `r̂ := 0` rule) where f64 sees
`ρ ~ 1e-16`: there `∇φ` and `Δφ` are singular anyway (the weights have a cone point at every data point).

## Measurements

Logs and `nvidia-smi dmon -s um` traces in `runs_tmp/pat-wp4/` (`timing-*.log`, `dmon-*.log`; every row records
`uptime` and `nvidia-smi` before and after). Clouds: Keenan Crane's spot (`data/models/spot`), sampled
area-uniformly with interpolated vertex normals; reference weights; 64-NN; default options. "Loaded" = WP7's CPU
training run (`train-… short_run`, 30–40 cores) was running (1-minute load 40–140 on 64 threads); "quiet" = load
3–25 (WP7 finished at ~03:28). Only the GPU was ours throughout (no other compute process in `nvidia-smi`).

### GPU utilization: why it was spiky, and what was done

**Profile.** No `nsys`/`ncu` on the machine (only `nvprof`, which does not support Ada), so: per-stage wall times with
device syncs (`BuildReport`, `PAT_GPU_PROFILE=1` for evaluation groups), `top -H` of the process, and `dmon` beside
every run.
- The host stages are not the gap: kNN (0.2–0.4 s per 1M), features (1–4 s of CPU per 1M, computed on the pool one
  group of 2,048 neighborhoods ahead of the GPU through a bounded channel) and the fit (0.1 s) are small and
  overlapped.
- The forward pass is **burn's per-op host dispatch**: three threads of the process (the caller and CubeCL's
  `DSU`/`DSD` stream threads) run at 80–90 % CPU each, which is the ~2.8 cores the orchestrator saw. Every batch
  records and dispatches a few hundred ops (fused by burn-fusion), then reads back 6 floats per neighborhood.
- With WP7 holding 36 cores, those three threads were starved and the GPU waited: SM 8–90 %, mean **40 %**
  (`dmon-dispatch.log`, Split3), and the throughput of the *same* call fluctuated 5k–42k/s between repetitions.
  With a quiet host the same code keeps the device fed: SM **96 %** mean over the 1M Split3 build, **91 %** over the
  1M TF32 build (`dmon-precompute-quiet.log`), **95 %** / **88 %** over the 4M builds (`dmon-4m-2.log`, load 19–26); the forward pass is then GPU bound (TF32 tensor-core matmuls and the
  fused element-wise kernels at ~9 TFLOP/s).

**What was tried** (`timing_dispatch`, now removed; `dispatch*.log`, `ab*.log`):

| change | effect |
|---|---|
| several batches queued before the read-back (depth 2, 4, 8) | no gain (TF32 39–43k/s at depth 1–8, quiet-ish); **with the first Split3 layout depth 2 changed the output bits** (burn-fusion segments the queued op stream differently, and fused kernels contract differently) → read back after every batch, which also pins the fusion segments |
| 2 or 4 host threads dispatching batches concurrently (outputs bitwise unchanged) | slower: TF32 30k → 8k/s; burn's dispatch serializes (SM 10–40 %) |
| Split3 as one K-concatenated matmul `[a_h, a_h, a_l]·[W_h; W_l; W_h]` (fewer ops) | Split3 steady at 10k/s under load but slower quiet (10k vs 15k/s): the concatenations cost more GPU time than the dispatch saved; reverted |
| q, k, v as one `E → 3E` matmul (kept) | ~2 % for Split3, TF32 unchanged; fewer ops per batch |
| batch shape 256 vs 512–2048 | 256 is best when the host is quiet (L2-resident activations); under load 1024 was better for TF32 (21k vs 14k/s) — not adopted, the quiet case is the one that matters |
| pinned staging buffers | not needed: 393 KB per batch up, 6 KB down |
| keep kNN → features → network → fit on the device | the host stages are ~2 % of the time and overlapped; the bound is inside burn's forward pass, which only a fused inference kernel of our own would remove (proposal 5) |

All variants kept the D4 properties: bitwise run to run, independent of the batch shape; the depth experiment is the
one that showed that overlap *can* change bits through burn-fusion, so it is not used.

### Network throughput (quiet, 50,000 neighborhoods, one network kept on the device, after warm-up)

| mode | neighborhoods/s | SM (dmon) |
|---|---|---|
| TF32 | 38,000–42,800 | 90–97 % |
| Split3 (default) | 14,700–15,100 | 92–97 % |

Without fusion (load 1–2): TF32 15.8k/s, Split3 5.9k/s (batch 512).

### Precompute

The largest cloud run is 4M points (1.8 GB of device memory: the features stream in groups, so device memory does not
grow with `n`; 7M would fit as well).

| n | device | kNN (CPU) | features | network | fit (CPU) | total | per 1M | load |
|---|---|---|---|---|---|---|---|---|
| 100k | GPU Split3 | 0.03 | 0.07, overlapped | 6.79 | 0.01 | **6.85** | 68.5 | quiet (7) |
| 100k | GPU TF32 | 0.03 | 0.06, overlapped | 2.51 | 0.01 | **2.58** | 25.8 | quiet (7) |
| 100k | CPU 64 threads | 0.02 | 0.77 (summed) | 21.9 wall | 0.01 | **21.9** | 219 | quiet → 23 |
| 1M | GPU Split3 | 0.40 | 1.06, overlapped | 67.0 | 0.10 | **67.7** | 67.7 | quiet (19–23) |
| 1M | GPU TF32 | 0.38 | 1.13, overlapped | 25.9 | 0.10 | **26.6** | 26.6 | quiet (19–24) |
| 1M | CPU 64 threads | 0.36 | 17.1 (summed) | 283.1 wall | 0.10 | **283.7** | 284 | 24 → 106 |
| 4M | GPU Split3 | 1.97 | 4.45, overlapped | 270.5 | 0.36 | **273.9** | 68.5 | 19–26 |
| 4M | GPU TF32 | 1.78 | 4.68, overlapped | 103.6 | 0.41 | **106.8** | 26.7 | 19–20 |
| 1M, loaded | GPU Split3 / TF32 / CPU | | | 156 / 65 / 410 | | 158 / 67 / 411 | | 50–140 |
| 4M, loaded | GPU Split3 / TF32 | 2.9 / 3.2 | 25 / 22 | 554 / 263 | 0.8 / 0.9 | 560 / 270 | 140 / 67 | 36–57 |

Stage comparisons at 1M (quiet): kd-tree 0.12 s (CPU); 64-NN fitting neighborhoods 0.20 s on 64 CPU threads vs
0.62 s on the GPU including transfers; fit 0.125 s CPU vs 0.133 s GPU including transfers. Hence
`gpu::USE_GPU_KNN = USE_GPU_FIT = false`: `PointsAsTori::build` keeps kNN and fit on the CPU
(`BuildReport::{knn_device, fit_device}` say so); the GPU functions stay public and are tested bitwise. The first
forward pass of a process JIT-compiles CubeCL's fused kernels (0.1–2.8 s, the warm-up rows).

**Against the floor** (§16.5: the paper's 40 s per 1M on a 3090; target ≤ 20 s on the 4080): **TF32 26.6 s per 1M
meets the paper's floor, Split3 (f32-accurate, default) 67.7 s does not**; neither meets 20 s. The CPU path takes
220–280 s per 1M. Proposal 5.

### Evaluation (spot 100k points, its tori, 32-NN, PerQuery C = 64; seconds; load 29–35)

| grid | device | φ | φ + ∇φ | φ + ∇φ + Δφ |
|---|---|---|---|---|
| 64³ (262k) | CPU f64, 64 threads | 0.150 | 0.144 | 0.147 |
| 64³ | GPU f64 | 0.061 | 0.068 | 0.063 |
| 64³ | GPU F32Stable | 0.059 | 0.063 | 0.058 |
| 256³ (16.8M) | CPU f64, 64 threads | 8.10 | 8.77 | — |
| 256³ | GPU f64 | 2.49 | 2.93 | — |
| 256³ | GPU F32Stable | 2.40 | 2.67 | — |

The device is fed during evaluation (SM 85–100 % in `dmon-eval.log`): per group of 2²⁰ queries 0.13–0.18 s are
upload and kernel, ~150 ns per query, about 10× the f64 arithmetic bound. The per-thread kNN (stack and insertion
list in local memory, divergent traversal) dominates, which is also why F32Stable (f32 blend, f64 kNN) is barely
faster than f64. Open item (a warp-cooperative or leaf-batched kNN).

**Per-query latency** (100k points, load 29): CPU sequential 32-NN 9.8e-6 s (φ) and 1.3e-5 s (φ + ∇φ) per query
(paper 1e-4–1e-3 s; target ≤ 5e-5 s: met). A one-query GPU batch takes 2.9e-4 s (launch, transfers), so single-point
queries stay on the CPU.

## SPEC proposals (for the orchestrator)

1. **§5.2/§5.4, D4.** burn-cuda's f32 matmul is TF32 on `sm_80+`; the default GPU inference uses the exact 3×TF32
   split (`PredictOptions::gpu_matmul = GpuMatmul::Split3`), `GpuMatmul::Tf32` is an opt-in. Burn's fusion is on for
   the `cuda` feature (through `burn-cubecl`). `GPU_BATCH = 256` (a constant). Measured tolerance against the CPU:
   raw outputs 4e-7, coefficients 4e-6 (`1e-5·(1 + |a|)` met). burn-cuda is run-to-run bitwise and independent of the
   batch shape; candle is not needed.
2. **Level 2 for the blend needs a crate-owned `exp`.** The CPU's `f64::exp` is the platform libm (glibc here; other
   platforms differ, so CPU results are already platform dependent at the ulp level). Proposal: one `exp` written with
   `+ − × /` only (Cody–Waite reduction and a fixed polynomial, ~0.5–1 ulp), used by `blend.rs` and `kernels.cu`. Then
   the blend is bitwise CPU = GPU and bitwise across platforms; the cost is ulp-level changes of the CPU results
   (far inside every parity tolerance). Until then: level 1 + the table above.
3. **§7.3 F32Stable.** Expand around the data point `p_i` (exact in f64 on the host) rather than `o = Q*(0,0)`: it
   needs no fit diagnostics (the field only has the tori), and with `K₂` and `d_p` from f64 it also removes the
   second cancellation `‖d‖ − |r|` (flat in both directions, `r ~ 1e10`). The Lean theorem `torus_sdf_stable` assumes
   `(o − c)·u = 0`; for `p_i` it needs the extra term `−2 (y·u)(w·u)` (drop the hypothesis) and a second identity
   `‖d‖² − r² = K₂ + 2 d_p·δ + |δ|²`. The cancellation-free blend arrangement (above) could get its own id
   (`pat-blend-gradient-stable`): it equals `pat-blend-gradient` / `pat-blend-laplacian` algebraically.
4. **§7.5.** GPU kNN for `k ≤ 128` (and `Radius` fallbacks `≤ 128`); larger `k < n` stays on the CPU (`Auto`) or is an
   error (`Cuda`). Single-point queries stay on the CPU in f64 whatever the device (latency, below).
5. **§16.5 floor; needs a user decision.** burn-cuda with fusion reaches the paper's floor only with TF32
   (26.6 s per 1M quiet, raw outputs 1.5e-4 off the CPU); f32-accurate Split3 takes 67.7 s per 1M quiet, with the GPU
   at 96 % SM (GPU bound on a quiet host; host-dispatch bound under outside CPU load). Options:
   (a) keep Split3 as the default and accept 1.6× the paper's time (current state); (b) make TF32 the default and
   document the 2e-4 tolerance (a saddle flips more often, cf. NOTES-WP1 item 1); (c) a fused f32 inference kernel of
   our own on NVRTC (one kernel per layer, deterministic, f32 SIMT or 3×TF32 `mma`), which would remove the host-side
   dispatch bound (estimate 50–100k/s) but departs from "a framework for the network" (D2) for GPU inference only,
   as §18.3 allowed for the CPU; (d) candle + cuBLAS (the §5.2 fallback; nvcc at build time; cuBLAS SGEMM is f32 by
   default), which was not tried because burn-cuda is repeatable. No candle dependency was added.

## Changes to files of other packages (minimal)

- `src/spatial/kdtree.rs` (WP1): `KdTree::flat()`, `pub(crate)`, `cfg(feature = "cuda")`.
- `src/network/mod.rs` (WP1): `GpuMatmul` and `PredictOptions::gpu_matmul`; `forward_features` and
  `predict_coefficients` dispatch to the GPU after their checks; `check_device` without the WP4 stub.
- `src/field.rs` (WP1): `EvalReport::device`, `BuildReport::{knn_device, fit_device}`, `TorusField::device()`,
  `TorusField::new` opens the GPU field (`open_gpu`), `evaluate` dispatches to it, `PointsAsTori::build` uses the GPU
  for kNN/fit when the precompute device is a GPU and `gpu::USE_GPU_{KNN,FIT}` say so.
- `Cargo.toml`: feature `cuda = ["burn/cuda", "dep:cudarc", "dep:burn-cubecl"]`, the optional `cudarc` and
  `burn-cubecl` (fusion) dependencies. The first WP4 commit also carried WP6's pending `ddg-core` line (their
  committed `mc.rs` needed it).
- `Cargo.lock`: once restored by `cargo metadata` after a mistaken `git checkout Cargo.lock` (WP4) dropped other
  packages' uncommitted lock entries; the regenerated entries are the same versions.

## How to resume

```
export CARGO_TARGET_DIR=/var/data/discrete-geometry/ddg-tools/target-pat4
cargo test --release -p ddg-points-as-tori --features cuda --test gpu -- --nocapture --test-threads 1
# timings (alone on the GPU; PAT_SIZES, PAT_CPU=0, PAT_EVAL_N)
cargo test --release -p ddg-points-as-tori --features cuda --test gpu -- --ignored timing_ --nocapture --test-threads 1
cargo test --release -p ddg-points-as-tori --features cuda --test gpu -- --ignored network_throughput --nocapture
```
`DDG_REQUIRE_GPU=1` turns a missing GPU into a failure; otherwise every GPU test prints `SKIP` and passes.
