# Points as Tori — contract (draft, 2026-11-04)

This file is the contract between:
- the Rust crate `crates/ddg-points-as-tori`,
- the Lean oracle (`lean/oracle`) and the Lean proofs (`lean/proofs`),
- the Python binding (`crates/ddg-py`) and the Blender add-on.

It implements:

> Nicole Feng, Ioannis Gkioulekas, Keenan Crane. *Points as Tori: Fast Pointwise Signed Distance
> for Point Clouds.* ACM Transactions on Graphics 45(4), Article 53, 2026.
> <https://doi.org/10.1145/3811385>
> (`../PointsAsTori/PointsAsTori.pdf`, "the paper"; equation, section and figure numbers are the
> paper's).

**Reference.** `../PointsAsTori/points-as-tori` at `253bf9f` (MIT, "Copyright (c) 2026 Nicole
Feng"; cited as `ref: <file>:<line>`).
- It may be read and ported. Ported code and the converted pretrained weights keep the notice in
  `LICENSES/points-as-tori-MIT.txt`. The crate credits all three authors.
- Its submodules (nanobind, libigl, fcpw, nanoflann) are empty here and are not ported.
- libigl's `copyleft/cgal/point_areas` and CGAL are GPL. Barill point areas are implemented from
  Barill et al. 2018, Sec. 5, not from that code.
- `NOTES-REFERENCE.md` is the survey: it gives the exact network, the clamps, the training
  pipeline, the data, and every paper/code discrepancy with line numbers. This file cites it as
  `N§k`.

**Conventions** (`CLAUDE.md`) apply:
- `f64` by default, except where §5 says otherwise. Stored indices are `u32`.
- Bad input gives an error enum naming the element, never a panic.
- Nothing is normalized, clamped or repaired silently: every clamp is an option and is counted in
  a report.
- The thread count is an option: an explicit pool, never the global rayon pool. Results are
  bitwise independent of the thread count (§5.4).
- Formula tags `[formula: pat-…]` (§13).

If an implementation or a test contradicts this file, stop and report it; do not quietly follow
either side.

User decisions of 2026-11-04 are in §18; they override any **[user]** item or recommendation
below that says otherwise (burn 0.22 on Rust 1.95.0, D4 accepted, own CPU kernel if burn is slow,
both 2D variants).

## 1. The problem

**Input.** An oriented point cloud `P = {(p_i, n_i)}`, with `n_i` unit normals.

**Output.** An analytic approximation `φ: R³ → R` of the signed distance to the surface the cloud
samples: negative inside, positive outside. It can be evaluated pointwise, with its gradient, at
any `x`.

**Construction.**
1. **Neighborhoods.** For each `i`, take the `k = 64` nearest points of `p_i`, `p_i` included and
   first (§3).
2. **Coefficients.** A fixed pretrained transformer (§6) maps the neighborhood, expressed in a
   local frame and scaled by the median neighbor distance `σ_i`, to six coefficients
   `a = (a00, a01, a10, a11, a02, a20)` of a height function
   `Q_i(s,t) = a00 + a10 s + a01 t + a11 st + a20 s² + a02 t²`
   over the tangent plane (Eq. 18).
3. **Tori.** A torus `T_i` is fitted analytically. Its equator passes through `p_i + a00 n_i`, and
   its principal curvatures and directions are those of `Q_i` (Eq. 19–22). Its SDF `g_i` is closed
   form (Eq. 23–24).
4. **Blend.** `φ(x) = Σ w_i(x) g_i(x) / Σ w_i(x)` with `w_i = exp(−λ_x(‖x − p_i‖ − σ_x))`
   (Eq. 25–26). The sum runs over a neighbor set of `x` (§7).

**Applications** (user, all):
- offsets, Booleans and morphology (Fig. 22–24),
- normal orientation (Fig. 27),
- sphere-traced rendering (Fig. 12),
- the 2D circles variant (Fig. 10; **not in the reference**, D15),
- baselines and the evaluation of Sec. 5.2 (Fig. 7, 17–20, 29),
- training the network in Rust from the authors' data (Sec. 4.3, 5.1),
- Python and Blender.

## 2. Data model

### 2.1 Point cloud

```rust
pub struct PointCloud { positions: Vec<[f64; 3]>, normals: Vec<[f64; 3]> }   // validated
impl PointCloud {
    pub fn new(positions, normals, policy: NormalPolicy) -> Result<(Self, CloudReport), PatError>;
}
pub enum NormalPolicy { RequireUnit { tol: f64 /* 1e-6 */ }, Normalize }
pub struct CloudReport { pub n: usize, pub normalized: usize, pub exact_duplicates: usize, pub bbox: [[f64; 3]; 2] }
```

**Errors** (`PatError`):
- `LengthMismatch`, `TooFewPoints { n, need }`, with `need = 2`.
- `NonFinite { index, what }`, `ZeroNormal { index }`.
- `NonUnitNormal { index, norm }` under `RequireUnit`.

**Normals.**
- The reference fits with the *raw* normal and builds features with the normalized one (N§17.4).
  Under our contract they are the same vector.
- `Normalize` is an explicit opt-in and is counted in the report. PLY input from photogrammetry
  usually needs it.

**Duplicates** are allowed. They are counted, and their order is fixed by index (§3.1).

`n < 64` is allowed. Neighborhoods then have `n` tokens (as in the reference, N§4), and the
`FitReport` says so.

### 2.2 Tori

```rust
pub struct Torus { pub center: [f64; 3], pub axis: [f64; 3], pub major: f64, pub minor: f64 }  // sign(T) = sign(minor)
pub struct Tori { pub tori: Vec<Torus>, pub coefficients: Vec<[f32; 6]>, pub sigma: Vec<f64> }
pub struct FitDiagnostics { pub kappa_min: f64, pub kappa_max: f64, pub v_min: [f64; 3], pub v_max: [f64; 3],
                            pub n_star: [f64; 3], pub clamps: ClampFlags }
```

- The axis is a unit vector. The reference stores spherical angles (D7). The angles are produced
  only for reference file export.
- `Tori` (de)serializes to:
  - our binary format: a little-endian header, f64 arrays, and the source cloud's hash;
  - the reference's pickle tuple layout, via the Python binding, for black-box tests only.

## 3. Neighbor search (D10)

### 3.1 kd-tree, our own (`spatial::KdTree`)

- Points are f64. The tree splits at the median of the widest axis and has buckets of 16
  (measured choice: 8–32).
- Construction is parallel over subtrees, deterministic (median by `select_nth_unstable_by` on
  `(coord, index)`).
- **Exact kNN**, with results ordered by `(distance², index)`, so ties are deterministic.
- Radius search sorted the same way.
- Batch queries on the explicit pool. Typical speed: 64-NN of 1M points in about 1 s on 32 cores
  (to be measured).

**Why our own.**
- The workspace has none for points: `ddg-repulsive-curves::Bvh` is an element BVH without kNN.
- We need exact results, a deterministic tie order (the network's first token must be `p_i`, and
  fixtures must be reproducible), f64, no global pool, and the same structure on the GPU (§7.5).
- `kiddo` (fast) has its own float and tie conventions and a fixed-size bucket API that would have
  to be wrapped anyway. The tree is about 400 lines.

### 3.2 Fitting neighborhoods

- `N_k(i) = [i] ++ (the k − 1 nearest j ≠ i by (‖p_j − p_i‖², j))`.
- The reference takes nanoflann's first k, which is `i` first unless an exact duplicate ties.
  Ours makes "self first" explicit (D11).

### 3.3 Evaluation neighborhoods

These are the same kNN or radius queries for arbitrary `x` (§7.2).

## 4. Features and output scaling (Sec. 4.3; N§2.2)

### 4.1 Features

For `N_k(i) = (j_0 = i, j_1, …, j_{k−1})`:

- **Unit normal.** `n̂_j = n_j / max(‖n_j‖, 1e-10)`. With validated input `n̂ = n`; the guard stays
  for the reference-order mode.
- **Scale.** `σ_i = median_m ‖p_{j_m} − p_i‖`, over all `k` values including `m = 0`. With `k` even,
  this is the mean of the `(k/2)`-th and `(k/2+1)`-th smallest. `[formula: pat-median-scale]`
- **Frame.** `ŝ_i, t̂_i` from `c = n̂_i`: `u = e_y` if `|c_x| > 0.9` else `e_x`;
  `ŝ = (c × u)/‖c × u‖`; `t̂ = c × ŝ`. `[formula: pat-local-frame]`
- **Feature vector.** `f_m = (d·ŝ, d·t̂, d·c, n̂_{j_m}·ŝ, n̂_{j_m}·t̂, n̂_{j_m}·c) ∈ R⁶`, with
  `d = (p_{j_m} − p_i)/σ_i`. `[formula: pat-features]`

**Precision (D3).**
- Default: compute `d` difference-first in f64, then round to f32.
- `FeatureOrder::Reference`: f32 throughout, with the reference's `p/σ − p_i/σ` on absolute
  coordinates (N§18). Used for parity tests against the JAX oracle.

**Invariance.**
- Features are invariant to translation and to uniform scaling of the cloud.
- They are equivariant to rotations only through the frame rule. The frame is not rotation
  invariant (the `|c_x| > 0.9` switch), and the network was trained with this fixed frame. No
  invariance claim beyond translation and scale.

`[formula: pat-feature-invariance]` (translation and scale).

### 4.2 Output scaling

The network output `y ∈ R⁶` maps to coefficients
`a = (σ y0, y1, y2, y3/σ, y4/σ, y5/σ)`.

`[formula: pat-output-scaling]`: with network outputs fixed, scaling the cloud by `c > 0` scales
the torus centers by `c` (about the origin of the cloud's scaling), radii by `c`, and the SDF by
`c`.

## 5. Precision, determinism, back ends (decisions D1–D5)

### 5.1 Numeric types

| stage | type | why |
|---|---|---|
| kd-tree, neighborhoods | f64 | exact neighbor sets |
| features | f64 → f32 (D3) | the network is f32 |
| network forward (inference) | **f32** | trained weights are f32; JAX trained in f32 (N§9.5); f64 forward of f32 weights changes outputs by ~1e-6 relative, below the training noise, and costs 2–64× |
| network f64 path | f64 (framework backend over f64, CPU only) | oracle comparisons and FD checks of the training gradient (§10.4) |
| coefficients → tori | f64 | the f32 coefficients are promoted exactly; fitting is cheap |
| torus SDF, blend, gradient (CPU) | **f64** | for flat regions `R, r ~ 1/|κ|` reach 10³–10⁸, and `‖d‖ − |r|` cancels: f32 error `ε·R` is up to 1e-4–10 (N§18) |
| GPU evaluation | f64 by default; f32 opt-in with the stable form (D8) | RTX 4080 f64 is 1/64 rate (~0.76 TFLOPS); 32-NN evaluation at 1080p × 40 steps is ~0.2 s/frame in f64, fine for offline images |
| training | f32 (framework), loss in f32 like the reference | matches the published weights' regime |

**Proposal (D1).** f32 network and f64 geometry.

Tradeoffs:
- f32 everywhere (the reference's demo) breaks in flat regions.
- f64 everywhere makes the network 2× slower on the CPU and 64× slower on this GPU, and does not
  match the trained weights better.
- The constant `C = 64` (Eq. 26) is kept as the default in both precisions. In f64 a larger `C`
  (up to ~700) is possible; whether it helps is a measurement (§16.3), not a default change.

### 5.2 Network back end (D2; user change 2026-11-04: use an existing framework)

The network (forward for inference, forward and backward for training) is written with an
existing Rust NN framework, not by hand.

**Candidates evaluated** (2026-11-04, prototypes in the scratchpad, numbers in §5.3):

| | **burn 0.21** (`burn`, `burn-flex`, `burn-cuda`, `burn-autodiff`) | **candle 0.10** (`candle-core`, `candle-nn`) |
|---|---|---|
| Rust stable / ed. 2024 | MSRV 1.92: builds on the installed stable 1.93.1. **0.22 needs 1.95** (toolchain bump **[user]**) | 0.10.2 builds on 1.93.1. **0.11.1 fails on 1.93** (`__cpuid_count` unsafe call in `quantized/repack_x86.rs`); needs a newer toolchain |
| system libraries | none for CPU; CUDA via CubeCL + `cudarc 0.19` with **NVRTC at run time** (no nvcc at build) | none for CPU; CUDA feature builds `candle-kernels` with **nvcc at build time** (bindgen_cuda), plus cuBLAS |
| `cudarc` | 0.19.0 (cubecl-cuda 0.10), compatible with the workspace's 0.19.10 | 0.19.4, compatible |
| MHA | `burn-nn` `MultiHeadAttention` with `mask_pad` (bias on; divides scores by `√d_k` after the product) | none; ~40 lines of matmul and softmax (as prototyped) |
| LayerNorm | ε configurable (default 1e-5, set 1e-6); two-pass variance | `candle_nn::LayerNorm`, ε configurable; two-pass |
| GELU (tanh) | `activation::gelu_approximate` | `Tensor::gelu` (tanh), `gelu_erf` |
| masking | `mask_pad` / `mask_fill` | `where_cond` / additive `−inf` |
| optimizer | `AdamW` (weight decay 0 = Adam) + `LinearLrScheduler`, `CosineAnnealing…`, `ComposedLrScheduler` | `AdamW`; schedule by `set_learning_rate` each step (trivial) |
| autodiff | `Autodiff<B>` decorator, first order; generic over CPU, CUDA, wgpu | `Var` + `backward()`, first order |
| CPU f64 | via backend element type (to be checked for Flex; `NdArray<f64>` exists, deprecated in 0.22) | native f64 |
| other GPUs | `wgpu` (Vulkan): a GPU path in Blender without CUDA | Metal only |
| CPU forward, prototype (10k neighborhoods, f32) | 86–95/s with the inner rayon pool; **335/s with 32 single-threaded chunk workers** (§5.3); max abs diff vs JAX 2.2e-6 (raw outputs) | 42/s (`broadcast_matmul`), 68/s (2D GEMMs); max abs diff 2.2e-5 (scaled coefficients) |
| CPU hot spots | 14 % GEMM (`gemm` crate); the rest is unfused element-wise ops (`powf`, `tanh`, strided copies, page faults) | similar (not profiled) |
| build | ~45 s for `burn-flex`; CubeCL + autodiff adds several minutes | ~35 s CPU |

**Both CPU paths are far below the reference's XLA CPU path** (~1,050 neighborhoods/s under the
same outside load). The hardware bound is ~5,000–15,000/s: 218 MFLOP per neighborhood against
2–3.6 TFLOPS.
- A CPU-only 1M-point precompute at 335/s takes ~50 min.
- The paper's floor (≈ 40 s for 1M on a 3090, §16.5) must come from the CUDA back end. The 4080
  needs ~5–10 s of GEMM time at 20–40 TFLOPS.

**Recommendation: burn 0.21** (pinned `=0.21.x`), with `burn-flex` on the CPU and `burn-cuda` on
the GPU, both behind `Autodiff` for training. Reasons:
- MHA with pad masks, LayerNorm, tanh GELU, Adam and the warmup + cosine schedulers are built in.
- CUDA kernels are compiled at run time by NVRTC through `cudarc 0.19`, like the workspace's own
  `cuda` feature. Building the default and `cuda` features needs no CUDA toolkit, which keeps
  wheels and Blender builds simple.
- One code path is generic over backends, including `wgpu`, a future non-NVIDIA path.

Risks:
- burn's API churn (`burn-ndarray` deprecated in 0.22) and MSRV jumps (0.22 needs 1.95). Mitigation:
  pin the version, keep the network module small (`network/burn.rs`), and keep the weights in
  framework-neutral safetensors.
- CPU throughput as measured. Mitigation:
  - run fixed-size chunks on single-threaded workers in our pool (§5.4);
  - use burn's fused `LayerNorm`/`gelu_approximate` instead of composed ops;
  - measure again in WP1.
  If WP1 cannot reach ≥ 1,000/s on the CPU, the fallback is **[user]**: either accept the CPU speed
  (the GPU path carries the floor), or write a CPU-only inference kernel on the `gemm` crate. That
  kernel would be cross-checked against the framework to ~1e-6. It would be a deviation from
  "framework for inference", and only for the CPU.

**candle is the fallback.** It is simpler, cuBLAS GEMMs are fast, and its determinism is
documented. But nvcc is needed at build time for `cuda`, it has no MHA or schedulers, it was slower
on the CPU, and the current release 0.11 needs a newer toolchain.

**Not used.**
- `tch-rs`: libtorch, ~2 GB of C++ shared libraries. That breaks "no system libraries beyond
  CUDA" and bloats wheels and Blender.
- `dfdx`: unmaintained since 2024, CUDA via its own cudarc pin.

**Second-order question (Eq. 27).** The eikonal term needs `∇_x φ` inside the loss, then
`∂/∂θ`. The reference writes `∇φ` in closed form (App. D, N§9.1) and lets autodiff differentiate
that closed form once. We do the same:
- `∇g_i` (§7.1) and the blend gradient (§7.4) are written as framework tensor ops.
- Only first-order autodiff is needed; neither framework offers double backward.
- The FD checks of §15.2 cover the composition.

### 5.3 Measurements behind §5.2 (2026-11-04)

- Machine: Threadripper 3970X, 32 cores / 64 threads, Zen 2, AVX2, no AVX-512. RTX 4080 16 GB,
  CUDA 12.4.
- Outside load: the reference-build agent's golden runs, ~27 cores busy (load 33–75).
- Input: 10,000 neighborhoods of `holey-bunny.ply`, the repository weights converted to
  safetensors, and features precomputed in f32 with the reference's order.

| run | neighborhoods/s | max abs diff vs JAX CPU |
|---|---|---|
| JAX 0.8.2 CPU (XLA), whole batch | ~1,050 | — |
| candle 0.10.2, `broadcast_matmul`, chunk 1000 | 42 | 2.2e-5 (scaled) |
| candle 0.10.2, 2D GEMMs, chunk 1000 | 68 | 2.2e-5 (scaled) |
| burn 0.21 Flex, rayon inside, chunk 1000 / 2500 | 86 / 95 | 2.3e-6 (raw) |
| burn 0.21 Flex, 32 workers × 1 thread, chunk 100 | 335 | 2.2e-6 (raw) |

**Bitwise equality.** burn-flex with single-threaded workers and fixed chunks of 100 gave the same
output hash (`40a18017c88c2c93`) on two runs with 32 workers and one run with 8.

**Scaling.** 8 workers gave 296/s, 32 gave 335–345/s. That is poor scaling, probably allocation and
page-fault bound (`kernel_init_pages` at 5 % in the profile); buffer reuse is a WP1 item.

The prototypes and logs are in the scratchpad (`proto-burn/`, `proto-candle/`, `jaxref/`); WP1
re-measures in the crate.

### 5.4 Determinism standard (D4) — **[user] confirm**

**CPU.** Bitwise reproducible run to run and **independent of the thread count**:
- the network runs on fixed chunks of `CHUNK = 128` neighborhoods (a constant, not an option);
- each chunk runs on one worker with framework-internal threading off;
- our pool distributes the chunks;
- all our kernels (features, fitting, blend, MC, baselines) are per-output loops in a fixed order.

The framework's GEMM (the `gemm` crate, FMA microkernels) is deterministic for a fixed shape. That
is a property to *test* (1, 3, 8, 32 threads), not to assume.

**GPU, our own kernels** (blend, fitting, eval, tracer, MC, baselines): the workspace rule from
`ddg-repulsive-surfaces` §8.5:
- one thread per output in the CPU's operation order, gathers not atomics;
- fixed-shape pairwise trees for reductions;
- `--fmad=false`.
Level 1: bitwise run to run. Level 2: bitwise equal to the CPU where the CPU uses the same tree.

**GPU, the network.** The framework's matmuls and reductions do not follow the rule: they use FMA,
tiling chosen per shape and GPU, and possibly split-K. burn-cuda (CubeCL) uses its own matmul
kernels; whether any reduction uses atomics must be tested. cuBLAS (the candle fallback) is
documented bitwise reproducible run to run on the same architecture and SM count when
`CUBLAS_WORKSPACE_CONFIG` is fixed and atomics mode is off.

The proposed standard:
- **level 1** (same device, same driver): bitwise run to run, tested;
- versus the CPU: **documented tolerance**, coefficients within `1e-5·(1 + |a|)` (measured in WP4);
- training on the GPU: run to run bitwise if the framework allows it; otherwise
  reproducibility of loss curves within tolerance, and the trained weights published with their
  seed and log.

This relaxes the workspace's GPU rule for the framework part only, and needs the user's
confirmation. If burn-cuda turns out non-deterministic run to run (atomics), use cuBLAS through
candle for GPU inference, or accept tolerance-only.

### 5.5 Threads

`threads: Option<usize>` in every options struct. It is an explicit `rayon::ThreadPool`, never the
global pool. The framework's own pool is pinned to 1 thread per worker (`burn-flex` is built without
its `rayon` feature, or with it scoped).

## 6. Network (Sec. 4.3; N§2)

### 6.1 Architecture (exact; the weights define it)

- `embed 128`, `layers 8`, `heads 8` (head dimension 16), `mlp 512`, `k 64`, `n_freq 0`, no
  dropout.
- Pre-LN blocks. LayerNorm ε = 1e-6 with "fast" variance in the reference. Our framework's
  two-pass variance differs in the last bits; this is a parity tolerance, not a mismatch.
- MHA with biases on q, k, v and out; query scaled by `1/4`; **no attention mask** with the
  pretrained weights. A pad mask is available for training on padded clouds (D19).
- FFN with tanh GELU.
- Final LayerNorm, center token, `Linear 128→6`.

`[formula: pat-network]` is cited: an architecture definition, not a theorem.

### 6.2 Weights (D6)

**Format.** safetensors: an 8-byte header length, a JSON header, then raw little-endian f32.
- It is framework neutral and parsed by the `safetensors` crate (MIT/Apache, no dependencies).
- Tensor names are the Flax paths of N§2.4. Shapes are as in Flax (`[128, 8, 16]` kernels etc.);
  our loader reshapes them.
- The `__metadata__` holds `config` (embed, layers, heads, mlp, k, n_freq), the source pkl's
  sha256, the reference commit, the MIT notice, and `source = "points-as-tori models/FundamentalFormPredictor.pkl"`.

**Committed.** `crates/ddg-points-as-tori/weights/pat-reference.safetensors`, ~6.37 MB.
- The MIT notice is in `LICENSES/points-as-tori-MIT.txt`, and `weights/README.md` says where the
  file comes from and that the notice covers it.
- It is embedded with `include_bytes!` behind the default feature `bundled-weights`, so the Python
  wheel and Blender get it without a path.
- The other pkl (`FundamentalFormPredictor_pretrained.pkl`, N§3) is converted by the same tool
  into `data/points-as-tori/` (gitignored), not committed.

**Converter.** `tools/convert_weights.py` (run with `uv run`, PEP 723 dependencies numpy and
safetensors).
- A restricted `pickle.Unpickler.find_class` allows exactly the 7 globals of N§2.4. Any other
  global is an error.
- Python runs `-I`. Never `pickle.load` unrestricted.
- The converter checks the config and the 136 names and shapes, and writes a sha256 manifest.

**Loader errors.** `PatError::Weights { reason }` for a wrong magic, a missing tensor, a wrong
shape, or a non-finite value.

### 6.3 Inference API

```rust
pub struct Model { /* config + framework module, one per device */ }
impl Model {
    pub fn reference() -> Result<Model, PatError>;                    // bundled weights
    pub fn from_safetensors(bytes: &[u8]) -> Result<Model, PatError>;
    pub fn config(&self) -> &ModelConfig;
}
pub enum Device { Cpu, Cuda { ordinal: usize }, Auto }
pub struct PredictOptions { pub device: Device, pub threads: Option<usize>, pub feature_order: FeatureOrder }
pub fn predict_coefficients(model, cloud, &Neighborhoods, &PredictOptions) -> Result<(Vec<[f32; 6]>, Vec<f64>, PredictReport), PatError>;
```

- `Cuda` without the feature, driver or device is `PatError::Gpu`. `Auto` falls back to the CPU.
- `PredictReport` holds the device, the seconds per stage, the chunk count, and non-finite outputs
  (error, with indices).

## 7. Torus fitting and SDF (Sec. 4.1, 4.2; N§6, N§7)

### 7.1 Fitting

All in f64; the formulas are those of N§6.1:
- `A`, `I`, `II` (Eq. 18, App. C): `[formula: pat-fundamental-forms]`
- `H`, `K` (Eq. 19): `[formula: pat-mean-gauss]`
- the eigenvector pair `(v1 ± disc, v2)`: generalized eigenvectors of `(II, I)` with eigenvalues
  `κ± = H ± √(H² − K)`. This is equivalent to Eq. 20, whose `w±` is the other row's solution.
  `[formula: pat-principal-pair]`
- Rayleigh quotient equals `κ±`: `[formula: pat-curvature-rayleigh]`
- 3D directions `v± = w_x(ŝ + a10 n) + w_y(t̂ + a01 n)`, tangent to `Q*` at 0:
  `[formula: pat-principal-directions]`
- `n*` (Eq. 21): unit, orthogonal to `v±`. `[formula: pat-surface-normal]`
- min/max by `|κ|`: `κmin = κ₊` iff `|κ₊| < |κ₋|`, as in the reference.
- `r = 1/|κmax|`, `R = 1/|κmin| − sgn(κ₊κ₋) r`, so `(R + sgn(κ₊κ₋) r)⁻¹ = |κmin|`:
  `[formula: pat-torus-radii]`
- sign (Eq. 22): `[formula: pat-torus-sign]`
- center `c = p + a00 n − (sign/|κmin|) n*` and axis `u = n* × v̂min`: `[formula: pat-torus-frame]`
- **contact.** The torus passes through `o = Q*(0,0)`, its tube normal there is `±n*`, and its
  principal curvatures at `o` are `κmin` and `κmax` with directions `v̂min, v̂max`:
  `[formula: pat-torus-contact]` (the paper's claim, proved).

**Clamps** are options and are counted:

```rust
pub struct FitOptions { pub eig_eps: f64 /* 1e-10 */, pub curvature_floor: f64 /* 1e-10 */, pub threads: Option<usize> }
pub struct FitReport { pub clamped_disc: usize, pub clamped_eigvec: usize, pub clamped_curvature: usize,
                       pub umbilic: usize, pub nonfinite: Vec<u32>, pub tokens_lt_k: usize }
```

- **Defaults are the JAX training rule (D9):** clamp `disc ≥ √1e-10`, `ev` components, `ev·I·ev`,
  and **`|κ| ≥ 1e-10`**.
- The C++ inference version has no `|κ|` floor and gives NaN on exact zeros. Since the network
  was trained with the JAX rule, the JAX rule is the right one.
- `sgn` is `jnp.sign` (0 at 0) for `sgn(κ₊κ₋)`. The eigenvector clamps use the C++
  `sign(0) = +1`, because JAX's 0 at exact zero is the degenerate case of N§6.3 and only changes
  `ev₋` at the exact umbilic, where the axis is irrelevant (`R = 0`).
- No clamp is silent: the report counts them.
- A non-finite torus is an error with its indices (it cannot occur with the floors above; the check
  stays).

```rust
pub fn fit_tori(cloud, coefficients: &[[f32; 6]], &FitOptions) -> Result<(Tori, FitReport), PatError>;
pub fn fit_diagnostics(cloud, coefficients, i) -> FitDiagnostics;   // κ, directions (orientation, App.)
```

### 7.2 Per-torus SDF

- `g(x) = sgn(r)(‖d‖ − |r|)`, `d = (‖(x − c) × u‖ − R, ⟨x − c, u⟩)` (Eq. 23–24).
  `[formula: pat-torus-sdf]`
  - The exact Euclidean SDF of the torus for `R ≥ |r|`.
  - For `R < |r|` (spindle) it is the distance to the generating circle's tube; it is cited as
    Quilez, and the identity "distance to a circle in its meridian plane" is proved.
- Gradient `∇g = sgn(r)(w₁ (I − uuᵀ)(x − c)/v + w₂ u)` with `w = d/‖d‖` and `v = ‖(x − c) × u‖`.
  `[formula: pat-torus-sdf-gradient]`
- **Ours: Laplacian.** `Δg = sgn(r)(1/‖d‖ + d₁/(‖d‖ v))`: the 2D Laplacian of `‖d‖` plus the
  axisymmetric term. This replaces the reference's placeholder. `[formula: pat-torus-sdf-laplacian]`

**Singular sets.** Where `v = 0` (axis) or `‖d‖ = 0` (core circle), `g` is continuous but not
differentiable. We return the one-sided limit along a fixed rule:
- On the axis: `∇g = sgn(r) w₂ u` with `w = (−R, d₂)/‖·‖`.
- On the core: `∇g = 0` (a measure-zero set), counted in `EvalReport::singular`.

No NaN (D12).

### 7.3 Stable form (ours, D8; f32 GPU only)

- With `o = Q*(0,0)` on the torus, evaluate `x − c = (x − o) + (o − c)`. `o − c` is stored
  exactly in f64.
- Compute `v − R = (v² − R²)/(v + R)`, with `v² − R²` expanded around `o`. This avoids `O(R)`
  cancellation.
- `[formula: pat-torus-sdf-stable]` proves equality with Eq. 23.
- Used only by `EvalPrecision::F32Stable` on the GPU. The CPU is f64 Eq. 23.

### 7.4 Blend and gradient (Eq. 25–26; App. D)

With a neighbor set `𝒩(x)`, take `σ_x = ½ max_{i∈𝒩} ‖x − p_i‖`, `λ_x` by rule (below),
`w_i = exp(−λ_x(ρ_i − σ_x))` with `ρ_i = ‖x − p_i‖`, and `φ = Σ w g / Σ w`.

- **Shift cancels** (Eq. 25 equals Eq. 1 / Eq. 17 discretized): `[formula: pat-blend-shift]`
- **Weight range.** For the "per-query" rule with `λ = C/σ`, `w_i ∈ [e^{−C}, e^{C}]`:
  `[formula: pat-blend-weight-range]`
- **Convexity.** `φ ∈ [min g_i, max g_i]`: `[formula: pat-blend-convex]`
- **Gradient with `λ`, `σ` frozen** (App. D, the reference's `∇φ`): `∇φ = A + φB`, with
  `A = Σ w(∇g − λ r̂ g)/W`, `B = λ Σ w r̂/W`. `[formula: pat-blend-gradient]`
- **Exact gradient** where the set and the argmax are locally constant (ours): add
  `(∂λ/∂x)·Σ w(σ − ρ)(g − φ)/W`, with `∂λ/∂x = −(C/σ²)·½ r̂_{argmax}`.
  `[formula: pat-blend-gradient-exact]`
- `r̂_i` at `x = p_i`: defined as 0. That is the limit of `w r̂ g` only in the sense of a fixed rule;
  `g_i` itself is smooth there. Counted in the report.

```rust
pub enum Neighborhood { Knn { k: usize /* 32 */ }, All, Radius { r: f64, fallback_k: usize /* 32 */ } }
pub enum LambdaRule { PerQuery { c: f64 /* 64 */ },               // reference (released)
                      Constant { lambda: f64 },
                      PaperHeuristic { c: f64 /* 64 */, factor: f64 /* 1e3 */, k_mean: usize /* 64 */ } } // λ = 10³/D, R_eval = 2C/λ
pub enum GradientMode { Frozen /* reference, App. D */, Exact }
pub struct EvalOptions { pub neighborhood: Neighborhood, pub lambda: LambdaRule, pub gradient: GradientMode,
                         pub area_weights: Option<Vec<f64>>, pub isovalue: f64, pub precision: EvalPrecision,
                         pub device: Device, pub threads: Option<usize> }
pub struct TorusField { /* cloud, tori, kd-tree, options */ }
impl TorusField {
    pub fn new(cloud, tori, EvalOptions) -> Result<Self, PatError>;
    pub fn distance(&self, x) -> f64; pub fn distance_gradient(&self, x) -> (f64, [f64; 3]);
    pub fn laplacian(&self, x) -> f64;                                    // §7.2 + App. D divergence, ours
    pub fn evaluate(&self, xs: &[[f64; 3]], out: EvalOut) -> Result<EvalReport, PatError>;  // batch
}
pub struct PointsAsTori;   // facade: build(cloud, model, opts) -> (TorusField, BuildReport): kNN + predict + fit
```

**Defaults (D13).** `Knn { k: 32 }`, `PerQuery { c: 64 }`, `Frozen`, no area weights. These are the
released reference, so the black-box comparisons and the pretrained weights' expectations hold.

`PaperHeuristic` implements Sec. 4.2 as printed: `λ` is constant, `σ_x` is still the shift, and
the 32-NN fallback applies when the radius set is empty. It is **without** the reference's dead
`λ·R/√3` factor (N§7.3).

**Edge cases.**
- An empty radius set falls back to kNN.
- `W` underflow cannot occur: the max weight is `e^{λσ} ≥ 1` for the nearest point.
- With `Exact` and `Radius`, the set changes discontinuously; documented.

### 7.5 Batch evaluation and GPU

**CPU.** Queries are chunked on the pool; the kNN per query uses the kd-tree.

**GPU** (feature `cuda`, our own NVRTC kernels in `src/gpu/*.cu`, following §5.4):
- The kd-tree is flattened: nodes and buckets in SoA arrays, with a stackless or fixed-stack
  traversal.
- Exact 32-NN per query is an insertion list of 32 sorted by `(d², index)`, the same order as the
  CPU.
- Then the blend in the CPU's order. f64 by default.
- Used by MC grids, the tracer and the evaluation. Level 2 (bitwise the CPU) is the target in f64
  with `--fmad=false`.

## 8. Applications (Sec. 5.3; user: all)

### 8.1 Offsets

`φ − c`: `isovalue = c`, no refit (as the reference). `[formula: pat-offset]` (the level set of
an exact SDF at `c` is the offset; cited for the approximation).

### 8.2 Booleans

`union = min(φ_A, φ_B)`, `intersection = max`, `difference = max(φ_A, −φ_B)`.

- Gradients are taken from the active argument; at ties, from `A`.
- Option: smooth min with radius `k` (Quilez polynomial smin), ours, off by default.
- `[formula: pat-boolean]`: min/max of SDFs is an exact SDF outside, and a bound inside (cited).

```rust
pub enum Csg { Field(TorusField), Union(Box<Csg>, Box<Csg>), Intersect(..), Subtract(..), Offset(Box<Csg>, f64) }
```

Every application takes `&dyn ImplicitField` (`distance`, `distance_gradient`, `bbox`).

### 8.3 Morphology (Fig. 24)

Erosion by `c`:
1. Sample the `−c` level set.
2. Take normals `∇φ/‖∇φ‖`.
3. Refit tori to the new cloud (network + fit).

Dilation is the same with `+c`.

**Sampling a level set (ours, D14).**
1. Marching cubes at resolution `res` (§8.6).
2. Area-uniform points on the triangles: target count = the input's count by default.
3. One to three Newton projections `x ← x − (φ(x) − ℓ)∇φ/‖∇φ‖²`, each with a residual check.

The reference's `sample_isosurface` ray marcher is buggy (N§17.6). The paper does not say how it
samples.

```rust
pub fn morphology(field, model, op: Morph { Erode(c), Dilate(c), Open(c), Close(c) }, &MorphOptions) -> Result<(TorusField, MorphReport), PatError>;
```

### 8.4 Normal orientation (Fig. 27; not in the reference, D16)

Iterate:
1. For each `i`, `ñ_i = Σ_{j∈N_k(i), j≠i} exp(−M_i(p_i − p_j)) ∇φ(p_j)`, with
   `M_i(r) = (κmin²⟨r, v̂min⟩² + κmax²⟨r, v̂max⟩² + ε²⟨r, n_i⟩²)^{1/2}`.
2. Flip `n_i` if `⟨n_i, ñ_i⟩ < 0`.
3. Refit and repeat until no flips or `max_iterations` (10, as Fig. 27).

Unspecified by the paper; our choices:
- the neighbors are `N_k(i)` with k = 64;
- `ε = 1/σ_i`, so the normal direction is weighted like a unit curvature at the neighborhood scale
  (the paper's `ε` is unspecified; measured in §16.3);
- simultaneous (Jacobi) flips;
- `∇φ(p_j)` with the full field (its own torus included).

`[formula: pat-orientation-metric]` (`M_i` is a norm when `κmin, κmax, ε ≠ 0`).

Report: flips per iteration, and the final count of `⟨n_i, ∇φ(p_i)⟩ < 0`.

### 8.5 Sphere tracing (Fig. 12)

`render::Tracer` on CPU (pool, per pixel) and GPU (`cuda`).

**Traversal.**
- Rays clip to the inflated bounding sphere (as the demo).
- Step `t += |φ − ℓ| / max(‖∇φ‖, g_min)·safety`, with `safety = 0.9` default (the demo uses 1.0
  and divides by `‖∇φ‖`; D17).
- Hit when `|φ − ℓ| < ε_hit` (1e-3 of the bbox diagonal); at most 128 steps (the demo uses 40).

**Shading.**
- Normal = normalized analytic gradient.
- Lambert plus the demo's IBL-GGX and Fresnel terms with its constants (N§14). The cube map is the
  demo's `ForestBlurred*.png` (the reference's MIT assets, copied with the notice) or a constant
  sky.
- Ground plane with a hard shadow; optional cut-plane slice with the demo's color map and isobands.

**Output.** PNG via the `png` crate (small, MIT/Apache). The camera is explicit
(`look_at`, `fov`, `resolution`).

The scene contains the CSG of §8.2, so offsets and Booleans render directly. Images are
deterministic (§5.4).

### 8.6 Marching cubes (D18)

- Our own, on a regular grid over a box, with sparse narrow-band evaluation: evaluate on a coarse
  grid; skip blocks with `|φ| > block diagonal·L` (Lipschitz bound `L = 1` from `‖∇φ‖ ≈ 1`, with
  a safety factor). This follows the paper's "masking voxels not within 8 voxels of the contour".
- Lorensen–Cline table with the **asymptotic decider on ambiguous faces** (Nielson–Hamann), so
  the output is crack-free.
- Edge vertices are shared, so the mesh is welded and watertight except at the box boundary.
- Vertices are placed by linear interpolation, then optionally projected by one Newton step.
- Output: `ddg_core` mesh plus per-vertex normals from `∇φ`.
- `[formula: pat-mc-interpolation]` (the root of the linear interpolant on an edge).

### 8.7 The 2D variant (Fig. 10; D15 **[user]**)

Not in the reference (N§13). Proposal: `pat2d` with circles, Eq. 25–26 in 2D, and per-point
circles from a *classical* local fit:
- weighted least-squares quadratic `Q(s) = a0 + a1 s + a2 s²` over the k nearest points in the
  tangent line;
- `κ = 2a2/(1 + a1²)^{3/2}`, center `Q*(0) − (sign/|κ|) n*`, SDF `sign·(‖x − c‖ − 1/|κ|)`;
- `|κ| < κ_floor` is a line (half-plane SDF).

Alternatives:
- (b) train a 2D network with the same pipeline (§10) on 2D shapes from the analytic SDFs of §11:
  new data, about a day of GPU;
- (c) a 2D slice of the 3D method (not Fig. 10).

Recommendation (a): it illustrates Fig. 10 with no learned part and is cheap. (b) is a research
add-on.

`[formula: pat2d-circle-fit]`, `pat2d-circle-sdf`, `pat2d-blend` (2D analogues).

## 9. Baselines (Sec. 5.2; user: all)

| method | definition | implementation |
|---|---|---|
| SHC (Eq. 28) | `φ = sgn(w)(−log|w|/λ + σ_x)`, `w = Σ A_i (λρ+1)⟨x−p_i, n_i⟩/(2πρ³) e^{−λ(ρ−σ)}` | same neighbor and λ rules as PAT (reference `signed_log_sum_exp`, N§7.6); gradient in closed form `[formula: pat-shc]`, `pat-shc-gradient` |
| SSPD (Eq. 29) | `Σ A⟨x−p, n⟩e / Σ A e` | `[formula: pat-sspd]`, gradient |
| point areas | Barill et al. 2018 §5: Delaunay triangulation of the k = 20 neighbors projected to the tangent plane, Voronoi area of the center | our own planar Delaunay (k ≤ 32, incremental with exact predicates via `robust`-style orientation in f64 with an adaptive fallback, or the workspace's if present) `[formula: pat-point-area]` (cited) |
| GWN, brute force | `Σ A_i ⟨n_i, p_i − x⟩/(4π‖p_i − x‖³)` (sign: inside = 1) | `[formula: pat-winding-dipole]` |
| fast winding number | Barill et al. 2018: octree, order-2 Taylor expansion, `β = 2` | our own, independent of libigl's code (MPL2 core is readable; we write from the paper) `[formula: pat-fwn-expansion]` (terms proved as Taylor coefficients; error bound cited) |
| regularized WN | Chen et al. 2024: `S(t) = erf t − (2t/√π)e^{−t²}`, `t = ρ/ε` | `[formula: pat-regularized-wn]` |
| LogSumExp (Madan–Levin 2022) | unsigned `−log Σ e^{−λρ}/λ` | `[formula: pat-lse]` |
| SPSR | Open3D `create_from_point_cloud_poisson` | black box through `uv run --with open3d` (Python 3.12), depth matched to the grid (Sec. 5.2); mesh SDF by our exact distance + GWN sign |
| SHM (Feng–Crane 2024) | grid signed heat method | **not planned** (the user did not list it). It would be a separate paper's port **[user]** |
| NN-VIPSS | — | not planned (not listed) |
| the authors' program | black box | golden outputs in `runs_tmp/pat-ref/` from the other agent; driver protocol in §16.7 |

**Ground truth for meshes:**
- exact unsigned distance via a triangle BVH (ddg-core / ddg-repulsive-curves `Obstacle` closest
  point), parallel;
- sign by GWN (Jacobson 2013: brute force for ≤ 20k faces, our FWN above);
- the JAX-compatible `gwn > π` threshold option for training parity.

## 10. Training (Sec. 4.3, 5.1; user: port to Rust, retrain, compare)

### 10.1 What is reproduced (N§9)

- The loss `L = mean|φ − φ_true| + mean|1 − ‖∇φ‖|` (Eq. 27) with the code's rules:
  - prediction by the blend over all valid tori of the cloud;
  - `σ` from the 64 nearest valid points to the query;
  - `λ = 64/σ`;
  - frozen-λ closed-form gradient.
  `[formula: pat-loss]`
- 120 queries per neighborhood (40 zero set / 40 narrow band ±0.2 / 40 cube √3·extent). Both the
  code's and the paper's box sizes are options; the code's are the default (D20).
- Adam (β 0.9/0.999, ε 1e-8, no weight decay), linear warmup, then cosine. The code's effective
  constant (`decay_steps = 1e8`) is the default; the paper's real cosine over the phase is an
  option (D21).
- The curriculum of N§9.3 including pretraining (phase table as data, `training/curriculum.toml`).
- Init: last layer zero kernel, bias `(0, 0, 0, 0, −½, −½)`; others lecun-normal. Seeded, D22.

### 10.2 What is fixed relative to the reference (D19, ours)

- **Padded clouds.**
  - Features and network run only on valid rows.
  - Attention uses a pad mask for `n < 64`.
  - No −1 → row 2047 indexing.
  - So the 3,588 padded simulated shapes contribute real gradients instead of NaN-zeroed steps
    (N§17.7).
  - A `reference_padding` option reproduces the reference's behavior for parity runs.
- **Normals of length 2** in the pretraining blob and validation files (N§10). They are
  normalized at load and counted; option `raw_normals` for parity.
- **NaN handling.**
  - Non-finite tori are masked, as in the reference.
  - NaN gradients are not zeroed silently. A step with a non-finite loss or gradient is skipped
    and counted in the `TrainReport`; the threshold is an option.
- **Step semantics.** All-padding query slices are skipped (the reference steps Adam on zero loss).
  Option to keep them.

### 10.3 Data (`train::data`)

**npz reader.** Our own, on the `zip` crate (pure Rust, `deflate` feature) plus a `.npy` header
parser for f64, i32, bool, C order. Members as in N§10. All files are validated:
- shapes agree;
- masks are prefixes;
- neighbor column 0 is self;
- indices are in range or −1.

Errors name the file, member and shape.

**Streaming.** One file at a time with a shuffle buffer. Phase 4 would otherwise need ~16 GB.
Deterministic order: files sorted by name (the reference uses unsorted `glob`, N§9.2), then
`permutation(seed, epoch)`.

**Ground truth.**
- Our exact point–triangle distance and GWN, in f64 on CPU or GPU (one thread per (row, query),
  gathers).
- Option: JAX's f32 `gwn > π` for parity.

### 10.4 Implementation with the framework

**Forward.**
- The network runs on all valid neighborhoods of the cloud: 2048 × 64 tokens × 8 layers.
- Fitting (§7.1, the JAX clamp rule) and the blend run as framework ops: f32, whole cloud,
  batched queries.

**Memory** (N training survey). Activations for 2048 × 64 tokens are about 10 GB, so the plan is
two passes per step:
1. Coefficients without the tape.
2. `∂L/∂a` from the blend with the tape on the 6 × 2048 coefficients only.
3. The network re-run in chunks of 256 neighborhoods with the tape, back-propagating
   `⟨∂L/∂a, a(θ)⟩`.

This is exact (the loss depends on θ only through `a`) and fits 16 GB.

**FD checks** (§15.2), in f64 on the CPU backend:
- `∂L/∂a` against central differences of `L(a)`;
- `∂L/∂θ` on a reduced network (2 layers, 2 heads, 16 tokens) against FD in f64;
- the blend's `∇φ` against FD of `φ` with frozen `λ`.

**Expected time** (RTX 4080).
- ≈ 1.34 TFLOP per step network forward + backward, plus ≈ 0.15 TFLOP blend, plus ground truth
  amortized: about 70–120 ms per step at 40–60 % of f32 peak with fused blend kernels.
- About 10–17 h for the ≈ 498k steps of the curriculum.
- With an unfused framework blend (materialized `1024 × 120 × 2048` tensors, ≈ 1 GB each) closer
  to the paper's 0.33 s per step, i.e. ≈ 45 h.
- The blend and ground truth may be written as custom CubeCL/NVRTC kernels if the framework's ops
  dominate (measured in WP10).

```rust
pub struct TrainOptions { pub curriculum: Curriculum, pub batch_rows: usize /* 1024 */, pub queries: QueryOptions,
                          pub optimizer: AdamOptions, pub schedule: Schedule, pub seed: u64, pub device: Device,
                          pub checkpoint_dir: PathBuf, pub reference_padding: bool, pub raw_normals: bool }
pub fn train(data_dir, init: Option<&Model>, &TrainOptions, progress: &mut dyn FnMut(&StepLog)) -> Result<(Model, TrainReport), PatError>;
```

Feature `train` (off by default) pulls in autodiff and the npz reader.

### 10.5 Retraining and comparison

1. Full curriculum from scratch, with seed and log committed.
2. Then evaluate on §16.2 against the authors' weights. Validation curves go to `TRAINING.md`.
3. Weights are published in `data/` (not committed unless the user asks).

## 11. Analytic test shapes (App. E; N§12)

- The 16 Quilez 2D SDFs (the 14 of App. E plus heart and isosceles), with exact gradients.
  - Revolution about y: `φ₂(‖(x,z)‖ − o, y)`.
  - Extrusion: `min(max(w), 0) + ‖max(w, 0)‖`.
  - `[formula: pat-revolution]`, `pat-extrusion`, `pat-extrusion-gradient`. The exactness of the 2D
    SDFs is cited (Quilez); the extrusion/revolution identities are proved.
- Parameters are drawn as App. E (`r ∈ [0.1, 0.9]`, `z ∈ [0.1, 0.9]`, `t ∈ [π/6, 5π/6]`,
  `o, h ∈ [0.1, 0.3]`). Perturbation directions are **uniform on the circle** (the reference's
  first-quadrant bias is a reference option).
- Each instance gets its own seed.
- The 70 instances are 14 shapes × 5, seeds 0..69, recorded in `data/points-as-tori/analytic.json`.
- Surface sampling:
  - area-uniform by construction: extrusion caps and sides are parametric; revolution by area
    weighting of the 2D boundary times `2π·radius`;
  - from 2D boundary samples, so no ray-casting bias (the reference's ray sampler is approximate
    and racy, N§12);
  - the ray sampler stays as a reference option.

## 12. Python and Blender

### 12.1 Python (`crates/ddg-py`, module `ddg_tools.points_as_tori`)

- `build(points, normals, device="cpu", threads=None) -> Field`.
- `Field.distance(q)`, `.gradient(q)`, `.distance_and_gradient(q)`.
- `Field.save(path)`, `load(path)`.
- `marching_cubes(field, res, bbox, level)`.
- `boolean(a, b, op)`, `offset`, `morphology`, `orient_normals`, `render(...) -> ndarray`.
- `baselines.*`.
- NumPy arrays in and out; errors as `ValueError` with the Rust message.
- The bundled weights are included. The wheel size grows by ~6.4 MB.
- burn without `cuda` adds about 3–6 MB to the extension module (to be measured; candle is similar).
  The `cuda` feature adds CubeCL (runtime NVRTC, no CUDA libraries in the wheel).

### 12.2 Blender (user decision 2026-11-04)

Operator "Points as Tori" on a point-cloud object:
- a mesh with vertices and normals (custom normals or the vertex normal attribute), or a PLY
  import.

**Cache.** The tori are built once and cached per object, keyed by the hash of positions and
normals. Stored in a sidecar file next to the `.blend`, or in memory for the session.

**Modes:**
- **Surface**: output a mesh at level `ℓ` and grid resolution (marching cubes).
- **Offset**: same, with `ℓ = c`.
- **Boolean** with a second point-cloud object (union, intersection, difference); both tori sets
  cached.
- **Orient normals**: writes the flipped normals back to the object's attribute, with a report of
  the flips.

**Progress.** Progress and cancel run through a modal operator, like the repulsive operators. GPU
(`cuda`) is used when available (`Device::Auto`).

**Files.**
- `blender/ddg_uv_unwrap/points_as_tori_operator.py`, `points_as_tori_io.py`;
- registration lines;
- `blender/tests/{test,ui}_points_as_tori.py`.

## 13. Formula ids

Rust functions carry `[formula: <id>]`. Each id is a Lean theorem in
`lean/proofs/DdgProofs/PointsAsTori/*.lean` ("proved") or a heading in `lean/proofs/CITED.md`
("cited").

| id | statement | Rust (planned) | Lean / CITED |
|---|---|---|---|
| `pat-median-scale` | median of an even list = mean of the middle two; `σ(c·P) = cσ(P)`, translation invariant | `features` | proved, `Features.lean` |
| `pat-local-frame` | `ŝ, t̂, c` orthonormal and right-handed when `‖c‖ = 1` (both branches of the switch) | `features` | proved, `Features.lean` |
| `pat-features` | definition; invariance under translation and uniform scaling of the cloud | `features` | proved, `Features.lean` |
| `pat-feature-invariance` | as above (for any frame rule depending only on `c`) | `features` | proved, `Features.lean` |
| `pat-network` | architecture (Vaswani 2017; Flax `MultiHeadAttention`, `LayerNorm`) | `network` | cited |
| `pat-output-scaling` | scaling the coefficients by `σ` makes tori and SDF scale covariant | `network` | proved, `Fitting.lean` |
| `pat-fundamental-forms` | Monge patch `Q*`: `I = [[1+a10², a10a01],[·, 1+a01²]]`, `II = (1/A)[[2a20, a11],[·, 2a02]]` at 0 (App. C) | `fit` | proved, `Fitting.lean` |
| `pat-mean-gauss` | `H = tr(I⁻¹II)/2`, `K = det II/det I` equal Eq. 19 | `fit` | proved, `Fitting.lean` |
| `pat-principal-pair` | `(II − κ±I)(v1 ± disc, v2) = 0` with `κ± = H ± √(H² − K)`, `disc` as N§6.1 (polynomial identity with `disc² = Δ`) | `fit` | proved, `Fitting.lean` |
| `pat-curvature-rayleigh` | Rayleigh quotient of a generalized eigenvector = its eigenvalue | `fit` | proved, `Fitting.lean` |
| `pat-principal-directions` | `v±` tangent to `Q*` at 0; `v₊ ⟂_I v₋` when `κ₊ ≠ κ₋` | `fit` | proved, `Fitting.lean` |
| `pat-surface-normal` | Eq. 21 is the unit normal of `Q*` at 0 | `fit` | proved, `Fitting.lean` |
| `pat-torus-radii` | `R + sgn(κ₊κ₋) r = 1/|κmin|` | `fit` | proved, `Fitting.lean` |
| `pat-torus-sign` | Eq. 22 equals "`κmax < 0 ∧ κmin ≠ 0`" | `fit` | proved, `Fitting.lean` |
| `pat-torus-frame` | center and axis (definitions); `u` unit, `u ⟂ n*` | `fit` | proved, `Fitting.lean` |
| `pat-torus-contact` | `g(Q*(0,0)) = 0`; the torus normal at `o` is `±n*`; principal curvatures of the torus at `o` are `κmin, κmax` (meridian `1/r`, parallel `1/(R + r)` with signs) | `fit` | proved (`Torus.lean`); curvature-of-torus facts proved for the parametrization |
| `pat-fit-contact` | `Geometric` rule (center `o + n*/κmin`): `g(o) = 0`, `∇g(o) = n*`. `Reference` rule: `g(o) = 0`, `∇g(o) = sgn(κmin κmax) n*` (−1 exactly at saddles; the torus there has curvatures `(−κmax, −κmin)`, the mirror image). Needs `κmin κmax ≠ 0` | `fit` (`SignRule`) | proved, `Torus.lean` (`PAT.fit_contact`) |
| `pat-torus-sdf` | Eq. 23 = distance to the torus surface for `R ≥ |r|` (the distance to a circle in its meridian half-plane) | `sdf` | proved for the meridian reduction, `Torus.lean`; exactness of the "circle twice" composition cited (Quilez) if the minimization is heavy |
| `pat-torus-sdf-gradient` | `∇g` (§7.2) | `sdf` | proved (`HasFDerivAt`), `Torus.lean` |
| `pat-torus-sdf-laplacian` | `Δg = s(1/‖d‖ + d₁/(‖d‖v))` | `sdf` | proved (2D identity + axisymmetric Laplacian), else FD-tested and cited |
| `pat-torus-sdf-stable` | stable form equals Eq. 23 | `gpu`, `sdf` | proved, `Torus.lean` |
| `pat-blend-shift` | the shift cancels | `blend` | proved, `Blend.lean` |
| `pat-blend-weight-range` | `λ = C/σ`, `0 ≤ ρ ≤ 2σ` ⇒ exponent ∈ [−C, C] | `blend` | proved, `Blend.lean` |
| `pat-blend-convex` | `φ` is a convex combination | `blend` | proved, `Blend.lean` |
| `pat-blend-gradient` | App. D with frozen `λ`, `σ` | `blend` | proved, `Blend.lean` |
| `pat-blend-gradient-exact` | adds the `∂λ` term | `blend` | proved, `Blend.lean` |
| `pat-blend-laplacian` | App. D divergence formula (with our `Δg`) | `blend` | proved (algebra of divergence), `Blend.lean` |
| `pat-laplace-method` | Eq. 15 and the self-normalized limit (Eq. 17) | docs | cited (Evans §4.5; Bender–Orszag §6.4) |
| `pat-signed-hopf-cole` | Eq. 8–14, App. A–B | docs | cited (the paper, App. A–B; Lipman 2021) |
| `pat-shc`, `pat-shc-gradient` | Eq. 28 and its gradient | `baselines` | proved (gradient), `Baselines.lean` |
| `pat-sspd`, `pat-sspd-gradient` | Eq. 29 | `baselines` | proved, `Baselines.lean` |
| `pat-winding-dipole` | point GWN term = solid-angle density of a dipole | `baselines` | cited (Barill 2018 Eq. 6) + definition |
| `pat-fwn-expansion` | order-0/1/2 Taylor terms of the dipole kernel | `baselines::fwn` | proved (derivatives), `Baselines.lean`; error bound cited |
| `pat-point-area` | Barill §5 point areas | `baselines` | cited |
| `pat-regularized-wn` | Chen et al. 2024 regularization | `baselines` | cited |
| `pat-lse` | Madan–Levin 2022 | `baselines` | cited |
| `pat-loss` | Eq. 27 as coded (definition) | `train` | proved (definition), `Loss.lean` |
| `pat-adam` | Adam update (Kingma–Ba 2015) | framework | cited |
| `pat-offset`, `pat-boolean` | offset and CSG of SDFs | `apps` | cited (Quilez; Ricci 1973) + min/max bound proved |
| `pat-orientation-metric` | `M_i` is a norm for nonzero weights | `apps::orient` | proved, `Apps.lean` |
| `pat-mc-interpolation` | linear edge root | `mc` | proved, `Apps.lean` |
| `pat-revolution`, `pat-extrusion`, `pat-extrusion-gradient` | Quilez revolution/extrusion and gradients | `analytic` | proved (gradients), cited (exactness) |
| `pat-quilez-2d` | the 16 2D SDFs | `analytic` | cited (Quilez 2019) |
| `pat2d-circle-fit`, `pat2d-circle-sdf`, `pat2d-blend` | §8.7 | `pat2d` | proved, `Pat2d.lean` |

## 14. Decisions: the paper as printed, the reference, and what we do

| D | topic | paper | reference | ours |
|---|---|---|---|---|
| D1 | precision | f32 constant C = 64 | network f32, fitting and eval f64 (C++), shader f32 | network f32; geometry f64; GPU eval f64, f32 stable opt-in (§5.1) |
| D2 | NN implementation | JAX | JAX/Flax | **burn 0.21** (framework, user change); candle fallback (§5.2) **[user: confirm framework and toolchain policy]** |
| D3 | feature arithmetic | — | f32, `p/σ − p₀/σ` | f64 difference-first then f32; reference order as option |
| D4 | determinism | — | OpenMP per query (deterministic), races in samplers | CPU bitwise, thread-independent; own GPU kernels by the workspace rule; framework GPU part run-to-run + tolerance **[user]** |
| D5 | threads | — | OpenMP global | explicit pool; framework single-threaded per chunk |
| D6 | weights | pkl | pickle of an NNX State | safetensors, committed 6.4 MB + `include_bytes!`, restricted-unpickler converter |
| D7 | axis storage | vector `u` | spherical angles | unit vector (angles only for export) |
| D8 | GPU SDF precision | — | f32 GLSL, unstable | f64 default; `F32Stable` (§7.3) |
| D9 | fitting clamps | none | JAX: 1e-10 on disc, eigvec, denominators, `|κ|`; C++: no `|κ|` floor (NaN on zero) | JAX rule (the trained regime), counted in `FitReport`; NaN impossible |
| D10 | kd-tree | nanoflann | nanoflann | own deterministic kd-tree (§3.1) |
| D11 | self in neighborhood | "`p_i` and its k nearest" (k+1 points?) | 64 including self, self first by nanoflann | 64 including self, self first by construction; `k` as the model's config |
| D12 | singular points of `g`, `r̂` at a data point | — | NaN (C++), 1e-10 guards (JAX) | fixed one-sided rules, counted (§7.2, §7.4) |
| D13 | evaluation neighborhood and λ | `R_eval = 2C/λ`, `λ = 10³/D`, 32-NN fallback | 32-NN, `λ = 64/σ_x` per query; radius mode dead and with a `R/√3` factor | default = reference (released); `PaperHeuristic` as printed without the dead factor; `All` (Eq. 26 exactly); **measured comparison in §16.3** |
| D14 | level-set sampling (morphology) | "sample the −c level set" | buggy ray marcher | MC + area sampling + Newton projection |
| D15 | 2D variant | Fig. 10 concept | absent | classical circle fit (§8.7) **[user: (a) classical, (b) trained 2D network, (c) skip]** |
| D16 | orientation | Fig. 27 procedure, `ε` unspecified | absent | §8.4: `ε = 1/σ_i`, Jacobi flips, k = 64, ≤ 10 iterations; `ε` measured |
| D17 | sphere tracing | Hart 1996 | step `|φ|/‖∇φ‖`, 40 steps, ε 1e-2 | step `0.9|φ|/max(‖∇φ‖, g_min)`, 128 steps, ε relative; demo constants as option |
| D18 | contouring | scikit-image MC, 1024³, 8-voxel mask | — | own MC + asymptotic decider, narrow band by Lipschitz bound |
| D19 | training padding | — | −1 → row 2047, no mask, NaN-zeroed steps | valid rows only + pad mask; `reference_padding` option |
| D20 | query boxes | cube side = extent; band width 0.2 | cube √3·extent, band ±0.2 | code values default (they trained the weights); paper values option |
| D21 | schedule | warmup + cosine | warmup + constant (decay over 1e8) | code default; real cosine option |
| D22 | data order and seeds | — | unsorted glob, fold_in keys | sorted files, our RNG (ChaCha8 seeded per (epoch, shape, row)); reference keys not reproduced (JAX threefry); documented |
| D23 | normals of length 2 in data | — | used raw in fitting | normalized at load, counted; `raw_normals` option |
| D24 | Laplacian | App. D | placeholder `Δg` | exact `Δg` (§7.2) |
| D25 | analytic shapes | 14, App. E ranges | 16, first-quadrant perturbations, racy sampler | 16 available, 14 × 5 = 70 evaluated with recorded seeds, uniform directions, area-uniform sampling |
| D26 | point areas | Barill via libigl (CGAL, GPL) | same | own from Barill 2018 §5 |
| D27 | baselines beyond the user's list | SHM, NN-VIPSS, GWN, SPSR | — | SHC, SSPD, GWN/FWN, regularized WN, LSE, SPSR (Open3D), authors' program; SHM and NN-VIPSS not planned **[user]** |
| D28 | area weights in the blend | omitted | option off; shader brute force on | option `area_weights`, default off |

## 15. Verification

### 15.1 Layers

1. **Analytic cases.**
   - Sphere coefficients give the tangent sphere.
   - Exact tori are fitted from the coefficients of an analytic torus's Monge patch, and the SDF
     reproduces it.
   - A plane gives a half-space up to `1e-10` curvature.
   - Cylinders.
   - Blend of identical tori = the torus.
2. **Finite differences** (§15.2).
3. **Invariances.** Translation and scaling of the cloud transform tori and SDF covariantly
   (network included); a rotation by a multiple of the frame switch is *not* invariant (documented).
4. **Parity with the reference.**
   - JAX CPU oracle coefficients (`tools/jax_oracle.py`, the scratchpad script made into a tool;
     Flax re-declaration, no reference import): coefficients within `2e-5` abs on the `data/`
     clouds.
   - The black-box golden files in `runs_tmp/pat-ref/` (other agent): coefficients within the
     JAX-GPU/TF32 tolerance; tori within `1e-6` relative given the same coefficients; SDF values
     within `1e-9` given the same tori (f64 vs f64).
5. **Lean oracle.** Features for one neighborhood, fitting from given coefficients, torus SDF and
   gradient, blend and gradient, SHC/SSPD.
6. **Lean proofs** (§13).
7. **Determinism tests** (§5.4): 1, 3, 8, 32 threads, two runs, bitwise; GPU level 1 and level 2
   where claimed.

### 15.2 Finite differences

Central differences, f64:
- `∇g` with step `1e-6·max(1, ‖x − c‖)` and relative error `< 1e-7`, away from the singular sets;
- `∇φ` frozen against FD with `λ`, `σ` held fixed; `Exact` against plain FD of `φ`;
- `Δg` and `Δφ` against second differences, `< 1e-5`;
- SHC and SSPD gradients;
- the loss: `∂L/∂a` and the reduced-network `∂L/∂θ` (§10.4) in f64, relative `< 1e-5`.
  The L1 kinks are avoided by query choice (`|φ − φ_true| > 10⁻³`).

### 15.3 Lean oracle (`lean/oracle`, `DdgOracle/PointsAsTori*.lean`)

- Float implementations of `pat-features` (given a neighborhood and frame rule), fitting, torus SDF
  and gradient, blend (all three λ rules), and SHC/SSPD.
- Fixtures in `lean/oracle/fixtures/pat/`, replayed by `tests/oracle.rs`.
- Rust equals the oracle to `1e-12` relative (f64), with the operation order documented.

### 15.4 Lean proofs and citations

`lean/proofs/DdgProofs/PointsAsTori/{Features,Fitting,Torus,Blend,Baselines,Loss,Apps,Pat2d}.lean`;
`CITED.md` section "Points as Tori".

## 16. Measurement plan

### 16.1 Parity

- Coefficients, tori and SDF on the reference's `data/` clouds and on 3 training clouds against
  the JAX oracle and the golden files (§15.1-4).
- Report the max/mean differences in `COMPARISON.md`.

### 16.2 Accuracy (Fig. 17–19)

**Datasets:**
- 100 ABC and 100 blob shapes not seen in training: from `ABC17_validation_3` and the training
  data's held-out split; the ABC chunk 2/17/55/58 test splits are not shipped, so we take shapes
  from the `nice/` files not used in phases 0–3 if any, otherwise state the overlap.
- 45 Thingi10K closed, manifold, oriented, non-self-intersecting meshes from
  `ddg-tools/data/thingi10k/` (the selection script records ids).
- 70 analytic shapes (§11).
- The reference's `data/` clouds.

**Protocol:**
- Meshes are centered (vertex mean, as `center_and_scale`) and scaled to [−1, 1]³.
- 512 points, uniform sampling.
- 64³ grid over [−1, 1]³, mean absolute error against ground truth (§9).

**Methods:** PAT (accelerated and all), SHC, SSPD, GWN-sign + distance, regularized WN, SPSR mesh
SDF, the authors' program, our retrained weights.

**Figures:** per-dataset bars (Fig. 18) and slices (Fig. 17).

### 16.3 Choices by measurement

- `LambdaRule`: reference vs paper heuristic vs `All`.
- `C`: 32, 64, 128 in f64.
- `k` evaluation: 16, 32, 64.
- Orientation `ε`.
- Tracer `safety`.

Each on the §16.2 subset (ABC + blobs, 50 each).

### 16.4 Density (Fig. 7, Fig. 29)

- One shape at 1024, 8192, 65,536 points, 128³ grid error (Fig. 7).
- Fig. 29: increasing density (110k → 29M) of one scan (David is not available; use the densest
  Thingi10K or a Keenan Crane model, recorded). Show the interior degradation, with SDF slices.

### 16.5 Performance (the floor)

Measured with no outside load, or reported with the measured load:
- precompute (kd-tree + network + fit) for 100k, 1M, 7M points on CPU (64 threads) and GPU;
- query time per point for 32-NN, sequential and batched;
- MC at 256³ and 1024³ narrow band;
- tracer at 1080p.

**Floor:** the paper's 40 s per 1M points (3090 + 3990X) and 1e-4–1e-3 s per query (sequential).
Targets: ≤ 20 s per 1M points on the 4080; query ≤ 5e-5 s sequential.

### 16.6 Training

Wall time per phase, loss curves against the reference's reported behavior, and a comparison of
retrained vs authors' weights on §16.2.

### 16.7 Black-box comparison with the authors' program

The other agent builds the reference and writes golden outputs to `runs_tmp/pat-ref/`
(`cpu_run1`, `cpu_run2`, …). WP14 reads them with a manifest:
- the inputs: cloud files and query grids;
- coefficients, tori and SDF values;
- the device and JAX precision mode.

Compare at the tolerances of §15.1-4.

## 17. Work packages

**General rules.**
- Each package is one agent in the main checkout on branch `points-as-tori`, with its own
  `CARGO_TARGET_DIR=target-pat<N>`, committing only its files.
- WP1 creates stubs of every module so parallel packages never edit the same file. `src/lib.rs` is
  edited only by WP1 and the orchestrator.
- Briefs repeat the standing rules (CLAUDE.md, parallel independent runs, Lean workspace rules,
  commit often, a notes file `NOTES-WP<N>.md`).

| WP | content | files | depends on | model |
|---|---|---|---|---|
| WP0 | toolchain bump to 1.95.0 (§18.1): `rust-toolchain.toml`, full workspace build/test/clippy, Python and Blender suites, minimal fixes for new lints | `rust-toolchain.toml`, lint fixes only | — | Sonnet |
| WP1 | **core**: deps (burn pinned, safetensors), `PatError`, `PointCloud` + validation, kd-tree (kNN/radius, deterministic, parallel), weights converter + committed safetensors + loader + `weights/README.md`, features (both orders), CPU f32 forward via burn-flex (chunked workers, §5.4), fitting with clamps and report, torus SDF/gradient, blend (3 λ rules, 3 neighborhoods, frozen/exact gradient), `PointsAsTori` facade, `tools/jax_oracle.py`, parity tests on the `data/` clouds; re-measure §5.3 in the crate; module stubs | `Cargo.toml`, `src/{lib,error,cloud,features,network/{mod,burn},fit,sdf,blend,field}.rs`, `src/spatial/*`, `tools/*`, `weights/*`, `tests/{parity,fit,sdf,determinism}.rs`, stubs | — | Opus |
| WP2 | Lean proofs + CITED (§13) | `lean/proofs/DdgProofs/PointsAsTori/*.lean`, the import line, `CITED.md` section | SPEC | Opus |
| WP3 | Lean oracle + fixtures | `lean/oracle/DdgOracle/PointsAsTori*.lean`, dispatch lines, `lean/oracle/fixtures/pat/`, oracle README row | SPEC (Rust replay test in WP15) | Sonnet |
| WP4 | CUDA: burn-cuda inference (chunks, tolerance, level-1 test), our NVRTC kernels for kNN + blend + gradient (f64, `F32Stable`), fitting on GPU, determinism tests, measurements | `src/gpu/*`, `src/network/cuda.rs`, `tests/gpu.rs`, `NOTES-GPU.md` | WP1 | Opus |
| WP5 | baselines: point areas (Barill, own planar Delaunay), GWN, FWN (Barill), regularized WN, SHC, SSPD, LSE; mesh ground truth (closest point BVH + GWN); SPSR via Open3D/uv script | `src/baselines/*`, `src/truth.rs`, `scripts/pat_spsr.py` (new file), `tests/baselines.rs` | WP1 | Opus |
| WP6 | marching cubes (asymptotic decider, narrow band), `ImplicitField`, CSG, offsets, level-set sampling, morphology | `src/{mc,csg,apps/mod,apps/morph}.rs`, `tests/{mc,apps}.rs` | WP1 | Sonnet |
| WP7 | training core: npz reader, data validation, query sampling (code/paper boxes), ground truth on CPU, loss as framework ops (fit + blend + closed-form gradient), two-pass memory plan, Adam + schedules, curriculum file, CPU training loop, FD checks (f64), short CPU run on the blob | `src/train/*`, `training/curriculum.toml`, `tests/train.rs` | WP1 | Opus |
| WP8 | analytic shapes (16 SDFs, revolve/extrude, gradients, area-uniform sampling, 70-instance set), blob generator (Perlin noise; own implementation), simulated scanner (cameras, grid rays vs triangle BVH), FPS | `src/analytic/*`, `src/generate/*`, `tests/{analytic,generate}.rs` | WP1 | Sonnet |
| WP9 | normal orientation; 2D variant per D15 | `src/apps/orient.rs`, `src/pat2d/*`, `tests/{orient,pat2d}.rs` | WP1, WP6 (sampling), D15 | Sonnet (orientation) / Opus (2D if trained) |
| WP10 | training on the GPU: burn-cuda autodiff, GPU ground truth and blend kernels if needed, the full curriculum run (long, background, logged), `TRAINING.md`, retrained weights in `data/` | `src/train/gpu.rs`, `TRAINING.md`, `runs_tmp/pat-train/` | WP4, WP7 | Opus |
| WP11 | renderer: CPU + CUDA sphere tracer, shading as the demo, cube map assets (with notice), cut-plane slices, PNG output, showcase images (Fig. 12 analogues) | `src/render/*`, `assets/*`, `examples/render.rs`, `showcase/points-as-tori/` | WP4, WP6 | Sonnet |
| WP12 | Python bindings | `crates/ddg-py/src/points_as_tori.rs`, registration lines, `.pyi`, README section, `python-tests/test_points_as_tori.py` | WP1, WP6 (WP4 for `device`) | Sonnet |
| WP13 | Blender operator (§12.2) | `blender/ddg_uv_unwrap/points_as_tori_{operator,io}.py`, registration lines, `blender/tests/*points_as_tori*.py`, README section | WP12, WP9 | Sonnet |
| WP14 | evaluation harness: datasets (ABC/blob validation, Thingi10K selection, analytic 70, `data/` clouds), 64³ errors for all methods (§16.2), Fig. 7 and 29 studies, black-box comparison (§16.7), `COMPARISON.md`, `EVALUATION.md` | `examples/{eval,density}.rs`, `scripts/pat_eval*.{sh,py}`, `COMPARISON.md`, `EVALUATION.md` | WP5, WP6, WP8 (and WP10 for the retrained row) | Sonnet |
| WP15 | oracle replay test, determinism audit, `check-formulas` pass | `tests/oracle.rs`, audit entries | WP1, WP3 | Sonnet |
| WP16 | benches (§16.5), timed runs alone, §16.3 sweeps that time things | `benches/*`, `[[bench]]` entries, results in this file | WP4, WP6 | Sonnet |

**Rounds** (at most 3 agents each; timed work alone at the end):

1. WP1, WP2, WP3.
2. WP4, WP5, WP7.
3. WP6, WP8, WP10. WP10 starts the long training run in the background as soon as WP4 + WP7 are
   in. The run continues through rounds 4–5, so the GPU is shared: other rounds avoid GPU timings
   while it runs.
4. WP9, WP11, WP12.
5. WP13, WP14, WP15.
6. Alone: WP16, then the §16.3/§16.5 timed parts of WP14 after the training run ends.

## 18. The user's decisions (2026-11-04) and open questions

**Decided by the user:**
- Training in Rust, retraining, comparison with the authors' weights.
- All applications.
- The baselines listed in §9.
- The datasets in §16.2.
- The Blender operator of §12.2.
- Standing rules.
- Framework for the network (change of 2026-11-04).

**Decided by the user on the draft (2026-11-04, second round):**
1. **D2 / toolchain: bump.** The workspace toolchain is pinned to Rust **1.95.0**
   (`rust-toolchain.toml`, WP0) and the crate uses **burn 0.22** (pinned `=0.22.x`). This
   supersedes "burn 0.21" wherever this file still says so; the 0.21 measurements of §5.3 are
   re-taken on 0.22 in WP1.
2. **D4: accepted.** Framework part on the GPU: bitwise run to run on the same device (tested),
   documented tolerance against the CPU. If burn-cuda is not repeatable run to run (atomics), GPU
   inference switches to cuBLAS through candle. Our own kernels keep the workspace rule.
3. **CPU inference speed: own kernel.** If WP1's tuning of burn stays below ~1,000
   neighborhoods/s on the CPU, WP1b writes a CPU-only inference path on the `gemm` crate (fused
   LayerNorm/GELU, reused buffers, fixed chunks, bitwise thread-independent), cross-checked
   against the burn path to ~1e-6. Training stays on burn.
4. **D15: both.** First the classical per-point circle fit (baseline), then a trained 2D network
   (same architecture on 2D, generated 2D training data: Quilez 2D SDFs, polygons, noise blobs,
   simulated 2D scans), compared against the fit. WP9 is split: WP9a orientation + classical 2D
   (Sonnet), WP9b the 2D network and its training (Opus, after WP7/WP10).

**Orchestrator defaults (the user did not pick these; they may override):**
5. **D27.** SHM and NN-VIPSS are not planned (not in the user's baseline list).
6. Retrained weights stay in `data/` until they match or beat the authors' on §16.2.

**Resolutions after the Lean oracle (WP3, 2026-11-04; details in `NOTES-WP3.md`,
format in `lean/oracle/README.md` "Points as Tori").** These override §7.4 where they differ.
1. Exponent shift: with the constant-λ rules (`PaperHeuristic` λ = 10³/D, fixed λ), the shift
   is `min ρ` over the neighborhood (λσ_x ≈ 1300 overflows `exp` in the 32-NN fallback). The
   `PerQuery` rule keeps σ_x = ½ max ρ and λ = C/σ_x.
2. A neighborhood with σ_x = 0 under `PerQuery` (a radius set containing only the query), and an
   SHC sum of exactly 0, return `DegenerateQuery` naming the query. Never NaN or inf.
3. The Laplacian at axis/core singular points of a torus is 0 and is counted. `pat-blend-laplacian`
   is the divergence of App. D, written out in the oracle README and checked by finite differences.
4. Fit counters: one per clamp (disc, eigenvector, denominators `ev·I·ev`, |κ| floor, and the
   remaining clamp as listed in NOTES-WP3). "Umbilic" means the disc clamp is active.
5. `sgn(0)` follows the C++ `sgn+` convention (it matters only at exact umbilics; JAX gives 0).
6. `radius auto` means R_eval = 2C/λ. The definition of D for the paper rule is the exact mean
   64-NN distance (the reference's C++ helper has a data race; `runs_tmp/pat-ref/README.md`).
7. The oracle agrees with the authors' program on all 34,222 tori of the golden set
   (relative ≤ 6e-15; axes 4e-10 from their angle round trip) and on the SDF at 1000 queries
   (≤ 1.3e-14). The CPU/GPU differences on near-flat tori are within the authors' program
   (CPU vs GPU), not between the formulas.

**Resolutions after the core (WP1, 2026-11-04; details in `NOTES-WP1.md` "SPEC change
proposals"; all accepted by the orchestrator).** These override §5, §7 and §15.1 where they differ.
1. **Saddle tori (new finding).** For κmin·κmax < 0, the paper's center (Eq. 22 sign with
   `c = Q* − sign/|κmin| n*`) gives the mirror image of the osculating torus: ∂g/∂n = −1 at the
   point. On the golden tori that is every saddle (holey-bunny 33/12,917, femur 4, snail 1).
   Switching those changes the 64³ SDF by up to 0.36. `FitOptions::sign_rule`:
   `SignRule::Reference` (default for parity: as trained, as the oracle) and
   `SignRule::Geometric` (center `o + n*/κmin`, signed; it differs only at saddles and exact
   zeros). WP14 picks the default by the §16.2 accuracy run. WP7/WP10 train with either rule
   (option), and the retrained comparison includes a `Geometric`-trained run. To be proved in
   Lean (WP2): the contact condition (∂g/∂n = +1 at the point) under `Geometric`, and its failure
   at saddles under `Reference`.
2. A plane gives a half-space (§15.1) only under `Geometric`. Under `Reference` an exact zero
   curvature gives a horn torus whose core passes through the point; this is counted
   (`zero_curvature`).
3. `CHUNK = 1` (4,000–4,650 neighborhoods/s at 32–64 threads; bitwise independent of chunk size
   and thread count). No `gemm` kernel: the §18.3 threshold is met.
4. Single-point evaluation returns `Result` (`DegenerateQuery`, non-finite query).
5. `PaperHeuristic` sets λ only; the radius is `Neighborhood::Radius { r, fallback_k }`,
   built by `TorusField::paper_options`. D = mean over points of the mean distance to the
   `k_mean` nearest, self included.
6. `All` and kNN with k ≥ n sum in index order (as the reference). The oracle sums in
   `(ρ², index)` order; WP15 compares at 1e-12.
7. Features normalize normals (`n / max(‖n‖, 1e-10)`) in both orders. Report fields as in
   NOTES-WP1 item 9.

**Resolutions after the Lean proofs (WP2, 2026-11-04; `NOTES-WP2.md`).**
1. New id `pat-fit-contact` (§13) proves the saddle finding of WP1 resolution 1.
2. **Degenerate eigenvector.** The reference pair `ev± = (v1 ± disc, v2)` is the zero vector for
   one sign when `v2 = 0` (off umbilics, e.g. `a10 = a01 = a11 = 0`, `a20 ≠ a02`). The 1e-10
   clamp then returns the other eigenvector, so both curvatures are equal. Near `v2 = 0`,
   `v1 − disc` cancels: relative error ≈ `ε·v1²/v2²`. In the f64 CPU fit this does not occur on
   the golden data: all 34,222 points match the closed form `H ± √(H²−K)` to 1.5e-13, and the
   minimum is `|v2/v1|` = 6e-6. In f32 it starts near `|v2/v1|` ≈ 3e-4, which affects the
   training loss. Hence `FitOptions::eigen_rule: Reference | Robust` (default `Reference`; switch
   counted). The training fit (f32) uses the robust form (WP7).
3. The weight range claim of §7.4 is replaced by: every weight lies in `[e^{−C}, e^{C}]` and
   `W ≥ |𝒩| e^{−C} > 0` (proved). "Max weight ≥ 1 for the nearest point" holds only when
   `ρ_min ≤ σ`.
4. `pat-output-scaling` (scale covariance) is exact only while no absolute 1e-10 clamp is active.
5. `pat-torus-sdf` is the exact distance for `R ≥ |r|`; spindle tori get the lower bound only.
   `pat-fwn-expansion` proves the gradient and Hessian of `1/(4π|y|)`; the order-2 term and the
   error bound are cited.
6. The ids `pat-shc-gradient`, `pat-sspd-gradient`, `pat-extrusion`, `pat-extrusion-gradient`,
   `pat2d-circle-sdf`, `pat2d-blend` each have their own Lean theorem. `pat2d-blend` covers shift
   and convexity only; there is no 2D gradient proof.

**Resolutions after the baselines (WP5, 2026-11-04; `NOTES-WP5.md`; proposals accepted).**
1. Point areas: the Voronoi cell comes from clipping half-planes in the tangent plane (same area
   as the planar Delaunay route, no robust predicates). Cells on the hull are clipped to the k-NN
   radius and counted.
2. New cited ids `pat-truth-closest-point`, `pat-solid-angle` (CITED.md, WP15 audit).
3. **"GWN-sign + distance" is defined** as unsigned distance to the nearest cloud point, signed
   by the (direct, fast or regularized) winding number at threshold ½. It is a baseline of ours,
   not in the paper. On the sanity run it is the strongest point method (512 points, 9 clouds,
   64³ mean abs. error: PAT 0.0183, GWN-sign 0.0169, LSE + GWN-sign 0.0165, SHC 0.0247,
   SSPD 0.0576, SPSR 0.0220 on 7 meshed clouds). WP14 reports it in every table next to the
   paper's methods, and also reports near-surface error (narrow band) and eikonality
   (‖∇φ‖ − 1), where nearest-point distance is weak, so that the grid mean is not the only
   criterion.
4. SHC with `w = 0` stays an error by default. The opt-in `shc_zero_value` (counted) exists
   because CAD grids put up to 167 grid points exactly in a face plane; WP14 uses it with value 0
   and reports the count.
5. FWN β: default 2 (SPEC); Barill's point-cloud value 2.3 is the option used in WP14's timing row.
6. `SignRule` default is still open: the network predicted saddles only on femur (4/3897;
   Geometric marginally better). WP14 adds clouds with saddles and meshes (e.g. Thingi10K
   samples, analytic saddle shapes from WP8).
