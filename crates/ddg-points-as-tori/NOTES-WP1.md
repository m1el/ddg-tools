# WP1 notes: core of the Points as Tori crate

Branch `points-as-tori`, build dir `target-pat1`, logs in `runs_tmp/pat-wp1/`. Toolchain 1.95.0, burn 0.22.1 (`=0.22.1`).

## What exists

| file | content |
|---|---|
| `Cargo.toml` | burn `=0.22.1` (`std`, `flex`, `simd`; **no** `rayon`), safetensors 0.6, rayon; features `bundled-weights` (default), `cuda`, `train` (empty, for WP4/WP7) |
| `src/error.rs` | `PatError` (adds `SizeMismatch`, `ThreadPool`, `Format`, `Backend`, `DegenerateQuery` to the SPEC list) |
| `src/cloud.rs` | `PointCloud`, `NormalPolicy`, `CloudReport`, `hash()` |
| `src/pool.rs` | `make_pool(threads)`: the explicit rayon pool of every stage |
| `src/spatial/{mod,kdtree}.rs` | kd-tree (median of the widest axis, bucket 16, tight node boxes, pre-order layout fixed by `n`), exact kNN and radius by `(d², index)`, `fitting_neighborhoods` (self first), `knn_batch`; unit tests against brute force with ties and duplicates, thread-independent construction |
| `src/features.rs` | frame, median, features in both orders (`FeatureOrder`), `neighborhood_features_f64` (the default order before rounding, for the oracle replay) |
| `src/network/{mod,burn}.rs` | weights loader (safetensors → `NetWeights`), `Model::reference()` (`include_bytes!`), `PatNet` burn module, `predict_coefficients`, `forward_features`, `scale_outputs{,_f64}` |
| `src/fit.rs` | fitting with all clamps counted, `SignRule`, `FitReport`, `fit_diagnostics`, `Tori::{to_bytes,from_bytes}` |
| `src/sdf.rs` | torus SDF, gradient, Laplacian, singular rules |
| `src/blend.rs` | blend, frozen/exact gradient, Laplacian, λ rules, shift rule |
| `src/field.rs` | `TorusField` (single and batch evaluation), `PointsAsTori::build`, `EvalOptions`, `BuildOptions`, `TorusField::paper_options` |
| `src/npy.rs` | `.npy` / uncompressed `.npz` reader (tests, tools; WP7 may reuse) |
| stubs | `analytic`, `apps/{mod,morph,orient}`, `baselines`, `csg`, `generate`, `gpu` (cfg cuda), `mc`, `network/cuda.rs` (cfg cuda), `pat2d`, `render`, `train/{mod,gpu}` (cfg train / cuda), `truth` |
| `tools/convert_weights.py` | restricted unpickler (the 7 globals), checks, safetensors + `.sha256` manifest |
| `tools/jax_oracle.py` | Flax re-declaration on our safetensors; reproduces the golden raw outputs, σ and coefficients **bitwise** (tooth) |
| `weights/` | `pat-reference.safetensors` (6,366,776 B, sha256 `fce4101a…`), its manifest, README; equal (max Δ = 0) to the parameters the authors' loader produces |
| `examples/throughput.rs` | CPU network throughput (§5.3) |
| `tests/{parity,fit,sdf,determinism}.rs`, `tests/common/` | see below |

## Parity with the authors' program (golden `runs_tmp/pat-ref/cpu_run1`, 12 inputs) — `tests/parity.rs`

`cargo test --release -p ddg-points-as-tori --test parity -- --nocapture --test-threads 1` (45 s; log `parity2.log`).

| stage | result (max over the 12 inputs) |
|---|---|
| kNN | identical rows except 16 rows with the same set in a different order (equal distances); no set differs |
| σ (reference order, f32) | relative ≤ 2.3e-7 (≤ 2 ulp) |
| raw network outputs (golden kNN, reference order) | ≤ 6.8e-6 abs (femur; others ≤ 2.2e-6) |
| coefficients `|Δa|/(1+|a|)` | ≤ 1.4e-5 (reference order), ≤ 2.0e-5 (default order, our kNN) |
| tori from the golden coefficients | centers, R, r ≤ 5.3e-15 relative; axis ≤ 4.3e-10 (their angle round trip) |
| SDF from the golden tori (32-NN and all points) | grid ≤ 1.8e-12, gradient ≤ 4.9e-10 (holey-bunny), random queries ≤ 5e-15 / 2e-13 |
| end to end (our kNN, network, fit, blend) vs golden grid SDF | median ≤ 8e-9, p99.9 ≤ 1.1e-7, max ≤ 2.6e-7; gradient max ≤ 6.4e-5 (both feature orders) |

Asserted tolerances: raw and coefficients ≤ 2e-5, tori ≤ 1e-6 relative, SDF ≤ 1e-9 and gradient ≤ 1e-7 given the tori, end-to-end
median ≤ 1e-6 and p99.9 ≤ 1e-4.

## CPU network throughput (§5.3 re-measured in the crate, burn 0.22.1 Flex)

holey-bunny, 12,917 neighborhoods, `examples/throughput`, outside load from `uptime` before each run (the load is mostly the run itself;
the machine was otherwise idle at 01:15). Log `tp-final.log`.

| threads | neighborhoods/s | load (1 min) |
|---|---|---|
| 1 | 147 | 11 |
| 8 | 1,105 | 13 |
| 32 | 3,960–4,010 | 4 |
| 64 | 4,580–4,680 | 9 |

Reference: JAX CPU ~1,050/s (survey, load ~27), the golden timing 12,917 in 9.56 s = 1,350/s (load 20–30).

**Chunk size is the tuning knob** (32 threads, 8,192 neighborhoods, load 30–35, `tp-burn-chunks*.log`):

| chunk | 1 | 2 | 4 | 8 | 16 | 32 | 64 | 128 | 256 |
|---|---|---|---|---|---|---|---|---|---|
| /s | 4,050 | 3,300 | 2,700–3,000 | 2,600–2,700 | 1,700–2,100 | 960 | 685 | 560 | 570 |

The activations of a chunk of 128 (128·64 rows × 512 floats = 16 MB for the FFN) leave the cache; at 1–8 they stay in L2/L3 and the
cores scale (147/s × 32 ≈ 4,700). **The outputs were bitwise identical for every chunk size from 1 to 512 and every thread count
from 1 to 64** (hash `c267515533d3391b` on holey-bunny; `fa03bf81fb79caec` on its first 8,192). So `CHUNK = 1`.

**Decision: no gemm kernel.** The tuned burn path is above the ~1,000/s threshold of SPEC §18.3 at 8 threads and 4× JAX at 32,
so the CPU-only kernel (WP1b) was not written; `CpuBackend` was removed. 147/s per core is ~32 GFLOPS of ~118 peak per core; a
hand kernel could gain perhaps 2–3× more if it is ever wanted.

burn 0.22 API facts (for WP4/WP7): `Tensor<D, K>` has no backend parameter (runtime dispatch); the CPU device is `Device::flex()`;
`burn-ndarray` is deprecated; `burn::nn::attention::MultiHeadAttention` takes `MhaInput::self_attn(x).mask_pad(m)` and divides the
scores by `√d_k` after the product (= 4 here, exact); `LayerNormConfig::with_epsilon(1e-6)`, computed by the backend op
`layer_norm` (fused in Flex, two-pass variance); `activation::gelu_approximate` is composed (`powf`, `tanh`), not fused;
`#[derive(Module)]` structs need no generic; `Param::from_tensor`; `TensorData::try_into_vec`. The module is `Sync` and is
shared by the workers.

## Determinism (`tests/determinism.rs`)

Whole pipeline (kNN, coefficients, σ, tori, SDF, gradient, Laplacian; both feature orders, 32-NN and radius) bitwise equal for
two runs on 1 thread and for 3, 8, 32 threads (2,500 points, 91 s in release). kd-tree construction thread-independent (unit test).

## Other tests

- `tests/fit.rs`: tangent sphere, torus at outer and inner equator (exact torus back), plane, cylinder, clamps and counters,
  errors, cloud validation, translation/scale covariance of features and coefficients (network included), tori round trip.
- `tests/sdf.rs`: `∇g` and `Δg` against differences (400 tori incl. spindles), singular rules, blend gradients against
  differences (constant λ frozen = exact; PerQuery exact), `Δφ` against Richardson second differences, convexity, identical tori,
  `DegenerateQuery`, no overflow with large constant λ, data-point queries.
- clippy `--all-targets -D warnings` clean; rustdoc clean.

## Findings that need a SPEC decision (proposed text; the orchestrator edits SPEC.md)

1. **The reference's tori are inside-out at saddles** (new; §7.1, D9, Eq. 22). For `κmin κmax < 0` the paper's center
   `c = Q*(0,0) − (sign(T)/|κmin|) n*` puts the center on the wrong side: the torus passes through `o` but is the mirror image of
   the osculating one, and `∂g/∂n = −1` at `o`. The center of the parallel circle is `o + n*/κmin` (signed); it equals the paper's
   formula exactly when `sign(T) = −sgn(κmin)`, i.e. away from saddles. Measured on the golden tori (numpy, scratch script
   `contact.py`): `∂g/∂n` at `o` is +0.5…+1.0 for every non-saddle torus and −0.94…−1.0 for every saddle torus (holey-bunny 33 of
   12,917, femur 4, snail 1). Switching only these tori changes the 64³ grid SDF by up to 0.36 (holey-bunny), 0.21 (snail), 0.12
   (femur); median 0. This is probably also the "isolated cells, 0.195 on holey-bunny, not investigated" of the golden CPU-vs-GPU
   comparison: a saddle flips with a tiny change of the coefficients. Implemented as `FitOptions::sign_rule`:
   `SignRule::Reference` (default, as trained and as the oracle) and `SignRule::Geometric` (ours; bitwise equal to `Reference`
   except at saddles and exact zeros). **Proposal:** keep `Reference` as the default for parity, add `Geometric` to §7.1 and D9,
   and decide the default by the §16.2 accuracy measurement (WP14, needs WP5's ground truth). Also relevant to WP7/WP10: the
   network was trained with the mirrored saddle tori; retraining with `Geometric` is a choice.
2. **Exact zero curvature** (§7.1, §15.1-1). With the JAX rule an exactly zero `κ` gives `s = −1`, `sgn = 0`: a plane gives
   `R = r = 1e10`, `c = o + 1e10 n*`, a horn torus whose *core* circle passes through `o` (`g(o) = 1e10`), not a half-space; a
   cylinder bending away from `n` is mirrored as above. SPEC §15.1 "a plane gives a half-space up to 1e-10" holds only under
   `Geometric` (which floors `κ` to `±floor`, `+floor` at 0). Tiny non-zero `κ` agree in both rules unless the round-off sign
   makes a saddle. Proposal: say so in §7.1 and §15.1.
3. **`CHUNK = 1`, not 128** (§5.4): measured above; the outputs do not depend on the chunk size.
4. **No `CpuBackend`/gemm kernel** (§18.3): threshold met by tuning.
5. **Single-point evaluation returns `Result`** (§7.4 API): `distance`, `distance_gradient`, `laplacian` return
   `Result<_, PatError>` for `DegenerateQuery` and non-finite queries (WP3 resolution 2: an error, not NaN). `blend_at` added.
6. **`PaperHeuristic` and its radius** (§7.4): the rule sets λ only; the evaluation radius is `Neighborhood::Radius { r, fallback_k }`.
   `TorusField::paper_options(cloud)` builds the paper's configuration (`R_eval = 2C/λ`, 32-NN fallback). `D` = mean over points of
   the mean distance to the `k_mean` nearest, self included (as the oracle).
7. **Sum order for `All` and `Knn { k ≥ n }`**: index order (as the reference, golden SDF to 1e-15); the oracle README sums over the
   `(ρ², index)`-sorted set. The difference is round-off (≤ 1e-15 relative); WP15 compares at 1e-12. Radius sets and kNN are sorted
   in both. The exact gradient's `r̂` is that of the member with the largest `(ρ², index)` in both (WP3 resolution 6).
8. **Features** normalize the normals (`n / max(‖n‖, 1e-10)`) in both orders (the SPEC says `n̂ = n` for the default; validated
   normals are unit within 1e-6, so this only follows the reference more closely). Default-order coefficients are scaled in f64 and
   then rounded to f32 (`scale_outputs_f64` exposes the f64 values).
9. **Report fields** beyond the SPEC: `FitReport { clamped_denominator, saddle, zero_curvature, seconds }` (`umbilic` = disc clamp, WP3
   resolution 4), `FitDiagnostics.clamps` has one flag per clamp; `EvalReport { singular_axis, singular_core, at_data_point,
   radius_fallbacks, terms }`; `BlendValue { sigma, shift }`; `PredictReport { features_seconds, forward_seconds }`.
10. **Laplacian at singular points** is 0 (WP3 resolution 3), replacing the SPEC §7.2 text ("one-sided limit" applies to the
    gradient only).
11. `pat-network` (cited in §13) is not yet in `lean/proofs/CITED.md`: `scripts/check-formulas.sh` reports it MISSING (WP2).
12. The crate's `rust-version` is `1.95` (burn 0.22's MSRV); the workspace still says 1.90.

## Not done / open

- GPU (`Device::Cuda`) returns `PatError::Gpu` (WP4). `EvalPrecision::F32Stable` likewise.
- `Tori` export to the reference's pickle layout (§2.2) is left to the Python binding (WP12); `Torus::axis_angles` exists.
- The second pkl (`FundamentalFormPredictor_pretrained.pkl`) is not converted yet (same tool; output to `data/points-as-tori/`).
- Facade timings (`BuildReport` has them) not tabulated; WP16.

## How to resume

```
export CARGO_TARGET_DIR=/var/data/discrete-geometry/ddg-tools/target-pat1
cargo test --release -p ddg-points-as-tori                      # all tests (parity skips without runs_tmp/pat-ref)
cargo test --release -p ddg-points-as-tori --test parity -- --nocapture --test-threads 1
cargo run --release -p ddg-points-as-tori --example throughput -- runs_tmp/pat-ref/inputs/holey-bunny/{points,normals}.npy 32 3
```
