# WP7 notes: training core

Branch `points-as-tori`, build dir `target-pat7`, logs in `runs_tmp/pat-wp7/`. burn 0.22.1 with
`autodiff` + `optim` behind the crate feature `train`. Everything runs on a burn `Device` passed by
the caller (`train::Ctx { dev, dtype }`): `Device::flex()` here; WP10 passes a CUDA device. Autodiff is
enabled per leaf tensor (`t.autodiff().require_grad()`), so plain and autodiff devices need no
separate code paths. f32 for training, f64 (`DType::F64` on Flex) for the finite-difference checks.

## What exists

| file | content |
|---|---|
| `src/train/npz.rs` | typed npz reader (deflate via `zip` 0.6, pure Rust; `<f8 <f4 <i4 <i8 |b1`, C order), member-by-member, errors name file + member; `write_npz` (stored) for fixtures and the JAX export |
| `src/train/data.rs` | `load_file`: shape agreement, mask and face-mask prefixes, neighbor column 0 = self, indices in `0..n` or a `−1` suffix, faces in range, finite positions/normals; normals normalized and counted (D23) unless `raw_normals`; `DataReport` |
| `src/train/truth.rs` | ground truth: BVH closest point (Ericson; degenerate triangles via edges), exact winding-number sign summed over all faces (`> π` ⇒ inside, the reference's test); `GwnPrecision::ReferenceF32` repeats the reference's f32 sign test |
| `src/train/queries.rs` | the reference's per-row sampling (40 zero set ∝ area in the regularized box, 40 narrow band, 40 cube), `QueryBoxes::{Code, Paper}` (D20), ChaCha8 per (seed, epoch, shape) key with the row as stream (D22) |
| `src/train/tfit.rs` | the torus fit as burn ops (both `SignRule`s, both `EigenRule`s, JAX clamps), any dtype/device |
| `src/train/loss.rs` | blend constants on the host (`build_plan`: σ from the 64 nearest valid tori, λ = 64/σ, ŵ = w/max(W, 1e-10), c = λ(Σŵ r̂ − r̂)), truncation `T`, the torus SDF + closed-form gradient + blend + Eq. 27 as burn ops, autodiff `∂L/∂a`; `host_blend` (f64 check through `crate::sdf`) |
| `src/train/net.rs` | parameters in the Flax order/names, WP1's `PatNet` built per pass (leaves tracked or not), lecun-normal init with the tangent-sphere head, safetensors writer in the weights layout (readable by `Model::from_safetensors`), `forward` with the pad mask |
| `src/train/step.rs` | the step in passes (plan, coefficients without tape, `∂L/∂a`, taped network chunks in fixed groups), padding modes (D19), `evaluate_shape` |
| `src/train/optim.rs` | burn's `Adam` per tensor (state keyed by parameter name), decoupled AdamW decay, the optax schedule (`Schedule::Code` default, `Cosine` option, D21) |
| `src/train/run.rs` | curriculum parser, file patterns (sorted), per-epoch buffered shuffle, the loop (`train`, `train_on`), JSON-lines log, checkpoints (`checkpoint.safetensors` + `checkpoint-adam.safetensors`, cursor in the metadata), resume, phase outputs, validation |
| `training/curriculum.toml` | the reference's six phases including pretraining (step counts 4,096 / 96,000 / 64,000 / 119,344 / 124,590 / 89,828 = 497,858, checked against the data in `tests/train.rs`) |
| `tests/train.rs` | 16 tests (below) + 6 ignored real-data runs |

**Edits to WP1's files** (minimal):
- `src/fit.rs`: `EigenRule { Reference (default), Robust }` in `FitOptions::eigen_rule`, `robust_pair`,
  `ClampFlags::eigen_switched`, `FitReport::eigen_switched` (all existing literals use `..Default`).
- `Cargo.toml`: feature `train = ["burn/autodiff", "burn/optim", "dep:zip", "dep:rand_chacha", "dep:toml"]`
  and the three optional dependency lines (all three were already in `Cargo.lock`; only the package's
  dependency list changed). Committed without WP4's uncommitted `cuda` lines (index built from HEAD +
  my lines).
- No change to `src/lib.rs` (the `train` module line was already there).

## Design

**Only `g_i` and `∇g_i` depend on the network.** The blend weights depend on the query and the
points only, so σ, λ, ŵ and `c_i = λ(Σ_j ŵ_j r̂_j − r̂_i)` are host constants (f64), and
`φ = Σ ŵ g`, `∇φ = Σ ŵ (∇g + c g)` is the reference's `A + φB` regrouped. The framework part per pair
is the fit (per torus), a gather, the torus SDF and its closed-form gradient (§5.2: first-order
autodiff only).

**Truncation.** Pairs with `λ(ρ_i − ρ_min) > T` (T = 40) are dropped: their weight is below
`e^{−40} = 4e−18` of the largest. The plan then uses the kd-tree (64-NN for σ, a radius query for the
kept set); the dropped mass is bounded and reported. On the blob at batch 1024: 3.66M pairs instead
of 251M dense. Check: the loss at T = 40 equals the untruncated (all-pairs, dense) loss to < 1e-12
(f64, `tensor_blend_matches_host_blend`).

**Chunks.** Queries are sorted by pair count and cut into chunks of ≤ 65,536 padded pairs, each on its
own autodiff graph; padding a chunk to its largest query was the first version's bottleneck (14 s →
0.37 s per step).

**Two passes** (§10.4) as planned, plus: rows whose `∂L/∂a` is exactly zero are not back-propagated
(exact), and the tori the plan never references are not even forwarded (exact; only matters when the
batch's queries cover part of the cloud).

**Determinism.** Taped chunks are grouped into a fixed number of groups (64) summed in order: the f32
step (loss and all gradients) is bitwise equal on 1, 4 and 13 threads (`f32_step_is_bitwise_independent_of_threads`).

## Tests (`cargo test --release -p ddg-points-as-tori --features train --test train`, 1.2 s)

| test | result |
|---|---|
| npz round trip, errors | ok |
| data validation (prefix masks, self column, norm-2 normals counted and normalized, `raw_normals`) | ok |
| ground truth vs brute force and the sphere | distance = brute force to 1e-14; sign right everywhere; f32 sign test agrees off the surface |
| tensor fit vs `fit::fit_torus` (f64, 200 random rows, 2 sign × 2 eigen rules) | ≤ 1e-11 relative, axes ≤ 1e-12 |
| tensor blend vs `host_blend` (f64) | loss terms ≤ 1e-11 |
| **FD ∂L/∂a** (f64, 132 entries, central h = 1e-6) | 2.1e-9 (Reference rules), 1.3e-9 (Geometric + Robust), relative to max·|g| |
| **FD ∇φ wrt x** (f64, frozen σ, 60 queries) | 9.0e-9 relative |
| **FD ∂L/∂θ** (reduced network: 2 layers, 2 heads, 16 tokens, E = 16, f64, 80 entries over all 40 tensors, through features → network → scaling → fit → blend → loss) | 5.9e-10 relative to max·|g| |
| robust rule continuous across `v1 = 0` (autodiff = FD on both sides) | ok |
| two passes independent of chunk size and groups (f64) | 1e-12 |
| f32 step bitwise thread-independent | ok |
| initial output = tangent sphere (raw output exactly `(0,0,0,0,−½,−½)`, `R = 0`, `r = σ`, center `p − σn`), 1,588,358 parameters, safetensors round trip | ok |
| Adam = optax's formula (two steps), schedules | ok |
| curriculum file and step counts against the data | ok |
| padded clouds (n < P, n < k, pad mask, reference padding) | ok |

## The saddle/eigenvector finding (orchestrator item, NOTES-WP2 item 1)

Measured on the f32 tensor fit with `a = (0, 0, 0, a11, −0.2, −0.5)` (`v2 = a11`, `v1 = −0.3`; ev₊
cancels), relative radius error against the f64 robust fit (`robust_eigen_rule`):

| a11 = v2 | 1e-1 | 1e-2 | 1e-3 | 1e-4 | 1e-5 | 1e-6 | 1e-7 | ≤ 1e-8 |
|---|---|---|---|---|---|---|---|---|
| Reference rule, f32 | 3e-8 | 1e-7 | 8e-8 | 1.0e-6 | 6.7e-5 | 6.6e-3 | 0.35 | 0.74–0.75 |
| Robust rule, f32 | 3e-8 | 1e-7 | 8e-8 | 7e-8 | 7e-10 | 7e-12 | 7e-14 | 0 |

So in f32 the reference form loses digits once `|v2/v1| ≲ 3e-4` and is wrong below `≈ 3e-6`
(somewhat later than the 1e-3 estimate). The robust rule (non-cancelling vector of the pair plus its
`I`-orthogonal complement `J·I·e`, exact for the symmetric pencil) has no clamp, no cancellation, and
is continuous with continuous derivatives at its switch `v1 = 0` (tested with autodiff = FD on both
sides). Available in both fits: `FitOptions::eigen_rule` (f64, counted in `FitReport::eigen_switched`)
and `TensorFitOptions::eigen_rule` (training). Effect on the loss with the authors' weights: below.

## The loss with the authors' weights vs the reference's (`reference_loss_export`, `runs_tmp/pat-wp7/refloss/`)

`ABC17_validation_3.npz` (10 shapes × 2 batches of 1024 rows × 120 queries), the bundled weights, our
queries and ground truth exported to npz; the reference's `SDF_loss` (copied verbatim into
`runs_tmp/pat-wp7/ref_loss.py`, gitignored, not part of the crate) with the authors' own network
(`FundamentalFormPredictor.load_saved_model` on their pkl), JAX 0.8.2 on the CPU, f32, dense over all
2048 tori. Ours: `FeatureOrder::Reference`, `raw_normals = true` (the reference fits with the stored
norm-2 normals), T = 40, f32.

| | distance | eikonal |
|---|---|---|
| max relative difference over the 20 batches | **7.1e-6** | **6.8e-6** |
| typical | 1e-7 – 4e-7 | 2e-7 – 1e-6 |

E.g. shape 0 batch 0: reference 0.8402562 / 19.35076, ours 0.8402562 / 19.35077. The difference is f32
round-off (theirs sums 2048 tori densely in f32, ours sums the kept pairs with f64 weights).

Ground truth (`runs_tmp/pat-wp7/gt_check.py`, 4 shapes × 30,720 queries against the reference's
`signed_distance_to_mesh`): unsigned distance ≤ 9e-8 apart (f32), **no sign difference off the
surface** (|d| > 1e-4); on the zero set the reference's f32 distance is often exactly 0 (sign 0).

**The reference validated with broken tori.** With the stored norm-2 normals its fit doubles the
center offset and `n*` (N§17.4), and the validation loss it prints is dominated by that:

| normals | mean distance | mean eikonal |
|---|---|---|
| raw (norm 2), as the reference | 0.988 | 10.58 |
| unit (ours, D23) | **0.0195** | **0.296** |

Its "validation" curve (and the pretraining phase, whose blob also has norm-2 normals) therefore
measured something else than the trained model's accuracy. Our default normalizes at load.

**Rule variants on the authors' weights** (max relative change of either term over the 20 batches):
`EigenRule::Robust` 8e-6 (unit normals), 4e-5 (raw normals); `FeatureOrder::DifferenceFirst` 1e-6;
`SignRule::Geometric` exactly 0, because the authors' network gives no saddle torus on these 10 clouds
(f64 fit: 0 saddles, 0 robust switches); `T = ∞` vs `T = 40` ≤ 4e-10. So training with the robust eigenvector rule does not
change the loss of the trained network measurably; it only removes the f32 failures near `v2 = 0`.

## CPU step time and memory (`step_timing`; Threadripper 3970X, 64 threads, authors' weights, blob)

| configuration | step | plan | coefficients (fwd, no tape) | loss + ∂L/∂a | taped network (pass 4) | peak RSS |
|---|---|---|---|---|---|---|
| first version, batch 1024 | 21 s | 0.25 | 0.5 | 14 (chunks padded to the largest query) | 6.6 | 19.5 GB |
| batch 1024 (phases 0–4) | **5.9 s** | 0.17 (kd-tree plan; 0.25 dense) | 0.48 | **0.37** | 4.8 | **3.4 GB** |
| batch 32 (pretraining) | **5.2 s** | 0.01 | 0.48 | 0.02 | 4.7 | 3.3 GB |
| batch 32, 1 thread | 50 s | 0.16 | 13.9 | 1.4 | 34 | 0.7 GB |
| batch 32, 8 threads | 8.5 s | | 1.9 | 1.0 | 5.6 | 1.3 GB |

Per shape the queries and ground truth take 0.16–0.18 s on 64 threads (5.5 s on one), shared by the
shape's steps. Pretraining's 32 rows are spread over the cloud (FPS order), so every torus is used and
nothing is skipped.

The taped network dominates: 2048 neighborhoods forward + backward ≈ 34 s on one core (2.5× the
forward), only 7× faster on 64 threads. The profile shows no lock contention (burn-autodiff has a
global graph-locator mutex, but it is not hot); the time is GEMM, element-wise zips and memmove/memset,
i.e. memory traffic: each chunk's backward materializes full weight-gradient tensors (6.4 MB) that are
then summed. Chunk 4 is the optimum (chunk 1: 7.6 s, 16: 8.5 s, 32: 16 s). glibc's mmap threshold
(`MALLOC_MMAP_THRESHOLD_=4G`) gains ~25 % on pass 4 (page faults); not set by the crate.

Full curriculum on this CPU: 498k steps × ~5.5 s ≈ 32 days. The CPU path is for tests and short runs.

## GPU estimate for WP10 (RTX 4080, not measured here: WP4 owns the GPU)

Per step (batch 1024, 2048 neighborhoods):
- network: pass 2 forward + pass 4 forward and backward = 4 × 218 MFLOP × 2048 ≈ **1.8 TFLOP**, plus
  ≈ 4 × 2048 × 7 MB ≈ 55 GB of activation traffic for LayerNorm/GELU/softmax (unfused). At 20–30
  TFLOPS effective and 600 GB/s: **≈ 150–190 ms**. Keeping pass 2's tape instead of re-running
  (single pass) would save a quarter but needs ≈ 14 GB of activations: not on 16 GB.
- loss: 3.7M pairs (T = 40; the reference's dense blend is 251M pairs) × ≈ 150 element-wise ops ≈
  4 GB of traffic plus ≈ 1,000 kernel launches: **≈ 10–20 ms**.
- host: queries + ground truth 0.08 s per step and the plan (now through the kd-tree) on the CPU, all
  independent of the parameters, so they can be computed for the next shape while the GPU trains
  (a prefetch thread; WP10).

So ≈ 0.17–0.21 s per step if the host work is overlapped: **≈ 24–29 h** for the 497,858 steps (the
paper: 45 h on a 3090). If burn-cuda reaches only 10 TFLOPS on these shapes, ≈ 0.3 s/step ≈ 42 h. The
largest lever is a fused LayerNorm/GELU (or burn's fusion feature); measure first.

GPU caveats for WP10: `select` (the gather of tori per pair) back-propagates as a scatter-add, which on
CUDA may use atomics (not run-to-run deterministic; D4 allows tolerance for the framework part, but
check); everything else in the loss is gathers and element-wise ops.

## Short CPU runs (pretraining setting: the blob, batch 32, Adam, from the tangent-sphere init)

`short_curve` (`runs_tmp/pat-wp7/curves/curve-*.tsv`): the loss on a fixed evaluation set (all
2048 rows of the blob, 245,760 fixed queries) every 16 steps; two runs side by side on 32 threads each.

| step | 0 | 16 | 32 | 48 | 64 | 80 | 96 | 112 | 128 |
|---|---|---|---|---|---|---|---|---|---|
| lr 1e-3 (the reference's) | 0.2630 | 0.2477 | 0.2480 | 0.2459 | 0.2455 | 0.2464 | 0.2470 | 0.2476 | 0.2449 |
| lr 1e-4 | 0.2630 | 0.2569 | 0.2503 | 0.2480 | 0.2483 | 0.2471 | 0.2466 | 0.2466 | 0.2462 |

The decrease comes from the eikonal term (0.228 → 0.202) while the distance term rises a little
(0.035 → 0.042): the L1 sum trades one for the other. For scale, the authors' final weights give
0.2694 on the same blob (unit normals): the tangent-sphere start is already near what this loss can
reach on a smooth blob; the eikonal floor of the blend is ≈ 0.2 here.

`short_run` (the `train()` loop with checkpoints and the JSON-lines log, `runs_tmp/pat-wp7/short/`):
256 steps (4 epochs of the pretraining phase), 1605 s, peak RSS 3.5 GB, 0 skipped steps; epoch means
0.2525, 0.2485, 0.2486, 0.2730 (per-step losses are over different rows and queries, std ≈ 0.02; the
last epoch has a spike at lr 1e-3, steps 225–256, distance up to 0.135).
`resume_is_exact`: 2 + 4 steps with a restart = 6 steps straight, bitwise (73 s).

## SPEC change proposals (for the orchestrator; nothing below is applied to SPEC.md)

1. **§10.4 blend as framework ops → host weights + framework SDF.** The weights `ŵ`, `σ`, `λ` and the
   vectors `c_i` do not depend on θ; computing them on the host (f64) and only `g_i`, `∇g_i` and the
   sums in the framework is exact and removes most of the framework work. With the truncation
   `T = 40` (weights below `e^{−40}` of the largest dropped, the kept set found by the kd-tree) the
   loss equals the dense one to < 1e-12 (f64) and the pair count drops 70× (3.7M vs 251M at batch
   1024). Proposed text: "`σ`, `λ`, `ŵ_i`, `c_i` are host constants; the framework evaluates the fit,
   the per-pair torus SDF and gradient and the weighted sums over the pairs with `λ(ρ_i − ρ_min) ≤ T`
   (`T = 40`, option, `∞` = dense)."
2. **Formula ids.** `[formula: pat-query-sampling]` (the per-row sampling as coded, a definition) and
   `[formula: pat-triangle-sampling]` (uniform barycentric sampling `(√u(1−v), √u v, 1 − √u)`;
   Osada, Funkhouser, Chazelle, Dobkin 2002, "Shape distributions", ACM TOG 21(4), §4.2) need CITED.md
   entries (WP2). `scripts/check-formulas.sh` reports them MISSING until then (the only two).
3. **§10.4 API.** As built: `train(data_dir, init, &TrainOptions, progress)` on the CPU and
   `train_on(data_dir, init, &TrainOptions, &Ctx, &ThreadPool, progress)` for any burn device (WP10);
   `TrainOptions { curriculum, step: StepOptions { loss, queries, feature_order, reference_padding,
   network_chunk, groups }, data, optimizer, schedule, seed, checkpoint_dir, checkpoint_every,
   max_steps, keep_empty_slices, files_in_buffer, resume, validate, threads }`. `batch_rows` lives in
   each curriculum phase (it differs between pretraining and the rest).
4. **Optimizer.** burn's `Adam` (`Optimizer::step` per tensor) with our own state keyed by parameter
   name, saved in safetensors. burn's `ModuleOptimizer` keys state by run-time `ParamId` and saves
   burnpack, which a restarted process cannot match to rebuilt parameters.
5. **Validation.** The reference validates (and pretrains) on files with norm-2 normals and fits with
   them (D23); its printed validation numbers are not the trained model's accuracy (0.99/10.6 vs
   0.0195/0.296 with unit normals, above). `TRAINING.md` (WP10) should report unit-normal validation,
   with the raw-normal number only for parity.
6. **EigenRule for training.** Keep `Reference` as the default for parity; propose that WP10's
   retraining uses `Robust` (it is what the f64 inference fit computes where the f32 reference form
   breaks, `|v2/v1| ≲ 3e-6`; on the authors' weights it changes the loss by ≤ 1e-6).
7. **Exact skips** (D19 addition): rows whose `∂L/∂a` is exactly zero are not back-propagated and
   tori no query pair references are not forwarded; both are exact.

## Open issues

- **GPU (WP10):** the backward of `select` (the per-pair gather) is a scatter-add; check whether
  burn-cuda uses atomics there (run-to-run determinism, D4). The host plan and ground truth must be
  prefetched on CPU threads to keep the GPU busy; GT is a brute-force winding number over all faces
  (≈ 0.16 s per shape on 64 threads), a hierarchical winding number would cut it if needed.
- **CPU pass 4** is memory bound (per-chunk weight-gradient tensors); a fused per-group
  accumulation would need burn internals. Not needed for WP10.
- **`reference_padding`** reproduces the reference's rules (all `P` rows through the network,
  `−1 → P − 1` in features and in the query boxes, no mask, NaN gradient entries zeroed); the NaN
  pattern is burn's, not bitwise JAX's. On a padded synthetic cloud 4,704 gradient entries come out
  NaN and are zeroed, as the reference would (`padded_clouds`); our default has none.
- **JAX keys** (threefry `fold_in`) are not reproduced (D22): our queries are different samples of the
  same distributions, so per-step losses are not comparable to a reference run, only statistics.
- The reference's file order in phases 1–4 came from unsorted `glob`; ours is sorted, and phase 1's
  `zip` truncation keeps `ABC17_train_mode0_{0..3}`.

## How to run

```
export CARGO_TARGET_DIR=/var/data/discrete-geometry/ddg-tools/target-pat7
cargo test --release -p ddg-points-as-tori --features train --test train                 # 16 tests, ~1 s
cargo test --release -p ddg-points-as-tori --features train --test train -- --ignored --nocapture step_timing
PAT_STEPS=256 cargo test ... -- --ignored --nocapture short_run      # train() loop, runs_tmp/pat-wp7/short/
PAT_LR=1e-3 PAT_STEPS=128 cargo test ... -- --ignored --nocapture short_curve
cargo test ... -- --ignored --nocapture reference_loss_export        # then (run from outside the reference repo):
JAX_PLATFORMS=cpu .../PointsAsTori-build/.venv/bin/python -I runs_tmp/pat-wp7/ref_loss.py runs_tmp/pat-wp7/refloss <reference repo>
cargo test ... -- --ignored --nocapture resume_is_exact
```
