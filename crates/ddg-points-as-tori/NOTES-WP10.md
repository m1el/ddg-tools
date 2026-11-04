# WP10 notes: GPU training and the full retraining run

Branch `points-as-tori`, build dir `target-pat10`, logs of the development runs in `runs_tmp/pat-wp10/`, the
long run in `runs_tmp/pat-train/` (how to run, monitor, resume and export: `TRAINING.md`). RTX 4080 16 GB,
Threadripper 3970X (32 cores / 64 threads), burn 0.22.1 + burn-cuda with fusion, TF32 matmuls (user decision
2026-11-04).

## 1. GPU inference default: `GpuMatmul::Tf32`

`PredictOptions::gpu_matmul` defaults to `GpuMatmul::Tf32` (user decision); `Split3` stays selectable
(`src/network/mod.rs`, docs in `src/network/cuda.rs`). `tests/gpu.rs::pipeline_gpu_vs_cpu` now runs both modes
with per-mode tolerances (5,000-point ellipsoid, reference weights, 4,000 SDF queries):

| mode | coefficients `max |Δa|/(1+|a|)` (tolerance) | SDF `max |Δφ|` (tolerance) | run to run |
|---|---|---|---|
| Tf32 (default) | **1.32e-3** (3e-3) | 4.0e-6 (2e-5); median 5.9e-7 | bitwise |
| Split3 | 3.94e-6 (1e-5) | 1.0e-8 (1e-6) | bitwise |

**SPEC proposal 1 (§18 WP4 item 1).** "Tf32 ≈ 2e-4 on coefficients" holds for the raw network outputs
(1.5e-4, measured again); the coefficients `a₃…a₅ = y/σ` amplify it by `1/σ`, so the coefficient tolerance
is 3e-3·(1 + |a|) (measured 1.3e-3), and the SDF moves by ≤ 4e-6 on this cloud.

## 2. GPU training (`src/train/gpu.rs`)

The step of WP7 (`step.rs`) is split into its host part and its device part:
`plan_slice` (blend constants on the kd-tree, the network rows it needs, and for the GPU the pair layout) and
`step_planned` (the passes); `step = plan_slice + step_planned`. On a CUDA `Ctx` (`gpu::cuda_ctx`,
`step::on_gpu`) `step_planned` calls three GPU passes instead of the CPU ones:

- **Pass 2, coefficients:** the network without the tape in batches of the fixed shape
  `StepOptions::network_chunk` (256 on the GPU) × 64 tokens, the last batch padded with copies of its first row,
  dispatched from one thread.
- **Pass 3, loss and ∂L/∂a:** the tori are fitted once (`fit_tensor`, the same code as the CPU); the pairs
  gather the *detached* tori into one leaf per dense query chunk (chunks as `loss::chunks`, budget
  `LossOptions::pair_chunk` = 2²⁰ padded pairs on the GPU); one backward gives `∂L/∂T` per pair; a **gather**
  over the transposed pair lists (`PairLayout`, built on the host ahead of time, counting sort) and a
  fixed-shape `sum_dim` give `∂L/∂T` per torus; a second backward of `⟨T(a), ∂L/∂T⟩` gives `∂L/∂a`. Exact
  (the chain rule regrouped), no `select`/`gather` backward, no atomics.
- **Pass 4, backprop:** the network with the tape on the rows with nonzero `∂L/∂y`, the same fixed batches,
  parameter gradients summed on the device in batch order.

**The training forward (`gpu::forward_split`).** The same network as `PatNet::forward`, the same operations
(checked in f64 on the CPU: outputs bitwise equal, parameter gradients 3e-16, with and without a pad mask,
`forward_split_equals_patnet`), but linears and LayerNorms on `[64, rows/64, ·]` views with the biases and
LayerNorm scales/shifts broadcast through `[64, 1, c]`. Reason (CubeCL kernel profile,
`CUBECL_DEBUG_LOG=… CUBECL_DEBUG_OPTION=profile`, `runs_tmp/pat-wp10/cubecl-prof*.log`): burn reduces a bias or
LayerNorm-parameter gradient over all 16,384 rows of a `[16384, c]` gradient in one reduction, which CubeCL runs
with its *unit* strategy on `c/4` threads walking every row: 0.24–0.9 ms per parameter, **46 % of the GPU time
of a step** (1,009 of 2,186 ms over 4 steps). Split into two short reductions it is 0.27 s, and the step's GPU
time falls from 2.19 to 1.64 s per 4 steps (−25 %). Other formulations measured: plain 3-D layout with
`[B, 1, c]` broadcasts (2.18 s: brings back slow weight-gradient matmul kernels), 2-D `matmul` for the weights
(2.54 s), groups 32/64/128 (1.82/1.64/1.79 s, within noise; 64 kept).

**Host side.** `build_plan`'s final merge is now parallel (same values and order); the GPU pair layout is
built in parallel. Per phase-0 shape on 60 threads, under load 20–90: queries + ground truth 0.34–0.41 s,
plan 0.22 s per slice (was 0.35 s: 0.10 blend constants + 0.07 GPU layout + merge).

**Prefetch (`run::train_on`).** A loader thread reads the next shuffle buffer's files (in parallel) one buffer
ahead; a preparer thread computes queries, ground truth and `plan_slice` up to 4 slices ahead; the training
loop only runs `step_planned` (on the calling thread on the GPU, inside the pool on the CPU) and the
checkpoints. Events carry the cursor, so checkpoints, resumption and the data order are unchanged (WP7's
`resume_is_exact` and all train tests pass). On the CPU the prefetching overlaps the host work with the steps
too. A burn-cuda failure surfaces as a panic at the next read-back; the loop catches it and returns
`PatError::Backend` (an OOM was provoked at batch 1024 and is reported, not a crash).

`TrainOptions::validate_every_steps` (new; default 0 = every epoch as WP7): validation at an epoch end only
after that many steps (always at a phase's last epoch). The run uses 1000, so pretraining (64 epochs of 64
steps) validates every 16 epochs, every other phase at every epoch.

## 3. Determinism (WP7 open issue)

**burn-cubecl's index backwards.**
- `select` backward = `float_select_assign(Add)` = `select_assign_kernel`: one thread per element of the
  *non-indexed* dims (9 threads for `[n, 9]` tori), each walking all indices serially. **No atomics,
  deterministic**, but serial: ~1 s per step at 3.7M pairs. Not used (pass 3 above).
- `gather` backward = `scatter_add_atomic` on CUDA (`fetch_add`): **not repeatable**. Not used.
- `float_sum` (full sum) = `SumStrategy::OneShot` with atomics when `autotune` is off (our build): not
  repeatable. Used only where its value is not read (the backward seed); all reductions whose values matter
  are `sum_dim` (fixed-shape reduce kernels) or host sums in f64.
- No autotune in the build (kernel choice is a function of the shape).

**Measured.** Two repetitions of 4–6 steps (loss, all gradients and the updated parameters hashed) from the
authors' weights and from the initialization: **bitwise identical** in the same process and across processes
(`examples/train.rs bench`, `runs_tmp/pat-wp10/bench*.log`), with every forward variant above. The GPU loss
unit test repeats bitwise. The long run itself: after a stop and resume, the 34 replayed steps were bitwise
equal to the first launch's. **Resumption** on the GPU: 128 straight steps vs 64 + resume + 64: all 128 losses,
all 136 weight tensors and all 272 Adam moments bitwise equal (`runs_tmp/pat-wp10/resume*`). So the D4
standard (bitwise run to run) holds for GPU training; no tolerance standard was needed. Not tested: bitwise
across different GPUs or driver versions (CubeCL compiles per device), across network batch sizes (the bits
depend on the batch shape: 256 vs 512 differ at round-off), and across burn versions.

## 4. GPU vs CPU

**One step from the same parameters** (`train compare`; same plan; CPU f32 Flex vs GPU TF32):

| shape, rows | weights | loss rel. diff | gradient `|Δg|/|g|` | cosine |
|---|---|---|---|---|
| blob (pretraining), 0..32 | authors | 1.9e-6 | 6.0e-2 (1.7e-2 before the grouped forward) | 0.9983 |
| blob, 0..32 | init | 1.4e-8 | 2.0e-3 | 0.999998 |
| phase-0 blob, 0..1024 | authors | 4.3e-7 | 5.3e-3 | 0.99999 |
| phase-0 blob, 0..1024 | init | 4.1e-9 | 3.3e-3 | 0.999995 |

The loss agrees at f32 level; the gradients at TF32 level (the network's matmuls; the grouped weight
gradients are 64 partial TF32 products, which moves the blob/authors case from 1.7e-2 to 6e-2). The loss part
alone (`gpu_loss_matches_cpu_and_repeats`, synthetic 800-point cloud, 121k pairs): loss 1.3e-7, `∂L/∂a`
**8.9e-6** with `EigenRule::Robust`; with `EigenRule::Reference` 2.2e-2 (the f32 cancellation of the
reference eigenvector form near `v2 = 0`, NOTES-WP7; CPU and GPU f32 round differently there). Training uses
Robust.

**Short run, same data and seed** (pretraining phase, 48 steps, the `train_on` loop on both devices,
`runs_tmp/pat-wp10/compare.log`): per-step losses agree to ≤ 1.07e-2 relative (median ~3e-3); steps 1–2 to
1e-8 / 4e-5, then the two trajectories separate at TF32 level while staying statistically identical (the
per-step loss varies by ~8 % between batches). GPU 48 steps in 26.6 s, CPU (64 threads, loaded) ~6 min.

## 5. Step time, memory, utilization

| condition | step (fwd + loss + bwd + update) | GPU memory | SM (dmon) |
|---|---|---|---|
| quiet host (load ~20), before the grouped forward | **0.56 s** (0.05 + 0.08 + 0.41 + 0.015) | 8.0 GB | ~60–100 % |
| GPU kernel time per step, kernel profile | 0.55 s → **0.41 s** with the grouped forward | | |
| loaded host (load 80–110, other agents' `nice 19` jobs), batch 256 | 0.54–0.96 s, typically **0.6 s** | 8.2 GB | 22–31 % |
| same, batch 512 | 0.49–0.74 s, typically 0.52 s | 12.5 GB | |
| batch 1024 | out of memory (with a 1 GB Blender process) | | |

**Utilization.** With the host loaded, the three dispatching threads (the caller, CubeCL's `DSU`/`DSD`) run at
92–96 % CPU and the GPU is busy only 22–31 % of the time: ~13,000 kernel launches per step, host-bound in burn's
per-op dispatch (NOTES-GPU's finding, now for training). The other agents' jobs already run at `nice 19`; they
still take the SMT siblings of our dispatch threads. On a quiet host the GPU time (0.41 s) bounds the step.
Batch 512 halves the launches of the network passes and is faster when host-bound, but needs 12.5 GB, which
leaves too little for the render agent's Blender processes (1–2.7 GB); the run uses 256 (8.2 GB) and can be
restarted with `--network-batch 512` when the GPU is ours alone (TRAINING.md).

**Host prefetch in phases 0–4** (measured 05:20, load ~105, one phase-0 shape, CPU only): queries + ground
truth 0.83 s wall on 16 threads, 0.45 s on 60; plan 0.149 / 0.134 s per slice (only partly parallel). A shape is
two steps ≈ 0.9 s of GPU work at 0.45 s/step, so the preparer needs ≈ 20–30 threads to stay ahead (pretraining
needs almost none: one shape per 64 steps). **Headroom for the run:** other agents ≤ ~24 threads in total
(machine load ≲ 32) gives the three dispatch threads idle SMT siblings and the preparer ~30 threads; sent to the
orchestrator at 05:21.

**Full curriculum** (497,858 steps): 0.6 s/step → **≈ 83 h** under the current load; ≈ 0.45 s/step → 62 h on a
quiet host; the 30–35 h estimate of SPEC §18 assumed WP4's inference throughput for the backward too, but the
backward is ~6× the forward in launches and time. The measured ETA of the run is in TRAINING.md.

## 5b. The run's pretraining: a jump at step ~810

In the live run (lr 1e-3, the reference's pretraining) the blob loss sat at 0.242 (distance 0.043, eikonal
0.200) for 12 epochs and then jumped within ~10 steps (809–818) to 0.28 (distance 0.137, eikonal 0.145), where it
stayed (epochs 13–16). WP7's CPU run had the same jump (distance up to 0.135 at its 4th epoch), so it is the
recipe, not the GPU path. `train eval` (CPU, unit normals) puts it in context:

| weights | blob distance / eikonal / sum | `ABC17_validation_3` distance / eikonal / sum |
|---|---|---|
| authors' final (`pat-reference`) | 0.0334 / 0.2322 / 0.2656 | 0.0195 / 0.2961 / 0.3156 |
| authors' `FundamentalFormPredictor_pretrained.pkl` (converted, not committed) | 0.0965 / 0.1930 / 0.2895 | 0.0572 / 0.2905 / 0.3477 |
| ours, step 1,088 (after the jump) | 0.1386 / 0.1409 / 0.2795 | 0.0895 / 0.2187 / 0.3082 |

The authors' `_pretrained` weights show the same high-distance/low-eikonal profile, which suggests (it does not
prove) that their pretraining ended in a similar state and phase 0 (lr 1e-5) moves away from it. Not
intervened; the phase-0 validation curve will tell. If phase 0 does not recover the distance term, a rerun of
pretraining at a lower lr (a deviation from the recipe) is the proposal.

## 5c. Optimization round (user decision "optimize first, then run", 05:40–05:57)

The run was stopped right after its step-2,880 checkpoint (05:39; checkpoint backed up in
`runs_tmp/pat-train/saved-step2880/` with SHA256SUMS, and verified by resuming it for one step).

**What changed** (commits `801b3c9` → `d9f1606`):
1. **CUDA graphs** for the network batches (`gpu.rs`, `burn::tensor::capture` / `Graph::replay`, burn 0.22.1 on
   CubeCL's hardware graphs). Pass 2 (forward, no tape) and pass 4 (forward with tape + backward + the 40 parameter
   gradients) of one fixed-shape batch are captured once per thread (warmup ×3, then capture; ~2 s) and replayed
   per batch: ~1,500 launches become one. The parameters live in uniquely owned stable buffers that are written in
   place (`slice_assign` on a refcount-1 tensor) from the current `NetParams` before each pass; features, output
   weights and the pad mask likewise. Graphs are keyed by device, architecture, batch shapes and mask presence
   (padded clouds get their own pair). `StepOptions::gpu_graphs` (default true), `--no-graphs`.
   Check: `graph_replay_matches_eager` (replays with changed inputs and parameters equal the eager closure
   bitwise).
2. **Flat-parameter training forward** (`forward_flat`, the same ops as `forward_split`; bitwise equal on the GPU to
   the previous binary with graphs off; `PatNet` equality in f64 kept).
3. One device-side finiteness test of the gradients (`Σ 0·g` per tensor, one read-back) instead of 40 read-backs.
4. Shapes re-chosen with graphs on: network batch 256 (512: bwd 0.255 vs 0.197 s), untaped forward batch 256
   (1024: 0.095 s vs 0.048 s; 2048: 0.112 s), loss pair chunk 2¹⁹ (2¹⁸: 0.08 s, 2²⁰: 0.08 s, 2¹⁹: **0.065 s**,
   2²²: 0.16 s).
5. `t_wait` in `log.jsonl`: seconds the loop waited for the preparer.

Tried and not kept: larger batches (above); more fusion (no knob beyond `max_ops = 64`); our own fused
attention/LayerNorm kernels (not reachable through burn's dispatch API without the device buffers; the attention
softmax is the largest remaining fused kernel, ~0.8 ms per layer per batch).

**Measured** (phase-0 shape, `bench`, load 24–30 with WP14 at ≤ 24 threads, nice 19):

| binary | fwd | loss | bwd | update | step | GPU memory |
|---|---|---|---|---|---|---|
| `801b3c9` (the run's, eager, forward 1024, chunk 2²⁰) | 0.095 | 0.080 | 0.240 | 0.016 | **0.44** | 7.1 GB |
| `d9f1606` (graphs, forward 256, chunk 2¹⁹) | 0.048 | 0.065 | 0.197 | 0.018 | **0.34** | 4.2 GB |

In the training loop from the same checkpoint (160 pretraining steps, load ~25): **0.495 → 0.301 s/step (−39 %)**,
mean `t_wait` 3 ms. Loss curves from the same checkpoint: step 1 equal to 1.7e-10 relative, then median 3.8e-3, max
2.0e-2 relative per step over 160 steps (TF32-level separation of the trajectories, as GPU vs CPU in §4); mean loss
0.2698 vs 0.2694. The new binary resumes the old binary's checkpoints. Run to run: bitwise (bench, 2 × 8 steps; the
loop, 80 steps from the same checkpoint). Tests: `tests/train.rs` (16, FD checks), the GPU unit tests,
`graph_replay_matches_eager`, `tests/gpu.rs` (12): pass.

The bits of `d9f1606` differ from `801b3c9`'s (burn-fusion segments a captured closure differently from the eager
stream), within the TF32 standard above.

**Phase 0 in the loop was slower than the bench** (the bench reuses one plan): with `d9f1606` the run's phase-0
steps took 0.60 s (`t_loss` 0.29–0.39 s, flat, so not compilation). Two causes, two fixes (commit `d7a8200`):

6. **Loss shapes vary every step**, so burn-fusion explores every step's stream anew. Chunk shapes (queries and
   pairs per query, tori and positions per torus) are now rounded up a fixed ladder (`gpu::ladder`: 1…8, then
   multiples of 8 growing ~1.25×); padded queries are masked out of the sums, padded torus rows read the zero row.
   In the loop: `t_loss` 0.328 → 0.257 s; loss curve vs `d9f1606` median 1.7e-5, max 1.3e-4 relative (only the
   zero padding of the sums changes).
7. **Our own preparer competes with the step for the host.** With 16 preparer threads the step drops to 0.40 s but
   the loop waits 0.35 s for data; 32 threads: 0.50 + 0.06; 60: 0.56 + 0.01. The ground truth is ~70 % of the
   preparer (≈ 13 core-s per phase-0 shape). **The winding-number sums now run on the GPU**
   (`gpu::solid_angle_sums`, an NVRTC f64 kernel, one thread per query over the faces in the CPU's operation order,
   no FMA, own stream; `QueryOptions::gpu_truth`, on by default for GPU runs in `examples/train.rs`,
   `--cpu-truth` to turn it off; distances stay on the CPU BVH; meshes below 512 faces stay on the CPU, where they
   are faster). Ground truth GPU vs CPU on 6 shapes (1.47M queries, blobs, ABC, simulated): **0 sign and 0 bit
   differences**. Kernel 0.186 s per 2048-face blob shape. In the loop from the same checkpoint the 120 losses with
   GPU and with CPU ground truth are **bitwise equal** to each other (and to the ladder run).

| binary, phase-0 loop from the step-4,096 checkpoint, 120 steps | load | fwd | loss | bwd | upd | wait | s/step |
|---|---|---|---|---|---|---|---|
| `801b3c9` (the original run binary, eager) | 37–50 | 0.124 | 0.260 | 0.507 | 0.010 | — | **0.90** |
| `d9f1606` (graphs) | 41–46 | 0.074 | 0.328 | 0.210 | 0.015 | 0.013 | 0.64 |
| ladder | 41–46 | 0.075 | 0.257 | 0.217 | 0.013 | 0.013 | 0.57 |
| ladder, CPU ground truth | 30–56 | 0.072 | 0.209 | 0.209 | 0.010 | 0.012 | 0.51 |
| **`d7a8200`** (ladder + GPU ground truth) | 34–41 | 0.150 | 0.108 | 0.198 | 0.008 | 0.012 | **0.48** |

So phase 0 went from 0.90 to **0.48 s/step (−47 %)** at load 30–55 (WP14 at ≤ 24 threads, nice 19). The run was
stopped once more (06:09, its last checkpoint the phase-0 start at step 4,096, backed up in `saved-step4096/`) and
resumed with `d7a8200`; steps 4,097–~4,700 of `d9f1606` appear twice in `log.jsonl` (later lines supersede).
Remaining bound: GPU time (network ≈ 0.25 s, loss ≈ 0.1 s, ground truth ≈ 0.09 s per step) with the host still
contended; an f32 ground-truth kernel (~30× cheaper) would change the sign of near-surface points and was not used.

## 6. SPEC proposals (nothing applied to SPEC.md)

1. §18 WP4 item 1: the Tf32 coefficient tolerance is 3e-3·(1 + |a|) (raw outputs 2e-4); above.
2. §18 WP7 item 9 / D4: answered: `select`'s backward is deterministic (serial, no atomics) but too slow;
   `gather`'s is atomic; full `sum` is atomic. The GPU loss uses detached gathers, transposed gather
   reductions and `sum_dim` only, and GPU training is bitwise run to run and resumes bitwise; D4's "bitwise run
   to run" therefore holds for training as well, with the caveat that the bits depend on the network batch
   shape (256 by default).
3. §10.4 "expected time": replace "≈ 30–35 h (TF32)" by the measured ≈ 0.48 s per phase-0 step with CUDA graphs,
   fixed loss shapes and GPU ground truth (≈ 66 h for the curriculum) when other agents keep to ≤ 24 threads; without graphs the bound was burn's
   per-op dispatch (13k launches per step). Further speed needs our own fused attention/LayerNorm kernels in the
   training graph, which departs from "the network in the framework" (D2); not done.
4. §10 / D19: `TrainOptions::validate_every_steps` (above), so that validating per epoch does not double the
   pretraining phase (64 epochs of 64 steps, ~40 s each).

## 7. Files

- `src/train/gpu.rs` (new): `is_cuda`, `cuda_ctx`, `PairLayout`/`pair_layout`, `loss_and_grad`,
  `coefficients`, `backprop`, `forward_split`; tests (pair layout transposition, GPU loss vs CPU and bitwise
  repeat, `forward_split` = `PatNet` in f64, an ignored profiling test).
- `src/train/step.rs` (WP7): `SlicePlan`, `plan_slice`, `step_planned`, `on_gpu`, `scale_raw`, `select_rows`
  `pub(crate)`; GPU dispatch at the three passes.
- `src/train/run.rs` (WP7): prefetching `train_on` (loader + preparer threads, event loop with the same
  cursor/checkpoint semantics), `validate_on`, `TrainOptions::validate_every_steps`, panics caught.
- `src/train/loss.rs` (WP7): parallel merge in `build_plan` (identical result).
- `src/network/mod.rs`, `src/network/cuda.rs` (WP4): the Tf32 default and docs; `tests/gpu.rs`: the per-mode
  pipeline test.
- `examples/train.rs` (new; `[[example]] required-features = ["train"]` in `Cargo.toml`): `run`, `bench`,
  `compare`, `export`.
- `Cargo.lock`: unchanged.

## 8. Open issues

- Host-bound under outside load (above). The run is at normal priority; the other agents' jobs are at nice 19.
- `forward_split` is used for GPU training only; GPU inference keeps WP4's `network/cuda.rs`.
- The loss chunks pad to 6.2M pairs for 3.8M real ones (1M-pair budget); padding is cheap on the GPU, a
  smaller budget costs launches.
- The GPU ground truth (brute-force winding number) stays on the CPU (0.34–0.41 s per shape, prefetched).
- A `Geometric`-rule run is deferred (SPEC §18 WP7 item 6, until WP14 shows the sign rule matters).
