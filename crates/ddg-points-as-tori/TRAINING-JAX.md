# Retraining with the authors' JAX trainer

The user decided (SPEC, 2026-11-04) to stop the burn retraining (`TRAINING.md`, `runs_tmp/pat-train/`, kept as is)
and retrain with **the authors' own training program**: `training/train.py` with `mesh_sampling_jax.py` and
`mesh_data.py` of points-as-tori (Feng, Gkioulekas, Crane, TOG 2026; code MIT, "Copyright (c) 2026 Nicole Feng",
`LICENSES/points-as-tori-MIT.txt`). It runs their curriculum, hyperparameters, data, seeds and f32 math, including
the known quirks: norm-2 normals in the pretraining blob and the validation file, a cosine schedule that is a
no-op, and the unsorted `glob` file order (NOTES-REFERENCE.md §9).

Live directory: `runs_tmp/pat-train-jax/` (gitignored). Versioned copies of the driver, patch and tests are in
`tools/train_jax/` and `tools/train_jax/geometric/`; the copies under `runs_tmp/` are the ones that run.

## Setup

| item | value |
|---|---|
| authors' code | clone `/var/data/discrete-geometry/PointsAsTori-build/points-as-tori`, commit `253bf9f`, unmodified (`git status` clean); C++ extension as built for the golden runs (`runs_tmp/pat-ref/README.md`); the installed `pointsastori` package equals `src/pointsastori` |
| sha256 | `train.py` b0972e06…, `mesh_sampling_jax.py` c1d56619…, `mesh_data.py` 6b2b02e3… |
| venv | `/var/data/discrete-geometry/PointsAsTori-build/.venv`, Python 3.11.15 |
| JAX stack | jax / jaxlib / jax-cuda12-plugin / jax-cuda12-pjrt 0.8.2, flax 0.12.2, optax 0.2.8, numpy 2.4.6, scipy 1.17.1; CUDA wheels: cudnn 9.27.0.42, cublas 12.9.3.13, nvcc/ptxas 12.9.86 |
| added for training | their README's list, `uv pip install --python <venv>/bin/python scipy numpy-stl matplotlib thingi10k py7zr pyvista trimesh pyfqmr noise`: numpy-stl 4.1.3, matplotlib 3.11.2, thingi10k 1.5.0, py7zr 1.1.4, pyvista 0.49.2 (vtk 9.7.2), trimesh 5.2.0, pyfqmr 0.5.0, noise 1.2.2, plus their dependencies (full list: `runs_tmp/pat-train-jax/venv-freeze.txt`; the only change to an existing package was fsspec 2026.9.0 → 2026.7.0). `mesh_data.py` imports all of them; training uses none of them. |
| data | `PAT_training_data/PAT_training_data/` as released (46 training npz, all padded to P = 2048, k = 64, V = F = 2048) |
| machine | RTX 4080 16 GB (driver 595.99.02, also drives the display), Threadripper 3970X, 94 GB RAM |

### Driver (`driver.py`)

Their code uses paths relative to the working directory (`datasets/…`, `models/`, `PlotTrainingHistory.png`).
The driver:
- chdirs into `runs_tmp/pat-train-jax/work/`, which holds `datasets → PAT_training_data/PAT_training_data`;
- puts the patched `train.py` (below) first on `sys.path`, then their `training/` directory;
- calls their `train()` unchanged.

It wraps `TrainFundamentalFormPredictor` (a module global that `train()` looks up) only to:
- skip phases already marked done in `work/checkpoints/<phase>.done`;
- move each phase's `PlotTrainingHistory.png` to `work/plots/<phase>.png`;
- stop before the first phase not in `PAT_PHASES` (used by the Geometric probe).

The driver also:
- turns on JAX's persistent compilation cache (`runs_tmp/pat-train-jax/jax-cache`, with
  `min_compile_time_secs 0` and `min_entry_size_bytes −1`; XLA's per-fusion autotune cache comes with it);
- counts compilations with a `jax.monitoring` listener. With `PAT_LOG_COMPILES=1` it logs each compile by name.

### Patch (`PATCHES.diff`, the only edit)

Their trainer saves only weights pickles once per epoch (`models/FundamentalFormPredictor_{epoch}.pkl`; the names
collide across phases), with no optimizer state and no way to resume. For a 54 h run that edit can't be avoided. A copy of
`train.py` (`runs_tmp/pat-train-jax/src/train.py`) adds **only new lines** (no original line is changed or removed):

- **Checkpoints.** `work/checkpoints/<phase>.npz` (npz, no pickle; previous one kept as `<phase>.prev.npz`) holds:
  - all model state (params and dropout RNG state) and all optimizer state (Adam count/mu/nu, schedule count, step);
  - a cursor (epoch, next shape in the epoch's permutation, phase step, batch count);
  - the logged loss lists.
- **When.** Every `PAT_CKPT_SECONDS` (900 s) at a shape boundary, and at every epoch end (after their save and
  validation).
- **Resume.** Restores everything. The epoch loop still draws each finished epoch's `np.random.permutation` before
  skipping it, so the data order and their JAX keys (`fold_in(epoch)`, `fold_in(global_idx)`) are unchanged.
  It then skips the finished shapes of the current epoch. Their own printouts (epoch averages, validation, plot)
  continue as if uninterrupted.
- **Progress line** every `PAT_LOG_EVERY` (100) steps (synchronizes on that step's loss):
  `[PAT] <UTC time> phase <tag> epoch e/E shape i/S step s/phase_total global g/497858 loss … d … e … lr …
  <s/step> eta_phase … h eta_total … h peak_mem … GiB compiles n (… s, cache hits …)`.
  `eta_total` extrapolates the current phase's rate, so it reads low during pretraining.
- `PAT_MAX_STEPS` (tests only): checkpoint and stop after that many steps.

## GPU memory: the one setting that differs

The phase-0…4 train step (the network on all 2048 points, the loss with second derivatives on 1024 × 120 queries
against all 2048 tori) needs **24.5 GB** before XLA rematerialization. On a 24 GB RTX 3090 (the paper's card) JAX's
default pool is 75 % = 18 GB, so XLA rematerializes there too. On 16 GB:

| `XLA_PYTHON_CLIENT_MEM_FRACTION` | pretraining | phase 0 |
|---|---|---|
| default 0.75 (12 GB) | OOM in `generate_shape_features` (brute-force ground truth, 2048 × 120 × 2048 query–face pairs) | — |
| 0.95 | fits (peak 14.57 GiB) | OOM: remat brings the step to 14.98 GiB, a single 14.97 GiB allocation that doesn't fit |
| **0.97** (used) | fits, peak 14.57 GiB | fits, peak 14.97–15.10 GiB of a 15.51 GiB pool; epoch-end validation at batch 1024 after a train step fits too (`val_memory_test.py`, peak 14.95 GiB) |

This is an allocator setting. The math is untouched; rematerialization only recomputes intermediates.
`--xla_gpu_deterministic_ops=true` needs 16.35 GiB in phase 0 and does not fit.

## Determinism and resume test

- **Default flags.** Their trainer is **not bitwise reproducible** on this GPU: two straight pretraining runs
  differ from step 16 on (XLA's default GPU reductions and autotuning). That holds for the authors' setup too.
  "Faithful" here means the same program, data, seeds and hyperparameters; the bits differ from run to run.
- **Deterministic ops (test only).** With `XLA_FLAGS=--xla_gpu_deterministic_ops=true`, a resumed run is
  **bitwise identical** to a straight one: every progress line, all of their epoch and validation printouts, and
  all 472 checkpoint arrays (weights, Adam moments, counters, RNG state).
  - Pretraining: 192 straight vs 128 + resume (`resume_test.sh`; the cursor sits at the end of an epoch's shapes).
  - Mid-epoch across shapes: their `TrainFundamentalFormPredictor` on one blob file at batch 256, 96 straight vs
    40 + resume (`resume_logic_test.py`; the cursor sits inside the shuffled permutation, at shape 5).

  Phase 0 itself at batch 1024 can't be tested this way (OOM, above).

So `run.sh` restarts after a failure. Up to 20 times, 120 s apart, it reruns the driver, which skips finished phases
and resumes the current one; at most 15 min of steps are repeated.

## Timing (RTX 4080, measured 2026-11-04, nothing else on the GPU)

| phase | s/step | notes |
|---|---|---|
| pretrain (batch 32) | **0.169** | plus ≈ 2.7 s per 64-step epoch (validation and checkpoint); ≈ 15 min in all |
| phase 0 (batch 1024) | **0.392–0.399** | 300-step test; phases 1–4 run the same executables (all files share the padded shapes) |

- **GPU.** SM 100 % and memory controller 99.98 % busy in steady state; 238–245 W of 340 W; clocks at boost
  (2.77 GHz SM, 10.8 GHz memory); frame buffer 15.4–15.8 GB (the pool).
- **Host.** About one core (94–150 % CPU), RSS 2.8 GB (pretraining) to 4 GB (phase 0). Phase 4 holds about 16 GB
  of training arrays in RAM.
- **Why slower than the paper's 0.33 s/step on a 3090.** The step is DRAM-bandwidth-bound, and the RTX 4080 has
  less memory bandwidth than the RTX 3090: 256-bit, ≈ 0.69–0.72 TB/s vs 384-bit, 0.94 TB/s. Bandwidth scaling
  predicts 0.33 × 0.94/0.69 ≈ 0.45 s/step; the measured 0.39 is better than that, presumably thanks to the 4080's
  larger L2. The 4080's higher FP32 rate doesn't help, and host RAM speed doesn't matter (one core, a few MB per shape).
- **Projection.** 4,096 × 0.17 s + 493,762 × 0.393 s + about 0.3 h of validation and data loading ≈ **54.5 h**.

  | phase | steps | expected |
  |---|---|---|
  | 0 | 96,000 | 10.5 h |
  | 1 | 64,000 | 7.0 h |
  | 2 | 119,344 | 13.0 h |
  | 3 | 124,590 | 13.6 h |
  | 4 | 89,828 | 9.8 h |

  Launched 2026-11-04 15:19 UTC, so the expected end is **2026-11-06 ≈ 22:00 UTC (± 1 h)**.

### Compilations (`PAT_LOG_COMPILES=1`)

- **Phase 0.** 51 compilations, all before step 10: `train_step` once (22 s, including XLA autotuning, the bulk of
  the ptxas calls), `generate_shape_features` once (7.6 s), and about 49 small eager ops (broadcasts, converts,
  parameter init, RNG). From step 10 to 300 (150 shapes, two of the four data files touched) the count stays at
  51. There is **no** compilation per step, per shape or per file.
- **Pretraining.** 48 at startup plus 5 at the first epoch end (validation: eager model call, a second
  `generate_shape_features` signature for the numpy validation inputs, `validation_loss`); none at later epoch ends.
- **Step-time windows.** The 100-step windows that contain an epoch end show 0.33 s/step in pretraining. That is
  validation and checkpoint time, not compiling.
- **Persistent cache.** A cold start spends 31 s compiling. With the cache, a restart (or a test, or the Geometric
  probe) spends 2 s (51 of 51 cache hits; cache 3.3 MB).

## Monitoring

All paths are under `/var/data/discrete-geometry/ddg-tools/runs_tmp/pat-train-jax/`.

```bash
grep "\[PAT\] 20" train.log | tail -1                       # step, phase, loss, s/step, ETA, memory, compiles
grep -E "^Epoch|validation:|PAT\] phase|PAT\] checkpoint|=== " train.log | tail
cat PIDS                                                    # wrapper, python (per start), dmon
kill -0 <python pid> && echo running
tail -3 dmon.log                                            # one line per minute (pwr, temp, sm, mem, clocks, fb)
ls work/checkpoints work/models work/plots
```

Their per-shape lines (`Epoch e, shape i/S: loss = … `; every shape in phase 0, every 100 shapes later) and their
epoch summaries with validation (norm-2 validation normals, so not comparable with TRAINING.md's unit-normal
numbers) are in `train.log` as well.

**Stopping.** Kill the wrapper first, then python (otherwise `run.sh` restarts it): `kill <wrapper>; kill <python>`.

**Resume by hand** (after a reboot, or after the restart budget is used up):

```bash
cd /var/data/discrete-geometry/ddg-tools/runs_tmp/pat-train-jax
nohup setsid ./run.sh > run.sh.out 2>&1 < /dev/null &
nohup setsid nvidia-smi dmon -s pucm -d 60 -c 5000 >> dmon.log 2>&1 < /dev/null &
```

`run.sh` appends the new PIDs to `PIDS`.

## Outputs and conversion

`work/models/` holds the authors' files:
- `FundamentalFormPredictor_pretrained.pkl` and `FundamentalFormPredictor_phase{0..4}.pkl`, one per phase end
  (each phase loads the previous phase's file);
- `FundamentalFormPredictor_{epoch}.pkl`, per epoch and overwritten across phases (their naming).

Convert with the restricted unpickler (never a plain `pickle.load`). This was tested on a pickle written by this
trainer: 136 tensors, 1,588,358 parameters.

```bash
cd /var/data/discrete-geometry/ddg-tools
uv run --no-project --with 'numpy>=2' --with 'safetensors>=0.4' -- python3 -I \
    crates/ddg-points-as-tori/tools/convert_weights.py \
    runs_tmp/pat-train-jax/work/models/FundamentalFormPredictor_phase4.pkl \
    data/points-as-tori/retrained-jax-reference.safetensors \
    --source "pat-train-jax: authors' train.py curriculum, end of phase 4 (RTX 4080, 2026-11)"
```

Intermediate phases convert the same way (e.g. `…_phase0.pkl` → `retrained-jax-reference-phase0.safetensors`).
`data/` is gitignored; evaluate as in `EVALUATION.md`.

## Data order on this machine

Phases 1–4 build their file lists with unsorted `glob`, so the order (and in phase 1, which 4 of the 6 ABC17 files
`zip` keeps) depends on the filesystem (ext4 here). The lists as their code computes them are in
`tools/train_jax/glob-order.txt`. Phase 1 gets blobs 1, 0, 3, 2 interleaved with ABC17 files 4, 0, 3, 2.
`train.log` also records every file loaded, in order.

## Geometric sign-rule probe (prepared, not launched)

`runs_tmp/pat-train-jax-geo/` (copies in `tools/train_jax/geometric/`) holds the probe of SADDLE-DEMO.md and
NOTES-WP1.md.

**The change.** A copy of the `pointsastori` package with one added line in the JAX torus fit
`torus_precompute_masked` (`network.py`; the training loss and validation both use it):

```python
center = jnp.where(4.0 * a20 * a02 - a11 * a11 < 0.0, 2.0 * (p_i + a00 * n_i) - center, center)
```

At a saddle (det II < 0) the centre moves to the osculating side, o + n*/κmin, mirrored through the code's contact
point o = p + a00 n. Radii, axis and sign are unchanged. `train.py` is the same patched copy as the main run.

**Check** (`check_geometric_fit.py`, `check_geometric_fit.log`: PASS):
- **Saddle demo exact forms (23 rows, 6 saddles).** Patched saddle centres equal `center_geometric` to ≤ 3.5e-8
  (original equals `center_reference`), with SDF(o) ≈ 0 and ∇g·n = **+1** (original −1). Non-saddle rows are
  unchanged bit for bit.
  - Rows 20–22 (a11 = 0) are off by 0.54 in both versions. That is the reference fit's degenerate eigenvector pair
    (`EigenRule::Reference`), as in the authors' C++.
- **4,096 random coefficient rows (2,631 saddles).** Patched centres equal the numpy one-liner
  `where(det II < 0, 2(p + a00 n) − c, c)` applied to the original output to 2.4e-7 (f32); axis and radii are
  bitwise unchanged.

**Launcher.** `runs_tmp/pat-train-jax-geo/run.sh` runs pretrain and phase 0 only (`PAT_PHASES`; 100,096 steps,
about 11 h), with the same driver, cache, memory setting, checkpoints and restarts. Output goes to
`runs_tmp/pat-train-jax-geo/{train.log,PIDS,dmon.log,work/}`. Start it only when the GPU is free:

```bash
nohup setsid /var/data/discrete-geometry/ddg-tools/runs_tmp/pat-train-jax-geo/run.sh > /dev/null 2>&1 < /dev/null &
```

Compare its `work/models/FundamentalFormPredictor_phase0.pkl` (converted as above) with the main run's phase-0
weights.
