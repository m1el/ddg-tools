# Retraining Points as Tori

The network of Points as Tori (Feng, Gkioulekas, Crane 2026) retrained from the authors' training data with this
crate's port of their training code (`src/train/`, MIT, Copyright (c) 2026 Nicole Feng; notice in
`LICENSES/points-as-tori-MIT.txt`), on the GPU (burn-cuda, TF32 matmuls). Design and measurements:
`NOTES-WP7.md` (training core) and `NOTES-WP10.md` (GPU training, determinism, timings). Contract: `SPEC.md` §10,
§18.

## The run

| | |
|---|---|
| data | `/var/data/discrete-geometry/PointsAsTori/PAT_training_data/PAT_training_data/` (read only) |
| curriculum | `training/curriculum.toml`: pretrain, phase0–phase4, 497,858 steps (table below) |
| fit rules | `SignRule::Reference` (as the authors trained), `EigenRule::Robust` (SPEC §18 WP7 item 6) |
| data rules | normals normalized at load (D23), the network on valid rows with an attention pad mask (D19), all-padding slices skipped |
| schedule | the code's (`Schedule::Code`: linear warmup, then constant; D21), Adam β 0.9/0.999, ε 1e-8 |
| seed | 42 (data order, queries, initialization; D22) |
| device | CUDA 0 (RTX 4080), TF32 matmuls (user decision 2026-11-04), network batch 256 neighborhoods (taped passes), 1024 (untaped forward) |
| checkpoints | every 2,000 steps (at the next shape boundary), at every epoch end and every phase end |
| validation | `ABC17_validation_3.npz` with unit normals, at epoch ends once 1,000 steps passed since the last one, and at every phase end |
| output | `/var/data/discrete-geometry/ddg-tools/runs_tmp/pat-train/` |
| code | commit `801b3c9` (binary frozen as `runs_tmp/pat-train/bin/train-801b3c9`); steps 1–576 ran with `7527474` (the same math without `--forward-batch`) |

| phase | files | epochs | steps | lr | rows per step |
|---|---|---|---|---|---|
| pretrain | `pretraining_blob.npz` | 64 | 4,096 | 1e-3 | 32 |
| phase0 | nice blobs | 24 | 96,000 | 1e-5 | 1024 |
| phase1 | blobs + ABC17, interleaved | 8 | 64,000 | 1e-5 | 1024 |
| phase2 | blobs, ABC17, ABC2 | 8 | 119,344 | 1e-6 | 1024 |
| phase3 | blobs, all nice ABC | 5 | 124,590 | 1e-6 | 1024 |
| phase4 | everything + simulated scans | 2 | 89,828 | 1e-7 | 1024 |

### Exact command

Started 2026-11-04 05:03 UTC; stopped twice for short A/B measurements of the batch shapes and restarted from
its checkpoints (05:06 from step 64; 05:16 from step 576 with the `--forward-batch 1024` binary); detached, at normal priority, through `runs_tmp/pat-train/run.sh`, which restarts the trainer after a failure
(up to 30 times, 2 minutes apart; every start resumes from the last checkpoint):

```bash
cd /var/data/discrete-geometry/ddg-tools
export CARGO_TARGET_DIR=/var/data/discrete-geometry/ddg-tools/target-pat10
cargo build --release -p ddg-points-as-tori --features train,cuda --example train
OUT=/var/data/discrete-geometry/ddg-tools/runs_tmp/pat-train
mkdir -p $OUT/bin && cp target-pat10/release/examples/train $OUT/bin/train-$(git rev-parse --short HEAD)
# run.sh loops over:
#   $OUT/bin/train-801b3c9 run --out $OUT --device cuda:0 --seed 42 --sign reference --eigen robust \
#       --schedule code --checkpoint-every 2000 --validate-every 1000 --forward-batch 1024 >> $OUT/train.log 2>&1
cd $OUT && setsid nohup ./run.sh > /dev/null 2>&1 < /dev/null &
setsid nohup timeout 5d nvidia-smi dmon -s um -d 10 -o T > $OUT/dmon.log 2>&1 < /dev/null &
```

PIDs are in `runs_tmp/pat-train/PIDS` (wrapper 144488, trainer 144491, dmon 7298 at the time of writing).
Before launching, the GPU was checked for other users (`nvidia-smi --query-compute-apps=pid,process_name,used_memory
--format=csv`): the run was started only when no `blender` process of the render agent was on the GPU.

## Monitoring

- `runs_tmp/pat-train/train.log`: human-readable progress every 30 s (step, phase, epoch, shape, loss and its
  two terms, learning rate, pairs, network rows, seconds of forward/loss/backward/update, a naive ETA), every
  validation, `=== start` / `=== exit` lines of the wrapper.
- `runs_tmp/pat-train/log.jsonl`: one JSON line per step (`phase, phase_name, epoch, step, shape, row, lr, loss,
  distance, eikonal, queries, pairs, rows_fwd, rows_bwd, t_prep, t_plan, t_fwd, t_loss, t_bwd, t_upd,
  skipped`) and one per validation (`val_distance`, `val_eikonal`). After a resume the steps after the
  checkpoint appear again (bitwise identical values, measured); later lines supersede earlier ones.
- `runs_tmp/pat-train/dmon.log`: `nvidia-smi dmon -s um` every 10 s (SM and memory-bandwidth %, frame-buffer MB).
- Quick looks:

```bash
cd /var/data/discrete-geometry/ddg-tools/runs_tmp/pat-train
tail -3 train.log; grep -E "validation|=== |SKIPPED" train.log | tail
kill -0 $(awk '{print $5}' PIDS | tr -d ,) && echo trainer alive     # or: pgrep -fx "/bin/bash ./run.sh"
python3 - <<'EOF'
import json
L=[json.loads(l) for l in open('log.jsonl')]
S=[r for r in L if 'val_distance' not in r]; V=[r for r in L if 'val_distance' in r]
last=S[-500:]; t=sum(r['t_fwd']+r['t_loss']+r['t_bwd']+r['t_upd'] for r in last)/len(last)
print(f"step {S[-1]['step']} {S[-1]['phase_name']} e{S[-1]['epoch']}; mean loss (last 500) "
      f"{sum(r['loss'] for r in last)/len(last):.5f}; {t:.3f} s/step; ETA {(497858-S[-1]['step'])*t/3600:.1f} h")
for v in V[-5:]: print(v['phase_name'], v['epoch'], v['step'], v['val_distance'], v['val_eikonal'])
EOF
```

A skipped step (non-finite loss or gradient) is logged with `"skipped":true` and counted; it does not stop the
run.

## Stopping and resuming

Stop the wrapper first, then the trainer (otherwise the wrapper restarts it):

```bash
kill <wrapper pid>; kill <trainer pid>        # see PIDS; wait with: while kill -0 <trainer pid>; do sleep 5; done
```

Resume by starting `run.sh` again (or the `train run …` command above with the same `--out`): the trainer reads
`checkpoint.safetensors` (weights, with the cursor `{phase, epoch, done, step, adam_count, adam_time}` in its
metadata) and `checkpoint-adam.safetensors` (Adam moments, keyed by parameter name), and continues at the next
shape after the cursor. The data order and the queries are keyed by (seed, phase, epoch, shape), so the resumed
run is bitwise identical to an uninterrupted one on the same GPU and binary (tested: 128 straight steps vs
64 + 64 resumed, all losses, weights and moments bitwise equal; and the live run's replayed steps after its
restart). At most the work since the last checkpoint (≤ 2,000 steps ≈ 20 min) is repeated.

`--network-batch 512` is ~15 % faster while the host is loaded (fewer kernel launches) but uses 12.5 GB of GPU
memory instead of 8.2 GB; switch to it (stop, restart with the option added to `run.sh`) only when no other
process needs the GPU. It changes the bits from then on (round-off), not the math.

## Results and export

Files in `runs_tmp/pat-train/`: `phase-<name>.safetensors` (weights at the end of each phase), `checkpoint*.safetensors`
(latest), `final.safetensors` (after the last phase). All are in the crate's weights layout
(`Model::from_safetensors` reads them). Export the final (or any phase's) weights to `data/` (gitignored; SPEC §18
item 6: they stay there until they match or beat the authors' weights on §16.2):

```bash
cd /var/data/discrete-geometry/ddg-tools
./runs_tmp/pat-train/bin/train-801b3c9 export --from runs_tmp/pat-train/final.safetensors \
    --to data/points-as-tori/retrained-2026-11-seed42.safetensors
# or a phase: --from runs_tmp/pat-train/phase-phase2.safetensors; or the latest checkpoint: --from runs_tmp/pat-train
```

The export keeps the training cursor and a licence line in the metadata and checks the round trip. Evaluate with
`GpuMatmul::Split3` or the CPU (f32-accurate inference) for the comparison tables against the authors' weights
(WP14, `EVALUATION.md`).

Validation numbers are with unit normals (D23). The authors' validation printout used the stored norm-2 normals
and is not comparable (NOTES-WP7: 0.99 / 10.6 vs 0.0195 / 0.296 for their weights with unit normals).

## Timing and expected end

Measured in the run: pretraining 0.6–0.8 s per step with the host at load 75–110 (other agents' jobs at nice 19;
the GPU at 22–31 % SM, bound by burn's per-op host dispatch, NOTES-WP10 §5). On a quiet host the GPU time is
0.41 s per step.

**Expected end** (497,858 steps): see the status line below, updated from the run's own step rate.
