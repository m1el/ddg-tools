#!/bin/bash
# The authors' Points-as-Tori curriculum (training/train.py) on the GPU, via driver.py + the checkpoint/resume
# patch (PATCHES.diff); crates/ddg-points-as-tori/TRAINING-JAX.md. Restarts after a failure: the driver skips
# finished phases (checkpoints/<phase>.done) and resumes the current phase from checkpoints/<phase>.npz.
D=/var/data/discrete-geometry/ddg-tools/runs_tmp/pat-train-jax
PY=/var/data/discrete-geometry/PointsAsTori-build/.venv/bin/python
export PAT_WORK=$D/work PAT_CKPT_SECONDS=900 PAT_LOG_EVERY=100 MPLBACKEND=Agg
# 16 GB RTX 4080: the phase-0..4 train step needs 14.97 GiB after XLA rematerialization; JAX's default pool
# (75 %) is 12 GB. 0.97 is the smallest tested fraction that fits (allocator setting only; no code or math change).
export XLA_PYTHON_CLIENT_MEM_FRACTION=0.97
echo "wrapper $$" >> $D/PIDS
for i in $(seq 1 20); do
  echo "=== start $i at $(date -Is)" >> $D/train.log
  $PY $D/driver.py >> $D/train.log 2>&1 &
  pid=$!
  echo "python $pid (start $i, $(date -Is))" >> $D/PIDS
  wait $pid
  rc=$?
  echo "=== exit $rc at $(date -Is)" >> $D/train.log
  [ $rc -eq 0 ] && break
  sleep 120
done
