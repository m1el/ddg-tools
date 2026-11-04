#!/bin/bash
# Geometric sign-rule probe (NOT launched yet; the GPU runs one job at a time): the authors' curriculum
# up to and including phase 0 (pretrain 4,096 + phase 0 96,000 steps, ~11 h on the RTX 4080), with the one-line
# Geometric centre fix in the JAX torus fit (PATCHES.diff). Same driver, resume and restart logic as
# ../pat-train-jax/run.sh; launch with
#   nohup setsid /var/data/discrete-geometry/ddg-tools/runs_tmp/pat-train-jax-geo/run.sh > /dev/null 2>&1 &
G=/var/data/discrete-geometry/ddg-tools/runs_tmp/pat-train-jax-geo
D=/var/data/discrete-geometry/ddg-tools/runs_tmp/pat-train-jax
PY=/var/data/discrete-geometry/PointsAsTori-build/.venv/bin/python
mkdir -p $G/work
ln -sfn /var/data/discrete-geometry/PointsAsTori/PAT_training_data/PAT_training_data $G/work/datasets
export PAT_WORK=$G/work PAT_SRC=$G/src PAT_PHASES=FundamentalFormPredictor_pretrained,FundamentalFormPredictor_phase0
export PAT_CKPT_SECONDS=900 PAT_LOG_EVERY=100 MPLBACKEND=Agg XLA_PYTHON_CLIENT_MEM_FRACTION=0.97
echo "wrapper $$" >> $G/PIDS
nvidia-smi dmon -s pucm -d 60 -c 1000 > $G/dmon.log 2>&1 &
echo "dmon $!" >> $G/PIDS
for i in $(seq 1 20); do
  echo "=== start $i at $(date -Is)" >> $G/train.log
  $PY $D/driver.py >> $G/train.log 2>&1 &
  pid=$!
  echo "python $pid (start $i, $(date -Is))" >> $G/PIDS
  wait $pid
  rc=$?
  echo "=== exit $rc at $(date -Is)" >> $G/train.log
  [ $rc -eq 0 ] && break
  sleep 120
done
