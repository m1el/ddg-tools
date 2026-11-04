#!/bin/bash
# Resume test with XLA's deterministic GPU ops (test only; the real run uses default flags, like the authors).
# phase 0: two straight 80-step runs (determinism) and 40 + resume to 80 (mid-epoch cursor);
# pretraining: straight 192 vs 128 + resume to 192 (cursor at the end of an epoch's shapes).
D=/var/data/discrete-geometry/ddg-tools/runs_tmp/pat-train-jax
T=$D/timing/resume
PY=/var/data/discrete-geometry/PointsAsTori-build/.venv/bin/python
export XLA_PYTHON_CLIENT_MEM_FRACTION=0.97 MPLBACKEND=Agg PAT_LOG_EVERY=2 XLA_FLAGS=--xla_gpu_deterministic_ops=true
DATA=/var/data/discrete-geometry/PointsAsTori/PAT_training_data/PAT_training_data
mkdir -p $T
for d in p0-a p0-b p0-c pre-a pre-c; do
  rm -rf $T/$d; mkdir -p $T/$d/models $T/$d/checkpoints; ln -s $DATA $T/$d/datasets
done
for d in p0-a p0-b p0-c; do
  cp $DATA/FundamentalFormPredictor_pretrained.pkl $T/$d/models/
  echo test > $T/$d/checkpoints/FundamentalFormPredictor_pretrained.done
done
run() { PAT_WORK=$T/$1 PAT_PHASES=$2 PAT_MAX_STEPS=$3 $PY $D/driver.py >> $T/$1/train.log 2>&1; echo "$1 $3 rc=$?"; }
run p0-a FundamentalFormPredictor_phase0 80
run p0-b FundamentalFormPredictor_phase0 80
run p0-c FundamentalFormPredictor_phase0 40
run p0-c FundamentalFormPredictor_phase0 80
run pre-a FundamentalFormPredictor_pretrained 192
run pre-c FundamentalFormPredictor_pretrained 128
run pre-c FundamentalFormPredictor_pretrained 192
C="$PY -I $D/compare_runs.py"
K0=checkpoints/FundamentalFormPredictor_phase0.npz
K1=checkpoints/FundamentalFormPredictor_pretrained.npz
echo "== phase0 straight a vs straight b"; $C $T/p0-a/train.log $T/p0-a/$K0 $T/p0-b/train.log $T/p0-b/$K0
echo "== phase0 straight a vs 40+resume c"; $C $T/p0-a/train.log $T/p0-a/$K0 $T/p0-c/train.log $T/p0-c/$K0
echo "== pretrain straight a vs 128+resume c"; $C $T/pre-a/train.log $T/pre-a/$K1 $T/pre-c/train.log $T/pre-c/$K1
