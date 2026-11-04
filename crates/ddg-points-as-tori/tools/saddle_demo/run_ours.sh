#!/usr/bin/env bash
# Our crate on every case, one after the other, 8 threads, nice 19, CPU only.
R=/var/data/discrete-geometry/ddg-tools
export RAYON_NUM_THREADS=8 CUDA_VISIBLE_DEVICES=
for c in $(ls $R/runs_tmp/pat-saddle/cases); do
  nice -n 19 $R/target-pat-saddle/release/examples/saddle_demo case $R/runs_tmp/pat-saddle $c --threads 8 || echo "FAILED $c"
done
