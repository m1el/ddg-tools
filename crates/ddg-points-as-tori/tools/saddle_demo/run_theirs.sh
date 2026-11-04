#!/usr/bin/env bash
# Run the authors' program (their venv) on CPU only, 8 threads, nice 19.
export JAX_PLATFORMS=cpu CUDA_VISIBLE_DEVICES= OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8 MKL_NUM_THREADS=8
export XLA_FLAGS="--xla_cpu_multi_thread_eigen=false intra_op_parallelism_threads=8"
exec nice -n 19 /var/data/discrete-geometry/PointsAsTori-build/.venv/bin/python -I \
  /var/data/discrete-geometry/ddg-tools/runs_tmp/pat-saddle/driver/theirs.py /var/data/discrete-geometry/ddg-tools/runs_tmp/pat-saddle "$@"
