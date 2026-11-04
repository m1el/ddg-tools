"""Run the authors' Points as Tori program (black box) on WP14's evaluation clouds (SPEC.md §16.7).

Not a uv project script: it runs in the venv of the reference build (see `runs_tmp/pat-ref/README.md`;
the venv was created by `uv venv` and has the authors' package `pointsastori` and its nanobind
extension). CPU JAX only (the GPU is reserved for training):

    JAX_PLATFORMS=cpu OMP_NUM_THREADS=4 \
    XLA_FLAGS="--xla_cpu_multi_thread_eigen=false intra_op_parallelism_threads=4" \
    uv run --no-project --python /var/data/discrete-geometry/PointsAsTori-build/.venv/bin/python -- \
        python -I crates/ddg-points-as-tori/scripts/pat_eval_ref.py <list-of-case-dirs.txt>

For every case directory `<case>` (with `points.npy`, `normals.npy`, float64 (n, 3)) it writes
`<case>/ref/`:

* `grid_sdf_accel.npy`, `grid_sdf_noaccel.npy` — `signed_distance` on the 64³ grid
  `meshgrid(linspace(-1, 1, 64), ..., indexing='ij')`, C order (`accelerate=True`: 32-NN; `False`: all
  points), as the golden driver `runs_tmp/pat-ref/driver/pat_golden.py`;
* `tori_centers.npy`, `tori_axes_angles.npy`, `tori_major_radius.npy`, `tori_minor_radius.npy` —
  from the authors' public constructor `PointsAsTori(points, normals)`;
* `net_coeffs.npy` (float32, (n, 6)) and `knn_indices.npy` (int32, (n, 64)) — `model(P[idx], N[idx])`
  with `idx = shape_3d.get_neighbors(P, 64, [])`, as the golden driver;
* `timing.json` — wall times of a warm second constructor call and of the two grid evaluations,
  with the thread settings (shared machine: times carry the outside load).

Cases that already have `ref/timing.json` are skipped.
"""

import json
import os
import sys
import time

import numpy as np
from flax import nnx

import pointsastori as pat
from pointsastori import infer, shape_3d
from pointsastori.network import FundamentalFormPredictor


def main() -> int:
    if len(sys.argv) != 2:
        print(__doc__)
        return 2
    cases = [l.strip() for l in open(sys.argv[1]) if l.strip()]
    ax = np.linspace(-1.0, 1.0, 64)
    G = np.stack(np.meshgrid(ax, ax, ax, indexing='ij'), -1).reshape(-1, 3)
    model, k_nb = FundamentalFormPredictor.load_saved_model(infer.DEFAULT_MODEL_PATH)
    nnx.eval_mode(model)
    for case in cases:
        out = os.path.join(case, 'ref')
        if os.path.isfile(os.path.join(out, 'timing.json')):
            continue
        try:
            P = np.load(os.path.join(case, 'points.npy'))
            N = np.load(os.path.join(case, 'normals.npy'))
            if P.ndim != 2 or P.shape[1] != 3 or P.shape != N.shape or not np.isfinite(P).all() or not np.isfinite(N).all():
                raise ValueError(f'bad cloud arrays {P.shape} {N.shape}')
            os.makedirs(out, exist_ok=True)
            sv = lambda k, v: np.save(os.path.join(out, f'{k}.npy'), np.asarray(v))
            t0 = time.perf_counter()
            obj = pat.PointsAsTori(P, N)
            first = time.perf_counter() - t0
            t0 = time.perf_counter()
            obj = pat.PointsAsTori(P, N)
            warm = time.perf_counter() - t0
            centers, axes, major, minor = obj._tdf.get_tori()
            sv('tori_centers', centers)
            sv('tori_axes_angles', axes)
            sv('tori_major_radius', major)
            sv('tori_minor_radius', minor)
            idx = shape_3d.get_neighbors(P, k_nb, np.array([]))
            sv('knn_indices', idx.astype(np.int32))
            sv('net_coeffs', np.array(model(P[idx], N[idx])))
            t0 = time.perf_counter()
            sv('grid_sdf_accel', obj.signed_distance(G, accelerate=True).reshape(64, 64, 64))
            t_acc = time.perf_counter() - t0
            t0 = time.perf_counter()
            sv('grid_sdf_noaccel', obj.signed_distance(G, accelerate=False).reshape(64, 64, 64))
            t_all = time.perf_counter() - t0
            timing = dict(n=int(len(P)), constructor_first_s=first, constructor_warm_s=warm, grid64_accel_s=t_acc, grid64_noaccel_s=t_all,
                          omp_threads=os.environ.get('OMP_NUM_THREADS', ''), xla_flags=os.environ.get('XLA_FLAGS', ''),
                          jax_platforms=os.environ.get('JAX_PLATFORMS', ''))
            with open(os.path.join(out, 'timing.json'), 'w') as fh:
                json.dump(timing, fh)
            print(f'{case}: n={len(P)} warm constructor {warm:.2f}s grid accel {t_acc:.2f}s noaccel {t_all:.2f}s', flush=True)
        except Exception as e:  # report and continue with the next case
            print(f'{case}: FAILED {type(e).__name__}: {e}', flush=True)
    return 0


if __name__ == '__main__':
    sys.exit(main())
