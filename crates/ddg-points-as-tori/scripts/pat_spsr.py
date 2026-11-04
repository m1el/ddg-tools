# /// script
# requires-python = ">=3.12,<3.13"
# dependencies = ["open3d", "numpy"]
# ///
"""Screened Poisson surface reconstruction (Kazhdan and Hoppe 2013) through Open3D, as a black box
baseline for Points as Tori (SPEC.md §9, paper Sec. 5.2).

The paper runs Open3D's wrapper single-threaded and limits the octree depth to the resolution of the
evaluation grid; for a 64^3 grid that is depth 6 (2^6 = 64), the default here (`--grid 64`). The
reconstructed mesh is written as `vertices.npy` (float64, (V, 3)) and `faces.npy` (int32, (F, 3)), with
`densities.npy` and `meta.json`; the crate's `truth::MeshTruth` (exact distance, winding-number sign,
`require_closed: false`) then gives its signed distance on the grid.

Run with uv (never pip):

    uv run crates/ddg-points-as-tori/scripts/pat_spsr.py --cloud runs_tmp/pat-ref/inputs/block --out runs_tmp/pat-wp5/spsr/block

`--cloud DIR` reads `DIR/points.npy` and `DIR/normals.npy` (float64, (n, 3), unit outward normals);
`--points/--normals` name the files directly. Nothing is trimmed (no density threshold).
"""

import argparse
import json
import math
import os
import sys
import time


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--cloud", help="directory with points.npy and normals.npy")
    ap.add_argument("--points", help="points.npy (n, 3)")
    ap.add_argument("--normals", help="normals.npy (n, 3)")
    ap.add_argument("--out", required=True, help="output directory")
    ap.add_argument("--grid", type=int, default=64, help="evaluation grid resolution; depth = ceil(log2(grid)) (default 64 -> 6)")
    ap.add_argument("--depth", type=int, default=None, help="octree depth (overrides --grid)")
    ap.add_argument("--scale", type=float, default=1.1, help="Open3D scale: reconstruction cube / bounding cube (default 1.1)")
    ap.add_argument("--threads", type=int, default=1, help="Open3D threads (paper: 1)")
    args = ap.parse_args()

    import numpy as np
    import open3d as o3d

    if args.cloud:
        pts_path = os.path.join(args.cloud, "points.npy")
        nrm_path = os.path.join(args.cloud, "normals.npy")
    elif args.points and args.normals:
        pts_path, nrm_path = args.points, args.normals
    else:
        print("error: give --cloud DIR or --points and --normals", file=sys.stderr)
        return 2
    p = np.load(pts_path).astype(np.float64)
    n = np.load(nrm_path).astype(np.float64)
    if p.ndim != 2 or p.shape[1] != 3 or n.shape != p.shape:
        print(f"error: points {p.shape} and normals {n.shape} must both be (n, 3)", file=sys.stderr)
        return 2
    if not (np.isfinite(p).all() and np.isfinite(n).all()):
        print("error: non-finite points or normals", file=sys.stderr)
        return 2
    depth = args.depth if args.depth is not None else int(math.ceil(math.log2(args.grid)))

    pcd = o3d.geometry.PointCloud()
    pcd.points = o3d.utility.Vector3dVector(p)
    pcd.normals = o3d.utility.Vector3dVector(n)
    t0 = time.perf_counter()
    mesh, dens = o3d.geometry.TriangleMesh.create_from_point_cloud_poisson(
        pcd, depth=depth, width=0, scale=args.scale, linear_fit=False, n_threads=args.threads
    )
    secs = time.perf_counter() - t0
    v = np.asarray(mesh.vertices, dtype=np.float64)
    f = np.asarray(mesh.triangles, dtype=np.int32)
    os.makedirs(args.out, exist_ok=True)
    np.save(os.path.join(args.out, "vertices.npy"), v)
    np.save(os.path.join(args.out, "faces.npy"), f)
    np.save(os.path.join(args.out, "densities.npy"), np.asarray(dens, dtype=np.float64))
    meta = dict(
        input_points=pts_path,
        input_normals=nrm_path,
        n_points=int(len(p)),
        depth=depth,
        width=0,
        scale=args.scale,
        linear_fit=False,
        threads=args.threads,
        vertices=int(len(v)),
        faces=int(len(f)),
        seconds=secs,
        open3d=o3d.__version__,
        python=sys.version.split()[0],
    )
    with open(os.path.join(args.out, "meta.json"), "w") as fh:
        json.dump(meta, fh, indent=1)
    print(json.dumps(meta))
    return 0


if __name__ == "__main__":
    sys.exit(main())
