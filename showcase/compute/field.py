"""Smoothest n-direction field (ddg_tools.smoothest_direction_field) on the intrinsic
Delaunay triangulation.

field.py <mesh> <n>"""
import sys, time
import numpy as np
import ddg_tools as d
from common import *

name, n = sys.argv[1], int(sys.argv[2])
pos, polys, _ = load_obj(model(name))
tri = triangles(polys)
t0 = time.time()
r = d.smoothest_direction_field(pos, tri, n=n, triangulation="delaunay")
print(name, n, "%.3fs" % (time.time() - t0), r.report)
idx = np.asarray(r.singularity_indices)
print("singularities", len(idx), "indices", sorted(set(idx.tolist())), "sum", idx.sum())
save(f"field_{name}_{n}", positions=pos, triangles=tri, n=n,
     vectors=np.asarray(r.vectors).reshape(-1, 3),
     sing_points=np.asarray(r.singularity_points).reshape(-1, 3), sing_index=idx)
