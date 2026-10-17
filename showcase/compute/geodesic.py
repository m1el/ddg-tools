"""Heat-method geodesic distance (ddg_tools.geodesic_distance) on the input mesh, its
intrinsic Delaunay triangulation and its intrinsic Delaunay refinement.

geodesic.py <mesh> [source-vertex]"""
import sys, time
import numpy as np
import ddg_tools as d
from common import *

name = sys.argv[1]
pos, polys, _ = load_obj(model(name))
tri = triangles(polys)
P = np.stack([pos[:, 0], -pos[:, 2], pos[:, 1]], 1)  # y-up -> z-up
src = int(sys.argv[2]) if len(sys.argv) > 2 else int(np.argmax(P[:, 2]))
out = dict(positions=pos, triangles=tri, source=src)
for t in ["input", "delaunay", "refined"]:
    t0 = time.time()
    r = d.geodesic_distance(pos, tri, np.array([src], dtype=np.int32), triangulation=t)
    out[t] = np.asarray(r.distances)
    print(t, "%.3fs" % (time.time() - t0), r.report, "max d %.4f" % out[t].max())
save(f"geodesic_{name}", **out)
