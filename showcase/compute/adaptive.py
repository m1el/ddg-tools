"""Adaptive refinement of an intrinsic triangulation (Sec. 5.4): model problems."""
import sys, time
import numpy as np
import ddg_tools as d
from common import *

name, problem = sys.argv[1], sys.argv[2]
maxv = int(sys.argv[3]) if len(sys.argv) > 3 else 20000
source = int(sys.argv[4]) if len(sys.argv) > 4 else 0
pos, polys, _ = load_obj(model(name))
tri = triangles(polys)
out = {}
for key, mv in (("coarse", len(pos)), ("adaptive", maxv)):
    t = time.time()
    r = d.adaptive_solve(pos, tri, problem, source=source, max_vertices=mv)
    out[key] = np.asarray(r.values).reshape(len(pos), r.n_components)
    out[key + "_nverts"] = r.report["n_triangulation_vertices"]
    print(name, problem, key, len(pos), r.report, "%.1fs" % (time.time() - t))
save(f"adaptive_{name}_{problem}", positions=pos, triangles=tri, source=source, **out)
