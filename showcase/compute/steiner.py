"""Steiner tree through all vertices (Sec. 5.1): intrinsic MST + Fermat points vs extrinsic MST."""
import sys
import numpy as np
import ddg_tools as d
from common import *

name = sys.argv[1]
pos, polys, _ = load_obj(model(name))
tri = triangles(polys)
out = {}
for key, kw in (("steiner", dict(intrinsic=True, steiner_points=True)),
                ("imst", dict(intrinsic=True, steiner_points=False)),
                ("emst", dict(intrinsic=False, steiner_points=False))):
    r = d.steiner_tree(pos, tri, **kw)
    out[f"{key}_points"] = np.asarray(r.points).reshape(-1, 3)
    out[f"{key}_offsets"] = np.asarray(r.offsets)
    out[f"{key}_length"] = r.report["length"]
    print(name, len(pos), key, r.report)
save(f"steiner_{name}", positions=pos, triangles=tri, **out)
