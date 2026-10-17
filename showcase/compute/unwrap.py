"""Conformal UV unwrap through the ddg-uv pipeline (optimal cones, cut, BFF), then packing."""
import sys, time
import numpy as np
import ddg_tools as d
from common import *

name = sys.argv[1]
lam = float(sys.argv[2]) if len(sys.argv) > 2 else 0.1
pos, polys, _ = load_obj(model(name))
cv, st, tr = corners(polys)
t = time.time()
r = d.unwrap(pos, cv, st, tr, lambda_=lam)
dt = time.time() - t
uv = np.frombuffer(r.uvs, dtype=np.float32).reshape(-1, 2).astype(np.float64)
charts = np.asarray(r.polygon_charts)
packed, prep = d.pack_charts(r.uvs, st, r.polygon_charts, margin=0.01)
packed = np.frombuffer(packed, dtype=np.float32).reshape(-1, 2).astype(np.float64)
fv, fs = flat_polys(polys)
rep = r.report
save(f"unwrap_{name}", positions=pos, poly_verts=fv, poly_sizes=fs, corner_vertices=cv,
     polygon_starts=st, uv=uv, packed=packed, charts=charts,
     cone_vertices=np.asarray(r.cone_vertices), cone_angles=np.asarray(r.cone_angles),
     cut_edges=np.asarray(r.cut_edges).reshape(-1, 2), seconds=dt)
print(name, "verts", len(pos), "polys", len(polys), "cones", len(r.cone_vertices),
      "charts", len(set(charts.tolist())), "cut edges", len(r.cut_edges) // 2, f"{dt:.2f}s")
for k, v in rep.items():
    if not isinstance(v, (list, dict)):
        print(" ", k, v)
print("pack", prep)
