"""Conformal map of spot to the sphere (BFF)."""
import numpy as np
import ddg_tools as d
from common import *

pos, polys, tex = load_obj(model("spot"))
tri = triangles(polys)
r = d.map_to_sphere(pos, tri)
S = np.asarray(r.positions).reshape(-1, 3)
uv, puv = tex
cuv = np.array([uv[i] for q in puv for i in q])  # per corner of the triangles (spot is triangles)
print({k: r.report[k] for k in ("is_injective", "quasi_conformal_mean", "quasi_conformal_max", "flipped_triangles", "seconds")})
save("sphere_spot", positions=pos, triangles=tri, sphere=S, corner_uv=cuv)
