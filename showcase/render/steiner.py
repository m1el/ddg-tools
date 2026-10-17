"""Steiner tree through all vertices (Navigating Intrinsic Triangulations, Sec. 5.1).

blender -b --python render/steiner.py -- <mesh> <steiner|imst|emst> <out.png> [azimuth]
"""

import sys

sys.path.insert(0, __file__.rsplit("/", 1)[0])
import blib as B  # noqa: E402
import numpy as np  # noqa: E402

a = B.args()
name, key, out = a[0], a[1], a[2]
az = float(a[3]) if len(a) > 3 else -30
D = B.load(f"steiner_{name}")
P, norm = B.normalize(B.rotz(B.y_up_to_z_up(D["positions"]), 180), size=2.2)
T = D["triangles"]
pts = B.apply_norm(B.rotz(B.y_up_to_z_up(D[f"{key}_points"]), 180), norm)
off = D[f"{key}_offsets"]
lines = [pts[off[k]:off[k + 1]] for k in range(len(off) - 1)]

# Original polygons for display (the control mesh is quads).
sys.path.insert(0, str(B.SHOW / "compute"))
from common import load_obj, model  # noqa: E402
_, polys, _ = load_obj(model(name))
flat = np.array([v for p in polys for v in p])
sizes = np.array([len(p) for p in polys])

B.reset()
B.setup(1200, 1200, samples=160)
surf = B.mesh_object("surface", P, flat, sizes, smooth=False,
                     material=B.mat_solid(B.lin((0.80, 0.82, 0.86)), roughness=0.55))
# The mesh's edges, thin and dark.
edges = set()
for p in polys:
    for u, v in zip(p, p[1:] + p[:1]):
        edges.add((min(u, v), max(u, v)))
B.tubes("edges", [P[[u, v]] for u, v in edges], 0.0035, B.mat_solid(B.lin((0.35, 0.37, 0.42)), 0.6))
color = (0.95, 0.60, 0.08) if key == "steiner" else (0.20, 0.45, 0.85)
B.tubes("tree", lines, 0.011, B.mat_solid(B.lin(color), roughness=0.3), resolution=5)
B.spheres("verts", P, 0.016, material=B.mat_solid(B.lin((0.12, 0.13, 0.16)), 0.4))
# Steiner points: polyline end points that are not mesh vertices.
ends = np.array([l[0] for l in lines] + [l[-1] for l in lines])
dmin = np.min(np.linalg.norm(ends[:, None, :] - P[None, :, :], axis=2), axis=1)
st = ends[dmin > 1e-6]
if len(st):
    st = np.unique(np.round(st, 9), axis=0)
    B.spheres("steiner", st, 0.022, material=B.mat_solid(B.lin((0.85, 0.15, 0.12)), 0.3))
    print("steiner points", len(st))

center = (P.min(0) + P.max(0)) / 2
radius = float(np.linalg.norm(P.max(0) - P.min(0))) / 2
B.studio(center, radius, ground_z=0.0)
B.fit_camera(P, azimuth=az, elevation=22, lens=55, margin=1.16, target=(P.min(0) + P.max(0)) / 2 - [0, 0, 0.08])
B.render_still(out)
