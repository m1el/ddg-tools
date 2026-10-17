"""Heat-method geodesic distance: the same mesh twice, distance on the input triangulation
(left) and on the intrinsic Delaunay refinement (right), as color bands with isolines.

blender -b --python render/geodesic.py -- <mesh> <input|delaunay|refined> still <out.png> <w> <h>
blender -b --python render/geodesic.py -- <mesh> <input|delaunay|refined> anim <w> <h>
"""

import json
import math
import sys

sys.path.insert(0, __file__.rsplit("/", 1)[0])
import blib as B  # noqa: E402
import numpy as np  # noqa: E402

a = B.args()
name, which, mode = a[0], a[1], a[2]
D = B.load(f"geodesic_{name}")
import os
ORIENT = os.environ.get("ORIENT", "yup")
P0 = B.y_up_to_z_up(D["positions"]) if ORIENT == "yup" else np.asarray(D["positions"]) * [1, 1, 1] if ORIENT == "raw" else -np.asarray(D["positions"])
ROT = int(os.environ.get("ROT", 0))
AZ, EL = int(os.environ.get("AZ", -14)), int(os.environ.get("EL", 22))
P0 = B.rotz(P0, ROT)
P, norm = B.normalize(P0, size=2.0)
T = D["triangles"]
src = int(D["source"])
DMAX = float(max(D["input"].max(), D["refined"].max()))

anim = mode == "anim"
B.reset()
W, H = (int(a[3]), int(a[4])) if anim else (int(a[4]), int(a[5]))
B.setup(W, H, samples=96 if anim else 160)

ext = P.max(0) - P.min(0)
gap = 0.12 * ext[0]
dx = ext[0] / 2 + gap / 2
period = DMAX / 18


def make(dist, x):
    Q = P + [x, 0, 0]
    o = B.mesh_object("surface", Q, T, smooth=True,
                      material=B.mat_attr("Col", roughness=0.45, lines="d", lines_period=period,
                                          lines_width=0.12, lines_strength=0.8))
    o.data.color_attributes.new("Col", "FLOAT_COLOR", "POINT")
    B.set_point_scalar(o, dist, "d")
    return o, Q


objs = [make(D[which], 0.0)]
dists = [D[which]]
src_mat = B.mat_solid(B.lin((0.95, 0.85, 0.2)), roughness=0.3)
for o, Q in objs:
    B.spheres("source", Q[src:src + 1], 0.022, material=src_mat)


def paint(front):
    """Color by distance up to `front`; beyond it the plain surface."""
    for (o, _), dist in zip(objs, dists):
        col = B.cmap(dist, "spectral", 0, DMAX)[:, ::1]
        col = B.cmap(DMAX - dist, "spectral", 0, DMAX)
        grey = np.array(B.lin((0.82, 0.82, 0.80)))
        w = np.clip((front - dist) / (0.02 * DMAX), 0, 1)[:, None]
        col = w * col + (1 - w) * grey
        rgba = np.ones((len(dist), 4), dtype=np.float32)
        rgba[:, :3] = col
        o.data.color_attributes["Col"].data.foreach_set("color", rgba.ravel())
        o.data.update()


allp = np.vstack([q for _, q in objs])
ground = allp[:, 2].min()
center = (allp.min(0) + allp.max(0)) / 2
radius = float(np.linalg.norm(allp.max(0) - allp.min(0))) / 2
B.studio(center, radius, ground_z=ground)
cam, ctr, dist = B.fit_camera(allp, azimuth=AZ, elevation=EL, lens=60, margin=1.25)
cam.data.shift_y = 0.03

if not anim:
    paint(10 * DMAX)
    B.render_still(a[3])
else:
    n = 150
    fronts = [DMAX * 1.02 * B.ease(i / (n - 30)) for i in range(n)]
    B.render_frames(f"geodesic_{name}_{which}", n, lambda i: paint(fronts[i] if fronts[i] > 0 else -1))
