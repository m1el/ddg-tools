"""Smoothest direction field: a glyph per vertex (cross for n = 4, line for n = 2, arrow for
n = 1) in the tangent plane, singularities as spheres (red: positive index, blue: negative).

blender -b --python render/field.py -- <mesh> <n> <out.png> [w h]
"""

import os
import sys

sys.path.insert(0, __file__.rsplit("/", 1)[0])
import blib as B  # noqa: E402
import numpy as np  # noqa: E402

a = B.args()
name, n, out = a[0], int(a[1]), a[2]
W, H = (int(a[3]), int(a[4])) if len(a) > 4 else (1800, 1200)
D = B.load(f"field_{name}_{n}")
P0 = B.y_up_to_z_up(D["positions"])
ROT = int(os.environ.get("ROT", {"spot": 0, "bob": 0, "blub": 0}.get(name, 0)))
P0 = B.rotz(P0, ROT)
P, norm = B.normalize(P0, size=2.2)
s = norm[0]
T = D["triangles"]
V = B.rotz(B.y_up_to_z_up(D["vectors"]), ROT)
SP = B.apply_norm(B.rotz(B.y_up_to_z_up(D["sing_points"]), ROT), norm)
SI = D["sing_index"]

B.reset()
B.setup(W, H, samples=160)
surf = B.mesh_object("surface", P, T, smooth=True,
                     material=B.mat_solid(B.lin((0.90, 0.88, 0.84)), roughness=0.5))

# Vertex normals.
fn = np.cross(P[T[:, 1]] - P[T[:, 0]], P[T[:, 2]] - P[T[:, 0]])
N = np.zeros_like(P)
for k in range(3):
    np.add.at(N, T[:, k], fn)
N /= np.linalg.norm(N, axis=1, keepdims=True)
# Field in the vertex tangent plane, unit length.
V = V - (V * N).sum(1, keepdims=True) * N
ok = np.linalg.norm(V, axis=1) > 1e-9
V[ok] /= np.linalg.norm(V[ok], axis=1, keepdims=True)

E = np.concatenate([T[:, [0, 1]], T[:, [1, 2]], T[:, [2, 0]]])
h = np.mean(np.linalg.norm(P[E[:, 0]] - P[E[:, 1]], axis=1))
L = float(os.environ.get("GLYPH", 0.42)) * h
r = 0.055 * h
lift = 1.2 * r
ink = B.mat_solid(B.lin((0.10, 0.16, 0.30)), roughness=0.4)
SPACING = float(os.environ.get("SPACING", 0))
if SPACING > 0:
    # Greedy subsample: glyph sites at least SPACING * h apart.
    from scipy.spatial import cKDTree
    kd = cKDTree(P)
    taken = np.zeros(len(P), bool)
    blocked = np.zeros(len(P), bool)
    for i in np.random.default_rng(1).permutation(len(P)):
        if blocked[i] or not ok[i]:
            continue
        taken[i] = True
        blocked[kd.query_ball_point(P[i], SPACING * h)] = True
    ok = ok & taken
C = P[ok] + lift * N[ok]
Vo, No = V[ok], N[ok]
segs = []
if n in (2, 4):
    segs += [np.stack([C - L * Vo, C + L * Vo], 1)]
if n == 4:
    Wv = np.cross(No, Vo)
    segs += [np.stack([C - L * Wv, C + L * Wv], 1)]
if segs:
    B.tubes("glyphs", list(np.concatenate(segs)), r, ink, resolution=1)
else:
    B.arrows("glyphs", C - 0.9 * L * Vo, Vo, 1.8 * L, r, ink, head=0.4)

colors = np.where((SI > 0)[:, None], [B.lin((0.86, 0.16, 0.12))], [B.lin((0.12, 0.40, 0.85))])
B.spheres("sing", SP, 0.9 * h * float(os.environ.get("SING", 1.0)), colors,
          material=B.mat_attr("Col", 0.25))

ground = P[:, 2].min()
center = (P.min(0) + P.max(0)) / 2
radius = float(np.linalg.norm(P.max(0) - P.min(0))) / 2
B.studio(center, radius, ground_z=ground)
AZ = int(os.environ.get("AZ", {"spot": -40, "bob": -40, "blub": -50}.get(name, -35)))
EL = int(os.environ.get("EL", 22))
cam, ctr, dist = B.fit_camera(P, azimuth=AZ, elevation=EL, lens=55,
                              margin=float(os.environ.get("MARGIN", 1.06)))
cam.data.shift_y = float(os.environ.get("SHIFT", 0.0))
B.render_still(out)
