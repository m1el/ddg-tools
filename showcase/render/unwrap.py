"""Conformal unwrap: the surface with a UV checker, seams (the cut) as tubes, cones as spheres,
and the morph from the surface to its flat layout.

blender -b --python render/unwrap.py -- <mesh> still3d <out.png>
blender -b --python render/unwrap.py -- <mesh> stilluv <out.png>
blender -b --python render/unwrap.py -- <mesh> morph
"""

import math
import os
import sys

sys.path.insert(0, __file__.rsplit("/", 1)[0])
import blib as B  # noqa: E402
import numpy as np  # noqa: E402

a = B.args()
name, mode = a[0], a[1]
D = B.load(f"unwrap_{name}")
P0 = B.y_up_to_z_up(D["positions"])
P, norm = B.normalize(P0, size=2.2)
fv, fs = D["poly_verts"], D["poly_sizes"]
UV = D["packed"]
cones_v, cones_a = D["cone_vertices"], D["cone_angles"]

AZ = float(os.environ.get("AZ", {"spot_quad": 135, "hand": -20, "bob_quad": -40}.get(name, -40)))
EL = float(os.environ.get("EL", 18))
CHECK = float(os.environ.get("CHECK", 48))

# ---- weld corners that share a vertex and a UV: vertices split only along the cut
key = np.round(UV * 1e6).astype(np.int64)
keys = np.stack([fv.astype(np.int64), key[:, 0], key[:, 1]], axis=1)
uniq, inv = np.unique(keys, axis=0, return_inverse=True)
inv = inv.ravel()
nw = len(uniq)
src = np.zeros(nw, dtype=np.int64)
src[inv] = fv
Uw = np.zeros((nw, 2))
Uw[inv] = UV
Pw = P[src]
polys_w = inv.astype(np.int32)

# Boundary edges of the welded mesh: the cut (and the mesh boundary).
starts = np.concatenate([[0], np.cumsum(fs)[:-1]])
nxt = np.arange(len(polys_w)) + 1
ends = starts + fs
nxt[ends - 1] = starts
e = np.stack([polys_w, polys_w[nxt]], axis=1)
es = np.sort(e, axis=1)
u_e, cnt = np.unique(es, axis=0, return_counts=True)
bnd = u_e[cnt == 1]

# ---- camera basis (same convention as blib.fit_camera)
az, el = math.radians(AZ), math.radians(EL)
dvec = np.array([math.sin(az) * math.cos(el), -math.cos(az) * math.cos(el), math.sin(el)])
fwd = -dvec
right = np.cross(fwd, [0, 0, 1.0])
right /= np.linalg.norm(right)
up = np.cross(right, fwd)


def poly_area3(X):
    tot = 0.0
    for s, n in zip(starts, fs):
        idx = polys_w[s:s + n]
        p0 = X[idx[0]]
        for k in range(1, n - 1):
            tot += 0.5 * np.linalg.norm(np.cross(X[idx[k]] - p0, X[idx[k + 1]] - p0))
    return tot


def poly_area2(Y):
    x, y = Y[:, 0], Y[:, 1]
    tot = 0.0
    for s, n in zip(starts, fs):
        idx = polys_w[s:s + n]
        tot += 0.5 * (x[idx] * np.roll(y[idx], -1) - np.roll(x[idx], -1) * y[idx]).sum()
    return tot


# Flat target: the layout in the plane facing the camera through the surface's centre, rotated
# in that plane to match the surface's projection (2D Procrustes, rotation only), scaled to the
# surface's area so checker squares keep their size.
cen = Pw.mean(0)
proj = np.stack([(Pw - cen) @ right, (Pw - cen) @ up], axis=1)
uvc = Uw - Uw.mean(0)
s = math.sqrt(poly_area3(Pw) / abs(poly_area2(Uw)))
H = (s * uvc).T @ proj
Uu, _, Vt = np.linalg.svd(H)
R = (Uu @ Vt).T
if np.linalg.det(R) < 0:
    Vt[-1] *= -1
    R = (Uu @ Vt).T
flat2 = (s * uvc) @ R.T
FLAT = cen + flat2[:, :1] * right + flat2[:, 1:2] * up
# Keep the flat layout above the ground.
ground = P[:, 2].min()
lift = max(0.0, ground + 0.02 - FLAT[:, 2].min())
FLAT[:, 2] += lift

morph = mode == "morph"
B.reset()
if morph:
    B.setup(1280, 800, samples=48)
elif mode == "stilluv":
    B.setup(900, 1200, samples=96)
else:
    B.setup(900, 1200, samples=160)

mat = B.mat_checker(scale=CHECK, c1=B.lin((0.96, 0.95, 0.93)), c2=B.lin((0.16, 0.42, 0.55)),
                    roughness=0.5)
obj = B.mesh_object("surface", Pw, polys_w, fs, smooth=True, material=mat)
B.set_uvs(obj, Uw[polys_w])
seam_mat = B.mat_solid(B.lin((0.85, 0.14, 0.10)), roughness=0.35)
cone_mat = B.mat_solid(B.lin((0.98, 0.72, 0.12)), roughness=0.3)
TUBE = 0.0065
CONE_R = 0.028
cone_w = np.array([np.where(src == v)[0][0] for v in cones_v], dtype=np.int64)
dyn = {"seam": None, "cones": None}


def boundary_lines(X):
    return [X[[i, j]] for i, j in bnd]


def place(t):
    t = B.ease(t)
    X = (1 - t) * Pw + t * FLAT
    obj.data.vertices.foreach_set("co", X.astype(np.float32).ravel())
    obj.data.update()
    for k in ("seam", "cones"):
        if dyn[k] is not None:
            B.bpy.data.objects.remove(dyn[k], do_unlink=True)
    dyn["seam"] = B.tubes("seam", boundary_lines(X), TUBE, seam_mat, resolution=2)
    dyn["cones"] = B.spheres("cones", X[cone_w], CONE_R, material=cone_mat) if len(cone_w) else None
    return X


center = (P.min(0) + P.max(0)) / 2
radius = float(np.linalg.norm(P.max(0) - P.min(0))) / 2

if mode == "still3d":
    place(0.0)
    B.studio(center, radius, ground_z=ground)
    B.fit_camera(Pw, azimuth=AZ, elevation=EL, lens=55, margin=1.05)
    B.render_still(a[2])
elif mode == "stilluv":
    place(1.0)
    B.studio(cen, radius, ground_z=None)
    ext = np.stack([(FLAT - cen) @ right, (FLAT - cen) @ up], axis=1)
    w, h = np.ptp(ext[:, 0]), np.ptp(ext[:, 1])
    # ortho_scale is the visible height (the longer side of a 900x1200 frame).
    scale = max(h, w / 0.75) * 1.08
    mid = cen + right * (ext[:, 0].max() + ext[:, 0].min()) / 2 + up * (ext[:, 1].max() + ext[:, 1].min()) / 2
    B.camera(tuple(mid + dvec * 10), tuple(mid), ortho_scale=scale)
    B.render_still(a[2])
else:
    B.studio(center, radius * 1.3, ground_z=ground)
    B.fit_camera(np.concatenate([Pw, FLAT]), azimuth=AZ, elevation=EL, lens=55, margin=1.12)
    hold0, go, hold1, back, hold2 = 20, 85, 40, 85, 20
    n = hold0 + go + hold1 + back + hold2

    def update(i):
        if i < hold0:
            t = 0.0
        elif i < hold0 + go:
            t = (i - hold0) / (go - 1)
        elif i < hold0 + go + hold1:
            t = 1.0
        elif i < hold0 + go + hold1 + back:
            t = 1 - (i - hold0 - go - hold1) / (back - 1)
        else:
            t = 0.0
        place(t)

    B.render_frames(f"unwrap_{name}", n, update)
