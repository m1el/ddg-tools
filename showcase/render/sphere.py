"""Conformal map to the sphere (Boundary First Flattening, ddg_tools.map_to_sphere).

blender -b --python render/sphere.py -- anim            spot <-> sphere, with spot's texture
blender -b --python render/sphere.py -- conformal <out> spot | sphere, Mercator checker of the
                                                         sphere pulled back onto spot
"""

import math
import sys

sys.path.insert(0, __file__.rsplit("/", 1)[0])
import blib as B  # noqa: E402
import numpy as np  # noqa: E402

a = B.args()
mode = a[0]
D = B.load("sphere_spot")
P, norm = B.normalize(B.rotz(B.y_up_to_z_up(D["positions"]), 180), size=2.2)
T = D["triangles"]
S = B.rotz(B.y_up_to_z_up(D["sphere"]), 180)
CUV = D["corner_uv"]

# Radius: the sphere of the same area as the cow; sitting on the ground under the cow's centre.
e1, e2 = P[T[:, 1]] - P[T[:, 0]], P[T[:, 2]] - P[T[:, 0]]
area = 0.5 * np.linalg.norm(np.cross(e1, e2), axis=1).sum()
R = math.sqrt(area / (4 * math.pi)) * 1.05
ctr = np.array([(P[:, 0].min() + P[:, 0].max()) / 2, (P[:, 1].min() + P[:, 1].max()) / 2, R])

# The map is defined up to a rotation: take the one that moves the cow least (Procrustes on
# the directions from the centre).
dirs = P - P.mean(0)
dirs /= np.linalg.norm(dirs, axis=1, keepdims=True)
U_, _, Vt = np.linalg.svd(S.T @ dirs)
Rm = (U_ @ Vt).T
if np.linalg.det(Rm) < 0:
    U_[:, -1] *= -1
    Rm = (U_ @ Vt).T
S = S @ Rm.T
Q = ctr + R * S


def mercator_uv(s, axis=2):
    """Per-corner conformal (Mercator) coordinates of unit vectors s, axis = pole."""
    ax = [(axis + 1) % 3, (axis + 2) % 3, axis]
    x, y, z = s[:, ax[0]], s[:, ax[1]], np.clip(s[:, ax[2]], -1 + 1e-12, 1 - 1e-12)
    lon = np.arctan2(y, x) / (2 * math.pi)
    v = np.arctanh(z) / (2 * math.pi)
    cu = lon[T].copy()  # (m, 3)
    # unwrap the longitude seam per triangle
    ref = cu[:, :1]
    pole_c = np.abs(z[T]) > 1 - 1e-6
    ref = np.where(pole_c[:, :1], cu[:, 1:2], ref)
    cu = np.where(cu - ref > 0.5, cu - 1, np.where(cu - ref < -0.5, cu + 1, cu))
    # A corner at the pole has no longitude: take the mean of the other two.
    others = (cu.sum(1, keepdims=True) - cu) / 2
    cu = np.where(pole_c, others, cu)
    return np.stack([cu.ravel(), v[T].ravel()], axis=1)


if mode == "anim":
    B.reset()
    B.setup(1280, 800, samples=64)
else:
    B.reset()
    B.setup(1800, 1200, samples=192)

tex_mat = B.mat_image(B.SHOW.parent / "data/models/spot/spot/spot_texture.png", roughness=0.5)

if mode == "anim":
    obj = B.mesh_object("spot", P, T, smooth=True, material=tex_mat)
    B.set_uvs(obj, CUV)
    allp = np.vstack([P, Q])
    rad = float(np.linalg.norm(allp.max(0) - allp.min(0))) / 2
    B.studio((allp.min(0) + allp.max(0)) / 2, rad, ground_z=0.0)
    # Frame the union over every turn angle: use a cylinder of the xy extent.
    rr = np.max(np.linalg.norm(allp[:, :2] - ctr[:2], axis=1))
    ring = np.array([[ctr[0] + rr * math.cos(t), ctr[1] + rr * math.sin(t), z]
                     for t in np.linspace(0, 2 * math.pi, 48) for z in (0.0, allp[:, 2].max())])
    B.fit_camera(ring, azimuth=-30, elevation=18, lens=50, margin=0.9)
    pivot = B.bpy.data.objects.new("pivot", None)
    B.bpy.context.collection.objects.link(pivot)
    pivot.location = (ctr[0], ctr[1], 0)
    obj.location = (-ctr[0], -ctr[1], 0)
    obj.parent = pivot
    hold0, morph, hold1 = 24, 66, 30
    seq = [0.0] * hold0 + [B.ease(i / (morph - 1)) for i in range(morph)] + [1.0] * hold1 + \
          [B.ease(1 - i / (morph - 1)) for i in range(morph)]
    n = len(seq)
    verts = obj.data.vertices

    def update(i):
        t = seq[i]
        # Blend along the straight line, pushed outwards a little in the middle so the cow
        # does not pass through itself.
        X = (1 - t) * P + t * Q
        bulge = 4 * t * (1 - t) * 0.12 * R
        d = X - ctr
        d /= np.maximum(np.linalg.norm(d, axis=1, keepdims=True), 1e-9)
        X = X + bulge * d
        verts.foreach_set("co", X.astype(np.float32).ravel())
        obj.data.update()
        pivot.rotation_euler[2] = 2 * math.pi * i / n

    B.render_frames("sphere_morph", n, update)
else:
    out = a[1]
    gap = 0.35
    dx = (P[:, 0].max() - ctr[0]) + gap + R
    cow = B.mesh_object("spot", P, T, smooth=True)
    sph = B.mesh_object("sphere", Q + np.array([dx + 0.15, 0, 0]), T, smooth=True)
    # Grid: pole along the cow's side axis (x) so the rings show on the flank.
    pole = a[2] if len(a) > 2 else "head"
    if pole == "head":
        # Pole at the sphere image of the snout tip: the rings centre on the face.
        p = S[np.argmin(P[:, 1])]
        z = np.array([0, 0, 1.0])
        v = np.cross(p, z)
        c = float(p @ z)
        vx = np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]])
        Rp = np.eye(3) + vx + vx @ vx / (1 + c)
        uv = mercator_uv(S @ Rp.T, axis=2)
    else:
        uv = mercator_uv(S, axis=int(pole))
    mat = B.mat_checker(scale=float(a[3]) if len(a) > 3 else 14.0, c1=B.lin((0.96, 0.95, 0.92)),
                        c2=B.lin((0.86, 0.33, 0.18)), roughness=0.45)
    for o in (cow, sph):
        B.set_uvs(o, uv)
        o.data.materials.append(mat)
    allp = np.vstack([P, Q + np.array([dx + 0.15, 0, 0])])
    rad = float(np.linalg.norm(allp.max(0) - allp.min(0))) / 2
    B.studio((allp.min(0) + allp.max(0)) / 2, rad, ground_z=0.0)
    B.fit_camera(allp, azimuth=float(a[4]) if len(a) > 4 else -22, elevation=14, lens=60, margin=1.04)
    B.render_still(out)
