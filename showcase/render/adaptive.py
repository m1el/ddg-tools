"""Adaptive refinement of an intrinsic triangulation (Navigating Intrinsic Triangulations,
Sec. 5.4), through ddg_tools.adaptive_solve.

blender -b --python render/adaptive.py -- green <out.png>   Green's function on spot, isolines
blender -b --python render/adaptive.py -- disk <out.png> <view>  face with the polar grid of its
                                                                  harmonic map to the unit disk
"""

import math
import sys

sys.path.insert(0, __file__.rsplit("/", 1)[0])
import blib as B  # noqa: E402
import numpy as np  # noqa: E402

a = B.args()
mode, out = a[0], a[1]
B.reset()
B.setup(1800, 1200, samples=192)

if mode == "green":
    D = B.load("adaptive_spot_green")
    P, _ = B.normalize(B.rotz(B.y_up_to_z_up(D["positions"]), 180), size=2.2)
    T = D["triangles"]
    g = D["adaptive"][:, 0]
    src = int(D["source"])
    lo, hi = np.percentile(g, 1), np.percentile(g, 99.7)
    obj = B.mesh_object("spot", P, T, smooth=True,
                        material=B.mat_attr("Col", roughness=0.45, lines="g", lines_period=(hi - lo) / 22,
                                            lines_width=0.05, lines_strength=0.5))
    B.set_point_colors(obj, B.cmap(g, "magma", lo - 0.9 * (hi - lo), hi + 0.05 * (hi - lo)))
    B.set_point_scalar(obj, g, "g")
    B.spheres("source", P[src][None], 0.03, material=B.mat_solid(B.lin((0.2, 0.75, 0.95)), 0.3))
    center = (P.min(0) + P.max(0)) / 2
    radius = float(np.linalg.norm(P.max(0) - P.min(0))) / 2
    B.studio(center, radius, ground_z=0.0)
    B.fit_camera(P, azimuth=-60, elevation=16, lens=55, margin=1.12,
                 target=center - np.array([0, 0, 0.06]))
else:
    D = B.load("adaptive_face_disk")
    P0 = D["positions"]
    T = D["triangles"]
    w = D["adaptive"]
    view = a[2] if len(a) > 2 else "a"
    # The face's own frame: z out of the face (the input's +z), y up.
    P, _ = B.normalize(B.y_up_to_z_up(P0), size=2.2)
    r = np.linalg.norm(w, axis=1)
    th = np.arctan2(w[:, 1], w[:, 0])
    # Polar grid: rings in log radius would be conformal; the map is harmonic, so plain rings.
    cu = th[T] / (2 * math.pi)
    ref = cu[:, :1]
    cu = np.where(cu - ref > 0.5, cu - 1, np.where(cu - ref < -0.5, cu + 1, cu))
    uv = np.stack([cu.ravel() * 3.0, r[T].ravel() * 1.0], axis=1)  # 24 sectors, 8 rings with scale 8
    face = B.mesh_object("face", P, T, smooth=True)
    B.set_uvs(face, uv)
    mat = B.mat_checker(scale=8.0, c1=B.lin((0.96, 0.95, 0.92)), c2=B.lin((0.14, 0.42, 0.62)), roughness=0.5)
    face.data.materials.append(mat)
    # The disk itself, standing next to the face.
    h = P[:, 2].max()
    R = 0.5 * h
    off = np.array([P[:, 0].max() + R + 0.25, 0.0, R])
    disk_pos = off + np.stack([w[:, 0] * R, np.zeros(len(w)), w[:, 1] * R], axis=1)
    disk = B.mesh_object("disk", disk_pos, T[:, ::-1], smooth=False)
    B.set_uvs(disk, uv.reshape(-1, 3, 2)[:, ::-1].reshape(-1, 2))
    disk.data.materials.append(mat)
    allp = np.vstack([P, disk_pos])
    center = (allp.min(0) + allp.max(0)) / 2
    radius = float(np.linalg.norm(allp.max(0) - allp.min(0))) / 2
    B.studio(center, radius, ground_z=0.0)
    B.fit_camera(allp, azimuth=float(view) if view != "a" else -15, elevation=8, lens=60, margin=1.06)
B.render_still(out)
