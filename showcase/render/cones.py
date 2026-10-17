"""Optimal cones: the surface colored by the log conformal scale factor u, cones as spheres
(red: positive angle, a cone tip; blue: negative, a saddle).

blender -b --python render/cones.py -- <mesh> still <run-index> <out.png>
blender -b --python render/cones.py -- <mesh> anim
"""

import json
import math
import sys

sys.path.insert(0, __file__.rsplit("/", 1)[0])
import blib as B  # noqa: E402
import numpy as np  # noqa: E402

a = B.args()
name, mode = a[0], a[1]
D = B.load(f"cones_{name}")
P0 = B.y_up_to_z_up(D["positions"])
P, norm = B.normalize(P0, size=2.2)
T = D["triangles"]
lams, U = D["lambdas"], D["u"]
off, CV, CA = D["cone_offsets"], D["cone_vertices"], D["cone_angles"]
n_runs = len(lams)

AZ = {"bob": float(__import__("os").environ.get("AZ", -40)), "camel": -60, "bunny": -20, "spot": -35}.get(name, -35)
EL = float(__import__("os").environ.get("EL", 28))

B._CMAPS["brbg"] = ["#01665e", "#35978f", "#80cdc1", "#c7eae5", "#e8e6e3", "#f6e8c3",
                    "#dfc27d", "#bf812d", "#8c510a"]
anim = mode == "anim"
B.reset()
if anim:
    B.setup(1280, 800, samples=48)
else:
    B.setup(1800, 1200, samples=160)

# Scale of the colormap: fixed over the sweep, from a mid run, symmetric about 0.
ref = n_runs // 2
UMAX = float(np.percentile(np.abs(U[ref] - np.mean(U[ref])), 98))

surf = B.mesh_object("surface", P, T, smooth=True,
                     material=B.mat_attr("Col", roughness=0.42, lines="u", lines_period=UMAX / 4,
                                         lines_width=0.07, lines_strength=0.35))
surf.data.color_attributes.new("Col", "FLOAT_COLOR", "POINT")
surf.data.attributes.new("u", "FLOAT", "POINT")

cone_mat = B.mat_attr("Col", roughness=0.25)
cone_obj = [None]
R0 = 0.022


def show_run(k, blend=None):
    """Show run k (or a blend of runs k and k+1 for u)."""
    u = U[k] if blend is None else (1 - blend) * U[k] + blend * U[min(k + 1, n_runs - 1)]
    u = u - np.mean(u)
    col = B.cmap(u, "brbg", -UMAX, UMAX)
    rgba = np.ones((len(P), 4), dtype=np.float32)
    rgba[:, :3] = col
    surf.data.color_attributes["Col"].data.foreach_set("color", rgba.ravel())
    surf.data.attributes["u"].data.foreach_set("value", u.astype(np.float32))
    surf.data.update()
    if cone_obj[0] is not None:
        B.bpy.data.objects.remove(cone_obj[0], do_unlink=True)
        cone_obj[0] = None
    v, ang = CV[off[k]:off[k + 1]], CA[off[k]:off[k + 1]]
    if len(v):
        r = R0 * np.clip(np.sqrt(np.abs(ang) / (math.pi / 2)), 0.45, 2.2)
        # Lift the sphere centre a little along the normal so it sits on the surface.
        c = P[v]
        colors = np.where((ang > 0)[:, None], [B.lin((0.86, 0.16, 0.12))], [B.lin((0.12, 0.32, 0.85))])
        cone_obj[0] = B.spheres("cones", c, r, colors, material=cone_mat)
    return len(v)


ground = P[:, 2].min()
center = (P.min(0) + P.max(0)) / 2
radius = float(np.linalg.norm(P.max(0) - P.min(0))) / 2
B.studio(center, radius, ground_z=ground)
if anim:
    # Frame every pose of the +-25 degree turntable, around the fixed centre.
    poses = np.concatenate([B.rotz(P - [center[0], center[1], 0], d) + [center[0], center[1], 0]
                            for d in range(-25, 26, 5)])
    cam, ctr, dist = B.fit_camera(poses, azimuth=AZ, elevation=EL, lens=55, margin=1.05,
                                  target=center)
else:
    cam, ctr, dist = B.fit_camera(P, azimuth=AZ, elevation=EL, lens=55, margin=1.06)

pivot = B.bpy.data.objects.new("pivot", None)
B.bpy.context.collection.objects.link(pivot)
pivot.location = (center[0], center[1], 0)
surf.parent = pivot
surf.location = (-center[0], -center[1], 0)

if not anim:
    k = int(a[2])
    UMAX = float(np.percentile(np.abs(U[k] - np.mean(U[k])), 98))
    n = show_run(k)
    if cone_obj[0] is not None:
        cone_obj[0].parent = pivot
        cone_obj[0].location = (-center[0], -center[1], 0)
    B.render_still(a[3])
    print(f"LABEL λ = {lams[k]:.4f} · {n} cones")
else:
    # Forward through the sweep, then back; two frames per run, u blended between runs.
    seq = [(k, 0.0) for k in range(n_runs)] + [(n_runs - 1, 0.0)] * 12 + \
          [(k, 0.0) for k in range(n_runs - 1, -1, -1)] + [(0, 0.0)] * 12
    frames = []
    for k, _ in seq:
        frames += [k, k]
    labels = []

    def update(i):
        k = frames[i]
        n = show_run(k)
        if cone_obj[0] is not None:
            cone_obj[0].parent = pivot
            cone_obj[0].location = (-center[0], -center[1], 0)
        pivot.rotation_euler[2] = math.radians(25 * math.sin(2 * math.pi * i / len(frames)))
        return n

    # Labels for the post-processing step.
    for i in range(len(frames)):
        k = frames[i]
        labels.append(f"λ = {lams[k]:.4f} · {off[k + 1] - off[k]} cones")
    (B.WORK / f"labels_cones_{name}.json").write_text(json.dumps(labels))
    B.render_frames(f"cones_{name}", len(frames), update)
