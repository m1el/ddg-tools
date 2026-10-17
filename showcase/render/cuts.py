"""Variational surface cutting: regions colored on the surface (crisp, from the per-corner
region functions: the argmin is taken per pixel of the interpolated functions, so the color
boundary is the piecewise-linear cut), the cut as dark tubes.

blender -b --python render/cuts.py -- still  <tag> <record> <out.png> [azimuth] [elevation]
blender -b --python render/cuts.py -- anim   <tag> <frames-name> [azimuth] [elevation] [turn-deg]
blender -b --python render/cuts.py -- unwrap <tag> <out.png> [azimuth] [elevation]
"""

import json
import os
import math
import sys

sys.path.insert(0, __file__.rsplit("/", 1)[0])
import blib as B  # noqa: E402
import numpy as np  # noqa: E402

REGION_COLORS = ["#f28e2b", "#4e79a7", "#59a14f", "#e15759", "#76b7b2", "#edc948", "#b07aa1",
                 "#ff9da7", "#9c755f", "#bab0ac"]
TWO_COLORS = ["#f28e2b", "#c8d2de"]


def region_colors(n):
    hexes = TWO_COLORS if n == 2 else REGION_COLORS[:n]
    return [B.lin(B._hex(h)) for h in hexes]


def mat_regions(n, roughness=0.5):
    """Base color = color of the region whose (interpolated) function f<k> is smallest."""
    mat = B.bpy.data.materials.new("regions")
    nt, bsdf = B._principled(mat)
    bsdf.inputs["Roughness"].default_value = roughness
    cols = region_colors(n)
    f0 = nt.nodes.new("ShaderNodeAttribute")
    f0.attribute_name = "f0"
    m = f0.outputs["Fac"]
    rgb = nt.nodes.new("ShaderNodeRGB")
    rgb.outputs[0].default_value = (*cols[0], 1)
    c = rgb.outputs[0]
    for k in range(1, n):
        fk = nt.nodes.new("ShaderNodeAttribute")
        fk.attribute_name = f"f{k}"
        lt = nt.nodes.new("ShaderNodeMath")
        lt.operation = "LESS_THAN"
        nt.links.new(fk.outputs["Fac"], lt.inputs[0])
        nt.links.new(m, lt.inputs[1])
        mn = nt.nodes.new("ShaderNodeMath")
        mn.operation = "MINIMUM"
        nt.links.new(fk.outputs["Fac"], mn.inputs[0])
        nt.links.new(m, mn.inputs[1])
        m = mn.outputs[0]
        ck = nt.nodes.new("ShaderNodeRGB")
        ck.outputs[0].default_value = (*cols[k], 1)
        mix = nt.nodes.new("ShaderNodeMix")
        mix.data_type = "RGBA"
        nt.links.new(lt.outputs[0], mix.inputs["Factor"])
        nt.links.new(c, mix.inputs["A"])
        nt.links.new(ck.outputs[0], mix.inputs["B"])
        c = mix.outputs["Result"]
    nt.links.new(c, bsdf.inputs["Base Color"])
    return mat


a = B.args()
mode, tag = a[0], a[1]
D = B.load(f"cuts_{tag}")
INFO = json.loads((B.WORK / f"cuts_{tag}.json").read_text())
N = INFO["n_regions"]
mesh_name = INFO["mesh"]

P0 = B.y_up_to_z_up(D["positions"])
P, NORM = B.normalize(P0, size=2.2)
CV, ST = D["corner_vertices"], D["polygon_starts"]
sizes = np.diff(np.concatenate([ST, [len(CV)]]))

DEFAULT_VIEW = {"spot": (-50, 18), "bob": (-40, 26)}.get(mesh_name, (-40, 20))


def cut_lines(rec):
    p = D["pts"][D["pt_off"][rec]:D["pt_off"][rec + 1]]
    s = D["starts"][D["st_off"][rec]:D["st_off"][rec + 1]]
    c = D["closed"][D["cl_off"][rec]:D["cl_off"][rec + 1]]
    p = B.apply_norm(B.y_up_to_z_up(p), NORM)
    lines = [p[s[i]:s[i + 1]] for i in range(len(s) - 1)]
    return lines, c


def scene(width, height, samples):
    B.reset()
    B.setup(width, height, samples=samples)


def add_surface(rec):
    obj = B.mesh_object("surface", P, CV, sizes, smooth=True, material=mat_regions(N))
    for k in range(N):
        obj.data.attributes.new(f"f{k}", "FLOAT", "CORNER")
    set_state(obj, rec)
    return obj


def set_state(obj, rec):
    st = D["states"][rec]
    for k in range(N):
        obj.data.attributes[f"f{k}"].data.foreach_set("value", st[k].astype(np.float32))
    obj.data.update()


INK = B.lin((0.07, 0.08, 0.10))
TUBE_R = 0.0075
CHECK = 56.0


def add_cut(rec, parent=None):
    lines, closed = cut_lines(rec)
    if not lines:
        return None
    obj = B.tubes("cut", lines, TUBE_R, B.mat_solid(INK, roughness=0.35), closed=closed)
    if parent is not None:
        obj.parent = parent
    return obj


def frame_scene(points, az, el, margin=1.05, lens=55):
    ground = points[:, 2].min()
    center = (points.min(0) + points.max(0)) / 2
    radius = float(np.linalg.norm(points.max(0) - points.min(0))) / 2
    B.studio(center, radius, ground_z=ground)
    return B.fit_camera(points, azimuth=az, elevation=el, lens=lens, margin=margin)


if mode == "still":
    rec = int(a[2]) if int(a[2]) >= 0 else len(D["states"]) - 1
    az = float(a[4]) if len(a) > 4 else DEFAULT_VIEW[0]
    el = float(a[5]) if len(a) > 5 else DEFAULT_VIEW[1]
    scene(1800, 1200, 160)
    add_surface(rec)
    add_cut(rec)
    frame_scene(P, az, el)
    B.render_still(a[3])

elif mode == "anim":
    frames_name = a[2]
    az = float(a[3]) if len(a) > 3 else DEFAULT_VIEW[0]
    el = float(a[4]) if len(a) > 4 else DEFAULT_VIEW[1]
    turn = float(a[5]) if len(a) > 5 else 60.0
    n_rec = len(D["states"])
    scene(1280, 800, 48)
    pivot = B.bpy.data.objects.new("pivot", None)
    B.bpy.context.collection.objects.link(pivot)
    surf = add_surface(0)
    surf.parent = pivot
    # Frame on the turning model: union of the rotated positions.
    rots = [B.rotz(P, t) for t in np.linspace(-turn / 2, turn / 2, 7)]
    frame_scene(np.concatenate(rots), az, el, margin=1.03)
    # Frame schedule: hold the start, then the steps with a slow-in (more frames on early,
    # fast steps), hold the end.
    hold0, hold1 = 18, 36
    n_flow = int(os.environ.get("CUTS_NFLOW", round(min(n_rec - 1, 200) * 1.6)))
    sched = [0] * hold0
    for i in range(n_flow):
        t = (i + 1) / n_flow
        sched.append(int(round((n_rec - 1) * (t ** 1.6))))
    sched += [n_rec - 1] * hold1
    n_frames = len(sched)
    labels = []
    E, L = D["energy"], D["cut_length"]
    for rec in sched:
        k = min(rec, len(E) - 1)
        labels.append(f"step {rec:3d}    cut length {L[k]:.2f}    energy {E[k]:.3g}")
    (B.WORK / f"labels_{frames_name}.json").write_text(json.dumps(labels))
    cut = [None]

    def update(i):
        rec = sched[i]
        set_state(surf, rec)
        if cut[0] is not None:
            B.bpy.data.objects.remove(cut[0], do_unlink=True)
        cut[0] = add_cut(rec, parent=pivot)
        pivot.rotation_euler[2] = math.radians(-turn / 2 + turn * B.ease(i / (n_frames - 1)))

    B.render_frames(frames_name, n_frames, update)

elif mode == "unwrap":
    az = float(a[3]) if len(a) > 3 else DEFAULT_VIEW[0]
    el = float(a[4]) if len(a) > 4 else DEFAULT_VIEW[1]
    UP = B.apply_norm(B.y_up_to_z_up(D["u_positions"]), NORM)
    ucv, ust = D["u_corner_vertices"], D["u_polygon_starts"]
    usz = np.diff(np.concatenate([ust, [len(ucv)]]))
    uv = np.load(B.WORK / f"cuts_{tag}_packed.npy")
    w = int(a[5]) if len(a) > 5 else 1800
    h = int(a[6]) if len(a) > 6 else 1200
    scene(w, h, 160)
    obj = B.mesh_object("unwrapped", UP, ucv, usz, smooth=True)
    B.set_uvs(obj, uv)
    reg = D["u_region"].astype(np.int64)
    cols = np.array(region_colors(N))
    pc = np.where(reg[:, None] >= 0, cols[np.clip(reg, 0, N - 1)], 0.8)
    pc = 0.35 + 0.65 * pc  # light squares: a pale tint of the region color
    B.set_face_colors(obj, pc, "Tint")
    obj.data.materials.append(B.mat_checker(scale=CHECK, c2=B.lin((0.16, 0.17, 0.2)),
                                            tint_attr="Tint", roughness=0.5))
    seams = D["u_seams"]
    lines = [UP[[i, j]] for i, j in seams]
    # Join seam edges into chains for clean tubes.
    sys.path.insert(0, str(B.SHOW / "compute"))
    from common import chain_edges  # noqa: E402
    chains = chain_edges(seams)
    B.tubes("seams", [UP[c] for c in chains], TUBE_R * 0.9, B.mat_solid(INK, roughness=0.35))
    frame_scene(UP, az, el, margin=1.12)
    B.render_still(a[2])

elif mode == "layout":
    # The packed UV layout as a flat mesh with the same checker, seen from above.
    size = int(a[3]) if len(a) > 3 else 1200
    uv = np.load(B.WORK / f"cuts_{tag}_packed.npy")
    ucv, ust = D["u_corner_vertices"], D["u_polygon_starts"]
    usz = np.diff(np.concatenate([ust, [len(ucv)]]))
    B.reset()
    B.setup(size, size, samples=64)
    # One vertex per corner: the layout is cut at the seams.
    pos = np.zeros((len(uv), 3))
    pos[:, :2] = uv
    obj = B.mesh_object("layout", pos, np.arange(len(uv)), usz, smooth=False)
    B.set_uvs(obj, uv)
    reg = D["u_region"].astype(np.int64)
    cols = np.array(region_colors(N))
    pc = np.where(reg[:, None] >= 0, cols[np.clip(reg, 0, N - 1)], 0.8)
    pc = 0.35 + 0.65 * pc
    B.set_face_colors(obj, pc, "Tint")
    obj.data.materials.append(B.mat_checker(scale=CHECK, c2=B.lin((0.16, 0.17, 0.2)),
                                            tint_attr="Tint", roughness=0.8))
    # Emission instead of the BSDF: the layout shows the colors exactly.
    mat = obj.data.materials[0]
    nt = mat.node_tree
    em = nt.nodes.new("ShaderNodeEmission")
    nt.links.new(nt.nodes["Checker Texture"].outputs["Color"], em.inputs["Color"])
    nt.links.new(em.outputs[0], nt.nodes["Material Output"].inputs["Surface"])
    # Chart outlines: sides whose two UV end points are not shared with another side.
    from collections import Counter
    q = np.round(uv * 2e5).astype(np.int64)
    keyof = lambda i, j: tuple(sorted([tuple(q[i]), tuple(q[j])]))
    cnt = Counter()
    sides = []
    for s0, n in zip(ust, usz):
        for k in range(n):
            i, j = s0 + k, s0 + (k + 1) % n
            key = keyof(i, j)
            cnt[key] += 1
            sides.append((key, i, j))
    lines = [pos[[i, j]] + [0, 0, 0.0005] for key, i, j in sides if cnt[key] == 1]
    B.tubes("outline", lines, 0.0012, B.mat_emission(INK, 1.0), resolution=1)
    lo, hi = uv.min(0), uv.max(0)
    cx, cy = (lo + hi) / 2
    B.camera((cx, cy, 5.0), (cx, cy, 0.0), ortho_scale=float(max(hi - lo)) * 1.06)
    sc = B.bpy.context.scene
    B.area_light("top", (0.5, 0.3, 4.0), (0.5, 0.5, 0), 6.0, 450)
    sc.world.node_tree.nodes["Background"].inputs["Strength"].default_value = 0.55
    B.render_still(a[2])
