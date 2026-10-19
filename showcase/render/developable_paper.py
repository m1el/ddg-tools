"""Developability scenes in the style of the paper's figures (Figs. 1, 22, 23): flat shading, so every
triangle is one tone and the developable pieces read as smooth sheets, and the seams (curves where the
surface turns) as thin tubes. Data: work/dp_<tag>.npz from compute/dev_paper.py (in_/out_ arrays and
the seam curves), work/dp_snap_<tag>.npz (states), work/dp_cut_<tag>.npz (charts).

blender -b --python render/developable_paper.py -- panel <tag> <in|out> <out.png> az el [color] [W H]
    color: plain | defect (angle defect, fixed scale) ; env DEV_NOSEAMS=1 draws no seam tubes
    The camera is fitted to the input and the result together, so the panels of a model line up.
blender -b --python render/developable_paper.py -- flow <snap-tag> <frames-name> az el [turn-deg] [per-state frames]
blender -b --python render/developable_paper.py -- chart3d <cut-tag> <out.png> az el [W H]
blender -b --python render/developable_paper.py -- layout <cut-tag> <out.png> [size]
"""
import json
import math
import os
import sys

sys.path.insert(0, __file__.rsplit("/", 1)[0])
import blib as B  # noqa: E402
import numpy as np  # noqa: E402
from developable import defect_colors, to_scene  # noqa: E402

sys.path.insert(0, str(B.SHOW / "compute"))
from dev_paper import seam_curves  # noqa: E402

bpy = B.bpy
SURF = (0.58, 0.65, 0.76)          # sRGB of the plain surface
SEAM = B.lin((0.86, 0.14, 0.10))
SEAM_R = float(os.environ.get("DEV_SEAM_R", 0.0048))   # tube radius, scene bounding-box diagonal 2
DEFECT_MAX = float(os.environ.get("DEV_DEFECT_MAX", 0.05))   # radians
CHART_COLORS = ["#4e79a7", "#f28e2b", "#59a14f", "#e15759", "#b07aa1", "#76b7b2", "#edc948", "#ff9da7",
                "#9c755f", "#86bc86", "#a0a8c4", "#d4a373"]


def flat_surface(name, P, F, color=None, rough=0.6):
    """Flat-shaded mesh; `color` (n, 3) linear point colors or None for the plain surface."""
    mat = B.mat_attr("Col", rough) if color is not None else B.mat_solid(B.lin(SURF), rough, name=name)
    obj = B.mesh_object(name, P, F, smooth=False, material=mat)
    if color is not None:
        B.set_point_colors(obj, color)
    return obj


def defect_surface_colors(defect):
    """Angle defect on a fixed scale (blue negative, red positive), fading to the plain surface
    color where it is small, so that the defect reads as spots and lines on the sheets."""
    c = defect_colors(defect, DEFECT_MAX)
    w = np.clip(np.abs(defect) / (0.25 * DEFECT_MAX), 0, 1)[:, None] ** 0.7
    return (1 - w) * np.array(B.lin(SURF)) + w * c


def seam_tubes(D, prefix, norm, name="seams", radius=None, key="", color=None):
    pts, sizes = D[f"{key}seam_pts"], D[f"{key}seam_sizes"]
    if not len(sizes):
        return None
    P = to_scene(pts, norm)
    lines, k = [], 0
    for s in sizes:
        lines.append(P[k:k + s])
        k += s
    return B.tubes(name, lines, radius or SEAM_R, B.mat_solid(color or SEAM, roughness=0.4), resolution=3)


def frame_scene(points, az, el, margin=1.06, lens=55):
    """Studio and camera as in the other scenes, with a smaller key light and a weak fill so that the
    flat triangles of a developable sheet differ in tone with their orientation."""
    ground = points[:, 2].min()
    center = (points.min(0) + points.max(0)) / 2
    radius = float(np.linalg.norm(points.max(0) - points.min(0))) / 2
    B.studio(center, radius, ground_z=ground, key=1.6)
    bpy.data.lights["key"].size *= 0.5
    bpy.data.lights["fill"].energy *= 0.3
    bpy.data.lights["rim"].energy *= 0.5
    bpy.context.scene.world.node_tree.nodes["Background"].inputs["Strength"].default_value = 0.12
    return B.fit_camera(points, azimuth=az, elevation=el, lens=lens, margin=margin)


def scene_norm(D):
    _, norm = B.normalize(B.y_up_to_z_up(D["in_P"]), size=2.0)
    return norm


if __name__ == "__main__":
    a = B.args()
    mode = a[0]
    if mode == "panel":
        tag, which, out = a[1], a[2], a[3]
        az, el = float(a[4]), float(a[5])
        kind = a[6] if len(a) > 6 else "plain"
        W = int(a[7]) if len(a) > 7 else 1000
        H = int(a[8]) if len(a) > 8 else 1000
        D = B.load(f"dp_{tag}")
        norm = scene_norm(D)
        B.reset()
        B.setup(W, H, samples=128)
        P = to_scene(D[f"{which}_P"], norm)
        col = defect_surface_colors(D["out_defect"]) if (kind == "defect" and which == "out") else None
        flat_surface("mesh", P, D[f"{which}_F"], col)
        if which == "out" and not os.environ.get("DEV_NOSEAMS"):
            seam_tubes(D, which, norm, color=B.lin((0.10, 0.11, 0.14)) if kind == "defect" else None,
                       radius=SEAM_R * (0.7 if kind == "defect" else 1))
        both = np.concatenate([to_scene(D["in_P"], norm), to_scene(D["out_P"], norm)])
        frame_scene(both, az, el, margin=1.04, lens=65)
        B.render_still(out)

    elif mode == "flow":
        snap, frames_name = a[1], a[2]
        az, el = float(a[3]), float(a[4])
        turn = float(a[5]) if len(a) > 5 else 40.0
        per = int(a[6]) if len(a) > 6 else 4
        S = dict(np.load(B.WORK / f"dp_snap_{snap}.npz"))
        meta = json.loads((B.WORK / f"dp_snap_{snap}.json").read_text())["states"]
        n = int(S["n_states"])
        # the input is state 0 of the schedule: loaded from the same npz of the full run
        D = B.load(f"dp_{os.environ.get('DEV_FULL', snap)}")
        _, norm = B.normalize(B.y_up_to_z_up(D["in_P"]), size=2.0)
        states = [dict(P=to_scene(D["in_P"], norm), F=D["in_F"], pts=np.zeros((0, 3)), sizes=np.zeros(0, int),
                       label="input", energy=None, round=-1, iteration=0)]
        for i in range(n):
            m = meta[i]
            states.append(dict(P=to_scene(S[f"s{i}_P"], norm), F=S[f"s{i}_F"], pts=S[f"s{i}_seam_pts"],
                               sizes=S[f"s{i}_seam_sizes"], energy=m["energy"], round=m["round"],
                               iteration=m["iteration"]))
        B.reset()
        B.setup(1100, 800, samples=48)
        pivot = bpy.data.objects.new("pivot", None)
        bpy.context.collection.objects.link(pivot)
        allP = np.concatenate([s["P"] for s in states])
        rots = [B.rotz(allP, t) for t in np.linspace(-turn / 2, turn / 2, 5)]
        frame_scene(np.concatenate(rots), az, el, margin=1.04, lens=60)
        sched = [(0, 0.0)] * 14
        for i in range(len(states) - 1):
            for j in range(per):
                sched.append((i, j / per))
        sched += [(len(states) - 1, 0.0)] * 40
        nfr = len(sched)
        labels = []
        for (i, t) in sched:
            s = states[i]
            if s["round"] < 0:
                labels.append(f"input    {len(s['F'])} faces")
            else:
                e = "E_P" if s["round"] == 0 else "E_lambda"
                labels.append(f"round {s['round'] + 1}/3 ({e})  iteration {s['iteration']:3d}"
                              f"    {len(s['F'])} faces    energy {s['energy']:.3g}")
        (B.WORK / f"labels_{frames_name}.json").write_text(json.dumps(labels))
        cur = []

        def update(f):
            i, t = sched[f]
            j = min(i + 1, len(states) - 1)
            s, s2 = states[i], states[j]
            same = t > 0 and s["F"].shape == s2["F"].shape and np.array_equal(s["F"], s2["F"])
            w = B.ease(t) if same else 0.0
            P = (1 - w) * s["P"] + w * s2["P"] if same else s["P"]
            for o in cur:
                bpy.data.objects.remove(o, do_unlink=True)
            cur.clear()
            obj = flat_surface("mesh", P, s["F"])
            obj.parent = pivot
            cur.append(obj)
            curves, _, _ = seam_curves(P, s["F"])
            if curves:
                tb = B.tubes("seams", curves, SEAM_R, B.mat_solid(SEAM, roughness=0.4), resolution=3)
                tb.parent = pivot
                cur.append(tb)
            pivot.rotation_euler[2] = math.radians(-turn / 2 + turn * B.ease(f / (nfr - 1)))

        B.render_frames(frames_name, nfr, update)

    elif mode in ("chart3d", "layout"):
        tag, out = a[1], a[2]
        C = dict(np.load(B.WORK / f"dp_cut_{tag}.npz"))
        info = json.loads((B.WORK / f"dp_cut_{tag}.json").read_text())
        nch = int(C["face_chart"].max()) + 1
        pal = np.array([B.lin(tuple(int(h[i:i + 2], 16) / 255 for i in (1, 3, 5))) for h in CHART_COLORS])
        fc = pal[C["face_chart"] % len(pal)]
        if mode == "chart3d":
            az, el = float(a[3]), float(a[4])
            W = int(a[5]) if len(a) > 5 else 1400
            H = int(a[6]) if len(a) > 6 else 1000
            _, norm = B.normalize(B.y_up_to_z_up(C["P"]), size=2.0)
            P = to_scene(C["P"], norm)
            B.reset()
            B.setup(W, H, samples=128)
            ob = B.mesh_object("mesh", P, C["F"], smooth=False, material=B.mat_attr("Col", 0.6))
            B.set_face_colors(ob, fc)
            for key, col, r in (("seam", SEAM, SEAM_R * 0.9), ("split", B.lin((0.05, 0.06, 0.08)), SEAM_R * 1.25)):
                sz = C[f"{key}_sizes"]
                if len(sz):
                    pts = to_scene(C[f"{key}_pts"], norm)
                    lines, k = [], 0
                    for s in sz:
                        lines.append(pts[k:k + s])
                        k += s
                    B.tubes(key, lines, r, B.mat_solid(col, roughness=0.4), resolution=3)
            frame_scene(P, az, el, margin=1.05, lens=65)
            B.render_still(out)
        else:
            size = int(a[3]) if len(a) > 3 else 1200
            uv = C["uv"]            # (3 F, 2) per triangle corner
            F3 = np.arange(len(uv)).reshape(-1, 3)
            pos = np.zeros((len(uv), 3))
            pos[:, :2] = uv
            B.reset()
            B.setup(size, size, samples=48)
            ob = B.mesh_object("layout", pos, F3, smooth=False)
            col = np.repeat(fc, 3, axis=0)
            # triangles covered more than once: darker
            if "face_overlap" in C:
                col = np.where(np.repeat(C["face_overlap"], 3)[:, None] > 0, col * 0.35 + np.array([0.35, 0, 0]), col)
            B.set_corner_colors(ob, col)
            mat = bpy.data.materials.new("flat")
            mat.use_nodes = True
            nt = mat.node_tree
            nt.nodes.remove(nt.nodes["Principled BSDF"])
            em = nt.nodes.new("ShaderNodeEmission")
            at = nt.nodes.new("ShaderNodeAttribute")
            at.attribute_name = "Col"
            nt.links.new(at.outputs["Color"], em.inputs["Color"])
            nt.links.new(em.outputs[0], nt.nodes["Material Output"].inputs["Surface"])
            ob.data.materials.append(mat)
            lo, hi = uv.min(0), uv.max(0)
            cx, cy = (lo + hi) / 2
            B.camera((cx, cy, 5.0), (cx, cy, 0.0), ortho_scale=float(max(hi - lo)) * 1.06)
            bpy.context.scene.world.node_tree.nodes["Background"].inputs["Strength"].default_value = 0.0
            B.render_still(out)
