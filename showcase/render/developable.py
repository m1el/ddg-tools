"""Developability scenes (Stein, Grinspun, Crane 2018): meshes shaded by the dihedral angle of
their edges (white = flat, dark = large angle, as in the paper's Figs. 14 and 26), colored by the
vertex energy or the angle defect, seams as tubes and terminals as spheres.

blender -b --python render/developable.py -- pair  <tag> <out.png> [az] [el] [color] [W H]
blender -b --python render/developable.py -- seams <tag> <out.png> [az] [el] [color] [W H]
(more modes below: flow, compare, cut)
color: energy (fixed log scale) | defect | plain
"""
import json
import math
import sys

sys.path.insert(0, __file__.rsplit("/", 1)[0])
import blib as B  # noqa: E402
import numpy as np  # noqa: E402

bpy = B.bpy
INK = B.lin((0.05, 0.06, 0.08))
SURF = (0.86, 0.86, 0.88)  # sRGB of the plain surface

# Fixed scales shared by every picture, so colors can be compared between pictures.
E_LO, E_HI = 1e-5, 1e-1       # vertex energy, log scale
DEFECT_MAX = 0.10             # radians, signed (blue negative, red positive)


# ----------------------------------------------------------------------------- data
def state(D, prefix):
    """Positions in the scene frame (z up, bbox diagonal 2.2 set by the first call's norm)."""
    return D[f"{prefix}_P"], D[f"{prefix}_F"]


B._CMAPS["energy"] = ["#f7f7f4", "#fff2b0", "#fdc761", "#f4823c", "#d03a2f", "#8a1a4a",
                      "#3b0f4f"]


def energy_colors(e, lo=E_LO, hi=E_HI):
    """Colors of log10(energy) on a fixed scale: near white at `lo` and below, through yellow and
    red to dark purple at `hi` and above."""
    v = np.log10(np.maximum(e, 1e-300))
    return B.cmap(v, "energy", np.log10(lo), np.log10(hi))


def defect_colors(x, vmax=DEFECT_MAX):
    return B.cmap(x, "rdbu", -vmax, vmax)  # positive defect red, negative blue


# ----------------------------------------------------------------------------- geometry
def face_vectors(P, F):
    """Per face the heights of its three corners over the opposite edges."""
    a, b, c = P[F[:, 0]], P[F[:, 1]], P[F[:, 2]]
    area2 = np.linalg.norm(np.cross(b - a, c - a), axis=1)
    le = np.stack([np.linalg.norm(c - b, axis=1), np.linalg.norm(a - c, axis=1),
                   np.linalg.norm(b - a, axis=1)], axis=1)
    return area2[:, None] / np.maximum(le, 1e-12)


def add_face_vector(obj, name, values):
    attr = obj.data.attributes.new(name, "FLOAT_VECTOR", "FACE")
    attr.data.foreach_set("vector", np.asarray(values, dtype=np.float32).ravel())


import os  # noqa: E402

EDGE_FLOOR = float(os.environ.get("DEV_FLOOR", 8.0))
EDGE_SCALE = float(os.environ.get("DEV_SCALE", 30.0))
EDGE_LW = float(os.environ.get("DEV_LW", 0.010))


def dihedral_strength(ang, scale_deg=None, floor_deg=None):
    """Edge darkness in [0, 1] from a dihedral angle in radians: 0 up to `floor_deg`, then a
    smooth rise to 1 at `scale_deg` (white = flat, dark = a large angle)."""
    scale_deg = EDGE_SCALE if scale_deg is None else scale_deg
    floor_deg = EDGE_FLOOR if floor_deg is None else floor_deg
    t = np.clip((np.degrees(ang) - floor_deg) / (scale_deg - floor_deg), 0, 1)
    return t * t * (3 - 2 * t) ** 0.7


def mat_dihedral(use_point_color, line_width, roughness=0.5, base=SURF, name="dihedral",
                 line_color=INK, bump=None):
    """Surface whose edges are drawn with a darkness given by the per-face vector 'dih'
    (strength of the edge opposite corner 0, 1, 2) and the per-face vector 'hgt' (corner heights)
    so the line has the constant world half-width `line_width`."""
    mat = bpy.data.materials.new(name)
    nt, bsdf = B._principled(mat)
    bsdf.inputs["Roughness"].default_value = roughness
    N = nt.nodes
    L = nt.links
    geo = N.new("ShaderNodeNewGeometry")
    sep_b = N.new("ShaderNodeSeparateXYZ")
    L.new(geo.outputs["Parametric"], sep_b.inputs[0])
    # barycentric of corner 0, 1, 2
    w0 = N.new("ShaderNodeMath")
    w0.operation = "SUBTRACT"
    w0.inputs[0].default_value = 1.0
    tmp = N.new("ShaderNodeMath")
    tmp.operation = "ADD"
    L.new(sep_b.outputs["X"], tmp.inputs[0])
    L.new(sep_b.outputs["Y"], tmp.inputs[1])
    L.new(tmp.outputs[0], w0.inputs[1])
    bary = [w0.outputs[0], sep_b.outputs["X"], sep_b.outputs["Y"]]
    dih = N.new("ShaderNodeAttribute")
    dih.attribute_name = "dih"
    hgt = N.new("ShaderNodeAttribute")
    hgt.attribute_name = "hgt"
    sd = N.new("ShaderNodeSeparateXYZ")
    sh = N.new("ShaderNodeSeparateXYZ")
    L.new(dih.outputs["Vector"], sd.inputs[0])
    L.new(hgt.outputs["Vector"], sh.inputs[0])
    strength = None
    for i, comp in enumerate(("X", "Y", "Z")):
        dist = N.new("ShaderNodeMath")
        dist.operation = "MULTIPLY"
        L.new(bary[i], dist.inputs[0])
        L.new(sh.outputs[comp], dist.inputs[1])
        ss = N.new("ShaderNodeMapRange")
        ss.interpolation_type = "SMOOTHSTEP"
        L.new(dist.outputs[0], ss.inputs["Value"])
        ss.inputs["From Min"].default_value = line_width * 0.55
        ss.inputs["From Max"].default_value = line_width * 1.25
        ss.inputs["To Min"].default_value = 1.0
        ss.inputs["To Max"].default_value = 0.0
        prod = N.new("ShaderNodeMath")
        prod.operation = "MULTIPLY"
        L.new(ss.outputs["Result"], prod.inputs[0])
        L.new(sd.outputs[comp], prod.inputs[1])
        if strength is None:
            strength = prod.outputs[0]
        else:
            mx = N.new("ShaderNodeMath")
            mx.operation = "MAXIMUM"
            L.new(strength, mx.inputs[0])
            L.new(prod.outputs[0], mx.inputs[1])
            strength = mx.outputs[0]
    if use_point_color:
        col = N.new("ShaderNodeAttribute")
        col.attribute_name = "Col"
        base_out = col.outputs["Color"]
    else:
        rgb = N.new("ShaderNodeRGB")
        rgb.outputs[0].default_value = (*B.lin(base), 1)
        base_out = rgb.outputs[0]
    mix = N.new("ShaderNodeMix")
    mix.data_type = "RGBA"
    L.new(strength, mix.inputs["Factor"])
    L.new(base_out, mix.inputs["A"])
    mix.inputs["B"].default_value = (*line_color, 1)
    L.new(mix.outputs["Result"], bsdf.inputs["Base Color"])
    return mat


def surface(name, P, F, ang_f, color=None, line_width=0.004, roughness=0.5, smooth=True,
            strength_scale=None):
    """A triangle mesh drawn with dihedral edges; `color` (n, 3) linear vertex colors or None."""
    mat = mat_dihedral(color is not None, line_width, roughness=roughness)
    obj = B.mesh_object(name, P, F, smooth=smooth, material=mat)
    if color is not None:
        B.set_point_colors(obj, color)
    set_edges(obj, P, F, ang_f, strength_scale)
    return obj


def set_edges(obj, P, F, ang_f, strength_scale=None):
    add_face_vector(obj, "dih", dihedral_strength(ang_f, strength_scale))
    add_face_vector(obj, "hgt", face_vectors(P, F))


def set_point_color(obj, rgb):
    me = obj.data
    rgba = np.ones((len(me.vertices), 4), dtype=np.float32)
    rgba[:, :3] = rgb
    me.color_attributes["Col"].data.foreach_set("color", rgba.ravel())


def seam_objects(P, edges, terminals, tube_r, sphere_r, root_color=None):
    from_dir = __import__("pathlib").Path(B.SHOW / "compute")
    sys.path.insert(0, str(from_dir))
    from common import chain_edges  # noqa: E402
    objs = []
    if len(edges):
        chains = chain_edges(np.asarray(edges))
        objs.append(B.tubes("seams", [P[c] for c in chains], tube_r,
                            B.mat_solid(B.lin((0.85, 0.12, 0.08)), roughness=0.35)))
    if len(terminals):
        objs.append(B.spheres("terminals", P[terminals], sphere_r,
                              colors=np.tile(B.lin((0.98, 0.78, 0.10)), (len(terminals), 1)),
                              material=B.mat_attr("Col", 0.3)))
    return objs


def frame_scene(points, az, el, margin=1.06, lens=55):
    ground = points[:, 2].min()
    center = (points.min(0) + points.max(0)) / 2
    radius = float(np.linalg.norm(points.max(0) - points.min(0))) / 2
    B.studio(center, radius, ground_z=ground)
    return B.fit_camera(points, azimuth=az, elevation=el, lens=lens, margin=margin)


def to_scene(p, norm):
    return B.apply_norm(B.y_up_to_z_up(p), norm)


def colors_for(kind, D, prefix):
    if kind == "energy":
        return energy_colors(D[f"{prefix}_energy"])
    if kind == "defect":
        return defect_colors(D[f"{prefix}_defect"])
    return None


def load(tag):
    D = B.load(f"developable_{tag}")
    info = json.loads((B.WORK / f"developable_{tag}.json").read_text())
    return D, info


DEFAULT_VIEW = {"spot": (-50, 16), "bob": (-40, 22), "blub": (-40, 20), "spot_control": (-50, 16), "bob_control": (-40, 22), "blub_control": (-40, 20),
                "crumpled": (-30, 35)}

if __name__ == "__main__":
    a = B.args()
    mode, tag = a[0], a[1]
    D, INFO = load(tag) if mode != "flow" else (None, {"mesh": ""})
    mesh = INFO["mesh"]
    if mode in ("pair", "seams", "single"):
        out = a[2]
        az = float(a[3]) if len(a) > 3 and a[3] != "-" else DEFAULT_VIEW.get(mesh, (-40, 20))[0]
        el = float(a[4]) if len(a) > 4 and a[4] != "-" else DEFAULT_VIEW.get(mesh, (-40, 20))[1]
        kind = a[5] if len(a) > 5 else "plain"
        W = int(a[6]) if len(a) > 6 else 2000
        H = int(a[7]) if len(a) > 7 else 1000
        B.reset()
        B.setup(W, H, samples=128)
        Pin, Fin = state(D, "in")
        Pout, Fout = state(D, "out")
        P0 = B.y_up_to_z_up(Pin)
        _, norm = B.normalize(P0, size=2.0)
        gap = float(os.environ.get('DEV_GAP', 1.6))
        ar = math.radians(az)
        # side by side across the view direction
        sx = np.array([math.cos(ar), math.sin(ar), 0.0]) * gap / 2
        Pi = to_scene(Pin, norm) - sx
        Po = to_scene(Pout, norm) + sx
        lw = EDGE_LW
        if mode == "single":
            prefix = a[8] if len(a) > 8 else "out"
            Ps = to_scene(D[f"{prefix}_P"], norm)
            surface("mesh", Ps, D[f"{prefix}_F"], D[f"{prefix}_dih"], colors_for(kind, D, prefix),
                    line_width=lw * (2 if prefix == "in" else 1))
            frame_scene(Ps, az, el, margin=1.04, lens=65)
            B.render_still(out)
            sys.exit(0)
        if mode == "pair":
            surface("input", Pi, Fin, D["in_dih"], colors_for(kind, D, "in"), line_width=lw * 2)
            surface("output", Po, Fout, D["out_dih"], colors_for(kind, D, "out"), line_width=lw)
        frame_scene(np.concatenate([Pi, Po]), az, el, margin=1.04, lens=65)
        B.render_still(out)

    # ------------------------------------------------------------------ seams / cut / layout
    elif mode in ("cut3d", "layout"):
        # cut3d <tag> <out> az el eps_index [W H]; layout <tag> <out> eps_index [size]
        C = dict(np.load(B.WORK / f"developable_cut_{tag}.npz"))
        CI = json.loads((B.WORK / f"developable_cut_{tag}.json").read_text())
        out = a[2]
        Pin = D["in_P"]
        _, norm = B.normalize(B.y_up_to_z_up(Pin), size=2.0)
        Pout, Fout = D["out_P"], D["out_F"]
        Ps = to_scene(Pout, norm)
        LS_MAX = float(os.environ.get("DEV_LS_MAX", 0.3))
        if mode == "cut3d":
            az, el, k = float(a[3]), float(a[4]), int(a[5])
            W = int(a[6]) if len(a) > 6 else 1600
            H = int(a[7]) if len(a) > 7 else 1200
            B.reset()
            B.setup(W, H, samples=128)
            if os.environ.get("DEV_COLOR") == "defect":
                surface("surface", Ps, Fout, D["out_dih"], defect_colors(D["out_defect"]),
                        line_width=EDGE_LW * 0.7)
            else:
                col = B.cmap(C[f"e{k}_logscale"], "rdbu", -LS_MAX, LS_MAX)
                ob = B.mesh_object("surface", Ps, Fout, smooth=True, material=B.mat_attr("Col", 0.5))
                B.set_corner_colors(ob, col)
            seam_objects(Ps, C[f"e{k}_seam_edges"], C[f"e{k}_terminals"], 0.0045, 0.011)
            frame_scene(Ps, az, el, margin=1.05, lens=65)
            B.render_still(out)
        else:
            k = int(a[3])
            size = int(a[4]) if len(a) > 4 else 1200
            uv = C[f"e{k}_uv"]
            ch = C[f"e{k}_chart"]
            # pack: charts side by side is the library's job (packed in the unit square)
            B.reset()
            B.setup(size, size, samples=64)
            pos = np.zeros((len(uv), 3))
            pos[:, :2] = uv
            nf = len(uv) // 3
            F3 = np.arange(len(uv)).reshape(-1, 3)
            col = B.cmap(C[f"e{k}_logscale"], "rdbu", -LS_MAX, LS_MAX)
            ob = B.mesh_object("layout", pos, F3, smooth=False)
            B.set_corner_colors(ob, col)
            mat = B.bpy.data.materials.new("flat")
            mat.use_nodes = True
            nt = mat.node_tree
            nt.nodes.remove(nt.nodes["Principled BSDF"])
            em = nt.nodes.new("ShaderNodeEmission")
            at = nt.nodes.new("ShaderNodeAttribute")
            at.attribute_name = "Col"
            nt.links.new(at.outputs["Color"], em.inputs["Color"])
            nt.links.new(em.outputs[0], nt.nodes["Material Output"].inputs["Surface"])
            ob.data.materials.append(mat)
            # chart outlines: triangle sides that no other triangle shares in UV
            from collections import Counter
            q = np.round(uv * 2e5).astype(np.int64)
            cnt, sides = Counter(), []
            for f in range(nf):
                for i in range(3):
                    i0, i1 = 3 * f + i, 3 * f + (i + 1) % 3
                    key = tuple(sorted([tuple(q[i0]), tuple(q[i1])]))
                    cnt[key] += 1
                    sides.append((key, i0, i1))
            lines = [pos[[i0, i1]] + [0, 0, 0.0005] for key, i0, i1 in sides if cnt[key] == 1]
            B.tubes("outline", lines, 0.0010, B.mat_emission(INK, 1.0), resolution=1)
            lo, hi = uv.min(0), uv.max(0)
            cx, cy = (lo + hi) / 2
            B.camera((cx, cy, 5.0), (cx, cy, 0.0), ortho_scale=float(max(hi - lo)) * 1.06)
            B.bpy.context.scene.world.node_tree.nodes["Background"].inputs["Strength"].default_value = 0.0
            B.render_still(out)

    # ------------------------------------------------------------------ flow animation
    elif mode == "flow":
        # flow <snap-tag> <frames-name> az el [turn-deg] [per-state frames]
        snap = a[1]
        S = dict(np.load(B.WORK / f"developable_snap_{snap}.npz"))
        SI = json.loads((B.WORK / f"developable_snap_{snap}.json").read_text())
        meta = SI["states"]
        n = int(S["n_states"])
        frames_name = a[2]
        az, el = float(a[3]), float(a[4])
        turn = float(a[5]) if len(a) > 5 else 40.0
        per = int(a[6]) if len(a) > 6 else 4
        kind = os.environ.get("DEV_KIND", "energy")
        _, norm = B.normalize(B.y_up_to_z_up(S["s0_P"]), size=2.0)
        allP = [to_scene(S[f"s{i}_P"], norm) for i in range(n)]
        B.reset()
        B.setup(1280, 800, samples=48)
        pivot = bpy.data.objects.new("pivot", None)
        bpy.context.collection.objects.link(pivot)
        rots = [B.rotz(np.concatenate(allP), t) for t in np.linspace(-turn / 2, turn / 2, 5)]
        frame_scene(np.concatenate(rots), az, el, margin=1.04, lens=60)
        # schedule: hold the input, then `per` frames between consecutive states, hold the end
        sched = [(0, 0.0)] * 14
        for i in range(n - 1):
            for j in range(per):
                sched.append((i, j / per))
        sched += [(n - 1, 0.0)] * 40
        nfr = len(sched)
        labels = []
        R = SI["subdivisions"] + 1
        e_first, e_last = meta[0]["energy"], meta[-1]["energy"]
        for (i, t) in sched:
            m = meta[i]
            if m["round"] < 0:
                labels.append(f"input    vertices {m['n_vertices']}    energy {m['energy']:.3g}")
            else:
                labels.append(f"round {m['round'] + 1}/{R}  iteration {m['iteration']:3d}"
                              f"    vertices {m['n_vertices']}    energy {m['energy']:.3g}")
        (B.WORK / f"labels_{frames_name}.json").write_text(json.dumps(labels))
        cur = [None]

        def update(f):
            i, t = sched[f]
            j = min(i + 1, n - 1)
            F = S[f"s{i}_F"]
            same = t > 0 and S[f"s{j}_F"].shape == F.shape and np.array_equal(S[f"s{j}_F"], F)
            w = B.ease(t) if same else 0.0
            P = (1 - w) * allP[i] + w * allP[j] if same else allP[i]
            e = S[f"s{i}_energy"]
            if same:
                e = np.exp((1 - w) * np.log(np.maximum(e, 1e-300)) + w * np.log(np.maximum(S[f"s{j}_energy"], 1e-300)))
            if cur[0] is not None:
                bpy.data.objects.remove(cur[0], do_unlink=True)
            col = energy_colors(e) if kind == "energy" else None
            obj = surface("mesh", P, F, S[f"s{i}_dih"], col, line_width=EDGE_LW * 0.8)
            obj.parent = pivot
            cur[0] = obj
            pivot.rotation_euler[2] = math.radians(-turn / 2 + turn * B.ease(f / (nfr - 1)))

        B.render_frames(frames_name, nfr, update)
