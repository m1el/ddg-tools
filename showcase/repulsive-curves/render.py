"""Repulsive Curves (Yu, Schumacher, Crane 2021): Blender renders of the frames written by
`examples/scenes.rs` (`frameNNNNN.obj`, curve OBJ files) as tubes.

blender -b --python showcase/repulsive-curves/render.py -- still <kind> <frame.obj> <out.png> [W H]
blender -b --python showcase/repulsive-curves/render.py -- anim  <kind> <frames-dir> <name> [W H] [N]
blender -b --python showcase/repulsive-curves/render.py -- sheet <kind> <out.png> <frame.obj>... (W H via RC_W, RC_H)

`anim` renders N (default 90) frames evenly chosen from the directory to
work/frames/<name>/NNNN.png with a camera fixed on all frames (SHOW_KEEP skips frames already
there). `sheet` puts several states of one run on one image (a row of panels). The kinds are in
KINDS; the reference's y axis is up, the planar scenes lie in y = 0 and are seen from above.
"""

import glob
import math
import os
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "render"))
import blib as B  # noqa: E402
import numpy as np  # noqa: E402

DATA = HERE.parent.parent / "data" / "repulsive-curves"
SCENES = DATA / "scenes"

INK = B.lin((0.07, 0.08, 0.10))
GOLD = B.lin((0.95, 0.60, 0.08))
BLUE = B.lin((0.12, 0.34, 0.80))
RED = B.lin((0.80, 0.12, 0.10))
GREY = B.lin((0.78, 0.79, 0.82))

# kind -> options
#   obstacle: OBJ path (triangle mesh drawn as a pale surface)
#   surface: ("torus", R, r) | ("sphere", r) | OBJ path   (the constraint surface)
#   top: planar scene seen from above (y = 0 plane)
#   az, el: camera angles (degrees)   r: tube radius as a fraction of the bounding diagonal
#   nodes: spheres at junctions / endpoints (vertices of degree != 2)
#   color: tube colour (or "arclength" for a gradient along each component)
#   glass: the obstacle is drawn transparent
KINDS = {
    "bunny": dict(obstacle=SCENES / "BunnyCurvePacking/bunny.obj", az=-35, el=12, r=0.0045,
                  color=GOLD, glass=True, frame_obstacle=True),
    "rings": dict(obstacle=SCENES / "RingsPacking/bunny.obj", az=-35, el=12, r=0.004,
                  color=GOLD, glass=True, frame_obstacle=True),
    "cup": dict(obstacle=SCENES / "SimpleCurvePacking/round-cup.obj", az=0, el=70, r=0.0035,
                color=GOLD, glass=True, frame_obstacle=True),
    "rect": dict(obstacle=SCENES / "SimpleCurvePacking/rectangle.obj", az=0, el=70, r=0.0035,
                 color=GOLD, glass=True, frame_obstacle=True),
    "torus": dict(surface=("torus", 1.0, 0.25), az=-20, el=38, r=0.006, color=GOLD),
    "sphere": dict(surface=("sphere", 1.0), az=-30, el=22, r=0.0045, color=GOLD),
    "doubletorus": dict(surface=SCENES / "ImplicitDoubleTorus/implicitRemeshed.obj", az=0, el=60,
                        r=0.005, color=GOLD),
    "bluenoise": dict(obstacle=SCENES / "BluenoiseCurves/bear-outline.obj", top=True, r=0.003,
                      color=GOLD, glass=False, flat_obstacle=True),
    "graph": dict(top=False, az=-20, el=35, r=0.007, color=BLUE, nodes=True),
    "graph-top": dict(top=True, r=0.006, color=BLUE, nodes=True),
    "capillary": dict(obstacle=SCENES / "CapillaryNetwork/hand-tri.obj", az=0, el=89.9, r=0.0018,
                      color=RED, glass=True, nodes=False),
    "spline": dict(top=True, r=0.008, color=GOLD, nodes=True),
    "unknot": dict(az=-30, el=20, r=0.011, color="arclength"),
}


def read_curve(path):
    V, E = [], []
    for line in open(path):
        p = line.split()
        if not p:
            continue
        if p[0] == "v":
            V.append([float(x) for x in p[1:4]])
        elif p[0] == "l":
            E.append([int(p[1]) - 1, int(p[2]) - 1])
    return np.array(V), np.array(E, dtype=np.int64)


def read_mesh(path):
    V, F = [], []
    for line in open(path):
        p = line.split()
        if not p:
            continue
        if p[0] == "v":
            V.append([float(x) for x in p[1:4]])
        elif p[0] == "f":
            idx = [int(t.split("/")[0]) - 1 for t in p[1:]]
            for k in range(1, len(idx) - 1):
                F.append([idx[0], idx[k], idx[k + 1]])
    return np.array(V), np.array(F, dtype=np.int64)


def chains(n, E):
    """Maximal paths of edges between vertices of degree != 2 (closed loops included)."""
    adj = [[] for _ in range(n)]
    for i, (a, b) in enumerate(E):
        adj[a].append((b, i))
        adj[b].append((a, i))
    used = np.zeros(len(E), dtype=bool)
    out, closed = [], []
    deg = np.array([len(a) for a in adj])

    def walk(start, nb, ei):
        path = [start]
        cur, e = nb, ei
        used[e] = True
        path.append(cur)
        while deg[cur] == 2:
            nxt = [(w, k) for w, k in adj[cur] if not used[k]]
            if not nxt:
                break
            w, k = nxt[0]
            used[k] = True
            path.append(w)
            cur = w
            if cur == start:
                break
        return path

    for v in range(n):
        if deg[v] != 2:
            for nb, ei in adj[v]:
                if not used[ei]:
                    out.append(walk(v, nb, ei))
                    closed.append(False)
    for v in range(n):
        for nb, ei in adj[v]:
            if not used[ei]:
                p = walk(v, nb, ei)
                out.append(p)
                closed.append(p[0] == p[-1])
    return out, closed, deg


def to_scene(p, top):
    # reference: y up. Blender: z up.
    return B.y_up_to_z_up(p)


def torus_mesh(R, r, nu=96, nv=48):
    u = np.linspace(0, 2 * np.pi, nu, endpoint=False)
    v = np.linspace(0, 2 * np.pi, nv, endpoint=False)
    U, Vv = np.meshgrid(u, v, indexing="ij")
    P = np.stack([(R + r * np.cos(Vv)) * np.cos(U), r * np.sin(Vv), (R + r * np.cos(Vv)) * np.sin(U)], -1)
    idx = np.arange(nu * nv).reshape(nu, nv)
    F = np.stack([idx, np.roll(idx, -1, 0), np.roll(np.roll(idx, -1, 0), -1, 1), np.roll(idx, -1, 1)], -1)
    return P.reshape(-1, 3), F.reshape(-1, 4)


def sphere_mesh(r, nu=64, nv=32):
    u = np.linspace(0, 2 * np.pi, nu, endpoint=False)
    v = np.linspace(0.02, np.pi - 0.02, nv)
    U, Vv = np.meshgrid(u, v, indexing="ij")
    P = np.stack([r * np.sin(Vv) * np.cos(U), r * np.cos(Vv), r * np.sin(Vv) * np.sin(U)], -1)
    idx = np.arange(nu * nv).reshape(nu, nv)
    F = np.stack([idx[:, :-1], np.roll(idx, -1, 0)[:, :-1], np.roll(idx, -1, 0)[:, 1:], idx[:, 1:]], -1)
    return P.reshape(-1, 3), F.reshape(-1, 4)


def glass_material(color=(0.75, 0.80, 0.90), alpha=0.18):
    mat = B.mat_solid(color, roughness=0.35, name="glass")
    bsdf = mat.node_tree.nodes["Principled BSDF"]
    bsdf.inputs["Alpha"].default_value = alpha
    return mat


class Dressing:
    """Everything of a kind besides the curve: the pale surface or obstacle."""

    def __init__(self, kind, cfg):
        self.cfg = cfg
        self.mesh = None  # (V, F) in Blender coordinates, quads/tris
        s = cfg.get("surface") or cfg.get("obstacle")
        if isinstance(s, tuple):
            if s[0] == "torus":
                V, F = torus_mesh(s[1], s[2])
            else:
                V, F = sphere_mesh(s[1])
            self.mesh = (B.y_up_to_z_up(V), F)
        elif s is not None and Path(s).exists():
            V, F = read_mesh(s)
            self.mesh = (B.y_up_to_z_up(V), F)
        if self.mesh is None and (cfg.get("surface") or cfg.get("obstacle")):
            print("warning: missing", s)

    def add(self):
        if self.mesh is None:
            return None
        V, F = self.mesh
        cfg = self.cfg
        if cfg.get("glass") or cfg.get("surface"):
            mat = glass_material() if cfg.get("glass") else B.mat_solid(B.lin((0.84, 0.85, 0.88)), 0.5)
        else:
            mat = B.mat_solid(B.lin((0.82, 0.83, 0.86)), 0.6)
        obj = B.mesh_object("dressing", V, F, smooth=not cfg.get("flat_obstacle"), material=mat)
        if cfg.get("glass") or cfg.get("surface"):
            obj.data.materials[0].blend_method = "BLEND" if hasattr(obj.data.materials[0], "blend_method") else None
            if cfg.get("surface"):
                bsdf = obj.data.materials[0].node_tree.nodes["Principled BSDF"]
                bsdf.inputs["Alpha"].default_value = 0.35
        return obj


def arclength_colors(path_pts):
    d = np.r_[0, np.cumsum(np.linalg.norm(np.diff(path_pts, axis=0), axis=1))]
    return d / max(d[-1], 1e-12)


def build_curve(kind, cfg, V, E, norm):
    P = apply(V, norm)
    n = len(P)
    ch, closed, deg = chains(n, E)
    R = cfg["r"] * norm[2]
    if cfg.get("color") == "arclength":
        # one object per segment group: colour each chain piece along its length
        pieces, cols = [], []
        for c in ch:
            pts = P[c]
            t = arclength_colors(pts)
            k = max(2, len(c) // 40)
            for s in range(0, len(c) - 1, k):
                seg = np.arange(s, min(s + k + 1, len(c)))
                pieces.append(pts[seg])
                cols.append(t[seg].mean())
        # a small palette from the colormap: one material per bin
        bins = 24
        groups = [[] for _ in range(bins)]
        for pc, cv in zip(pieces, cols):
            groups[min(bins - 1, int(cv * bins))].append(pc)
        for b, g in enumerate(groups):
            if not g:
                continue
            col = B.cmap(np.array([(b + 0.5) / bins]), "spectral")[0]
            B.tubes(f"curve{b}", g, R, B.mat_solid(B.lin(tuple(col)), roughness=0.3), resolution=6)
    else:
        B.tubes("curve", [P[c] for c in ch], R, B.mat_solid(cfg["color"], roughness=0.3),
                closed=closed, resolution=6)
    if cfg.get("nodes"):
        special = np.where(deg != 2)[0]
        if len(special):
            B.spheres("nodes", P[special], 1.9 * R, material=B.mat_solid(INK, roughness=0.35))
    return P


def apply(p, norm):
    s, off, _ = norm
    return (B.y_up_to_z_up(p) + off) * s


def make_norm(allpts, size=2.0):
    """Scale so that the bounding diagonal of the union of `allpts` (in Blender coordinates) is
    `size`; the third entry is the scaled diagonal (for the tube radius)."""
    p = B.y_up_to_z_up(np.asarray(allpts))
    lo, hi = p.min(0), p.max(0)
    s = size / np.linalg.norm(hi - lo)
    off = -np.array([(lo[0] + hi[0]) / 2, (lo[1] + hi[1]) / 2, lo[2]])
    return (s, off, size)


def scene_setup(kind, cfg, frames_pts, W, H, samples=96, dressing=True):
    """Union bounding box of all frames (and the dressing when `frame_obstacle`)."""
    dr = Dressing(kind, cfg)
    allpts = [np.asarray(f) for f in frames_pts]
    cat = np.concatenate(allpts)
    if dr.mesh is not None and cfg.get("frame_obstacle"):
        # obstacle in reference coordinates for the normalisation
        V0 = np.stack([dr.mesh[0][:, 0], dr.mesh[0][:, 2], -dr.mesh[0][:, 1]], -1)
        cat = np.concatenate([cat, V0])
    norm = make_norm(cat)
    return dr, norm, cat


def frame_camera(cfg, norm, cat, W, H):
    P = apply(cat, norm)
    c = (P.min(0) + P.max(0)) / 2
    rad = float(np.linalg.norm(P.max(0) - P.min(0))) / 2
    B.studio(c, rad, ground_z=float(P[:, 2].min()) - 0.02 * rad)
    if cfg.get("top"):
        lo, hi = P.min(0), P.max(0)
        ext = max(hi[0] - lo[0], (hi[1] - lo[1]) * W / H * (H / W) if False else hi[1] - lo[1])
        asp = W / H
        sx, sy = (hi[0] - lo[0]), (hi[1] - lo[1])
        ortho = max(sx, sy * asp) * 1.12
        B.camera((c[0], c[1], c[2] + 10), (c[0], c[1], c[2]), ortho_scale=ortho)
    else:
        B.fit_camera(P, azimuth=cfg.get("az", -30), elevation=cfg.get("el", 25), lens=50,
                     margin=1.1)


def add_dressing(dr, norm):
    obj = dr.add()
    if obj is None:
        return
    s, off, _ = norm
    # dressing coordinates are already z-up; apply the same offset and scale
    me = obj.data
    co = np.empty(len(me.vertices) * 3, dtype=np.float32)
    me.vertices.foreach_get("co", co)
    co = (co.reshape(-1, 3) + off) * s
    me.vertices.foreach_set("co", co.astype(np.float32).ravel())
    me.update()


def clear_curve_objects():
    for o in list(B.bpy.data.objects):
        if o.name.startswith(("curve", "nodes")):
            B.bpy.data.objects.remove(o, do_unlink=True)


def pick_frames(files, n):
    if len(files) <= n:
        return list(range(len(files)))
    idx = np.unique(np.round(np.linspace(0, len(files) - 1, n)).astype(int))
    return [int(i) for i in idx]


def main():
    a = B.args()
    mode, kind = a[0], a[1]
    cfg = KINDS[kind]
    if mode == "still":
        frame, out = a[2], a[3]
        W = int(a[4]) if len(a) > 4 else 1600
        H = int(a[5]) if len(a) > 5 else 1000
        V, E = read_curve(frame)
        B.reset()
        B.setup(W, H, samples=int(os.environ.get("RC_SAMPLES", 128)))
        dr, norm, cat = scene_setup(kind, cfg, [V], W, H)
        add_dressing(dr, norm)
        build_curve(kind, cfg, V, E, norm)
        frame_camera(cfg, norm, cat, W, H)
        B.render_still(out)
    elif mode == "anim":
        d, name = a[2], a[3]
        W = int(a[4]) if len(a) > 4 else 960
        H = int(a[5]) if len(a) > 5 else 720
        n = int(a[6]) if len(a) > 6 else 90
        files = sorted(glob.glob(os.path.join(d, "frame*.obj")))
        sel = pick_frames(files, n)
        curves = [read_curve(files[i]) for i in sel]
        B.reset()
        B.setup(W, H, samples=int(os.environ.get("RC_SAMPLES", 64)))
        dr, norm, cat = scene_setup(kind, cfg, [c[0] for c in curves], W, H)
        add_dressing(dr, norm)
        frame_camera(cfg, norm, cat, W, H)

        def update(i):
            clear_curve_objects()
            V, E = curves[i]
            build_curve(kind, cfg, V, E, norm)

        B.render_frames(name, len(curves), update)
        # frame -> step labels for post.py
        import json
        json.dump([os.path.basename(files[i]) for i in sel], open(B.FRAMES / name / "frames.json", "w"))
    else:
        raise SystemExit("mode: still | anim")


main()
