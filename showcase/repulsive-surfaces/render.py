"""Repulsive Surfaces (Yu, Brakensiek, Schumacher, Crane 2021): Blender renders of the OBJ frames
written by `examples/apps` and `examples/ishii`, in the studio look of the Repulsive Curves renders
(`showcase/render/blib.py`: soft area lights, shadow catcher, Cycles).

blender -b --python showcase/repulsive-surfaces/render.py -- still <spec.json> <out.png> [W H]
blender -b --python showcase/repulsive-surfaces/render.py -- sheet <spec.json> <outdir> [W H]
blender -b --python showcase/repulsive-surfaces/render.py -- anim  <spec.json> <name> [W H]

A *panel* is {"layers": [layer, ...], "az": -35, "el": 22, "margin": 1.1}.
  still: spec is one panel.            -> out.png (transparent; post.py flattens it)
  sheet: spec is {"panels": [panel]}   -> outdir/pNN.png, each panel framed on its own; with
         "common_scale": true every panel is framed with the largest bounding diagonal, so the
         sizes of the shapes can be compared.
  anim:  spec is {"panels": [panel], "az": .., "el": ..}: one frame per panel, one camera and one
         normalisation (union of the framed layers of all panels) -> showcase/work/frames/<name>/NNNN.png.
A *layer*:
  {"obj": path}                 triangle mesh, smooth shaded
     "color": [r, g, b] linear-ish sRGB triple | "components" (one palette colour per component)
     "alpha": 1.0 (<1: translucent glass)   "roughness": 0.35   "flat": false (flat shading)
     "wire": false (also draw the mesh edges)   "frame": true (include in the camera framing)
  {"points": path}              the vertices of an OBJ as small spheres ("radius" as a fraction of the
                                bounding diagonal, "index": [subset], "color")
  {"cylinder": {"r":, "z0":, "z1":, "center": [x, y]}}   an axis-parallel (z) cylinder
  {"box": {"center": [x, y, z], "half": [hx, hy, hz]}}   a box
  {"torus": {"R":, "r":, "center": [x, y, z]}}           a torus with the z axis as its axis
anim specs may set "track": true: every frame is normalised to the same size on its own (a zoom that
follows the shape; the camera is fitted to the union of the normalised frames).
The scenes are z up (the reference's applications are written that way); nothing is rotated.
"""

import json
import os
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "render"))
import blib as B  # noqa: E402
import numpy as np  # noqa: E402

GOLD = B.lin((0.95, 0.60, 0.08))
BLUE = B.lin((0.12, 0.34, 0.80))
RED = B.lin((0.80, 0.12, 0.10))
GREEN = B.lin((0.16, 0.58, 0.32))
PURPLE = B.lin((0.50, 0.25, 0.70))
TEAL = B.lin((0.10, 0.62, 0.66))
PINK = B.lin((0.85, 0.35, 0.55))
GREY = B.lin((0.55, 0.57, 0.62))
INK = B.lin((0.07, 0.08, 0.10))
PALETTE = [GOLD, BLUE, RED, GREEN, PURPLE, TEAL, PINK, GREY]
PALE = B.lin((0.78, 0.82, 0.90))


def read_obj(path):
    V, F = [], []
    with open(path) as f:
        for line in f:
            if line.startswith("v "):
                V.append(line.split()[1:4])
            elif line.startswith("f "):
                idx = [int(t.split("/")[0]) - 1 for t in line.split()[1:]]
                for k in range(1, len(idx) - 1):
                    F.append([idx[0], idx[k], idx[k + 1]])
    return np.array(V, dtype=np.float64), np.array(F, dtype=np.int64).reshape(-1, 3)


def components(n, F):
    """Connected-component label of every vertex (union-find)."""
    parent = list(range(n))

    def find(a):
        while parent[a] != a:
            parent[a] = parent[parent[a]]
            a = parent[a]
        return a

    for a, b, c in F.tolist():
        ra = find(a)
        parent[find(b)] = ra
        parent[find(c)] = ra
    roots = np.array([find(i) for i in range(n)])
    _, lab = np.unique(roots, return_inverse=True)
    return lab


def cylinder_mesh(r, z0, z1, cx, cy, n=48):
    t = np.linspace(0, 2 * np.pi, n, endpoint=False)
    ring = np.stack([cx + r * np.cos(t), cy + r * np.sin(t)], -1)
    V = np.concatenate([np.c_[ring, np.full(n, z0)], np.c_[ring, np.full(n, z1)],
                        [[cx, cy, z0], [cx, cy, z1]]])
    F = []
    for i in range(n):
        j = (i + 1) % n
        F += [[i, j, n + j], [i, n + j, n + i], [2 * n, j, i], [2 * n + 1, n + i, n + j]]
    return V, np.array(F)


def torus_mesh(R, r, cx=0.0, cy=0.0, cz=0.0, n=96, m=32):
    u = np.linspace(0, 2 * np.pi, n, endpoint=False)
    v = np.linspace(0, 2 * np.pi, m, endpoint=False)
    U, Vv = np.meshgrid(u, v, indexing="ij")
    X = (R + r * np.cos(Vv)) * np.cos(U) + cx
    Y = (R + r * np.cos(Vv)) * np.sin(U) + cy
    Z = r * np.sin(Vv) + cz
    V = np.stack([X, Y, Z], -1).reshape(-1, 3)
    F = []
    for i in range(n):
        for j in range(m):
            a, b = i * m + j, ((i + 1) % n) * m + j
            c, d = ((i + 1) % n) * m + (j + 1) % m, i * m + (j + 1) % m
            F += [[a, b, c], [a, c, d]]
    return V, np.array(F)


def box_mesh(c, h):
    c, h = np.asarray(c, float), np.asarray(h, float)
    V = np.array([[sx, sy, sz] for sz in (-1, 1) for sy in (-1, 1) for sx in (-1, 1)], float) * h + c
    F = [[0, 2, 1], [1, 2, 3], [4, 5, 6], [5, 7, 6], [0, 1, 4], [1, 5, 4], [2, 6, 3], [3, 6, 7],
         [0, 4, 2], [2, 4, 6], [1, 3, 5], [3, 7, 5]]
    return V, np.array(F)


class Layer:
    """Geometry of one layer in scene coordinates (before normalisation)."""

    def __init__(self, d, cache):
        self.d = d
        self.frame = d.get("frame", True)
        if "obj" in d:
            if d["obj"] not in cache:
                cache[d["obj"]] = read_obj(d["obj"])
            self.V, self.F = cache[d["obj"]]
            self.kind = "mesh"
        elif "points" in d:
            if d["points"] not in cache:
                cache[d["points"]] = read_obj(d["points"])
            V = cache[d["points"]][0]
            self.V = V[d["index"]] if "index" in d else V
            self.F = None
            self.kind = "points"
        elif "cylinder" in d:
            c = d["cylinder"]
            self.V, self.F = cylinder_mesh(c["r"], c["z0"], c["z1"], *c.get("center", [0, 0]))
            self.kind = "mesh"
        elif "torus" in d:
            t = d["torus"]
            self.V, self.F = torus_mesh(t["R"], t["r"], *t.get("center", [0, 0, 0]))
            self.kind = "mesh"
        elif "box" in d:
            self.V, self.F = box_mesh(d["box"]["center"], d["box"]["half"])
            self.kind = "mesh"
            self.d = dict(d, flat=True)
        else:
            raise ValueError(f"bad layer {d}")

    def add(self, norm, idx):
        s, off = norm
        diag = 2.0
        P = (self.V + off) * s
        d = self.d
        alpha = d.get("alpha", 1.0)
        if self.kind == "points":
            col = B.lin(tuple(d["color"])) if "color" in d else INK
            B.spheres(f"pts{idx}", P, d.get("radius", 0.01) * diag, subdiv=d.get("subdiv", 2),
                      material=B.mat_solid(col, roughness=0.35))
            return
        rough = d.get("roughness", 0.35)
        parts = []
        if d.get("color") == "components":
            lab = components(len(P), self.F)
            for c in range(lab.max() + 1):
                vm = lab == c
                fm = vm[self.F[:, 0]]
                remap = -np.ones(len(P), dtype=np.int64)
                remap[vm] = np.arange(vm.sum())
                parts.append((P[vm], remap[self.F[fm]], PALETTE[c % len(PALETTE)]))
        else:
            col = B.lin(tuple(d["color"])) if "color" in d else (PALE if alpha < 1 else GOLD)
            parts.append((P, self.F, col))
        for k, (pv, fv, col) in enumerate(parts):
            mat = B.mat_solid(col, roughness=rough, name=f"m{idx}_{k}", coat=0.0 if alpha < 1 else 0.15)
            if alpha < 1:
                mat.node_tree.nodes["Principled BSDF"].inputs["Alpha"].default_value = alpha
            obj = B.mesh_object(f"surf{idx}_{k}", pv, fv, smooth=not d.get("flat", False), material=mat)
            if d.get("wire"):
                obj.data.materials.append(B.mat_solid(INK, roughness=0.6, name="wire"))
                mod = obj.modifiers.new("wire", "WIREFRAME")
                mod.thickness = d.get("wire_thickness", 0.0025) * diag
                mod.use_replace = False
                mod.material_offset = 1


def clear_geometry():
    for o in list(B.bpy.data.objects):
        if o.type == "MESH" and o.name != "ground":
            B.bpy.data.objects.remove(o, do_unlink=True)


def framed_points(layers):
    return np.concatenate([l.V for l in layers if l.frame])


def make_norm(points, size=2.0):
    lo, hi = points.min(0), points.max(0)
    s = size / np.linalg.norm(hi - lo)
    off = -np.array([(lo[0] + hi[0]) / 2, (lo[1] + hi[1]) / 2, lo[2]])
    return (s, off)


def studio_and_camera(P, panel, margin, frame_pts=None):
    """Lights, ground and camera for the (normalised) points P; the camera frames `frame_pts`."""
    c = (P.min(0) + P.max(0)) / 2
    rad = float(np.linalg.norm(P.max(0) - P.min(0))) / 2
    B.studio(c, rad, ground_z=float(P[:, 2].min()) - 0.02 * rad)
    B.fit_camera(P if frame_pts is None else frame_pts, azimuth=panel.get("az", -35),
                 elevation=panel.get("el", 22), lens=50, margin=margin)


def new_scene(W, H, samples):
    B.reset()
    B.setup(W, H, samples=samples)


def load_panel(panel, cache):
    return [Layer(d, cache) for d in panel["layers"]]


def main():
    a = B.args()
    mode, spec = a[0], json.load(open(a[1]))
    cache = {}
    if mode == "still":
        W = int(a[3]) if len(a) > 3 else 1600
        H = int(a[4]) if len(a) > 4 else 1000
        layers = load_panel(spec, cache)
        pts = framed_points(layers)
        norm = make_norm(pts)
        new_scene(W, H, int(os.environ.get("RS_SAMPLES", 128)))
        for i, l in enumerate(layers):
            l.add(norm, i)
        studio_and_camera((pts + norm[1]) * norm[0], spec, spec.get("margin", 1.1))
        B.render_still(a[2])
    elif mode == "sheet":
        outdir = Path(a[2])
        outdir.mkdir(parents=True, exist_ok=True)
        W = int(a[3]) if len(a) > 3 else 600
        H = int(a[4]) if len(a) > 4 else 600
        panels = [load_panel(p, cache) for p in spec["panels"]]
        diag = max(np.linalg.norm(np.ptp(framed_points(ls), axis=0)) for ls in panels)
        for i, (p, layers) in enumerate(zip(spec["panels"], panels)):
            pts = framed_points(layers)
            norm = make_norm(pts)
            P = (pts + norm[1]) * norm[0]
            fp = None
            if spec.get("common_scale"):
                # frame a cube of the common diagonal around this shape: equal scale in all panels
                norm = (2.0 / diag, make_norm(pts)[1])
                P = (pts + norm[1]) * norm[0]
                c = (P.min(0) + P.max(0)) / 2
                e = 1.0 / np.sqrt(3.0)
                fp = np.array([c - e, c + e])
            new_scene(W, H, int(os.environ.get("RS_SAMPLES", 96)))
            for k, l in enumerate(layers):
                l.add(norm, k)
            studio_and_camera(P, p, p.get("margin", spec.get("margin", 1.1)), fp)
            B.render_still(str(outdir / f"p{i:02d}.png"))
    elif mode == "anim":
        name = a[2]
        W = int(a[3]) if len(a) > 3 else 960
        H = int(a[4]) if len(a) > 4 else 720
        panels = [load_panel(p, cache) for p in spec["panels"]]
        allpts = np.concatenate([framed_points(ls) for ls in panels])
        if spec.get("track"):
            norms = [make_norm(framed_points(ls)) for ls in panels]
            allpts = np.concatenate([(framed_points(ls) + n[1]) * n[0] for ls, n in zip(panels, norms)])
            norm = (1.0, np.zeros(3))
        else:
            norms = None
            norm = make_norm(allpts)
        new_scene(W, H, int(os.environ.get("RS_SAMPLES", 64)))
        studio_and_camera((allpts + norm[1]) * norm[0], spec, spec.get("margin", 1.1))

        def update(i):
            clear_geometry()
            for k, l in enumerate(panels[i]):
                l.add(norms[i] if norms else norm, k)

        B.render_frames(name, len(panels), update)
    else:
        raise SystemExit("mode: still | sheet | anim")


main()
