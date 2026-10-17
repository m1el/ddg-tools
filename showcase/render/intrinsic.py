"""Intrinsic triangulations: input edges (thin, dark) and the edges of the intrinsic
triangulation traced over the input surface (colored), inserted vertices as dots.

blender -b --python render/intrinsic.py -- <mesh> still <stage> <out.png> [w h]
blender -b --python render/intrinsic.py -- <mesh> anim
"""

import json
import math
import sys

sys.path.insert(0, __file__.rsplit("/", 1)[0])
import blib as B  # noqa: E402
import numpy as np  # noqa: E402

a = B.args()
name, mode = a[0], a[1]
D = B.load(f"intr_{name}")
zup = name == "torus"
P0 = D["positions"] if zup else B.y_up_to_z_up(D["positions"])
P, norm = B.normalize(P0, size=2.2)
T = D["triangles"]
NS = int(D["n_stages"])
labels = [str(x) for x in D["labels"]]
nin = len(P)


def stage(k):
    pts = B.apply_norm(D[f"s{k}_pts"] if zup else B.y_up_to_z_up(D[f"s{k}_pts"]), norm)
    off = D[f"s{k}_off"]
    lines = [pts[off[i]:off[i + 1]] for i in range(len(off) - 1)]
    verts = B.apply_norm(D[f"s{k}_verts"] if zup else B.y_up_to_z_up(D[f"s{k}_verts"]), norm)
    return lines, D[f"s{k}_straight"], verts


anim = mode == "anim"
W, H = (1280, 800) if anim else (int(a[4]), int(a[5])) if len(a) > 5 else (1800, 1200)
B.reset()
B.setup(W, H, samples=128)

surf_mat = B.mat_solid(B.lin((0.86, 0.85, 0.83)), roughness=0.55)
surf = B.mesh_object("surface", P, T, smooth=False, material=surf_mat)

# Vertex normals, to lift lines off the surface so tubes are not half buried.
from scipy.spatial import cKDTree
fn = np.cross(P[T[:, 1]] - P[T[:, 0]], P[T[:, 2]] - P[T[:, 0]])
VN = np.zeros_like(P)
for k in range(3):
    np.add.at(VN, T[:, k], fn)
VN /= np.linalg.norm(VN, axis=1, keepdims=True)
KD = cKDTree(P)


def lift(pts, h):
    _, idx = KD.query(pts, k=3)
    n = VN[idx].sum(1)
    n /= np.linalg.norm(n, axis=1, keepdims=True)
    return pts + h * n


# Input edges.
E = set()
for t in T:
    for i in range(3):
        u, v = int(t[i]), int(t[(i + 1) % 3])
        E.add((min(u, v), max(u, v)))
E = np.array(sorted(E))
size = 2.2
r_in = 0.0011 * size
r_int = 0.0026 * size
ink = B.mat_solid(B.lin((0.30, 0.31, 0.34)), roughness=0.6)
in_obj = B.tubes("input_edges", [lift(P[e], 0.6 * r_in) for e in E], r_in, ink, resolution=2)

int_mat = B.mat_solid(B.lin((0.93, 0.36, 0.08)), roughness=0.35)
dot_mat = B.mat_solid(B.lin((0.10, 0.45, 0.80)), roughness=0.3)
dyn = []


def show(k):
    for o in dyn:
        B.bpy.data.objects.remove(o, do_unlink=True)
    dyn.clear()
    lines, straight, verts = stage(k)
    moved = [l for l, s in zip(lines, straight) if not s]
    if moved:
        dyn.append(B.tubes("intrinsic", [lift(l, 0.9 * r_int) for l in moved], r_int, int_mat, resolution=3))
    if len(verts) > nin:
        dyn.append(B.spheres("inserted", lift(verts[nin:], r_int), 0.0075 * size, material=dot_mat, subdiv=2))
    for o in dyn:
        o.parent = pivot


ground = P[:, 2].min()
center = (P.min(0) + P.max(0)) / 2
radius = float(np.linalg.norm(P.max(0) - P.min(0))) / 2
B.studio(center, radius, ground_z=ground)
AZ, EL = {"torus": (-20, 42)}.get(name, (-35, 25))
# Room for a caption under the subject in the animation: frame the points plus a strip below.
Pf = P if not anim else np.vstack([P, P - [0, 0, 0.0]])
cam, ctr, dist = B.fit_camera(P, azimuth=AZ, elevation=EL, lens=55, margin=1.04 if not anim else 1.16)
panel = (not anim) and len(a) > 5
if panel:
    cam, ctr, dist = B.fit_camera(P, azimuth=AZ, elevation=EL, lens=55, margin=1.2)
    cam.data.shift_y = float(__import__("os").environ.get("SHIFT", 0.035))
if anim:
    cam.data.shift_y = -0.05
pivot = B.bpy.data.objects.new("pivot", None)
B.bpy.context.collection.objects.link(pivot)
surf.parent = pivot
in_obj.parent = pivot

if not anim:
    k = int(a[2])
    show(k)
    B.render_still(a[3])
else:
    # Every stage held a few frames; holds at input, Delaunay, refined, ODT.
    n_flip = sum(1 for l in labels if l.startswith("flip"))
    seq = []
    for k in range(NS):
        hold = 3
        if k == 0 or k == n_flip - 1 or k == NS - 2:
            hold = 40
        if k == NS - 1:
            hold = 50
        seq += [k] * hold
    mina = D["min_angle"]
    caps = []
    for k in seq:
        l = labels[k]
        if l.startswith("flip"):
            done, tot = l.split()[1].split("/")
            text = "input mesh" if done == "0" else f"flipping to intrinsic Delaunay   {done} / {tot} flips"
        elif l.startswith("refine"):
            done, tot = l.split()[1].split("/")
            text = f"intrinsic Delaunay refinement   +{done} vertices"
        else:
            text = f"intrinsic optimal Delaunay   +{len(D[f's{k}_verts']) - nin} vertices"
        caps.append(f"{text}      min angle {mina[k]:.1f}°")
    (B.WORK / f"labels_intrinsic_{name}.json").write_text(json.dumps(caps))
    # Render every stage once, then lay out the frames (holds are copies).
    import shutil
    uniq = B.FRAMES / f"intrinsic_{name}_stages"
    uniq.mkdir(parents=True, exist_ok=True)
    for k in range(NS):
        f = uniq / f"{k:04d}.png"
        if not f.exists():
            show(k)
            B.render_still(f)
    out = B.FRAMES / f"intrinsic_{name}"
    shutil.rmtree(out, ignore_errors=True)
    out.mkdir(parents=True)
    for i, k in enumerate(seq):
        shutil.copy(uniq / f"{k:04d}.png", out / f"{i:04d}.png")
