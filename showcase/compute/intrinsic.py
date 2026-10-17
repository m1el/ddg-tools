"""Intrinsic triangulation stages (Rust helper `intrinsic`): flips to Delaunay, refinement,
optimal Delaunay, every edge traced over the input mesh.

intrinsic.py <mesh> <flip-steps> <refine-steps>
<mesh>: a dataset/model name, or 'torus' (a sheared torus grid, generated here)."""
import json, subprocess, sys
import numpy as np
from common import *


def sheared_torus(nu=44, nv=10, shear=4, R=1.0, r=0.42):
    """A torus grid whose quads are split along a diagonal that skips `shear` columns: the
    surface is smooth, the triangles are long slivers with obtuse angles."""
    pos = []
    for j in range(nv):
        for i in range(nu):
            u, v = 2 * np.pi * i / nu, 2 * np.pi * j / nv
            pos.append([(R + r * np.cos(v)) * np.cos(u), (R + r * np.cos(v)) * np.sin(u), r * np.sin(v)])
    idx = lambda i, j: (j % nv) * nu + (i % nu)
    tris = []
    for j in range(nv):
        # ring j to ring j+1, ring j+1 shifted by `shear` columns: a strip of slivers
        for i in range(nu):
            a, b = idx(i, j), idx(i + 1, j)
            c, d = idx(i + shear, j + 1), idx(i + 1 + shear, j + 1)
            tris += [[a, b, d], [a, d, c]]
    return np.array(pos), np.array(tris, dtype=np.int32)


name = sys.argv[1]
fs, rs = sys.argv[2], sys.argv[3]
extra = sys.argv[4:5]
if name == "torus":
    pos, tri = sheared_torus()
else:
    pos, polys, _ = load_obj(model(name))
    tri = triangles(polys)
mj = WORK / f"intr_{name}.mesh.json"
json.dump({"positions": pos.tolist(), "triangles": tri.tolist()}, open(mj, "w"))
out = WORK / f"intr_{name}.json"
subprocess.run([str(WORK / "target/release/showcase-data"), "intrinsic", str(mj), str(out), fs, rs, *extra], check=True)
st = json.load(open(out))["stages"]
arrs = dict(positions=pos, triangles=tri, n_stages=len(st))
for k, s in enumerate(st):
    lines = s["lines"]
    arrs[f"s{k}_pts"] = np.concatenate([np.array(l) for l in lines])
    arrs[f"s{k}_off"] = np.cumsum([0] + [len(l) for l in lines])
    arrs[f"s{k}_straight"] = np.array(s["straight"])
    arrs[f"s{k}_verts"] = np.array(s["vertices"])
    print(k, s["label"], "edges", len(lines), "verts", len(s["vertices"]),
          "min %.2f max %.2f" % (s["min_angle"], s["max_angle"]))
arrs["labels"] = np.array([s["label"] for s in st])
arrs["min_angle"] = np.array([s["min_angle"] for s in st])
arrs["max_angle"] = np.array([s["max_angle"] for s in st])
save(f"intr_{name}", **arrs)
