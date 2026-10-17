"""Optimal cones: raw optimum of the canonical problem over a sweep of lambda (Rust helper),
plus the pipeline's cones for the hero still."""
import json, subprocess, sys
import numpy as np
from common import *

name = sys.argv[1] if len(sys.argv) > 1 else "bunny"
pos, polys, _ = load_obj(model(name))
tri = triangles(polys)
mj = WORK / f"{name}.mesh.json"
json.dump({"positions": pos.tolist(), "triangles": tri.tolist()}, open(mj, "w"))
lams = np.geomspace(float(sys.argv[2]) if len(sys.argv) > 2 else 0.003, float(sys.argv[3]) if len(sys.argv) > 3 else 0.6, 72)
out = WORK / f"cones_raw_{name}.json"
subprocess.run([str(WORK / "target/release/showcase-data"), "cones", str(mj), str(out),
                *[repr(float(l)) for l in lams]], check=True)
runs = json.load(open(out))["runs"]
u = np.array([r["u"] for r in runs])
cv = [np.array(r["cone_vertices"], dtype=np.int32) for r in runs]
ca = [np.array(r["cone_angles"]) for r in runs]
off = np.cumsum([0] + [len(c) for c in cv])
save(f"cones_{name}", positions=pos, triangles=tri, lambdas=lams, u=u,
     cone_offsets=off, cone_vertices=np.concatenate(cv), cone_angles=np.concatenate(ca),
     distortion=np.array([r["distortion"] for r in runs]))
for l, c in zip(lams, cv):
    print(f"{l:.4f} {len(c)}")
