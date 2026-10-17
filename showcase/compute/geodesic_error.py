"""Error of the showcase's heat-method distances against the exact polyhedral distance
(Rust helper `exact`) for the same source vertex: mean |d - d*| / d* over the input vertices."""
import json, subprocess, sys
import numpy as np
from common import *

name = sys.argv[1]
D = np.load(WORK / f"geodesic_{name}.npz")
mj, out = WORK / f"{name}.mesh.json", WORK / f"geodesic_exact_{name}.json"
json.dump({"positions": D["positions"].tolist(), "triangles": D["triangles"].tolist()}, open(mj, "w"))
subprocess.run([str(WORK / "target/release/showcase-data"), "exact", str(mj), str(out), str(int(D["source"]))], check=True)
ex = np.array(json.load(open(out))["distances"])
m = ex > 0
res = {}
for t in ["input", "delaunay", "refined"]:
    res[t] = float(np.mean(np.abs(D[t][m] - ex[m]) / ex[m]))
    print(f"{t}: mean relative error {100 * res[t]:.2f} %")
json.dump(res, open(WORK / f"geodesic_error_{name}.json", "w"))
