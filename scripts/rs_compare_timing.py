#!/usr/bin/env python3
"""ms/step table of the timed pass (each run alone, 1 thread): rs_compare_timing.py RUNDIR/timed"""
import csv
import os
import statistics
import sys

d = sys.argv[1]
SCENES = ["EDT", "Hopf", "Handcuffs", "Walnut", "Lung", "Bunny", "Recon", "kbb"]


def col(path, key):
    if not os.path.exists(path):
        return None
    rows = list(csv.DictReader(open(path)))[1:]
    v = [float(r[key]) for r in rows]
    return v or None


print("| scene | faces | reference | ours compat | ours default | ours noedge | steps timed |")
print("|---|---|---|---|---|---|---|")
for s in SCENES:
    ref = col(f"{d}/ref_{s}/steps.csv", "wall_ms")
    cells = []
    for m in ["compat", "default", "noedge"]:
        v = col(f"{d}/ours_{s}_{m}.csv", "step_ms")
        cells.append(f"{statistics.mean(v):.0f} (med {statistics.median(v):.0f})" if v else "-")
    faces = ""
    try:
        faces = list(csv.DictReader(open(f"{d}/ref_{s}/steps.csv")))[0]["nF"]
    except Exception:
        pass
    r = f"{statistics.mean(ref):.0f} (med {statistics.median(ref):.0f})" if ref else "-"
    print(f"| {s} | {faces} | {r} | {cells[0]} | {cells[1]} | {cells[2]} | {len(ref) if ref else 0} |")
