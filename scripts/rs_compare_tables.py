#!/usr/bin/env python3
"""Markdown tables of the black-box comparison (crates/ddg-repulsive-surfaces/COMPARISON.md).

usage: rs_compare_tables.py RUNDIR [scene ...]

Reads, per scene `S` and run `R` in RUNDIR (produced by `compare_reference run` and the reference's
`rsurfaces_headless`, see scripts/rs-compare-reference.sh): `ref_S/steps.csv`, `ours_S_R.csv`,
`final_ref_S.csv`, `final_ours_S_R.csv` (`compare_reference stats` of the final meshes) and
`stats_ref_S.csv`, `stats_ours_S_R.csv` (`stats` of the dumped meshes: self-intersection timeline).
"""
import csv
import os
import sys

run = sys.argv[1]
SC = sys.argv[2:] or ["EDT", "Hopf", "Borromean", "Bunny", "Handcuffs", "Walnut", "Lung", "Recon", "kbb", "kb2", "kb1"]
LIMITS = {x.split(":")[0]: int(x.split(":")[1]) for x in SC if ":" in x}
scenes = [x.split(":")[0] for x in SC]
CUR = {"limit": None}
MODES = ["compat", "default", "compat-noedge", "noedge"]


def rd(p):
    p = os.path.join(run, p)
    if not os.path.exists(p):
        return []
    try:
        rows = list(csv.DictReader(open(p)))
        lim = CUR["limit"]
        if lim is not None and rows and "iter" in rows[0] and os.path.basename(p) == "steps.csv":
            rows = [r for r in rows if int(r["iter"]) <= lim]
        if lim is not None and rows and "step" in rows[0] and os.path.basename(p).startswith("ours_"):
            rows = [r for r in rows if int(r["step"]) <= lim]
        return rows
    except Exception:
        return []


def last_row(rows):
    rows = [r for r in rows if r and (r.get("V") or "").isdigit()]
    return rows[-1] if rows else None


def inter_summary(rows):
    n = [r for r in rows if r.get("intersects") in ("0", "1", "0.0", "1.0")]
    bad = [r for r in n if r["intersects"].startswith("1")]
    return f"{len(bad)}/{len(n)}" if n else "n/a"


def curve_table(s, steps=(0, 1, 5, 10, 25, 50, 100, 150, 200, 300, 350, 600)):
    ref = rd(f"ref_{s}/steps.csv")
    runs = {"reference": {int(r["iter"]): float(r["total_objective"]) for r in ref}}
    for m in MODES:
        o = rd(f"ours_{s}_{m}.csv")
        if o:
            runs[f"ours {m}"] = {int(r["step"]): float(r["F"]) for r in o}
    last = max((max(v) for v in runs.values() if v), default=0)
    cols = [k for k in steps if k <= last]
    print(f"\nObjective F against step (each run's own Barnes-Hut value, penalties included):\n")
    print("| run | " + " | ".join(str(k) for k in cols) + " |")
    print("|---|" + "---|" * len(cols))
    for name, d in runs.items():
        print(f"| {name} | " + " | ".join(f"{d[k]:.5g}" if k in d else "-" for k in cols) + " |")


for s in scenes:
    CUR["limit"] = LIMITS.get(s)
    ref = rd(f"ref_{s}/steps.csv")
    oursets = {m: rd(f"ours_{s}_{m}.csv") for m in MODES}
    if not ref and not any(oursets.values()):
        continue
    print(f"\n#### {s}\n")
    print("| run | steps | ms/step | F (own, BH) | exact TP energy of the final mesh | V / F | area | volume | extents | sd1,sd2,sd3 | area drift, volume drift | dumped states intersecting |")
    print("|---|---|---|---|---|---|---|---|---|---|---|---|")
    base = next((o[0] for o in oursets.values() if o), None)
    a0 = float(base["area"]) if base else float("nan")
    v0 = float(base["volume"]) if base else float("nan")

    def drift(r):
        try:
            return f"{float(r['area'])/a0-1:+.1e}, {float(r['volume'])/v0-1:+.1e}"
        except Exception:
            return "-"

    fr = last_row(rd(f"final_ref_{s}.csv"))
    if ref:
        ms = [float(r["wall_ms"]) for r in ref[1:]] or [0.0]
        tot = float(ref[-1]["total_objective"])
        row = [f"reference", ref[-1]["iter"], f"{sum(ms)/len(ms):.0f}", f"{tot:.5g}"]
        if fr:
            row += [f"{float(fr['energy']):.5g}", f"{fr['V']} / {fr['F']}", f"{float(fr['area']):.5g}", f"{float(fr['volume']):.5g}",
                    f"{float(fr['bbox_x']):.2f}x{float(fr['bbox_y']):.2f}x{float(fr['bbox_z']):.2f}",
                    f"{float(fr['sd1']):.3f},{float(fr['sd2']):.3f},{float(fr['sd3']):.3f}"]
        else:
            row += [f"-", f"{ref[-1]['nV']} / {ref[-1]['nF']}", "-", "-", "-", "-"]
        row.append(drift(fr) if fr else "-")
        row.append(inter_summary(rd(f"stats_ref_{s}.csv")))
        print("| " + " | ".join(row) + " |")
    for m in MODES:
        o = oursets[m]
        if not o:
            continue
        ms = [float(r["step_ms"]) for r in o[1:]] or [0.0]
        step_inter = [r["intersects"] for r in o if r["intersects"] not in ("NaN", "")]
        fo = last_row(rd(f"final_ours_{s}_{m}.csv"))
        row = [f"ours {m}", o[-1]["step"], f"{sum(ms)/len(ms):.0f}", f"{float(o[-1]['F']):.5g}"]
        if fo:
            row += [f"{float(fo['energy']):.5g}", f"{fo['V']} / {fo['F']}", f"{float(fo['area']):.5g}", f"{float(fo['volume']):.5g}",
                    f"{float(fo['bbox_x']):.2f}x{float(fo['bbox_y']):.2f}x{float(fo['bbox_z']):.2f}",
                    f"{float(fo['sd1']):.3f},{float(fo['sd2']):.3f},{float(fo['sd3']):.3f}"]
        else:
            row += ["-", f"{o[-1]['V']} / {o[-1]['Fc']}", f"{float(o[-1]['area']):.5g}", f"{float(o[-1]['volume']):.5g}", "-", "-"]
        row.append(drift(fo) if fo else drift(o[-1]))
        if step_inter:
            bad = sum(1 for x in step_inter if x.startswith("1"))
            row.append(f"{bad}/{len(step_inter)} step ends")
        else:
            row.append(inter_summary(rd(f"stats_ours_{s}_{m}.csv")))
        print("| " + " | ".join(row) + " |")
    curve_table(s)
