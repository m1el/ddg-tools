#!/usr/bin/env python3
"""Tables of the remeshing study (crates/ddg-uv/examples/remesh_study.rs).

Usage: uv run scripts/remesh_tables.py runs_tmp/remesh/real/*.log

Reads the `ROW` lines (one run and cut variant each, `key=value`) and prints, per mesh and over
all meshes: convergence (steps, stop reasons, rises, motion at the end), quality on the common
yardstick relative to the baseline of the same case (median and range of the relative
difference, number of cases where the treatment is lower), the effect of the bend merge, and
cost. Only the standard library is used.
"""

import statistics
import sys
from collections import defaultdict

TREATS = ["input", "flip", "flipS", "u2", "u2S", "u4", "u4S", "ka2", "ka2S", "ka4", "ka4S",
          "odt2", "odt4"]
COARSE = {"spot_control", "bob_control", "blub_control"}


def group(m):
    """coarse (control meshes), camel (the camel family, base length TARGET=area), fine."""
    if m in COARSE:
        return "coarse"
    if m.startswith("camel"):
        return "camel"
    return "fine"


GROUPS = ["coarse", "fine", "camel"]


def parse(paths):
    rows = []
    for path in paths:
        with open(path) as f:
            for line in f:
                if not line.startswith("ROW "):
                    continue
                d = {}
                for tok in line.split()[1:]:
                    if "=" in tok:
                        k, v = tok.split("=", 1)
                        d[k] = v
                rows.append(d)
    return rows


def num(d, k):
    try:
        return float(d[k])
    except (KeyError, ValueError):
        return float("nan")


def med(xs):
    xs = [x for x in xs if x == x]
    return statistics.median(xs) if xs else float("nan")


def fmt(x, p=3):
    if x != x:
        return "–"
    return f"{x:+.{p}f}"


def rel_table(rows, key, base, title, meshes):
    """Relative difference of `key` against treatment `base` of the same mesh and case."""
    by = {}
    for r in rows:
        if r.get("variant") != "cut" or "error" in r:
            continue
        by[(r["mesh"], r["case"], r["treat"])] = num(r, key)
    print(f"\n### {title}\n")
    print(f"relative difference of `{key}` against `{base}` of the same case: median "
          f"[min, max] in %, and cases lower than `{base}` / cases")
    print()
    header = "| treatment | " + " | ".join(meshes) + " | all |"
    print(header)
    print("|" + "---|" * (len(meshes) + 2))
    for t in TREATS:
        if t == base:
            continue
        cells = []
        allv = []
        for m in meshes + [None]:
            vals = []
            for (mm, case, tt), v in by.items():
                if tt != t or (m is not None and mm != m):
                    continue
                b = by.get((mm, case, base))
                # A run whose cut vanished has no energy (a closed surface without a cut is left
                # out of the energy): not compared, counted in the convergence table.
                if b is None or b != b or v != v or b == 0 or v == 0:
                    continue
                vals.append((v - b) / abs(b))
            if not vals:
                cells.append("")
                continue
            wins = sum(1 for x in vals if x < 0)
            cells.append(f"{100 * med(vals):+.1f} [{100 * min(vals):+.1f}, "
                         f"{100 * max(vals):+.1f}] {wins}/{len(vals)}")
        if any(cells):
            print(f"| {t} | " + " | ".join(cells) + " |")


def convergence(rows, meshes):
    print("\n### Convergence (variant cut; per treatment over all cases of a mesh group)\n")
    print("| treatment | group | runs | steps (median) |"
          " Converged / Stalled / Rose / Max / Vanished |"
          " rises (median) | rises / steps | motion at the end (input edges, median) |")
    print("|---|---|---|---|---|---|---|---|")
    for t in TREATS:
        for gname in GROUPS:
            rs = [r for r in rows if r.get("treat") == t and r.get("variant") == "cut"
                  and "error" not in r and group(r["mesh"]) == gname]
            if not rs:
                continue
            stops = defaultdict(int)
            for r in rs:
                stops[r["stop"]] += 1
            steps = [num(r, "steps") for r in rs]
            rises = [num(r, "rises") for r in rs]
            frac = [num(r, "rises") / max(num(r, "steps"), 1) for r in rs]
            print(f"| {t} | {gname} | {len(rs)} | {med(steps):.0f} | "
                  f"{stops['Converged']} / {stops['EnergyStalled']} / {stops['EnergyRose']} / "
                  f"{stops['MaxSteps']} / {stops['CutVanished']} | {med(rises):.0f} | {med(frac):.2f} | "
                  f"{med([num(r, 'motion_end') for r in rs]):.3f} |")


def by_alpha(rows, base, key):
    print(f"\n### By length weight: `{key}` against `{base}`, median in % (all meshes, both "
          "energies, both starts)\n")
    factors = sorted({r.get("alpha_factor") for r in rows if "alpha_factor" in r},
                     key=lambda x: -float(x))
    by = {}
    for r in rows:
        if r.get("variant") == "cut" and "error" not in r:
            by[(r["mesh"], r["case"], r["treat"])] = (num(r, key), r["alpha_factor"],
                                                       group(r["mesh"]))
    print("| treatment | " + " | ".join(f"α×{f} {g}" for f in factors for g in GROUPS) + " |")
    print("|" + "---|" * (len(GROUPS) * len(factors) + 1))
    for t in TREATS:
        if t == base:
            continue
        cells = []
        for f in factors:
            for coarse in GROUPS:
                vals = []
                for (m, case, tt), (v, ff, c) in by.items():
                    if tt != t or ff != f or c != coarse:
                        continue
                    b = by.get((m, case, base))
                    if b and b[0] == b[0] and v == v and b[0] != 0 and v != 0:
                        vals.append((v - b[0]) / abs(b[0]))
                cells.append(f"{100 * med(vals):+.1f}" if vals else "")
        if any(cells):
            print(f"| {t} | " + " | ".join(cells) + " |")


def merge_table(rows, meshes):
    print("\n### The bend merge (Session treatments): energy after the merge against the cut as "
          "it came back, relative, median [max] in 1e-4; bends removed (median); largest "
          "displacement (fraction of an edge, max)\n")
    by = {}
    for r in rows:
        if "error" not in r:
            by[(r["mesh"], r["case"], r["treat"], r["variant"])] = r
    print("| treatment | tolerance | spread yardstick | polyhedral yardstick | length |"
          " bends removed | displacement |")
    print("|---|---|---|---|---|---|---|")
    for t in ["flip", "odt2", "odt4"]:
        for v in ["m0.02", "m0.1"]:
            ds, dp, dl, nm, disp = [], [], [], [], []
            for (m, case, tt, vv), r in by.items():
                if tt != t or vv != v:
                    continue
                c = by.get((m, case, t, "cut"))
                if not c:
                    continue
                for key, out in (("ref_total_exact_len", ds), ("poly_total", dp),
                                 ("exact_length", dl)):
                    a, b = num(r, key), num(c, key)
                    if a == a and b == b and b != 0:
                        out.append((a - b) / abs(b))
                nm.append(num(r, "merged"))
                disp.append(num(r, "disp_rel"))
            if not ds:
                continue
            print(f"| {t} | {v[1:]} | {1e4 * med(ds):+.1f} [{1e4 * max(ds, key=abs):+.1f}] | "
                  f"{1e4 * med(dp):+.1f} [{1e4 * max(dp, key=abs):+.1f}] | "
                  f"{1e4 * med(dl):+.1f} [{1e4 * max(dl, key=abs):+.1f}] | {med(nm):.0f} | "
                  f"{max(disp):.3f} |")


def cost(rows, meshes):
    print("\n### Cost (from the parallel runs; times are not single-thread): vertices, build s, "
          "ms per step (median), run s (median)\n")
    print("| treatment | " + " | ".join(meshes) + " |")
    print("|" + "---|" * (len(meshes) + 1))
    for t in TREATS:
        cells = []
        for m in meshes:
            rs = [r for r in rows if r.get("treat") == t and r.get("mesh") == m
                  and r.get("variant") == "cut" and "error" not in r]
            if not rs:
                cells.append("")
                continue
            cells.append(f"{int(num(rs[0], 'V'))}, {med([num(r, 'build_s') for r in rs]):.2f}, "
                         f"{med([num(r, 'ms_step') for r in rs]):.1f}, "
                         f"{med([num(r, 'run_s') for r in rs]):.1f}")
        print(f"| {t} | " + " | ".join(cells) + " |")


def checks(rows):
    bad = [r for r in rows if "error" in r]
    conf = [r for r in rows if num(r, "conflicts") > 0 or num(r, "wrong_sides") > 0]
    lens = [abs(num(r, "ref_length") / num(r, "exact_length") - 1) for r in rows
            if "error" not in r and num(r, "exact_length") > 0]
    print(f"\nrows {len(rows)}, errors {len(bad)}, yardstick conflicts {len(conf)}, "
          f"reference length against exact length: median {med(lens):.2e}, max "
          f"{max(lens) if lens else float('nan'):.2e}")
    for r in bad[:20]:
        print("  error:", r.get("mesh"), r.get("case"), r.get("treat"), r.get("error")[:200])


def main():
    rows = parse(sys.argv[1:])
    meshes = sorted({r["mesh"] for r in rows}, key=lambda m: (GROUPS.index(group(m)), m))
    checks(rows)
    convergence(rows, meshes)
    for base in ["flip", "flipS"]:
        rel_table(rows, "ref_total_exact_len", base,
                  f"Quality, spread yardstick, against {base}", meshes)
        rel_table(rows, "poly_total", base, f"Quality, polyhedral yardstick, against {base}",
                  meshes)
    for base in ["flip", "flipS"]:
        rel_table(rows, "on_input", base, f"Energy of the cut restricted to the input vertices, "
                  f"on the input mesh, against {base}", meshes)
    rel_table(rows, "e_final", "flip", "The flow's own final energy against flip (not "
              "comparable across triangulations; for orientation)", meshes)
    by_alpha(rows, "flip", "ref_total_exact_len")
    by_alpha(rows, "flipS", "ref_total_exact_len")
    by_alpha(rows, "flip", "poly_total")
    merge_table(rows, meshes)
    cost(rows, meshes)


if __name__ == "__main__":
    main()
