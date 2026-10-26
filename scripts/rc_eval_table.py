# /// script
# requires-python = ">=3.10"
# dependencies = ["numpy", "matplotlib"]
# ///
"""Sec. 7 comparison tables and scatter plots (crates/ddg-repulsive-curves/SPEC.md 14.2).

  uv run scripts/rc_eval_table.py --iters-dir runs_tmp/rc-eval/runA [--times-dir runs_tmp/rc-eval/runB]
                                  --out runs_tmp/rc-eval/report [--factor 1.1]

Reads results/*.jsonl of eval_knots. Per knot and group (hard, soft) the reference method's minimal
energy E_ref (Hs projected gradient; Hs NCG); per method the iterations and wall time (from its own
history) to the first iterate with E <= factor * E_ref. A method *fails* on a knot when it never
gets there: `stuck` (the run ended with the line search stuck), `nonconvergent` (cap reached or
converged above the target), `error`. Knots on which every method of the group ends stuck are
excluded and counted (as in Fig. 19 / Supp. Figs. 2, 5). Writes failure_table.md, summary.json,
and the four scatter plots {hard,soft}_{iterations,time}.png (x: the reference method, y: the
method; failures on the dashed margins; below the diagonal the method is faster).
"""
import argparse
import glob
import json
import os
from collections import defaultdict

import numpy as np

GROUPS = {
    "hard": ("hs-pg", ["hs-pg", "l2-pg", "h1-pg", "h2-pg", "hs-armijo"]),
    "soft": ("soft-hs-ncg", ["soft-hs-ncg", "soft-hs-lbfgs", "soft-h1-ncg", "soft-h1-lbfgs"]),
}
LABEL = {
    "hs-pg": "H^s projected gradient (ref.)", "l2-pg": "L2 projected gradient",
    "h1-pg": "H1 projected gradient", "h2-pg": "H2 projected gradient",
    "hs-armijo": "H^s, normalized Armijo",
    "soft-hs-ncg": "H^s NCG (ref.)", "soft-hs-lbfgs": "H^s L-BFGS",
    "soft-h1-ncg": "H1 NCG", "soft-h1-lbfgs": "H1 L-BFGS",
}
# Paper, Sec. 7 / Fig. 19: failure percentages (stuck for Hs, failed for the others), where given.
PAPER = {
    "hs-pg": "6.5 % stuck", "h1-pg": "64.5 % failed", "h2-pg": "92.1 % failed",
    "l2-pg": ">= 97.7 % failed", "hs-armijo": "(not in paper)",
    "soft-hs-ncg": "3.3 % stuck", "soft-hs-lbfgs": "14.8 % stuck",
    "soft-h1-ncg": "(not in paper)", "soft-h1-lbfgs": "(not in paper)",
}


def load(d):
    out = {}
    for f in glob.glob(os.path.join(d, "results", "*.jsonl")):
        for line in open(f):
            line = line.strip()
            if line:
                r = json.loads(line)
                out[(r["knot"], r["method"])] = r
    return out


def reach(r, target):
    """(iteration, seconds) of the first history entry with E-hat <= target, else None."""
    for it, t, f, e, g in r["history"]:
        if e is not None and e <= target:
            return it, t
    return None


def analyse(A, T, factor):
    knots = sorted({k for k, _ in A})
    res = {}
    for gname, (ref, methods) in GROUPS.items():
        rows, excluded = [], []
        for k in knots:
            runs = {m: A[(k, m)] for m in methods if (k, m) in A}
            if ref not in runs or len(runs) < len(methods):
                continue  # incomplete knot
            if all(r["status"] == "stuck" for r in runs.values()):
                excluded.append(k)
                continue
            eref = min(e for _, _, _, e, _ in runs[ref]["history"] if e is not None)
            target = factor * eref
            row = {"knot": k, "eref": eref, "vertices": runs[ref]["vertices"], "m": {}}
            for m, r in runs.items():
                hit = reach(r, target)
                th = None
                if T is not None and (k, m) in T:
                    th = reach(T[(k, m)], factor * min(e for _, _, _, e, _ in T[(k, ref)]["history"]) if (k, ref) in T else target)
                row["m"][m] = {
                    "status": r["status"], "iters": r["iterations"], "seconds": r["seconds"],
                    "hit": hit, "hit_timed": th,
                    "violations": r["audit"]["violations"],
                    "final_ratio": min(e for _, _, _, e, _ in r["history"] if e is not None) / eref,
                }
            rows.append(row)
        res[gname] = {"rows": rows, "excluded": excluded}
    return res


def failure_class(c):
    if c["hit"] is not None:
        return "ok"
    return {"stuck": "stuck", "error": "error"}.get(c["status"], "nonconvergent")


def table(res):
    L = []
    for g, (ref, methods) in GROUPS.items():
        rows, ex = res[g]["rows"], res[g]["excluded"]
        L.append(f"### {g}: {len(rows)} knots compared, {len(ex)} excluded (every method stuck)\n")
        L.append("| method | knots | reached 1.1 E_ref | stuck | nonconvergent | error | failed total | "
                 "ended stuck (any) | audit violations | paper |")
        L.append("|---|---:|---:|---:|---:|---:|---:|---:|---:|---|")
        for m in methods:
            cs = [r["m"][m] for r in rows]
            n = max(len(cs), 1)
            cnt = defaultdict(int)
            for c in cs:
                cnt[failure_class(c)] += 1
            anystuck = sum(c["status"] == "stuck" for c in cs)
            viol = sum(c["violations"] for c in cs)
            pc = lambda x: f"{x} ({100 * x / n:.1f} %)"
            fail = cnt["stuck"] + cnt["nonconvergent"] + cnt["error"]
            L.append(f"| {LABEL[m]} | {len(cs)} | {pc(cnt['ok'])} | {pc(cnt['stuck'])} | "
                     f"{pc(cnt['nonconvergent'])} | {pc(cnt['error'])} | {pc(fail)} | {pc(anystuck)} | "
                     f"{viol} | {PAPER[m]} |")
        L.append("")
    L.append("`stuck`/`nonconvergent`/`error`: runs that never reach 1.1 E_ref, by how they ended. The reference "
             "method always reaches 1.1 E_ref of its own minimum, so its stuck rate (paper: 'stuck') is the "
             "column `ended stuck (any)`; for the other methods the paper's 'failed' is `failed total`.\n")
    return "\n".join(L)


def scatter(res, g, kind, outdir, have_times):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    ref, methods = GROUPS[g]
    others = [m for m in methods if m != ref]
    rows = res[g]["rows"]
    key = "hit_timed" if (kind == "time" and have_times) else "hit"
    idx = 0 if kind == "iterations" else 1
    pts = {}
    allv = []
    for m in others:
        P = []
        for r in rows:
            x = r["m"][ref][key]
            y = r["m"][m][key]
            P.append((None if x is None else max(x[idx], 1e-3), None if y is None else max(y[idx], 1e-3),
                      failure_class(r["m"][m])))
            allv += [v for v in P[-1][:2] if v is not None]
        pts[m] = P
    if not allv:
        return None
    lo, hi = min(allv) / 1.5, max(allv) * 1.5
    ylim = hi * 4  # failure margin
    fig, axs = plt.subplots(1, len(others), figsize=(4.2 * len(others), 4.4), squeeze=False)
    unit = "iterations" if kind == "iterations" else "seconds"
    for ax, m in zip(axs[0], others):
        ok = [(x, y) for x, y, c in pts[m] if x is not None and y is not None]
        if ok:
            ax.scatter(*zip(*ok), s=14, c="#1f77b4", label="reached")
        for cls, col, mk in (("stuck", "#d62728", "x"), ("nonconvergent", "#ff7f0e", "^"), ("error", "#7f7f7f", "s")):
            f = [(x, ylim) for x, y, c in pts[m] if x is not None and y is None and c == cls]
            if f:
                ax.scatter(*zip(*f), s=22, c=col, marker=mk, label=f"{cls} ({len(f)})")
        refail = [(lo, y) for x, y, c in pts[m] if x is None and y is not None]
        if refail:
            ax.scatter(*zip(*refail), s=14, c="#2ca02c", marker="<", label="reference failed")
        ax.plot([lo, hi], [lo, hi], "k:", lw=0.8)
        ax.axhline(ylim, color="k", ls="--", lw=0.6)
        ax.set_xscale("log"); ax.set_yscale("log")
        ax.set_xlim(lo, hi); ax.set_ylim(lo, ylim * 1.6)
        ax.set_xlabel(f"{LABEL[ref]}: {unit} to 1.1 E_ref")
        ax.set_ylabel(f"{LABEL[m]}: {unit}")
        ax.set_title(LABEL[m], fontsize=10)
        ax.legend(fontsize=7, loc="lower right")
    fig.tight_layout()
    p = os.path.join(outdir, f"{g}_{kind}.png")
    fig.savefig(p, dpi=130)
    plt.close(fig)
    return p


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--iters-dir", required=True)
    ap.add_argument("--times-dir")
    ap.add_argument("--out", required=True)
    ap.add_argument("--factor", type=float, default=1.1)
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    A = load(a.iters_dir)
    T = load(a.times_dir) if a.times_dir else None
    res = analyse(A, T, a.factor)
    md = table(res)
    open(os.path.join(a.out, "failure_table.md"), "w").write(md + "\n")
    print(md)
    plots = []
    for g in GROUPS:
        for kind in ("iterations", "time"):
            p = scatter(res, g, kind, a.out, T is not None)
            if p:
                plots.append(p)
    json.dump({g: {"compared": [r["knot"] for r in v["rows"]], "excluded": v["excluded"]} for g, v in res.items()},
              open(os.path.join(a.out, "summary.json"), "w"), indent=1)
    print("plots:", *plots)


main()
