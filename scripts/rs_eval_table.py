#!/usr/bin/env python3
"""WP16B: tables and plots from runs_tmp/wp16b/{fig12,fig5}_*.{csv,log}.

uv run --with matplotlib scripts/rs_eval_table.py [DIR]

Energies are the unweighted tangent-point energy (alpha, beta = 6, 12), which scales like 1/area, so the
soft-problem methods (area penalty, area may grow) are reported as E * A / A0 (A0 = initial area); `raw` columns
are not normalised. Fig. 12 values at fixed times: Barnes-Hut energy of the last step finished by that time;
`end` = exact energy at the end (summary line) when available.
"""
import csv, glob, os, re, sys

D = sys.argv[1] if len(sys.argv) > 1 else "runs_tmp/wp16b"
METHODS = ["hs", "l2", "h1", "h2", "aqp", "h1-lbfgs", "bcqn"]
NAMES = {"hs": "Hs", "l2": "L2", "h1": "H1", "h2": "H2", "aqp": "AQP", "h1-lbfgs": "H1-LBFGS", "bcqn": "BCQN"}


def load(prefix, m):
    p = f"{D}/{prefix}_{m}"
    if not os.path.exists(p + ".csv"):
        return None
    rows = list(csv.DictReader(open(p + ".csv")))
    summ = {}
    if os.path.exists(p + ".log"):
        for line in open(p + ".log"):
            if line.startswith("summary"):
                summ = dict(re.findall(r"(\w+)=(\S+)", line))
    return rows, summ


def norm(r, a0):
    return float(r["E_bh"]) * float(r["area"]) / a0


def at_time(rows, a0, t):
    best = None
    for r in rows:
        if float(r["wall"]) <= t:
            best = r
    return None if best is None else norm(best, a0)


def at_step(rows, a0, s):
    for r in rows:
        if int(r["step"]) == s:
            return norm(r, a0), float(r["E_bh"])
    return None, None


def f(x, w=9):
    return " " * (w - 1) + "-" if x is None else f"{x:{w}.2f}"


out = []
out.append("## Fig. 12 (E*A/A0 at wall time; end = exact energy at the end, normalised)\n")
for i in range(1, 6):
    out.append(f"### mesh {i}\n")
    out.append("| method | faces | steps | secs | 60 s | 300 s | 1200 s | end (BH) | end (exact) | ending | GMRES/solve |")
    out.append("|---|---|---|---|---|---|---|---|---|---|---|")
    for m in METHODS:
        r = load(f"fig12_{i}", m)
        if not r or not r[0]:
            out.append(f"| {NAMES[m]} | no data |")
            continue
        rows, s = r
        a0 = float(s.get("area0", rows[0]["area"]))
        ex = float(s["E_exact_end"]) * float(rows[-1]["area"]) / a0 if s and s.get("E_exact_end", "nan") != "nan" else None
        out.append(
            f"| {NAMES[m]} | {s.get('faces', '?')} | {s.get('steps', len(rows))} | {s.get('seconds', '?')} | "
            f"{f(at_time(rows, a0, 60))} | {f(at_time(rows, a0, 300))} | {f(at_time(rows, a0, 1200))} | "
            f"{f(norm(rows[-1], a0))} | {f(ex)} | {s.get('ending', 'running/unknown')} | {s.get('gmres_iters_per_solve', '-')} |"
        )
    out.append("")

out.append("## Fig. 5 (E*A/A0 after n steps, Barnes-Hut; same surface at three resolutions)\n")
out.append("| method | faces | steps done | secs | n=10 | n=50 | n=100 | n=300 | ratio E(300)/E(10) | s/step | ending | GMRES/solve |")
out.append("|---|---|---|---|---|---|---|---|---|---|---|---|")
for fcount in (11000, 40000, 135000):
    for m in METHODS:
        r = load(f"fig5_{fcount}", m)
        if not r or not r[0]:
            out.append(f"| {NAMES[m]} | {fcount} | no data |")
            continue
        rows, s = r
        a0 = float(s.get("area0", rows[0]["area"]))
        v = {n: at_step(rows, a0, n)[0] for n in (10, 50, 100, 300)}
        last = norm(rows[-1], a0)
        e10 = v[10]
        final = v[300] if v[300] is not None else last
        ratio = None if e10 is None else final / e10
        n = len(rows)
        out.append(
            f"| {NAMES[m]} | {s.get('faces0', fcount)} | {n} | {s.get('seconds', '?')} | {f(v[10])} | {f(v[50])} | {f(v[100])} | "
            f"{f(v[300])} | {f(ratio, 7)} | {float(rows[-1]['wall']) / n:.2f} | {s.get('ending', 'running/unknown')} | {s.get('gmres_iters_per_solve', '-')} |"
        )
text = "\n".join(out)
print(text)
open(f"{D}/tables.md", "w").write(text + "\n")

try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
except ImportError:
    sys.exit(0)

fig, axs = plt.subplots(1, 5, figsize=(24, 4.5))
for i in range(1, 6):
    for m in METHODS:
        r = load(f"fig12_{i}", m)
        if r and r[0]:
            rows, s = r
            a0 = float(s.get("area0", rows[0]["area"]))
            axs[i - 1].plot([float(x["wall"]) for x in rows], [norm(x, a0) for x in rows], label=NAMES[m])
    axs[i - 1].set_xscale("log"); axs[i - 1].set_yscale("log"); axs[i - 1].set_title(f"Fig. 12 mesh {i}")
    axs[i - 1].set_xlabel("wall seconds (1 thread)")
axs[0].set_ylabel("E * A / A0"); axs[0].legend()
fig.tight_layout(); fig.savefig(f"{D}/fig12.png", dpi=110)

fig, axs = plt.subplots(1, 3, figsize=(15, 4.5))
for k, fcount in enumerate((11000, 40000, 135000)):
    for m in METHODS:
        r = load(f"fig5_{fcount}", m)
        if r and r[0]:
            rows, s = r
            a0 = float(s.get("area0", rows[0]["area"]))
            axs[k].plot([int(x["step"]) for x in rows], [norm(x, a0) for x in rows], label=NAMES[m])
    axs[k].set_yscale("log"); axs[k].set_title(f"Fig. 5, {fcount} faces"); axs[k].set_xlabel("step")
axs[0].set_ylabel("E * A / A0"); axs[0].legend()
fig.tight_layout(); fig.savefig(f"{D}/fig5.png", dpi=110)
