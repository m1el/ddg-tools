# /// script
# requires-python = ">=3.10"
# dependencies = ["numpy"]
# ///
"""One-thread check (SPEC 17n): the two reference methods, 1 thread, run alone (12 jobs, idle machine)
against the same runs in phase A (1 thread, 60 jobs in parallel) and phase B (64 threads).
  uv run scripts/rc_eval_check1.py runs_tmp/rc-eval/runA runs_tmp/rc-eval/runB runs_tmp/rc-eval/runB1"""
import glob, json, sys
import numpy as np


def load(d):
    out = {}
    for f in glob.glob(d + "/results/*.jsonl"):
        r = json.loads(open(f).readline())
        out[(r["knot"], r["method"])] = r
    return out


A, B, C = (load(d) for d in sys.argv[1:4])
em = lambda r: min(h[3] for h in r["history"])


def tt(r, t):
    for h in r["history"]:
        if h[3] <= t:
            return h[0], h[1]


print("| knot | method | target E | alone, 1 thr: it / s | phase A (60 parallel): it / s | phase B (64 thr): it / s | A / alone | alone / 64 thr | alone ended |")
print("|---|---|---:|---|---|---|---:|---:|---|")
rs = []
f = lambda x: f"{x[0]} / {x[1]:.0f}" if x else "-"
for (k, m), c in sorted(C.items()):
    t = 1.1 * em(A[(k, m)])
    ta, tb, tc = tt(A[(k, m)], t), tt(B[(k, m)], t), tt(c, t)
    ra = ta[1] / tc[1] if ta and tc else None
    rb = tc[1] / tb[1] if tb and tc else None
    print(f"| {k} | {m} | {t:.2f} | {f(tc)} | {f(ta)} | {f(tb)} | {ra:.2f} | {rb:.1f} | {c['end_reason']} {c['seconds']:.0f} s |" if ra and rb else f"| {k} | {m} | {t:.2f} | {f(tc)} | {f(ta)} | {f(tb)} | - | - | {c['end_reason']} |")
    if ra and rb:
        rs.append((ra, rb, m))
for m in ("hs-pg", "soft-hs-ncg"):
    x = [r for r in rs if r[2] == m]
    if x:
        print(f"\n{m}: median (phase A time / alone time) {np.median([r[0] for r in x]):.2f}, median (alone / 64 threads) {np.median([r[1] for r in x]):.2f}, n = {len(x)}")
