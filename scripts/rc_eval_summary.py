# /// script
# requires-python = ">=3.10"
# dependencies = ["numpy"]
# ///
"""Phase A per-method summary (SPEC 17n): end reasons, iterations and 1-thread seconds to 1.1 E_ref
(reached runs), final energy over E_ref, wall hours, back ends.  uv run scripts/rc_eval_summary.py runs_tmp/rc-eval/runA"""
import collections, glob, json, sys
import numpy as np

REF = {"hs-pg": "hs-pg", "l2-pg": "hs-pg", "h1-pg": "hs-pg", "h2-pg": "hs-pg", "hs-armijo": "hs-pg",
       "soft-hs-ncg": "soft-hs-ncg", "soft-hs-lbfgs": "soft-hs-ncg", "soft-h1-ncg": "soft-hs-ncg", "soft-h1-lbfgs": "soft-hs-ncg"}
R = {}
for f in glob.glob(sys.argv[1] + "/results/*.jsonl"):
    r = json.loads(open(f).readline())
    R[(r["knot"], r["method"])] = r
em = lambda r: min(h[3] for h in r["history"])
print("| method | runs | energy_rate | contact_floor | stuck (library) | time cap | reached 1.1 E_ref | median it | median s (1 thr, 60 parallel) | median Ê_final / E_ref | total core-hours | back ends |")
print("|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|")
for m, ref in REF.items():
    rows = [(k, r) for (k, mm), r in R.items() if mm == m]
    end = collections.Counter(r["end_reason"] for _, r in rows)
    hits, ratio = [], []
    for k, r in rows:
        t = 1.1 * em(R[(k, ref)])
        h = next(((x[0], x[1]) for x in r["history"] if x[3] <= t), None)
        if h:
            hits.append(h)
        ratio.append(r["history"][-1][3] / em(R[(k, ref)]))
    be = collections.Counter(r["energy_backend"] + "+" + r["solver_backend"] for _, r in rows)
    print(f"| {m} | {len(rows)} | {end['energy_rate']} | {end['contact_floor']} | {end['stuck']} | {end['time']} | {len(hits)} | "
          f"{np.median([h[0] for h in hits]) if hits else float('nan'):.0f} | {np.median([h[1] for h in hits]) if hits else float('nan'):.0f} | "
          f"{np.median(ratio):.3g} | {sum(r['seconds'] for _, r in rows) / 3600:.0f} | {dict(be)} |")
allr = list(R.values())
print("\naudit violations:", sum(r["audit"]["violations"] for r in allr), "| invalid final curves:", sum(bool(r.get("final_curve_invalid")) for r in allr),
      "| errors:", sum(r["status"] == "error" for r in allr), "| total core-hours:", round(sum(r["seconds"] for r in allr) / 3600))
