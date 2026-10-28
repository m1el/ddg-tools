# /// script
# requires-python = ">=3.10"
# dependencies = ["numpy"]
# ///
"""Phase B summary (SPEC 17n): per method on the 24-knot subset, 64 threads, run alone.

  uv run scripts/rc_eval_phaseb.py runs_tmp/rc-eval/runB [runs_tmp/rc-eval/runA]

Time to 1.1 x E_ref (E_ref: the Hs reference's minimal corrected energy of the same knot; from phase B
itself, and with a second E_ref from the 1-thread phase A run when the directory is given), wall
seconds, speed relative to the reference, and end-reason counts.
"""
import glob, json, sys, collections
import numpy as np

REF = {"hard": "hs-pg", "soft": "soft-hs-ncg"}
GROUP = {"hard": ["hs-pg", "l2-pg", "h1-pg", "h2-pg", "hs-armijo"],
         "soft": ["soft-hs-ncg", "soft-hs-lbfgs", "soft-h1-ncg", "soft-h1-lbfgs"]}


def load(d):
    out = {}
    for f in glob.glob(d + "/results/*.jsonl"):
        r = json.loads(open(f).readline())
        out[(r["knot"], r["method"])] = r
    return out


def emin(r):
    return min(h[3] for h in r["history"] if h[3] is not None)


def tt(r, target):
    for h in r["history"]:
        if h[3] is not None and h[3] <= target:
            return h[0], h[1]
    return None


B = load(sys.argv[1])
A = load(sys.argv[2]) if len(sys.argv) > 2 else None
knots = sorted({k for k, _ in B})
print(f"{len(knots)} knots, {len(B)} runs")
for g, ms in GROUP.items():
    for ref_src, label in ((B, "E_ref from phase B reference"), (A, "E_ref from phase A (24-min, 1-thread) reference")):
        if ref_src is None:
            continue
        print(f"\n### {g}: target 1.1 x {label}\n")
        print("| method | reached | median s to target | median s (reached, Hs-ref also reached) | median iterations | ended: " + " / ".join(["energy_rate", "contact_floor", "stuck", "time"]) + " |")
        print("|---|---:|---:|---:|---:|---|")
        for m in ms:
            rows = []
            end = collections.Counter()
            for k in knots:
                if (k, m) not in B or (k, REF[g]) not in ref_src:
                    continue
                target = 1.1 * emin(ref_src[(k, REF[g])])
                hit = tt(B[(k, m)], target)
                refhit = tt(B[(k, REF[g])], target)
                rows.append((hit, refhit))
                end[B[(k, m)]["end_reason"]] += 1
            ok = [h for h, _ in rows if h]
            both = [(h, r) for h, r in rows if h and r]
            med = lambda v: f"{np.median(v):.1f}" if len(v) else "-"
            print(f"| {m} | {len(ok)}/{len(rows)} | {med([h[1] for h in ok])} | {med([h[1] for h, r in both])} | {med([h[0] for h in ok])} | "
                  + " / ".join(str(end[x]) for x in ["energy_rate", "contact_floor", "stuck", "time"]) + " |")
    # reference total
print("\nWall time by method (sum over knots, hours):", {m: round(sum(r["seconds"] for (k, mm), r in B.items() if mm == m) / 3600, 2) for m in sum(GROUP.values(), [])})
print("audit violations:", sum(r["audit"]["violations"] for r in B.values()), "invalid curves:", sum(bool(r.get("final_curve_invalid")) for r in B.values()))
