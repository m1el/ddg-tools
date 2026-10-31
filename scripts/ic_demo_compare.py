"""WP10, SPEC 10.3: our runs against the authors' demo (runs_tmp/wp10/demo/ from scripts/ic-demo-eval.sh).
Appends the section 'Against the reference' to runs_tmp/wp10/summary.md (also printed)."""
import csv, glob, json, os, statistics, sys, collections, math
root = sys.argv[1]
W = os.path.join(root, "runs_tmp/wp10")
meta = {r["id"]: r for r in csv.DictReader(open(os.path.join(root, "data/thingi10k/converted.tsv")), delimiter="\t")}
kept = [i for i, r in meta.items() if r["status"] == "kept"]

def stats(b, i, op):
    p = os.path.join(W, "demo", b, f"{i}_{op}_stats.tsv")
    if not os.path.exists(p): return None
    lines = open(p).read().splitlines()
    if len(lines) < 2: return {"state": "crash"}
    h = [x.strip() for x in lines[0].split("\t")]; v = lines[1].split("\t")
    d = dict(zip(h, v))
    try:
        d["cs"] = int(float(d["commonSubdivisionVertices"])); d["V"] = int(d["outputVertices"])
    except Exception:
        return {"state": "crash"}
    d["state"] = "ok" if d["cs"] >= 0 else "crash"
    return d

def demo_state(b, i, op):
    s = stats(b, i, op)
    if s is None:
        return "missing" if not os.path.exists(os.path.join(W, "demo", b, f"{i}_{op}_done")) else "crash"
    if s["state"] == "ok": return "ok"
    log = os.path.join(W, "demo", b, f"{i}_{op}_out.log")
    txt = open(log).read() if os.path.exists(log) else ""
    if "exit=137" in txt or "exit=124" in txt: return "timeout"
    return "crash"

def ours(b, i):
    p = os.path.join(W, "jobs", b, f"{i}.json")
    return json.load(open(p)) if os.path.exists(p) else None

L = []; w = L.append
w("\n## Against the reference (the authors' demo, `--noGUI`, headless build)\n")
w("Demo runs as in `benchmark/run_benchmark.py` (`--refineMaxInsertions=0`, one process per mesh, 600 s cap). The demo's integer backend mollifies by 1e-5 itself; its signpost backend does not mollify (geometry-central's signpost constructor takes the input as it is), unlike our signposts which use the same mollified lengths as our integers.\n")
w("| run | meshes | success IDT | success refinement |")
w("|---|---|---|---|")
dem = {}
for b in ("integer", "signpost"):
    ids = [i for i in kept if os.path.exists(os.path.join(W, "demo", b, f"{i}_flip_done")) or os.path.exists(os.path.join(W, "demo", b, f"{i}_refine_done"))]
    fl = collections.Counter(demo_state(b, i, "flip") for i in ids); rf = collections.Counter(demo_state(b, i, "refine") for i in ids)
    dem[b] = (ids, fl, rf)
    n = len(ids)
    w(f"| demo, {b} | {n} | {100*fl['ok']/max(n,1):.2f} % ({fl['ok']}; {dict(fl)}) | {100*rf['ok']/max(n,1):.2f} % ({rf['ok']}; {dict(rf)}) |")
    for label, sel in (("no boundary", lambda i: meta[i]["n_boundary_edges"] == "0"), ("with boundary", lambda i: meta[i]["n_boundary_edges"] != "0")):
        sub = [i for i in ids if sel(i)]
        a = sum(demo_state(b, i, "flip") == "ok" for i in sub); c = sum(demo_state(b, i, "refine") == "ok" for i in sub)
        w(f"| demo, {b}, {label} | {len(sub)} | {100*a/max(len(sub),1):.2f} % ({a}) | {100*c/max(len(sub),1):.2f} % ({c}) |")
w("")
w("Demo failures, integer backend: " + ", ".join(f"{i} ({op}: {demo_state('integer', i, op)})" for i in dem["integer"][0] for op in ("flip", "refine") if demo_state("integer", i, op) not in ("ok", "missing"))[:3000] + "\n")
w("Demo failures, signpost backend (first 60): " + ", ".join(f"{i} ({op}: {demo_state('signpost', i, op)})" for i in dem["signpost"][0] for op in ("flip", "refine") if demo_state("signpost", i, op) not in ("ok", "missing"))[:1500] + "\n")

# IDT: integer vs integer demo, separately for meshes without Delaunay ties and with
cmpd = os.path.join(W, "cmp")
from collections import Counter
res = {0: Counter(), 1: Counter()}; worst = {0: 0.0, 1: 0.0}; bad = []
for i in kept:
    c = os.path.join(cmpd, f"{i}.json")
    if not os.path.exists(c): continue
    try: o = json.load(open(c))
    except Exception: continue
    s = stats("integer", i, "flip")
    if o.get("status") != "ok" or not s or s["state"] != "ok": continue
    t = 1 if o["idt_ties"] > 0 else 0
    R = res[t]; R["n"] += 1
    R["V1 equal"] += s["V"] == o["v0"]
    R["V^S equal"] += s["cs"] == o["idt_sub_v"]
    vp = os.path.join(W, "demo/integer", f"{i}_flip_vef.txt")
    if os.path.exists(vp):
        dv = tuple(int(x) for x in open(vp).read().split())
        be = int(meta[i]["n_boundary_edges"])
        mine = (o["idt_sub_v"], (3 * o["idt_sub_tf"] + be) // 2, o["idt_sub_tf"])
        R["vef n"] += 1
        if dv == mine: R["vef equal"] += 1
        elif t == 0: bad.append((i, "vef", dv, mine))
    lp = os.path.join(W, "demo/integer", f"{i}_flip_idt_lengths.txt"); op = os.path.join(cmpd, f"{i}.idt")
    if os.path.exists(lp) and os.path.exists(op):
        def distinct(xs):
            # distinct values up to 1e-11 relative (the demo writes 16 digits; equal lengths of a
            # symmetric mesh may differ in the last place)
            xs = sorted(xs); out = []
            for x in xs:
                if out and abs(x - out[-1]) <= 1e-11 * abs(x): continue
                out.append(x)
            return out
        da = distinct(float(x) for x in open(lp).read().split())
        db = distinct(float(x) for x in open(op).read().split())
        R["lengths n"] += 1
        # Two-pointer matching of the sorted sets within a relative tolerance: the number of
        # lengths of either set that have no partner (a different flip of one tie shifts a
        # positional comparison, so the sets are matched, not zipped).
        def unmatched(tol):
            i1 = j1 = 0; miss = 0
            while i1 < len(da) and j1 < len(db):
                x, y = da[i1], db[j1]
                if abs(x - y) <= tol * max(abs(x), abs(y)): i1 += 1; j1 += 1
                elif x < y: miss += 1; i1 += 1
                else: miss += 1; j1 += 1
            return miss + (len(da) - i1) + (len(db) - j1)
        u9 = unmatched(1e-9); u5 = unmatched(1e-5)
        R["lengths equal (rel 1e-9)"] += u9 == 0
        R["lengths equal (rel 1e-5)"] += u5 == 0
        R["unmatched total"] += u9; R["distinct total"] += len(db)
        worst[t] = max(worst[t], u9 / max(len(db), 1))
        if u9 and t == 0: bad.append((i, "unmatched lengths", u9, len(db)))
w("### IDT, integer coordinates (ours against the demo's integer backend)\n")
w("The IDT is unique up to flips of ties (cocircular quads, common in CAD meshes), and the normal coordinates, `Sum n+` and the subdivision depend on which way a tie is flipped. `ties` below is the number of interior edges within the Delaunay criterion's tolerance of a tie at the end of our flipping. The sets of distinct edge lengths are compared (the demo writes each interior edge twice).\n")
w("| subset | meshes | `V1` equal | `V^S` equal | subdivision V/E/F (triangulated) equal | lengths equal 1e-9 | lengths equal 1e-5 | unmatched distinct lengths (1e-9) |")
w("|---|---|---|---|---|---|---|---|")
for t, name in ((0, "no ties"), (1, "with ties")):
    R = res[t]
    w(f"| {name} | {R['n']} | {R['V1 equal']} | {R['V^S equal']} | {R['vef equal']} of {R['vef n']} | {R['lengths equal (rel 1e-9)']} of {R['lengths n']} | {R['lengths equal (rel 1e-5)']} | {R['unmatched total']} of {R['distinct total']} ({100*worst[t]:.1f} % worst mesh) |")
w("")
if bad: w("Differences among the meshes without ties (first 20): " + "; ".join(str(x) for x in bad[:20]) + "\n")

# Refinement statistics, for both of our Delaunay tests
for ob, label in (("integer-ref", "authors' Delaunay test (cotan weight < -1e-6)"), ("integer", "default Delaunay test (angle sum, 1e-10)")):
    rows = []
    for i in kept:
        s = stats("integer", i, "refine"); o = ours(ob, i)
        if s and s["state"] == "ok" and o and o.get("status") == "ok" and o.get("refine_ok") == 1:
            rows.append((i, s, o))
    r_v1 = [o["v1"] / s["V"] for _, s, o in rows]; r_cs = [o["refine_sub_v"] / s["cs"] for _, s, o in rows]
    w(f"### Refinement, integer coordinates, ours with the {label} (statistics, not equality, I23)\n")
    w(f"Meshes compared: {len(rows)}. `|V1|` ours/demo: mean {statistics.mean(r_v1):.4f}, median {statistics.median(r_v1):.4f}, min {min(r_v1):.3f}, max {max(r_v1):.3f}, equal on {sum(1 for _, s, o in rows if s['V'] == o['v1'])}, within 1 % on {sum(1 for x in r_v1 if abs(x-1) <= 0.01)}, within 5 % on {sum(1 for x in r_v1 if abs(x-1) <= 0.05)}; `|V^S|` ours/demo: mean {statistics.mean(r_cs):.4f}, median {statistics.median(r_cs):.4f}, min {min(r_cs):.3f}, max {max(r_cs):.3f}, equal on {sum(1 for _, s, o in rows if s['cs'] == o['refine_sub_v'])}.")
    tot_o = sum(o["v1"] for _, s, o in rows); tot_d = sum(s["V"] for _, s, o in rows)
    w(f"Totals: `Sum |V1|` ours {tot_o}, demo {tot_d} ({100*(tot_o/tot_d-1):+.2f} %); `Sum |V^S|` ours {sum(o['refine_sub_v'] for _,s,o in rows)}, demo {sum(s['cs'] for _,s,o in rows)}.")
    amin = [(float(s['outputMinValidAngleDeg']), o['min_angle_valid_deg'], i) for i, s, o in rows]
    w(f"Minimum angle over valid faces: demo below 25 degrees on {sum(1 for a,_,_ in amin if a < 25-1e-6)} meshes ({[i for a,_,i in amin if a < 25-1e-6][:12]}), ours on {sum(1 for _,b,_ in amin if b < 25-1e-6)} ({[i for _,b,i in amin if b < 25-1e-6][:12]}).")
    big = sorted(rows, key=lambda r: -abs(math.log(r[2]['v1'] / r[1]['V'])))[:8]
    w("Largest `|V1|` ratios: " + ", ".join(f"{i} ({o['v1']}/{s['V']})" for i, s, o in big) + "\n")
# signposts demo vs ours
for b, ob in (("signpost", "sp-ruppert"),):
    rs = []
    for i in kept:
        s = stats(b, i, "refine"); o = ours(ob, i)
        if s and s["state"] == "ok" and o and o.get("status") == "ok" and o.get("refine_ok") == 1: rs.append((i, s, o))
    if rs:
        r = [o["v1"] / s["V"] for _, s, o in rs]
        w(f"Signposts, demo against our Ruppert-rule signposts ({len(rs)} meshes both completed): `|V1|` ratio mean {statistics.mean(r):.4f}, equal on {sum(1 for _,s,o in rs if s['V']==o['v1'])}.\n")
# Disagreements in success
dis = []
for i in kept:
    o = ours("integer-ref", i)
    ds = demo_state("integer", i, "refine")
    if o is None or ds == "missing": continue
    if (ds == "ok") != (o.get("status") == "ok" and o.get("refine_ok") == 1): dis.append((i, "demo " + ds, "ours " + (o.get("status") if o.get("status") != "ok" else ("ok" if o.get("refine_ok") == 1 else "invalid"))))
w(f"Meshes where the demo's integer refinement and ours disagree on success: {len(dis)}: {dis[:40]}\n")
# Timings (loaded)
rows = [(i, stats("integer", i, "refine"), ours("integer-ref", i)) for i in kept]
rows = [r for r in rows if r[1] and r[1]["state"] == "ok" and r[2] and r[2].get("status") == "ok"]
tr = [(float(s["refinementDuration"]), o["refine_seconds"]) for _, s, o in rows if "refinementDuration" in s]
if tr:
    w(f"Refinement seconds under load (40-56 parallel jobs, not timings): demo median {statistics.median(a for a,_ in tr):.4f}, ours median {statistics.median(b for _,b in tr):.4f}; Sum demo {sum(a for a,_ in tr):.1f} s, Sum ours {sum(b for _,b in tr):.1f} s.\n")
out = "\n".join(L) + "\n"
p = os.path.join(W, "summary.md")
base = open(p).read().split("\n## Against the reference")[0] if os.path.exists(p) else ""
open(p, "w").write(base.rstrip("\n") + "\n" + out)
print(out)
