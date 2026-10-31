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
w("")
w("Demo failures, integer backend: " + ", ".join(f"{i} ({op}: {demo_state('integer', i, op)})" for i in dem["integer"][0] for op in ("flip", "refine") if demo_state("integer", i, op) not in ("ok", "missing"))[:3000] + "\n")
w("Demo failures, signpost backend (first 60): " + ", ".join(f"{i} ({op}: {demo_state('signpost', i, op)})" for i in dem["signpost"][0] for op in ("flip", "refine") if demo_state("signpost", i, op) not in ("ok", "missing"))[:1500] + "\n")

# IDT: integer vs integer demo
cmpd = os.path.join(W, "cmp")
n = same_v = same_cs = same_vef = lens_ok = lens_n = 0; worst = 0.0; bad = []
for i in kept:
    c = os.path.join(cmpd, f"{i}.json")
    if not os.path.exists(c): continue
    try: o = json.load(open(c))
    except Exception: continue
    s = stats("integer", i, "flip")
    if o.get("status") != "ok" or not s or s["state"] != "ok": continue
    n += 1
    same_v += s["V"] == o["v0"]
    same_cs += s["cs"] == o["idt_sub_v"]
    vp = os.path.join(W, "demo/integer", f"{i}_flip_vef.txt")
    if os.path.exists(vp):
        dv = tuple(int(x) for x in open(vp).read().split())
        mine = (o["idt_sub_v"], o["idt_sub_tf"] + (o["idt_sub_e"] - o["idt_sub_f"] + o["idt_sub_f"]) * 0, o["idt_sub_tf"])
        # triangulated: F = tf, E = (3 tf + boundary)/2
        be = int(meta[i]["n_boundary_edges"])
        mine = (o["idt_sub_v"], (3 * o["idt_sub_tf"] + be) // 2, o["idt_sub_tf"])
        if dv == mine: same_vef += 1
        else: bad.append((i, "vef", dv, mine))
    lp = os.path.join(W, "demo/integer", f"{i}_flip_idt_lengths.txt"); op = os.path.join(cmpd, f"{i}.idt")
    if os.path.exists(lp) and os.path.exists(op):
        a = [float(x) for x in open(lp).read().split()]
        # the demo lists each edge twice (once per face), ours once per edge
        a = sorted(set(a)) if False else a
        b2 = [float(x) for x in open(op).read().split()]
        # demo: F x 3 lengths = each interior edge twice, boundary edges once
        from collections import Counter
        ca = Counter(a); cb = Counter(b2)
        da = sorted(ca); db = sorted(cb)
        lens_n += 1
        if len(da) == len(db):
            m = max(abs(x - y) / max(abs(y), 1e-300) for x, y in zip(da, db))
            worst = max(worst, m); lens_ok += m < 1e-9
            if m >= 1e-9: bad.append((i, "lengths", m))
        else:
            bad.append((i, "n_distinct_lengths", len(da), len(db)))
w("### IDT, integer coordinates (ours against the demo's integer backend)\n")
w(f"Meshes compared: {n} (both completed, up to 20000 faces for the edge lengths). `|V1|` (= `|V0|`) equal: {same_v}; `|V^S|` = `|V1| + Sum n+` equal: {same_cs}; subdivision V, E, F of the triangulated subdivision equal: {same_vef}; sorted sets of distinct edge lengths equal within 1e-9 relative: {lens_ok} of {lens_n} (worst relative difference {worst:.2e}).\n")
if bad: w("Differences (first 15): " + "; ".join(str(x) for x in bad[:15]) + "\n")

# Refinement statistics
rows = []
for i in kept:
    s = stats("integer", i, "refine"); o = ours("integer", i)
    if s and s["state"] == "ok" and o and o.get("status") == "ok" and o.get("refine_ok") == 1:
        rows.append((i, s, o))
def ratio(f): return [f(s, o) for _, s, o in rows]
r_v1 = ratio(lambda s, o: o["v1"] / s["V"]); r_cs = ratio(lambda s, o: o["refine_sub_v"] / s["cs"])
w("### Refinement, integer coordinates (statistics, not equality, I23)\n")
w(f"Meshes compared: {len(rows)}. `|V1|` ours/demo: mean {statistics.mean(r_v1):.4f}, min {min(r_v1):.3f}, max {max(r_v1):.3f}, equal on {sum(1 for _, s, o in rows if s['V'] == o['v1'])}; `|V^S|` ours/demo: mean {statistics.mean(r_cs):.4f}, min {min(r_cs):.3f}, max {max(r_cs):.3f}, equal on {sum(1 for _, s, o in rows if s['cs'] == o['refine_sub_v'])}.")
tot_o = sum(o["v1"] for _, s, o in rows); tot_d = sum(s["V"] for _, s, o in rows)
w(f"Totals: `Sum |V1|` ours {tot_o}, demo {tot_d} ({100*(tot_o/tot_d-1):+.2f} %); `Sum |V^S|` ours {sum(o['refine_sub_v'] for _,s,o in rows)}, demo {sum(s['cs'] for _,s,o in rows)}.")
amin = [(float(s['outputMinValidAngleDeg']), o['min_angle_valid_deg'], i) for i, s, o in rows]
w(f"Minimum angle over valid faces: demo below 25 degrees on {sum(1 for a,_,_ in amin if a < 25-1e-6)} meshes, ours on {sum(1 for _,b,_ in amin if b < 25-1e-6)}; both below on {sum(1 for a,b,_ in amin if a < 25-1e-6 and b < 25-1e-6)}.")
big = sorted(rows, key=lambda r: -abs(math.log(r[2]['v1'] / r[1]['V'])))[:10]
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
    o = ours("integer", i)
    ds = demo_state("integer", i, "refine")
    if o is None or ds == "missing": continue
    if (ds == "ok") != (o.get("status") == "ok" and o.get("refine_ok") == 1): dis.append((i, "demo " + ds, "ours " + (o.get("status") if o.get("status") != "ok" else ("ok" if o.get("refine_ok") == 1 else "invalid"))))
w(f"Meshes where the demo's integer refinement and ours disagree on success: {len(dis)}: {dis[:40]}\n")
# Timings (loaded)
tr = [(float(s["refinementDuration"]), o["refine_seconds"]) for _, s, o in rows if "refinementDuration" in s]
if tr:
    w(f"Refinement seconds under load (40-56 parallel jobs, not timings): demo median {statistics.median(a for a,_ in tr):.4f}, ours median {statistics.median(b for _,b in tr):.4f}; Sum demo {sum(a for a,_ in tr):.1f} s, Sum ours {sum(b for _,b in tr):.1f} s.\n")
out = "\n".join(L) + "\n"
p = os.path.join(W, "summary.md")
base = open(p).read().split("\n## Against the reference")[0] if os.path.exists(p) else ""
open(p, "w").write(base.rstrip("\n") + "\n" + out)
print(out)
