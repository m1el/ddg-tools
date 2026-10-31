"""WP10: summarizes runs_tmp/wp10/jobs/*/*.json into runs_tmp/wp10/summary.md (Table 1 etc.)."""
import csv, json, os, sys, glob, statistics, collections

root = sys.argv[1]
out_path = os.path.join(root, "runs_tmp/wp10/summary.md")
meta = {r["id"]: r for r in csv.DictReader(open(os.path.join(root, "data/thingi10k/converted.tsv")), delimiter="\t")}
kept = [i for i, r in meta.items() if r["status"] == "kept"]
BACKENDS = [("integer", "integer coordinates (Chew)"), ("sp-ruppert", "signposts (Ruppert rules)"), ("sp-chew", "signposts (Chew rules)")]
data = {}
for b, _ in BACKENDS:
    d = {}
    for f in glob.glob(os.path.join(root, "runs_tmp/wp10/jobs", b, "*.json")):
        try:
            r = json.load(open(f))
        except Exception:
            r = {"id": os.path.basename(f)[:-5], "status": "crash", "error": "unreadable result"}
        d[r["id"]] = r
    data[b] = d

def closed(i): return meta[i]["n_boundary_edges"] == "0"
def cat(r):
    if r is None: return "missing"
    if r["status"] != "ok": return r["status"]
    return "ok"
def idt_ok(r): return r is not None and r["status"] == "ok" and r.get("idt_ok") == 1
def ref_ok(r): return r is not None and r["status"] == "ok" and r.get("refine_ok") == 1
def pct(n, d): return f"{100.0 * n / d:.2f} %" if d else "-"
def q(xs, p):
    xs = sorted(xs)
    return xs[min(len(xs) - 1, int(round(p * (len(xs) - 1))))] if xs else float("nan")

L = []
w = L.append
w("# WP10 summary: Thingi10k evaluation\n")
w(f"Meshes: {len(kept)} kept by `thingi_convert` ({sum(closed(i) for i in kept)} without boundary, {sum(not closed(i) for i in kept)} with). Mollified, `eps = 1e-5 * mean edge`; IDT then refinement to 25 degrees; success = the common subdivision is recovered and passes the N12-style checks (valid fan-triangulated `TriMesh`, Euler characteristic, components and boundary loops of the input).\n")
w("## Table 1 as reproduced\n")
w("Paper (7696 meshes): explicit overlay 100 % / -, signposts 96.0 % / 69.1 %, integer coordinates 100 % / 100 %.\n")
w("| representation | subset | meshes | IDT success | refinement success |")
w("|---|---|---|---|---|")
for b, name in BACKENDS:
    for label, sel in (("all", lambda i: True), ("no boundary", closed), ("with boundary", lambda i: not closed(i))):
        ids = [i for i in kept if sel(i)]
        a = sum(idt_ok(data[b].get(i)) for i in ids); c = sum(ref_ok(data[b].get(i)) for i in ids)
        w(f"| {name} | {label} | {len(ids)} | {pct(a, len(ids))} ({a}) | {pct(c, len(ids))} ({c}) |")
w("")
w("## Outcomes per representation\n")
for b, name in BACKENDS:
    c = collections.Counter(cat(data[b].get(i)) for i in kept)
    w(f"* {name}: " + ", ".join(f"{k} {v}" for k, v in sorted(c.items())))
w("")
w("## Failures (every one, with cause)\n")
for b, name in BACKENDS:
    bad = [(i, data[b].get(i)) for i in kept if not (idt_ok(data[b].get(i)) and ref_ok(data[b].get(i)))]
    w(f"### {name}: {len(bad)} meshes\n")
    groups = collections.defaultdict(list)
    for i, r in bad:
        if r is None: key = "missing result"
        elif r["status"] != "ok": key = f"{r['status']} (exit {r.get('exit')})"
        else:
            e = r.get("idt_error") or r.get("refine_error") or ""
            stage = "idt" if r.get("idt_ok") != 1 else "refine"
            import re
            e = re.sub(r"\d+(\.\d+)?(e-?\d+)?", "N", e)[:140]
            key = f"{stage}: {e}"
        groups[key].append(i)
    for k, ids in sorted(groups.items(), key=lambda kv: -len(kv[1])):
        shown = ", ".join(sorted(ids, key=int)[:40]) + (" ..." if len(ids) > 40 else "")
        w(f"* **{k}** ({len(ids)}): {shown}")
    w("")
w("## Sizes (successful refinements)\n")
w("Paper: `|V^S|/|V0|` after refinement mean 20, 95th percentile 45; refined `|V1|/|V0|` 3.7 and 7.8.\n")
w("| representation | meshes | `|V^S|/|V0|` mean | p95 | `|V1|/|V0|` mean | p95 | IDT `|V^S|/|V0|` mean | p95 |")
w("|---|---|---|---|---|---|---|---|")
for b, name in BACKENDS:
    rs = [data[b][i] for i in kept if ref_ok(data[b].get(i)) and idt_ok(data[b].get(i))]
    s = [r["refine_sub_v"] / r["v0"] for r in rs]; v = [r["v1"] / r["v0"] for r in rs]; si = [r["idt_sub_v"] / r["v0"] for r in rs]
    if rs:
        w(f"| {name} | {len(rs)} | {statistics.mean(s):.2f} | {q(s, .95):.2f} | {statistics.mean(v):.2f} | {q(v, .95):.2f} | {statistics.mean(si):.2f} | {q(si, .95):.2f} |")
w("")
w("## Counters (integer coordinates, over all completed meshes)\n")
rs = [data["integer"][i] for i in kept if data["integer"].get(i, {}).get("status") == "ok"]
def tot(k): return sum(r.get(k, 0) or 0 for r in rs)
def cnt(k): return sum(1 for r in rs if (r.get(k, 0) or 0) > 0)
w("| counter | total | meshes affected |")
w("|---|---|---|")
for k, label in (("inserted", "vertices inserted"), ("removals", "vertices removed"), ("removals_refused", "removals refused (I8, I27)"),
                 ("clamped_crossings", "clamped crossings in refinement (I5)"), ("idt_clamps", "clamped crossings, IDT subdivision"),
                 ("refine_sub_clamps", "clamped crossings, refined subdivision"), ("trace_failures", "trace failures"),
                 ("face_split_refused", "face splits refused (App. A recovery, I6; includes the I28 case)"), ("places_moved", "places moved into their face (I16)"),
                 ("near_arc", "face splits within 1e-9 of an arc"), ("n_below_bound", "corners below 25 degrees at the end"), ("insertion_limit", "insertion limit reached"),
                 ("n_not_delaunay", "interior edges not Delaunay at the end"), ("n_encroached", "encroached boundary edges at the end")):
    w(f"| {label} | {tot(k)} | {cnt(k)} |")
w(f"| certificate failures (IDT) | {sum(1 for r in rs if r.get('idt_certificate_ok') != 1)} | |")
w(f"| certificate failures (refinement) | {sum(1 for r in rs if r.get('refine_certificate_ok') != 1)} | |")
mv = [r["min_angle_valid_deg"] for r in rs if r.get("min_angle_valid_deg") is not None]
w(f"\nMinimum angle over valid faces (integers): min {min(mv):.3f}, 1st percentile {q(mv, .01):.3f}, meshes below 25 degrees: {sum(1 for x in mv if x < 25 - 1e-9)} of {len(mv)}.\n")
w("## Counters (signposts, Chew and Ruppert rules)\n")
for b in ("sp-ruppert", "sp-chew"):
    rs2 = [data[b][i] for i in kept if data[b].get(i, {}).get("status") == "ok"]
    mv2 = [r["min_angle_valid_deg"] for r in rs2 if r.get("min_angle_valid_deg") is not None]
    w(f"* {b}: removals refused {sum(r.get('removals_refused', 0) or 0 for r in rs2)}, trace failures {sum(r.get('trace_failures', 0) or 0 for r in rs2)}, snapped {sum(r.get('snapped', 0) or 0 for r in rs2)}, meshes with min valid angle below 25: {sum(1 for x in mv2 if x < 25 - 1e-9)}, insertion limit {sum(1 for r in rs2 if r.get('insertion_limit') == 1)}")
w("")
w("## Cost (jobs ran 56 at a time: the seconds are for orientation, not timings)\n")
w("| representation | median s | p95 s | max s (mesh) | max peak RSS MB (mesh) |")
w("|---|---|---|---|---|")
for b, name in BACKENDS:
    rs3 = [r for r in data[b].values() if r.get("status") == "ok"]
    if not rs3: continue
    t = [r["seconds_total"] for r in rs3]
    mt = max(rs3, key=lambda r: r["seconds_total"]); mm = max(rs3, key=lambda r: r["peak_rss_mb"])
    w(f"| {name} | {statistics.median(t):.3f} | {q(t, .95):.2f} | {mt['seconds_total']:.1f} ({mt['id']}) | {mm['peak_rss_mb']} ({mm['id']}) |")
w("")
w("## Integers against signposts (same operations?)\n")
same = diff = 0; ex = []
for i in kept:
    a = data["integer"].get(i); b2 = data["sp-chew"].get(i)
    if ref_ok(a) and ref_ok(b2):
        if a["v1"] == b2["v1"] and a["refine_sub_v"] == b2["refine_sub_v"]: same += 1
        else:
            diff += 1
            if len(ex) < 15: ex.append(f"{i} ({a['v1']}/{b2['v1']})")
w(f"Both Chew: refined `|V1|` and `|V^S|` equal on {same} meshes, differ on {diff}; first differing (integer/signpost `|V1|`): {', '.join(ex)}\n")
idts = [(i, data['integer'][i], data['sp-chew'][i]) for i in kept if idt_ok(data['integer'].get(i)) and idt_ok(data['sp-chew'].get(i))]
bad = [(i, a['idt_sum_nplus'], b2['idt_sum_nplus']) for i, a, b2 in idts if a['idt_sum_nplus'] != b2['idt_sum_nplus'] or a['idt_flips'] != b2['idt_flips']]
w(f"IDT: flips and `Sum n+` equal on {len(idts) - len(bad)} of {len(idts)} meshes; differences: {bad[:15]}\n")
open(out_path, "w").write("\n".join(L) + "\n")
print("\n".join(L))
