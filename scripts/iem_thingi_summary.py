"""WP9b: summarizes $OUT/jobs/<mode>/*.json (scripts/iem-thingi.sh) into $OUT/summary.md.

usage: iem_thingi_summary.py ROOT OUT
"""
import csv, json, os, sys, glob, collections

root, out = sys.argv[1], sys.argv[2]
meta = {r["id"]: r for r in csv.DictReader(open(os.path.join(root, "data/thingi10k/converted.tsv")), delimiter="\t")}
kept = [i for i, r in meta.items() if r["status"] == "kept"]
MODES = [("refine", "pre-refinement (integer, Chew, 25 deg)"), ("noref", "no pre-refinement"), ("refine-sp", "pre-refinement (signposts)")]

data = {}
for m, _ in MODES:
    d = {}
    for f in glob.glob(os.path.join(out, "jobs", m, "*.json")):
        try:
            r = json.load(open(f))
        except Exception:
            r = {"id": os.path.basename(f)[:-5], "status": "crash", "error": "unreadable result"}
        d[r["id"]] = r
    if d:
        data[m] = d


def pct(n, d):
    return f"{100.0 * n / d:.2f} %" if d else "-"


def q(xs, p):
    xs = sorted(xs)
    return xs[min(len(xs) - 1, int(round(p * (len(xs) - 1))))] if xs else float("nan")


def fmt(x):
    return "-" if x != x else f"{x:.3g}"


L = []
w = L.append
w("# WP9b summary: Thingi10k benchmark of intrinsic simplification\n")
w(f"Meshes: {len(kept)} kept by `thingi_convert`. Target: as far as possible (`AllRemovable`), default ICE metric. "
  "`ratio = |V_final| / |V_input|` (the paper normalizes both counts by the number of components, which cancels); "
  "`ratio_mean_component` is the mean over components of the component's final / input vertices. "
  "Timings are indicative only (24 parallel jobs on a shared machine).\n")

for m, name in MODES:
    if m not in data:
        continue
    d = data[m]
    rs = [d[i] for i in kept if i in d]
    ok = [r for r in rs if r.get("status") == "ok"]
    w(f"## {name} (`{m}`)\n")
    c = collections.Counter(d[i]["status"] if i in d else "missing" for i in kept)
    w("Outcomes: " + ", ".join(f"{k} {v}" for k, v in sorted(c.items())) + "\n")
    n = len(kept)
    # a run "reaches" a ratio if it finished ok and its ratio is below the bound
    for lim in (0.1, 0.01):
        k = sum(1 for r in ok if r.get("ratio") is not None and r["ratio"] <= lim)
        w(f"* ratio <= {lim:.0%}: {pct(k, n)} of all meshes ({k}/{n}), {pct(k, len(ok))} of the {len(ok)} completed")
    for lim in (0.1, 0.01):
        k = sum(1 for r in ok if r.get("ratio_mean_component") is not None and r["ratio_mean_component"] <= lim)
        w(f"* mean per-component ratio <= {lim:.0%}: {pct(k, n)} ({k}/{n})")
    w("  (paper, Fig. 17, without pre-refinement: 98 % and 84 %)\n")

    # failures
    w("### Failures\n")
    bad = [r for r in rs if r.get("status") != "ok"]
    g = collections.defaultdict(list)
    for r in bad:
        key = r["status"] + (": " + (r.get("error") or "")[:90] if r.get("error") else "")
        g[key].append(r["id"])
    w(f"Job failures (no result): {len(bad)}.\n")
    for k, ids in sorted(g.items(), key=lambda kv: -len(kv[1])):
        w(f"* {k}: {len(ids)} (e.g. {', '.join(ids[:5])})")
    w("")
    cf = sum(1 for r in ok if r.get("commit_failures_total", 0) > 0)
    ff = sum(1 for r in ok if r.get("flatten_failures", 0) > 0)
    cert = [r for r in ok if r.get("certify_ok") == 0]
    clamp = sum(1 for r in ok if r.get("n_clamped", 0) > 0)
    w(f"Completed meshes with commit failures: {cf} ({pct(cf, len(ok))}); with flattening failures: {ff} ({pct(ff, len(ok))}); "
      f"certificate failures: {len(cert)}; with point clamps: {clamp}.\n")
    reasons = collections.Counter()
    for r in ok:
        for k, v in (r.get("commit_failures") or {}).items():
            reasons[("commit", k)] += v
    skips = collections.Counter()
    skip_meshes = collections.Counter()
    for r in ok:
        for k, v in (r.get("skips") or {}).items():
            skips[k] += v
            skip_meshes[k] += 1
    if reasons:
        w("Commit failures by reason (events): " + ", ".join(f"{k[1]} {v}" for k, v in reasons.most_common()) + "\n")
    w("Vertices left unremovable (skips) by reason, summed / meshes: " +
      (", ".join(f"{k} {v} / {skip_meshes[k]}" for k, v in skips.most_common()) or "none") + "\n")
    for r in cert[:10]:
        w(f"* certify failed on {r['id']}: {r.get('certify_error')}")
    w("")

    # ratio distribution
    w("### Ratio distribution\n")
    rat = [r["ratio"] for r in ok if r.get("ratio") is not None]
    w("quantiles of `ratio`: " + ", ".join(f"p{int(p*100)} {fmt(q(rat, p))}" for p in (0.5, 0.75, 0.9, 0.95, 0.98, 0.99)) + f", max {fmt(max(rat) if rat else float('nan'))}\n")
    stops = collections.Counter(r.get("stop") for r in ok)
    w("stop reasons: " + ", ".join(f"{k} {v}" for k, v in stops.most_common()) + "\n")

    # time
    w("### Time\n")
    spr = [r["sec_per_removal"] for r in ok if r.get("sec_per_removal") is not None and r.get("removals", 0) >= 1000]
    w(f"seconds per removal (meshes with >= 1000 removals, simplify phase only, {len(spr)} meshes): median {fmt(q(spr, .5))}, p90 {fmt(q(spr, .9))}, p99 {fmt(q(spr, .99))}, max {fmt(max(spr) if spr else float('nan'))}")
    if spr:
        w(f"removals per second: median {fmt(1/q(spr, .5))}, p10 {fmt(1/q(spr, .9))}  (paper: about 10,000 per second)")
    tt = [r["t_total"] for r in ok]
    w(f"total seconds per mesh: median {fmt(q(tt, .5))}, p90 {fmt(q(tt, .9))}, p99 {fmt(q(tt, .99))}, max {fmt(max(tt) if tt else float('nan'))}; sum {fmt(sum(tt)/3600)} CPU hours")
    if m != "noref":
        tr = [r.get("t_refine", 0) + r.get("t_locate", 0) for r in ok]
        ts = [r.get("t_simplify", 0) for r in ok]
        w(f"share of time: pre-refinement {pct(sum(tr), sum(tt))}, simplify {pct(sum(ts), sum(tt))}")
        ad = [r.get("refine_added", 0) / max(1, r["v0"]) for r in ok if "v0" in r]
        w(f"vertices added by pre-refinement / |V|: median {fmt(q(ad, .5))}, p90 {fmt(q(ad, .9))}, max {fmt(max(ad) if ad else float('nan'))}")
    rss = [r["peak_rss_mb"] for r in rs if "peak_rss_mb" in r]
    w(f"peak RSS (MB): median {fmt(q(rss, .5))}, p99 {fmt(q(rss, .99))}, max {fmt(max(rss) if rss else float('nan'))}\n")

    # outliers
    w("### Outliers (ratio > 10 %) and genus\n")
    out_r = [r for r in ok if r.get("ratio") is not None and r["ratio"] > 0.1]
    def genus_note(r):
        c, chi, b = r.get("components", 1), r.get("chi", 0), r.get("boundary_edges", 0)
        # closed: chi = sum (2 - 2 g_k): the total genus is (2c - chi)/2
        return (2 * c - chi) / 2 if b == 0 else None
    gen = [genus_note(r) for r in out_r]
    hi = sum(1 for g in gen if g is not None and g >= 10)
    w(f"{len(out_r)} meshes above 10 % ({pct(len(out_r), n)}); of these {hi} closed with total genus >= 10, "
      f"{sum(1 for r in out_r if r.get('boundary_edges', 0) > 0)} with boundary.\n")
    allg = [(genus_note(r), r["ratio"]) for r in ok if genus_note(r) is not None]
    for lo, hi_ in ((0, 1), (1, 10), (10, 100), (100, 1e9)):
        sel = [rt for g, rt in allg if lo <= g < hi_]
        if sel:
            w(f"* total genus in [{lo}, {hi_:g}): {len(sel)} meshes, median ratio {fmt(q(sel, .5))}, share <= 10 %: {pct(sum(1 for x in sel if x <= .1), len(sel))}")
    w("")

    # worst table
    w("### Worst cases (by ratio)\n")
    w("| id | |V| | comps | chi | bdry edges | final | ratio | skips | commit fails | status |")
    w("|---|---|---|---|---|---|---|---|---|---|---|")
    worst = sorted(ok, key=lambda r: -(r.get("ratio") or 0))[:25]
    for r in worst:
        sk = ",".join(f"{k}:{v}" for k, v in (r.get("skips") or {}).items())
        w(f"| {r['id']} | {r.get('v0')} | {r.get('components')} | {r.get('chi')} | {r.get('boundary_edges')} | {r.get('v_final')} | {fmt(r.get('ratio'))} | {sk} | {r.get('commit_failures_total')} | {r['status']} |")
    w("")
    w("Slowest meshes:\n")
    w("| id | |V| | t_total (s) | removals | RSS MB |")
    w("|---|---|---|---|---|")
    for r in sorted(ok, key=lambda r: -r["t_total"])[:10]:
        w(f"| {r['id']} | {r.get('v0')} | {fmt(r['t_total'])} | {r.get('removals')} | {r.get('peak_rss_mb')} |")
    w("")

if "refine" in data and "noref" in data:
    w("## With against without pre-refinement\n")
    both = [i for i in kept if data["refine"].get(i, {}).get("status") == "ok" and data["noref"].get(i, {}).get("status") == "ok"]
    better = sum(1 for i in both if data["refine"][i]["ratio"] < data["noref"][i]["ratio"] - 1e-12)
    worse = sum(1 for i in both if data["refine"][i]["ratio"] > data["noref"][i]["ratio"] + 1e-12)
    w(f"{len(both)} meshes ok in both: refinement gives a smaller final ratio on {better}, a larger on {worse}, equal on {len(both)-better-worse}.\n")

open(os.path.join(out, "summary.md"), "w").write("\n".join(L) + "\n")
print("\n".join(L))
