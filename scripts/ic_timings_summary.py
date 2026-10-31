"""WP10: tables of runs_tmp/wp10/timings/ (scripts/ic-timings.sh) -> stdout and runs_tmp/wp10/timings.md."""
import json, os, sys, statistics
root = sys.argv[1]; T = os.path.join(root, "runs_tmp/wp10/timings")
def ours(i, b):
    p = os.path.join(T, f"{i}_ours_{b}.json")
    try: return json.load(open(p))
    except Exception: return None
def demo(i, b, op):
    p = os.path.join(T, f"{i}_demo_{b}_{op}_stats.tsv")
    try:
        l = open(p).read().splitlines(); h = [x.strip() for x in l[0].split("\t")]; return dict(zip(h, l[1].split("\t")))
    except Exception: return None
def f(x, d=3):
    try: return f"{float(x):.{d}f}"
    except Exception: return "-"
L = []; w = L.append
def row(i):
    oi, os_ = ours(i, "integer"), ours(i, "sp-ruppert")
    di_f, di_r, ds_f, ds_r = demo(i, "integer", "flip"), demo(i, "integer", "refine"), demo(i, "signpost", "flip"), demo(i, "signpost", "refine")
    g = lambda d, k: d.get(k) if d else None
    return (i, oi, os_, di_f, di_r, ds_f, ds_r)
hdr = "| mesh | `|V0|` | ours int IDT | ours int refine | ours int subdiv | ours sp IDT | ours sp refine | ours sp subdiv | demo int IDT | demo int refine | demo int subdiv | demo sp IDT | demo sp refine |"
sep = "|---|---|---|---|---|---|---|---|---|---|---|---|---|"
def line(r):
    i, oi, os_, di_f, di_r, ds_f, ds_r = r
    v = lambda d, k, dg=3: f(d[k], dg) if d and k in d else "-"
    return f"| {i} | {oi['v0'] if oi else '-'} | {v(oi,'idt_seconds')} | {v(oi,'refine_seconds')} | {v(oi,'refine_sub_seconds')} | {v(os_,'idt_seconds')} | {v(os_,'refine_seconds')} | {v(os_,'refine_sub_seconds')} | {v(di_f,'flippingDuration')} | {v(di_r,'refinementDuration')} | {v(di_r,'commonSubdivisionTracingDuration')} | {v(ds_f,'flippingDuration')} | {v(ds_r,'refinementDuration')} |"
paper = [l.strip() for l in open(os.path.join(T, "paper.txt")) if l.strip()]
w("Seconds, one pinned core, run alone. 'refine' of ours is `refine` including the certificate (integers); 'subdiv' is the construction of the common subdivision after refinement (ours: `common_subdivision` / traces + overlay; demo: `commonSubdivisionTracingDuration`, demo integer only: its signpost run stops at the assertion or is not timed for the subdivision).\n")
w("### The paper's meshes (44395 Fig. 15, 54674 Fig. 16, 719791 Fig. 19)\n")
w("Paper: 0.2 s (Fig. 15), 0.6 s (Fig. 16 top), 719791: 4 min with integers, 1.5 min with signposts (i9-9980XE).\n")
w(hdr); w(sep)
for i in paper: w(line(row(i)))
ids = [l.strip() for l in open(os.path.join(T, "subset50.txt")) if l.strip()]
rows = [row(i) for i in ids]
w("\n### The 50-mesh subset (every 155th kept mesh by id)\n")
w("Ids: " + ", ".join(ids) + "\n")
w(hdr); w(sep)
for r in rows: w(line(r))
def tot(sel):
    s = 0.0
    for r in rows:
        x = sel(r)
        if x is not None: s += float(x)
    return s
S = lambda k, n: tot(lambda r: r[n][k] if r[n] and k in r[n] else None)
w("\n**Totals over the subset (s):** " + "; ".join([
    f"ours integers IDT {S('idt_seconds',1):.2f}, refine {S('refine_seconds',1):.2f}, subdivision {S('refine_sub_seconds',1):.2f}",
    f"ours signposts IDT {S('idt_seconds',2):.2f}, refine {S('refine_seconds',2):.2f}, subdivision {S('refine_sub_seconds',2):.2f}",
    f"demo integers IDT {S('flippingDuration',3):.2f}, refine {S('refinementDuration',4):.2f}, subdivision {S('commonSubdivisionTracingDuration',4):.2f}",
    f"demo signposts IDT {S('flippingDuration',5):.2f}, refine {S('refinementDuration',6):.2f}"]) + "\n")
ok = [r for r in rows if r[1] and r[4] and "refine_seconds" in r[1] and "refinementDuration" in r[4]]
ratios = [float(r[4]["refinementDuration"]) / r[1]["refine_seconds"] for r in ok if r[1]["refine_seconds"] > 0.01]
if ratios: w(f"Refinement time, demo / ours, on the {len(ratios)} meshes where ours takes over 10 ms: median {statistics.median(ratios):.2f}, min {min(ratios):.2f}, max {max(ratios):.2f}.\n")
try:
    lo = [float(l.split()[1].rstrip(",")) for l in open(os.path.join(T, "uptime.log"))]
    w(f"Load average (1 min), sampled every 15 s during the runs: median {statistics.median(lo):.2f}, 90th percentile {sorted(lo)[int(0.9*len(lo))]:.2f}, max {max(lo):.2f} ({len(lo)} samples); one job at a time on core 7, the other agent's work on the machine is the rest. `load.log` has the load before and after every mesh.\n")
except Exception as e: w(f"(no uptime.log: {e})")
out = "\n".join(L) + "\n"; open(os.path.join(root, "runs_tmp/wp10/timings.md"), "w").write(out); print(out)
