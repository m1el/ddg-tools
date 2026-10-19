# /// script
# requires-python = ">=3.9"
# ///
"""Aggregates the logs of cut_study.sh into CSV + markdown tables.

uv run crates/ddg-developability/examples/cut_study.py runs_tmp/dev-cut-study
"""
import csv, glob, os, re, statistics as st, sys

root = sys.argv[1]
KV = re.compile(r"(\w+)=(\S+)")


def result(path):
    try:
        for line in reversed(open(path).read().splitlines()):
            if line.startswith("RESULT "):
                return dict(KV.findall(line))
    except OSError:
        pass
    return None


rows = []
for p in sorted(glob.glob(f"{root}/cut/*.log")):
    r = result(p)
    base = os.path.basename(p)[:-4]
    m = re.match(r"(.+)-(covariance|intrinsic|combinatorial)-(\d+)_(\w+)_([\de.-]+)_(off|[\d.]+)$", base)
    if not m or not r or r.get("cut") != "ok":
        print("skipped", base, None if not r else r.get("cut"), file=sys.stderr)
        continue
    mesh, fe, budget, ce, eps, tol = m.groups()
    nf = int(r["faces_cut"])
    rows.append(
        dict(
            mesh=mesh, flow=fe, budget=int(budget), cut_energy=ce, eps=eps, tol=tol, faces=nf,
            flips=int(r["flipped_uv"]), flip_pct=100.0 * int(r["flipped_uv"]) / nf, mean_u=float(r["mean_abs_u"]),
            terminals=int(r["terminals"]) + int(r["added_terminals"]), cut_len=float(r["cut_len_frac"]),
            charts=int(r["charts"]), om_non=float(r["om_nonterminals"]), nonemb=int(r["nonembedded"]),
            nonemb_non=int(r["nonembedded_nonterminal"]),
        )
    )
with open(f"{root}/all.csv", "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0]))
    w.writeheader()
    w.writerows(rows)

CONTROL = ("spot_control_mesh", "bob_controlmesh", "blub_control_mesh")


def cls(r):
    if r["mesh"].startswith("Crumpled"):
        return "crumpled"
    return "control" if r["mesh"] in CONTROL else "dense"


def variant(r):
    return (r["cut_energy"], r["tol"])


out = []
# table 1: per (class, flow energy, budget) x (cut energy, tol) at each eps: worst flips %, max mean|u|, mean cut length
for c in ("control", "dense", "crumpled"):
    for fe in ("covariance", "intrinsic", "combinatorial"):
        for budget in sorted({r["budget"] for r in rows if cls(r) == c and r["flow"] == fe}):
            sel = [r for r in rows if cls(r) == c and r["flow"] == fe and r["budget"] == budget]
            if not sel:
                continue
            out.append(f"\n#### {c} meshes, flow {fe}, {budget} it/round ({len({r['mesh'] for r in sel})} meshes)\n")
            out.append("| terminal energy | defect tol | eps | worst flips % | mean flips % | worst mean abs u | median mean abs u | mean cut length % | mean terminals |")
            out.append("|---|---|---|---|---|---|---|---|---|")
            keys = sorted({(r["cut_energy"], r["tol"], r["eps"]) for r in sel}, key=lambda k: (k[0], k[1], float(k[2])))
            for k in keys:
                g = [r for r in sel if (r["cut_energy"], r["tol"], r["eps"]) == k]
                out.append(
                    f"| {k[0]} | {k[1]} | {k[2]} | {max(r['flip_pct'] for r in g):.2f} | {st.mean(r['flip_pct'] for r in g):.2f} | "
                    f"{max(r['mean_u'] for r in g):.2e} | {st.median(r['mean_u'] for r in g):.2e} | "
                    f"{100*st.mean(r['cut_len'] for r in g):.2f} | {st.mean(r['terminals'] for r in g):.0f} |"
                )
open(f"{root}/summary.md", "w").write("\n".join(out) + "\n")
print("\n".join(out))
