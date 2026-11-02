# /// script
# dependencies = ["matplotlib"]
# ///
"""WP16a: tables, convergence fits and (if matplotlib is present) the Fig. 10 plots.

Reads the logs of `examples/consistency.rs` in `runs_tmp/wp16a/` (see NOTES-WP16A.md).
Run: uv run scripts/rs_consistency_plot.py [logdir]
"""
import glob
import math
import re
import sys

LOG = sys.argv[1] if len(sys.argv) > 1 else "runs_tmp/wp16a"
# the smooth reference values: the Duffy rule at the finest level (the cell rule is the check)
REF = {
    ("pert", 6, 12): "smooth_pert_6_12_trap:72_duffy:12:24:24.log",
    ("pert", 4, 8): "smooth_pert_4_8_trap:96_duffy:10:22:22.log",
    ("torus", 6, 12): "smooth_torus_6_12_traps:128_duffy:12:24:24.log",
    ("torus", 4, 8): "smooth_torus_4_8_traps:128_duffy:12:24:24.log",
}


def kv(line):
    return dict(p.split("=", 1) for p in line.split()[1:] if "=" in p)


def smooth(key):
    return float(kv(open(f"{LOG}/{REF[key]}").read().strip())["value"])


def fit(xs, ys):
    lx, ly = [math.log(x) for x in xs], [math.log(y) for y in ys]
    n = len(xs)
    mx, my = sum(lx) / n, sum(ly) / n
    return sum((a - mx) * (b - my) for a, b in zip(lx, ly)) / sum((a - mx) ** 2 for a in lx)


def discrete(surf):
    rows = []
    for f in glob.glob(f"{LOG}/disc_{surf}_*.log"):
        for line in open(f):
            if line.startswith("DISC"):
                rows.append(kv(line))
    return sorted(rows, key=lambda r: int(r["faces"]))


def main():
    series = {}
    for surf in ("pert", "torus"):
        rows = discrete(surf)
        for a, b in ((6, 12), (4, 8)):
            E = smooth((surf, a, b))
            print(f"\n{surf} (alpha, beta) = ({a}, {b}), smooth energy {E:.10f}")
            print("   faces        h         E_h          rel.err")
            hs, es = [], []
            for r in rows:
                eh = float(r[f"E_{a}_{b}"])
                h = float(r["h"])
                err = abs(eh - E) / E
                hs.append(h)
                es.append(err)
                print(f"{int(r['faces']):8d}  {h:.5f}  {eh:12.6f}  {err:.4e}  (signed {(eh - E) / E:+.3e})")
            print("  fitted rate, all points:", round(fit(hs, es), 3), " last four:", round(fit(hs[-4:], es[-4:]), 3))
            series[(surf, a, b)] = (hs, es)
    print("\nBarnes-Hut error vs theta (signed, relative to the exact discrete energy)")
    thetas = {}
    for f in sorted(glob.glob(f"{LOG}/theta_pert_*.log")) + sorted(glob.glob(f"{LOG}/disc_pert_*.log")):
        for line in open(f):
            if line.startswith("THETA"):
                d = kv(line)
                t = float(d["theta"])
                if t == 0:
                    continue
                signed = (float(d["bh"]) - float(d["exact"])) / float(d["exact"])
                thetas.setdefault((d["faces"], d["alpha"]), {})[t] = signed
    for (faces, alpha), dct in sorted(thetas.items(), key=lambda kv_: (kv_[0][1], int(kv_[0][0]))):
        ts = sorted(dct)
        print(f"faces {faces} alpha {alpha}: " + "  ".join(f"{t:g}:{dct[t]:+.2e}" for t in ts))
        hi = [t for t in ts if t >= 0.35]
        if len(hi) >= 3 and all(dct[t] > 0 for t in hi):
            print("   fitted order, theta >= 0.35:", round(fit(hi, [dct[t] for t in hi]), 2))
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        return
    fig, ax = plt.subplots(1, 2, figsize=(10, 4))
    for (surf, a, b), (hs, es) in series.items():
        ax[0].loglog(hs, es, "o-", label=f"{surf} ({a},{b})")
    hs = [0.03, 0.3]
    ax[0].loglog(hs, [0.05 * (h / 0.03) for h in hs], "k:", label="O(h)")
    ax[0].loglog(hs, [0.02 * (h / 0.03) ** 2 for h in hs], "k--", label="O(h^2)")
    ax[0].set_xlabel("h")
    ax[0].set_ylabel("relative error vs smooth")
    ax[0].legend()
    for (faces, alpha), dct in sorted(thetas.items()):
        if alpha == "6":
            ts = sorted(dct)
            ax[1].loglog(ts, [abs(dct[t]) for t in ts], "o-", label=f"{faces} faces")
    ts = [0.1, 1.0]
    ax[1].loglog(ts, [1e-4 * (t / 0.1) ** 2 for t in ts], "k--", label="O(theta^2)")
    ax[1].set_xlabel("theta")
    ax[1].set_ylabel("|relative error| vs exact discrete")
    ax[1].legend()
    fig.tight_layout()
    fig.savefig(f"{LOG}/consistency.png", dpi=130)
    print("wrote", f"{LOG}/consistency.png")


main()
