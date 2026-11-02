#!/usr/bin/env python3
"""Stills and animations of the Repulsive Surfaces showcase, rendered with render.py (Blender
Cycles) from the OBJ frames of `run_all.sh` (WP15), `examples/ishii` and `examples/fidelity`
(handcuffs unlinking). Written to showcase/media/repulsive_surfaces_* (in the main checkout,
untracked).

  uv run showcase/repulsive-surfaces/make_media.py [genus ishii handcuffs pole walnut shrinkwrap
                                                   envelopes cobblestones h klein]   (default: all)

Inputs (absolute, all in the main checkout):
  MAIN/showcase/media/repulsive-surfaces/<run>/{frame_NNNNN,final}.obj, steps.csv   (WP15)
  MAIN/runs_tmp/ishii/<name>_final.obj, <name>.log                                  (Fig. 16)
  MAIN/runs_tmp/handcuffs-unlink/{start.obj,frames/,run.csv}                        (unlinking run)
  <worktree>/data/repulsive-surfaces/scenes (scripts/rs-fetch-scenes.sh)            (obstacles)
"""
import csv
import json
import os
import re
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

HERE = Path(__file__).resolve().parent
SHOW = HERE.parent
ROOT = SHOW.parent
MAIN = Path(os.environ.get("DDG_MAIN", "/var/data/discrete-geometry/ddg-tools"))
RUNS = MAIN / "showcase" / "media" / "repulsive-surfaces"
MEDIA = MAIN / "showcase" / "media"
ISHII = MAIN / "runs_tmp" / "ishii"
UNLINK = MAIN / "runs_tmp" / "handcuffs-unlink"
LOGS = MAIN / "runs_tmp" / "rs-showcase"
SCENES = ROOT / "data" / "repulsive-surfaces" / "scenes"
WORK = SHOW / "work" / "rs"
FRAMES = SHOW / "work" / "frames"
sys.path.insert(0, str(SHOW))
import post  # noqa: E402

JOBS = int(os.environ.get("RS_JOBS", 6))
BG = post.BG


def blender(args, log):
    cmd = ["blender", "-b", "--python", str(HERE / "render.py"), "--"] + [str(a) for a in args]
    LOGS.mkdir(parents=True, exist_ok=True)
    for attempt in range(6):
        with open(LOGS / log, "w") as f:
            r = subprocess.run(cmd, stdout=f, stderr=subprocess.STDOUT)
        if r.returncode == 0 and "Traceback" not in open(LOGS / log).read():
            return
    raise RuntimeError(f"blender failed: {LOGS / log}")


def spec_file(name, spec):
    WORK.mkdir(parents=True, exist_ok=True)
    p = WORK / f"{name}.json"
    json.dump(spec, open(p, "w"))
    return p


def frames(run):
    return sorted((RUNS / run).glob("frame_*.obj"))


def steps(run):
    """{iteration: (energy_before, energy_after)} of a WP15 steps.csv."""
    out = {}
    p = RUNS / run / "steps.csv"
    for row in csv.DictReader(open(p)):
        out[int(row["iteration"])] = (float(row["energy_before"]), float(row["energy_after"]))
    return out


def step_of(path):
    return int(re.search(r"(\d+)\.obj$", str(path)).group(1))


def energy_at(st, k):
    if k == 0:
        return st[min(st)][0]
    return st[k][1] if k in st else st[max(st)][1]


def nverts(path):
    return sum(1 for l in open(path) if l.startswith("v "))


def sub(name):
    """'6_16' -> '6' + subscript digits."""
    a, b = name.split("_")
    return a + "".join("₀₁₂₃₄₅₆₇₈₉"[int(c)] for c in b)


def panel(layers, az=-35, el=25, margin=1.1):
    return {"layers": layers, "az": az, "el": el, "margin": margin}


def sheet(name, panels, labels, cols, W, tile, caption=None, common_scale=False, ps=22, margin=1.1):
    """Render the panels, label them, tile them into MEDIA/repulsive_surfaces_<name>.png."""
    d = FRAMES / f"rs_{name}"
    d.mkdir(parents=True, exist_ok=True)
    sp = spec_file(f"sheet_{name}", {"panels": panels, "common_scale": common_scale, "margin": margin})
    blender(["sheet", sp, d, W, W], f"sheet_{name}.log")
    comp = []
    for i, lab in enumerate(labels):
        c = d / f"c{i:02d}.png"
        post.composite(d / f"p{i:02d}.png", c, lab, pointsize=ps)
        comp.append(str(c))
    out = MEDIA / f"repulsive_surfaces_{name}.png"
    rows = (len(comp) + cols - 1) // cols
    subprocess.run(["montage", *comp, "-tile", f"{cols}x{rows}", "-geometry", f"{tile}x{tile}+4+4",
                    "-background", BG, str(out)], check=True)
    print("wrote", out)


def anim(name, panels, labels, fps, W=960, H=720, az=-35, el=25, margin=1.1, gif_width=640):
    sp = spec_file(f"anim_{name}", {"panels": panels, "az": az, "el": el, "margin": margin})
    fdir = FRAMES / f"rs_{name}"
    if fdir.exists():
        for f in fdir.glob("*.png"):
            f.unlink()
    blender(["anim", sp, f"rs_{name}", W, H], f"anim_{name}.log")
    json.dump(labels, open(fdir / "labels.json", "w"))
    post.video(fdir, MEDIA / f"repulsive_surfaces_{name}", fdir / "labels.json", fps=fps, gif_width=gif_width)


# ---------------------------------------------------------------------------------------- genus
# Exact energies and status after 800 steps (NOTES-FIDELITY, with the fixed eigen-solver).
GENUS_E = [2.234, 31.58, 38.59, 52.21, 55.92, 54.31, 55.87, 65.60, 71.91, 72.44, 93.72, 90.76]
GENUS_NOTE = {
    0: "round", 1: "torus, 1.6 %", 2: "same shape at conv.", 3: "not converged",
    4: "not checked", 5: "not checked", 6: "sym. to 1.9 %", 7: "not converged", 8: "not converged",
    9: "not converged", 10: "not converged", 11: "not converged",
}


def do_genus():
    panels = [panel([{"obj": str(RUNS / f"genus{g:02d}" / "final.obj")}], az=-35, el=30) for g in range(12)]
    labels = [f"g = {g}   E {GENUS_E[g]:.4g}\n{GENUS_NOTE[g]}" for g in range(12)]
    sheet("genus", panels, labels, 4, 800, 400, ps=26)


# ---------------------------------------------------------------------------------------- ishii
def ishii_results():
    res = {}
    for l in sorted(ISHII.glob("*.log")):
        for line in open(l):
            if line.startswith("RESULT"):
                kv = dict(t.split("=", 1) for t in line.split()[1:] if "=" in t)
                res[kv["name"]] = kv
    return res


def do_ishii():
    names = ["4_1", "5_1", "5_2", "5_3", "5_4"] + [f"6_{i}" for i in range(1, 17)]
    res = ishii_results()
    panels, labels = [], []
    for n in names:
        f = ISHII / n / f"{n}_final.obj"
        if not f.exists():
            print("ishii: missing", n)
            continue
        panels.append(panel([{"obj": str(f)}], az=-25, el=55, margin=1.05))
        st = res.get(n, {})
        mark = "" if st.get("status") == "Converged" else " *"
        labels.append(sub(n) + mark)
    sheet("ishii", panels, labels, 7, 480, 240, ps=30, margin=1.05)


# ------------------------------------------------------------------------------------- handcuffs
def do_handcuffs():
    fs = [UNLINK / "start.obj"] + sorted((UNLINK / "frames").glob("frame_*.obj"))
    rows = list(csv.DictReader(open(UNLINK / "run.csv")))
    # run.csv of the fidelity tool: iteration, energy_after
    en = {int(r["iteration"]): float(r["energy_after"]) for r in rows}
    e0 = float(rows[0]["energy_before"])
    pans, labs = [], []
    for f in fs:
        k = 0 if f.name == "start.obj" else step_of(f)
        pans.append(panel([{"obj": str(f)}], az=-30, el=40))
        e = e0 if k == 0 else en.get(k, en[max(en)])
        labs.append(f"Linked handcuffs, no plane attractor   step {k}   energy {e:.4g}")
    anim("handcuffs", pans, labs, fps=15, el=40, az=-30)


# ------------------------------------------------------------------------------------------ pole
def do_pole():
    cyl = {"cylinder": {"r": 0.075, "z0": -1.0, "z1": 1.0}, "color": [0.12, 0.34, 0.80], "alpha": 1.0,
           "roughness": 0.4}
    pans, labs = [], []
    for stage, run, tag in ((1, "pole_stage1", "stage 1, cylinder weight 0.1"),
                            (2, "pole_stage2", "stage 2, cylinder weight 0.000125")):
        st = steps(run)
        fs = frames(run)
        for f in fs[::2] + [fs[-1]]:
            k = step_of(f)
            pans.append(panel([{"obj": str(f)}, cyl], az=-40, el=22))
            labs.append(f"Pole scene, {tag}   step {k}   energy {energy_at(st, k):.4g}")
    anim("pole", pans, labs, fps=18, az=-40, el=22)


# ---------------------------------------------------------------------------------------- walnut
def do_walnut():
    shell = {"obj": str(SCENES / "Walnut" / "shell2.obj"), "alpha": 0.14, "color": [0.55, 0.65, 0.85],
             "roughness": 0.2}
    st = steps("walnut")
    fs = frames("walnut")
    pans, labs = [], []
    for f in fs + [fs[-1]] * 8:
        k = step_of(f)
        pans.append(panel([{"obj": str(f)}, shell], az=-35, el=20))
        labs.append(f"Walnut growing in its shell, area x8   step {k}   vertices {nverts(f)}   energy {energy_at(st, k):.4g}")
    anim("walnut", pans, labs, fps=8, az=-35, el=20)


# ----------------------------------------------------------------------------------- the stills
def still(name, p, W=1600, H=1000, label=None):
    d = FRAMES / f"rs_{name}"
    d.mkdir(parents=True, exist_ok=True)
    sp = spec_file(f"still_{name}", p)
    raw = d / "raw.png"
    blender(["still", sp, raw, W, H], f"still_{name}.log")
    post.composite(raw, MEDIA / f"repulsive_surfaces_{name}.png", label)
    print("wrote", MEDIA / f"repulsive_surfaces_{name}.png")


def do_shrinkwrap():
    cloud = str(SCENES / "Reconstruction" / "bunny-scaled.obj")
    n = nverts(cloud)
    idx = list(range(0, n, 3))
    fs = frames("shrinkwrap_bunny_points")
    picks = [fs[0], fs[len(fs) // 3], fs[-1]]
    pts = {"points": cloud, "index": idx, "radius": 0.007, "color": [0.12, 0.34, 0.80]}
    # the surface is translucent: at 220 steps it still encloses the whole cloud
    panels = [panel([{"obj": str(f), "alpha": 0.3, "color": [0.95, 0.60, 0.08], "roughness": 0.3}, pts], az=-35, el=20)
              for f in picks]
    labs = [f"step {step_of(f)}" for f in picks]
    sheet("shrinkwrap", panels, labs, 3, 800, 520, ps=26)


def do_envelopes():
    d = RUNS
    spot = {"obj": str(SCENES / "Reconstruction" / "spot_triangulated.obj"), "color": [0.95, 0.60, 0.08]}
    lv = [([0.85, 0.35, 0.55], 0.5), ([0.16, 0.58, 0.32], 0.35), ([0.10, 0.62, 0.66], 0.25), ([0.12, 0.34, 0.80], 0.16)]
    layers = [spot] + [{"obj": str(d / f"envelope_level{k}" / "final.obj"), "alpha": a, "color": c, "roughness": 0.2}
                       for k, (c, a) in zip(range(1, 5), lv)]
    still("envelopes", panel(layers, az=-35, el=18), 1600, 1100)


def do_cobblestones():
    layers = [{"obj": str(MAIN / "runs_tmp" / "rs-cobble" / "cobblestones" / "final.obj"), "color": "components"},
              {"box": {"center": [0, 0, 0], "half": [3.2, 3.2, 0.9]}, "alpha": 0.12, "color": [0.55, 0.65, 0.85],
               "roughness": 0.2}]
    still("cobblestones", panel(layers, az=-30, el=32), 1600, 1000)


def do_h():
    pins = {"points": str(RUNS / "h_sparse_volume" / "final.obj"), "index": list(range(24)), "radius": 0.012,
            "color": [0.07, 0.08, 0.10]}
    def pins_of(run):
        return dict(pins, points=str(RUNS / run / "final.obj"))
    box = {"box": {"center": [1.5, 2.0, 0.5], "half": [1.9, 2.4, 0.9]}, "alpha": 0.1, "color": [0.55, 0.65, 0.85],
           "roughness": 0.2}
    def surf(run):
        return {"obj": str(RUNS / run / "final.obj"), "wire": True, "wire_thickness": 0.0012}
    panels = [
        panel([{"obj": str(RUNS / "h_sparse_volume" / "frame_00000.obj"), "wire": True, "wire_thickness": 0.0012}, pins_of("h_sparse_volume")], az=-40, el=35),
        panel([surf("h_sparse_volume"), pins_of("h_sparse_volume")], az=-40, el=35),
        panel([surf("h_sparse_area"), pins_of("h_sparse_area")], az=-40, el=35),
        panel([surf("h_sparse_volume_box"), pins_of("h_sparse_volume_box"), box], az=-40, el=35),
    ]
    labs = ["start", "volume x2", "area x1.5", "volume x2, in a box"]
    sheet("h_pins", panels, labs, 4, 800, 400, ps=26, common_scale=True)


def do_klein():
    panels = [panel([{"obj": str(RUNS / "klein_subcritical" / "frame_00000.obj")}], az=-35, el=25),
              panel([{"obj": str(RUNS / "klein_subcritical" / "final.obj")}], az=-35, el=25)]
    sheet("klein", panels, ["start (70 intersecting pairs)", "step 700, (3.5, 7), 49 pairs"], 2, 900, 700, ps=28)


ALL = {"genus": do_genus, "ishii": do_ishii, "handcuffs": do_handcuffs, "pole": do_pole, "walnut": do_walnut,
       "shrinkwrap": do_shrinkwrap, "envelopes": do_envelopes, "cobblestones": do_cobblestones, "h": do_h,
       "klein": do_klein}


def main():
    sel = sys.argv[1:] or list(ALL)
    MEDIA.mkdir(parents=True, exist_ok=True)
    with ThreadPoolExecutor(JOBS) as ex:
        futs = {n: ex.submit(ALL[n]) for n in sel}
        for n, f in futs.items():
            try:
                f.result()
            except Exception as e:  # keep going: report at the end
                print("FAILED", n, e)


if __name__ == "__main__":
    main()
