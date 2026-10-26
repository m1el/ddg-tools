#!/usr/bin/env python3
"""Stills and animations of the Repulsive Curves showcase from the runs of `run_all.sh`
(showcase/work/rc/<run>/frameNNNNN.obj + report.json), rendered with render.py (Blender Cycles),
labelled with the numbers of the reports and written to showcase/media/repulsive_*.

  uv run showcase/repulsive-curves/make_media.py [stills|anims|all] [name ...]

Renders run in parallel (RC_JOBS, default 6 Blender processes; they share one GPU).
"""
import json
import math
import os
import subprocess
import textwrap
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

HERE = Path(__file__).resolve().parent
SHOW = HERE.parent
RUNS = SHOW / "work" / "rc"
MEDIA = SHOW / "media"
FRAMES = SHOW / "work" / "frames"
sys.path.insert(0, str(SHOW))
import post  # noqa: E402

JOBS = int(os.environ.get("RC_JOBS", 6))
# control points to mark (vertex indices of the start curve), see render.py
PINS = {"spline_points": "62,94,126,158,222,254,286,318", "spline_tangents": "0,50,99"}


def report(run):
    return json.load(open(RUNS / run / "report.json"))


def frames(run):
    return sorted((RUNS / run).glob("frame*.obj"))


def curve_length(path):
    V, L = [], 0.0
    E = []
    for line in open(path):
        p = line.split()
        if p and p[0] == "v":
            V.append([float(x) for x in p[1:4]])
        elif p and p[0] == "l":
            E.append((int(p[1]) - 1, int(p[2]) - 1))
    return sum(math.dist(V[a], V[b]) for a, b in E), len(V)


def blender(args, log, pins=""):
    cmd = ["blender", "-b", "--python", str(HERE / "render.py"), "--"] + [str(a) for a in args]
    for attempt in range(8):  # several Blenders share one GPU: retry on out-of-memory
        with open(log, "w") as f:
            r = subprocess.run(cmd, stdout=f, stderr=subprocess.STDOUT, env={**os.environ, "RC_PINS": pins})
        if r.returncode == 0 and "Traceback" not in open(log).read():
            return
    raise RuntimeError(f"blender failed: {log}")


# name -> (run, kind, frame picks as fractions, caption function)
def stills():
    def pk(run, fr):
        f = frames(run)
        return [f[min(len(f) - 1, int(round(x * (len(f) - 1))))] for x in fr]

    def cap_grow(r):
        return (f"{r['steps']} steps, {r['vertices'][0]} to {r['vertices'][1]} vertices, "
                f"length x{r['length'][1] / r['length'][0]:.1f}, "
                f"min distance between edges {r['audit_end']['min_edge_distance']:.3g}, {r['seconds']:.0f} s")

    def cap_graph(r):
        return (f"{r['steps']} steps, {r['vertices'][1]} vertices, energy {r['energy'][0]:.3g} to "
                f"{r['energy'][1]:.3g}, min distance {r['audit_start']['min_edge_distance']:.3g} to "
                f"{r['audit_end']['min_edge_distance']:.3g}, {r['seconds']:.0f} s")

    return {
        "packing_domain": ("bunny", "bunny", [0, 0.5, 1], cap_grow, (640, 640)),
        "packing_cup": ("cup", "cup", [0, 0.5, 1], cap_grow, (640, 640)),
        "packing_surface": ("torus-multi", "torus", [0, 0.5, 1], cap_grow, (640, 640)),
        "packing_sphere": ("sphere", "sphere", [0, 0.5, 1], cap_grow, (640, 640)),
        "packing_boundary": ("bluenoise", "bluenoise", [0, 0.5, 1], cap_grow, (640, 640)),
        "graph_drawing": ("planar-graph", "graph-top", [0, 0.25, 1], cap_graph, (640, 640)),
        "graph_isometric_k5": ("k5", "graph", [0, 0.5, 1], cap_graph, (640, 640)),
        "graph_isometric_k33": ("k33", "graph", [0, 0.5, 1], cap_graph, (640, 640)),
        "vascular": ("capillary-low-bh", "capillary", [0, 1], cap_graph, (960, 1100)),
        "spline_points": ("interp", "spline", [0, 0.5, 1], cap_grow, (640, 640)),
        "spline_tangents": ("b-spline-tangents", "spline", [0, 1], cap_grow, (640, 640)),
    }, pk


def make_still(name, spec, pk):
    run, kind, fr, cap, (W, H) = spec
    if not (RUNS / run / "report.json").exists():
        print("skip", name, "(no report)")
        return
    r = report(run)
    files = pk(run, fr)
    tmp = FRAMES / f"rc_{name}"
    tmp.mkdir(parents=True, exist_ok=True)
    panels = []
    for i, f in enumerate(files):
        out = tmp / f"p{i}.png"
        blender(["still", kind, f, out, W, H], HERE.parent.parent / "runs_tmp" / f"rc-wp13/render_{name}_{i}.log", PINS.get(name, ""))
        comp = tmp / f"c{i}.png"
        post.composite(out, comp)
        panels.append(comp)
    row = tmp / "row.png"
    subprocess.run(["convert", *map(str, panels), "+append", "-background", post.BG, "-gravity", "south",
                    "-splice", "0x110", str(row)], check=True)
    post.composite(row, MEDIA / f"repulsive_{name}.png", "\n".join(textwrap.wrap(cap(r), 30 * len(panels) * W // 640)), pointsize=34)
    print("wrote", MEDIA / f"repulsive_{name}.png")


def make_anim(name, run, kind, W, H, n, title):
    d = RUNS / run
    r = report(run)
    log_by_step = {l[0]: l for l in r["log"]}
    blender(["anim", kind, d, f"rc_{name}", W, H, n], HERE.parent.parent / f"runs_tmp/rc-wp13/render_anim_{name}.log")
    fdir = FRAMES / f"rc_{name}"
    names = json.load(open(fdir / "frames.json"))
    files = frames(run)
    # frame k of the run is step k * every, the last frame is the last step
    ev = int(round((r["steps"]) / max(1, len(files) - 1))) if len(files) > 1 else 1
    labels = []
    for nm in names:
        k = int(nm[5:10])
        step = r["steps"] if k == len(files) - 1 else k * ev
        L, nv = curve_length(d / nm)
        e = log_by_step.get(step, [0, 0, r["energy"][0]])[2] if step else r["energy"][0]
        labels.append(f"{title}   step {step}   vertices {nv}   length {L:.2f}   energy {e:.4g}")
    json.dump(labels, open(fdir / "labels.json", "w"))
    post.video(fdir, MEDIA / f"repulsive_{name}", fdir / "labels.json", fps=24, gif_width=640)


def main():
    a = sys.argv[1:]
    mode = a[0] if a else "all"
    sel = set(a[1:])
    MEDIA.mkdir(exist_ok=True)
    jobs = []
    with ThreadPoolExecutor(JOBS) as ex:
        if mode in ("stills", "all"):
            specs, pk = stills()
            for name, spec in specs.items():
                if not sel or name in sel:
                    jobs.append(ex.submit(make_still, name, spec, pk))
        if mode in ("anims", "all"):
            if not sel or "unknot" in sel:
                jobs.append(ex.submit(make_anim, "unknot", "freedman600", "unknot", 960, 720, 120, "Freedman-He-Wang unknot"))
            if not sel or "packing" in sel:
                jobs.append(ex.submit(make_anim, "packing", "bunny", "bunny", 960, 720, 90, "Packing in the bunny"))
        for j in jobs:
            j.result()


if __name__ == "__main__":
    main()
