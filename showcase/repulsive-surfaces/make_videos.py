#!/usr/bin/env python3
"""Videos of the Repulsive Surfaces results, rendered with render.py (Blender Cycles) from the OBJ
frames of the runs. Same studio, camera conventions and label style as make_media.py (whose helpers
it reuses). Written to showcase/media/repulsive_surfaces_<name>.mp4/.gif.

  uv run showcase/repulsive-surfaces/make_videos.py [dodeca dodeca_perturbed genus genus_plane envelopes
                                                    shrinkwrap h_pins slideoff cobblestones klein
                                                    spindles eversion]   (default: all)

Multi-panel videos: every panel is its own Blender animation (own framing; `track` zooms to the shape of
each frame), the labelled panels are tiled with ImageMagick. Panels share one time axis (the step
number); a run that has ended holds its last frame. Identical consecutive frames are rendered once.
Env: RS_JOBS (parallel blender processes, default 4), RS_SAMPLES (Cycles samples, default 64).

Inputs (all in the main checkout): showcase/media/repulsive-surfaces/<run>/ (WP15 frames),
runs_tmp/genus2/g<G>/ (converged genus series), runs_tmp/dodeca2/{ours,ref}_*/ (dodecahedron study),
the reference's scenes (RepulsiveSurfaces/repulsive-surfaces/scenes, MIT) for obstacles.
"""
import csv
import json
import math
import os
import re
import shutil
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import make_media as M  # noqa: E402
import post  # noqa: E402

os.environ.setdefault("RS_SAMPLES", "64")
RUNS, MEDIA, MAIN, FRAMES = M.RUNS, M.MEDIA, M.MAIN, M.FRAMES
DODECA = MAIN / "runs_tmp" / "dodeca2"
GENUS = M.GENUS
SCENES = Path(os.environ.get("RS_SCENES", "/var/data/discrete-geometry/RepulsiveSurfaces/repulsive-surfaces/scenes"))
if not SCENES.exists():
    SCENES = M.SCENES
JOBS = int(os.environ.get("RS_JOBS", 4))
M.LOGS = MAIN / "runs_tmp" / "rs-videos"
GOLD = [0.95, 0.60, 0.08]


# ------------------------------------------------------------------------------------------ series
def csv_energy(d):
    """{iteration: energy_after} of a steps.csv (NaN rows skipped)."""
    out = {}
    p = Path(d) / "steps.csv"
    if p.exists():
        for r in csv.DictReader(open(p)):
            try:
                e = float(r["energy_after"])
            except ValueError:
                continue
            if math.isfinite(e):
                out[int(r["iteration"])] = e
    return out


def series(d, pattern="frame_*.obj", final=True):
    """[(step, path)] of a run: the frames, and final.obj as the last step of steps.csv."""
    d = Path(d)
    fs = sorted((step_num(f), f) for f in d.glob(pattern))
    if final and (d / "final.obj").exists():
        en = csv_energy(d)
        last = max(en) if en else fs[-1][0]
        if last > fs[-1][0]:
            fs.append((last, d / "final.obj"))
    return fs


def step_num(path):
    return int(re.search(r"(\d+)\.obj$", str(path)).group(1))


def at(ser, t):
    """The last frame with step <= t (the run holds its last frame once it has ended)."""
    best = ser[0]
    for s in ser:
        if s[0] <= t:
            best = s
    return best


def energy_near(en, k):
    if not en:
        return None
    ks = [i for i in en if i <= k]
    return en[max(ks)] if ks else en[min(en)]


def times_union(sers, hold=8, maxn=None):
    ts = sorted({s for ser in sers for s, _ in ser})
    if maxn and len(ts) > maxn:
        ts = ts[:: math.ceil(len(ts) / maxn)] + [ts[-1]]
    return ts + [ts[-1]] * hold


# ------------------------------------------------------------------------------------------- grid
def cell(layers_fn, label_fn, ts, az=-35, el=25, track=False, margin=1.1):
    """A panel of the video: per time t, its layers and its label."""
    return {"panels": [M.panel(layers_fn(t), az=az, el=el, margin=margin) for t in ts],
            "labels": [label_fn(t) for t in ts], "track": track, "az": az, "el": el, "margin": margin}


def render_cell(name, i, c, cw, ch):
    d = FRAMES / f"rsv_{name}_c{i}"
    if d.exists():
        shutil.rmtree(d)
    uniq, idx = [], []
    for p in c["panels"]:
        k = json.dumps(p, sort_keys=True)
        if not uniq or uniq[-1][0] != k:
            uniq.append((k, p))
        idx.append(len(uniq) - 1)
    sp = M.spec_file(f"vid_{name}_c{i}", {"panels": [p for _, p in uniq], "az": c["az"], "el": c["el"],
                                        "margin": c["margin"], "track": c["track"]})
    M.blender(["anim", sp, f"rsv_{name}_c{i}", cw, ch], f"vid_{name}_c{i}.log")
    return d, idx


def video(name, cells, cols, cw, ch, fps, title=None, title_color="#1f2329", ps=None, gif_width=800):
    """Render the cells (parallel blender jobs), label and tile them, encode MP4 + GIF."""
    T = len(cells[0]["panels"])
    ps = ps or max(14, cw // 26)
    with ThreadPoolExecutor(JOBS) as ex:
        rend = list(ex.map(lambda a: render_cell(name, a[0], a[1], cw, ch), enumerate(cells)))
    out = FRAMES / f"rsv_{name}"
    if out.exists():
        shutil.rmtree(out)
    out.mkdir(parents=True)
    rows = (len(cells) + cols - 1) // cols

    def compose(t):
        tiles = []
        for i, c in enumerate(cells):
            d, idx = rend[i]
            tile = out / f"t{t:04d}_{i:02d}.png"
            post.composite(d / f"{idx[t]:04d}.png", tile, c["labels"][t], pointsize=ps)
            tiles.append(str(tile))
        f = out / f"{t:04d}.png"
        if len(tiles) == 1:
            shutil.move(tiles[0], f)
        else:
            subprocess.run(["montage", *tiles, "-tile", f"{cols}x{rows}", "-geometry", f"{cw}x{ch}+0+0",
                            "-background", post.BG, str(f)], check=True)
            for x in tiles:
                os.remove(x)
        if title:
            subprocess.run(["convert", str(f), "-font", post.FONT, "-pointsize", str(int(ps * 1.25)),
                            "-fill", title_color, "-gravity", "north", "-annotate", "+0+10", title, str(f)], check=True)

    with ThreadPoolExecutor(12) as ex:
        list(ex.map(compose, range(T)))
    post.video(out, MEDIA / f"repulsive_surfaces_{name}", None, fps=fps, gif_width=gif_width)


def fmt_e(e):
    return "" if e is None else f"   E {e:.4g}"


# ------------------------------------------------------------------------------------- dodecahedron
def do_dodeca_pair(name, tag, ours_dir, ref_dir, ours_note, ref_note, ours_stop, ref_stop):
    ours = series(DODECA / ours_dir)
    en_o = csv_energy(DODECA / ours_dir)
    ref = sorted((step_num(f), f) for f in (DODECA / ref_dir / "out").glob("step_*.obj"))
    ts = list(range(0, max(ours[-1][0], ref[-1][0]) + 1, 100))
    ts = ts + [ts[-1]] * 8
    # ours' last frame is final.obj at its stop step; shown from the first t at or after it
    def pick_ours(t):
        return ours[-1] if t >= ours[-1][0] else at(ours, t)
    def lo(t):
        s, f = pick_ours(t)
        end = "   stopped" if f.name == "final.obj" else ""
        return f"Ours (ddg-tools)   step {s}{fmt_e(energy_near(en_o, s))}{end}\n{ours_note}"
    def lr(t):
        s, f = at(ref, t)
        return f"Reference program (authors' demo, MIT)   step {s}\n{ref_note}"
    kw = dict(az=-25, el=30, track=True)
    cells = [cell(lambda t: [{"obj": str(pick_ours(t)[1])}], lo, ts, **kw),
             cell(lambda t: [{"obj": str(at(ref, t)[1])}], lr, ts, **kw)]
    video(name, cells, 2, 960, 720, fps=5,
          title=f"Genus 11 from a chain of 11 fused tori ({tag}): area and barycenter fixed", ps=24)


def do_dodeca():
    do_dodeca_pair("dodecahedron", "hole size 0.25", "ours_h25", "ref_h25",
                   "reaches the dodecahedral cage (E 76.5)", "reaches the dodecahedral cage (E 75.8)", 2601, 3000)


def do_dodeca_perturbed():
    do_dodeca_pair("dodecahedron_perturbed", "hole size 0.30, start perturbed by 0.05 edge lengths, seed 3",
                   "ours_h30_p005_s3", "ref_h30_p005_s3",
                   "ends elongated (E 80.8), no dodecahedron", "reaches the dodecahedron (E 75.0)", 3455, 3000)


# --------------------------------------------------------------------------------------- genus
def genus_series(g):
    d = GENUS / f"g{g}"
    return series(d), csv_energy(d)


def do_genus():
    data = [genus_series(g) for g in range(12)]
    ts = list(range(0, max(s[-1][0] for s, _ in data) + 1, 50))
    ts = ts + [ts[-1]] * 10
    cells = []
    for g, (ser, en) in enumerate(data):
        end = ser[-1][0]
        def pick(t, ser=ser, end=end):
            return ser[-1] if t >= end else at(ser, t)
        def lab(t, g=g, pick=pick, en=en, end=end):
            s, _ = pick(t)
            tail = "   converged" if t >= end else ""
            return f"g = {g}   step {s}{fmt_e(energy_near(en, s))}{tail}"
        cells.append(cell(lambda t, pick=pick: [{"obj": str(pick(t)[1])}], lab, ts, az=-35, el=30, track=True))
    video("genus", cells, 4, 480, 360, fps=8, ps=17, gif_width=960)


def plane_layer():
    return {"box": {"center": [0, 0, -0.01], "half": [2.6, 2.6, 0.01]}, "alpha": 0.16,
            "color": [0.55, 0.65, 0.85], "roughness": 0.2, "frame": False}


def do_genus_plane():
    runs = [(g, w) for w in (5, 100) for g in range(2, 7)]
    sers = {r: series(RUNS / f"genus{r[0]:02d}_plane{r[1]}") for r in runs}
    ens = {r: M.steps(f"genus{r[0]:02d}_plane{r[1]}") for r in runs}
    ts = sorted({s for ser in sers.values() for s, _ in ser})
    ts = ts + [ts[-1]] * 10
    cells = []
    for r in runs:
        ser, en = sers[r], ens[r]
        lab = lambda t, r=r, ser=ser, en=en: (
            f"g = {r[0]}, plane weight {r[1]}   step {at(ser, t)[0]}   E {M.energy_at(en, at(ser, t)[0]):.4g}")
        cells.append(cell(lambda t, ser=ser: [{"obj": str(at(ser, t)[1])}, plane_layer()], lab, ts,
                          az=-35, el=18, track=False))
    video("genus_plane", cells, 5, 384, 360, fps=12, ps=14, gif_width=960,
          title="Genus 2..6 with an attracting plane: weight 5 (top), 100 (bottom)")


# ------------------------------------------------------------------------------------ envelopes
def do_envelopes():
    spot = {"obj": str(SCENES / "Reconstruction" / "spot_triangulated.obj"), "color": GOLD}
    lv = [([0.85, 0.35, 0.55], 0.5), ([0.16, 0.58, 0.32], 0.35), ([0.10, 0.62, 0.66], 0.25), ([0.12, 0.34, 0.80], 0.16)]
    sers = [series(RUNS / f"envelope_level{k}") for k in range(1, 5)]
    frames = []  # (level index, step, path)
    for k, ser in enumerate(sers):
        for s, f in ser:
            frames.append((k, s, f))
    frames += [frames[-1]] * 8
    def layers(j):
        k, s, f = frames[j]
        out = [spot]
        for q in range(k):  # earlier envelopes stay as translucent shells
            out.append({"obj": str(sers[q][-1][1]), "alpha": lv[q][1], "color": lv[q][0], "roughness": 0.2})
        out.append({"obj": str(f), "alpha": 0.6, "color": lv[k][0], "roughness": 0.2})
        return out
    ts = list(range(len(frames)))
    c = cell(layers, lambda j: f"Nested envelopes of the spot: level {frames[j][0] + 1}   step {frames[j][1]}", ts,
             az=-35, el=18)
    video("envelopes", [c], 1, 1280, 720, fps=9, ps=26, gif_width=800)


# ------------------------------------------------------------------------------------ shrinkwrap
class Rng:
    """splitmix64 as in ddg_repulsive_surfaces::generate::Rng."""
    M64 = (1 << 64) - 1

    def __init__(self, seed):
        self.s = seed

    def unit(self):
        self.s = (self.s + 0x9E3779B97F4A7C15) & self.M64
        z = self.s
        z = ((z ^ (z >> 30)) * 0xBF58476D1CE4E5B9) & self.M64
        z = ((z ^ (z >> 27)) * 0x94D049BB133111EB) & self.M64
        z ^= z >> 31
        return (z >> 11) / float(1 << 53)


def write_soup(path):
    """The polygon soup of shrinkwrap_bunny_soup: 35 % of the bunny's faces deleted (seed 7)."""
    V, F = read_obj(SCENES / "Reconstruction" / "bunny-scaled.obj")
    rng = Rng(7)
    keep = [f for f in F if rng.unit() >= 0.35]
    with open(path, "w") as o:
        for f in keep:
            for i in f:
                o.write("v %.9g %.9g %.9g\n" % tuple(V[i]))
        for k in range(len(keep)):
            o.write(f"f {3 * k + 1} {3 * k + 2} {3 * k + 3}\n")
    return len(keep), len(F)


def read_obj(path):
    V, F = [], []
    for line in open(path):
        if line.startswith("v "):
            V.append([float(x) for x in line.split()[1:4]])
        elif line.startswith("f "):
            idx = [int(t.split("/")[0]) - 1 for t in line.split()[1:]]
            for k in range(1, len(idx) - 1):
                F.append([idx[0], idx[k], idx[k + 1]])
    return np.array(V), F


def do_shrinkwrap():
    M.WORK.mkdir(parents=True, exist_ok=True)
    soup = M.WORK / "bunny_soup.obj"
    kept, tot = write_soup(soup)
    bunny = str(SCENES / "Reconstruction" / "bunny-scaled.obj")
    spot = str(SCENES / "Reconstruction" / "spot_triangulated.obj")
    runs = [("shrinkwrap_bunny_points", "Bunny point cloud", {"points": bunny, "radius": 0.006, "color": [0.12, 0.34, 0.80],
                                                              "index": list(range(0, M.nverts(bunny), 3))}),
            ("shrinkwrap_bunny_soup", f"Bunny polygon soup ({kept} of {tot} faces)",
             {"obj": str(soup), "color": [0.12, 0.34, 0.80], "flat": True, "alpha": 1.0}),
            ("shrinkwrap_spot_points", "Spot point cloud", {"points": spot, "radius": 0.006, "color": [0.12, 0.34, 0.80],
                                                          "index": list(range(0, M.nverts(spot), 3))})]
    sers = [series(RUNS / r[0]) for r in runs]
    ts = times_union(sers, hold=8)
    cells = []
    for (run, name, obst), ser in zip(runs, sers):
        en = csv_energy(RUNS / run)
        cells.append(cell(lambda t, ser=ser, obst=obst: [
            {"obj": str(at(ser, t)[1]), "alpha": 0.4, "color": GOLD, "roughness": 0.3}, obst],
            lambda t, ser=ser, name=name, en=en: f"{name}   step {at(ser, t)[0]}{fmt_e(energy_near(en, at(ser, t)[0]))}",
            ts, az=-35, el=20))
    video("shrinkwrap", cells, 3, 640, 600, fps=5, ps=18, gif_width=960,
          title="Shrink-wrapping: a sphere falls onto the data (translucent: the cloud stays visible)")


# ------------------------------------------------------------------------------------------- h pins
def do_h_pins():
    runs = [("h_sparse_volume", "sparse pins, volume x2"), ("h_sparse_area", "sparse pins, area x1.5"),
            ("h_sparse_sqerr_volume", "squared-error penalty, volume x2"), ("h_nopins_volume", "no pins, volume x2"),
            ("h_dense_volume", "dense pins, volume x2"), ("h_sparse_volume_box", "sparse pins, volume x2, box"),
            ("h_dense_volume_box", "dense pins, volume x2, box")]
    box = {"box": {"center": [1.5, 2.0, 0.5], "half": [1.9, 2.4, 0.9]}, "alpha": 0.1, "color": [0.55, 0.65, 0.85],
           "roughness": 0.2}
    cells = []
    sers = [series(RUNS / r) for r, _ in runs]
    ts = times_union(sers, hold=8)
    for (run, lab), ser in zip(runs, sers):
        first = ser[0][1]
        V0, _ = read_obj(first)
        V1, _ = read_obj(ser[-1][1])
        mv = np.linalg.norm(V1 - V0, axis=1)
        pins = list(range(24)) if "sparse" in run and "sqerr" not in run else (
            [int(i) for i in np.where(mv < 1e-9)[0]] if "dense" in run else [])
        def layers(t, ser=ser, run=run, pins=pins):
            f = str(at(ser, t)[1])
            out = [{"obj": f, "wire": True, "wire_thickness": 0.0012}]
            if pins:
                out.append({"points": f, "index": pins, "radius": 0.012, "color": [0.07, 0.08, 0.10]})
            if "box" in run:
                out.append(box)
            return out
        cells.append(cell(layers, lambda t, ser=ser, lab=lab: f"{lab}   step {at(ser, t)[0]}", ts, az=-40, el=35))
    video("h_pins", cells, 4, 480, 400, fps=6, ps=15, gif_width=960,
          title="The H of Fig. 22: pins, growth, and a confining box (black dots: pinned vertices)")


# ------------------------------------------------------------------------------------------ slideoff
def do_slideoff():
    torus = {"torus": {"R": 0.8, "r": 0.1, "center": [0, 0, 0]}, "color": [0.12, 0.34, 0.80], "roughness": 0.4}
    plane = {"box": {"center": [0, 0, 0.19], "half": [1.6, 1.6, 0.005]}, "alpha": 0.14, "color": [0.55, 0.65, 0.85],
             "roughness": 0.2, "frame": False}
    frames = []
    for stage, run, tag in ((1, "slideoff_stage1", "stage 1, torus bar weight 1"),
                            (2, "slideoff_stage2", "stage 2, bar weight 0.001, attracting plane")):
        st = M.steps(run)
        fs = M.frames(run)
        for f in fs[::2] + [fs[-1]]:
            k = M.step_of(f)
            frames.append((f, tag, k, M.energy_at(st, k), stage))
    frames += [frames[-1]] * 8
    def layers(j):
        out = [{"obj": str(frames[j][0])}, torus]
        if frames[j][4] == 2:
            out.append(plane)
        return out
    c = cell(layers, lambda j: f"Handle sliding off a bar, {frames[j][1]}   step {frames[j][2]}   energy {frames[j][3]:.4g}",
             list(range(len(frames))), az=-40, el=30)
    video("slideoff", [c], 1, 1280, 720, fps=14, ps=24, gif_width=800)


# --------------------------------------------------------------------------------------- cobblestones
def do_cobblestones():
    slab = {"box": {"center": [0, 0, 0], "half": [3.2, 3.2, 0.9]}, "alpha": 0.12, "color": [0.55, 0.65, 0.85],
            "roughness": 0.2}
    runs = [("cobblestones", "slab weight 1"), ("cobblestones_w30", "slab weight 30"),
            ("cobblestones_w300", "slab weight 300"), ("cobblestones_noremesh", "slab weight 1, no remeshing")]
    sers = [series(RUNS / r) for r, _ in runs]
    ts = times_union(sers, hold=8)
    cells = []
    for (run, lab), ser in zip(runs, sers):
        cells.append(cell(lambda t, ser=ser: [{"obj": str(at(ser, t)[1]), "color": "components"}, slab],
                          lambda t, ser=ser, lab=lab: f"14 growing spheres in a slab, {lab}   step {at(ser, t)[0]}",
                          ts, az=-30, el=32))
    video("cobblestones", cells, 2, 960, 540, fps=6, ps=19, gif_width=960)


# --------------------------------------------------------------------------------------------- klein
def do_klein():
    ser = series(RUNS / "klein_subcritical")
    en = csv_energy(RUNS / "klein_subcritical")
    ts = [s for s, _ in ser] + [ser[-1][0]] * 8
    c = cell(lambda t: [{"obj": str(at(ser, t)[1])}],
             lambda t: f"Klein bottle, exponents (3.5, 7) below the critical ones   step {at(ser, t)[0]}{fmt_e(energy_near(en, at(ser, t)[0]))}",
             ts, az=-35, el=25)
    video("klein", [c], 1, 1280, 720, fps=10, ps=24, gif_width=800)


def do_spindles():
    runs = [("subcritical_spindle_p3.5", "(3.5, 7), r = 1.02"), ("subcritical_spindle_p3.5_r1.01", "(3.5, 7), r = 1.01"),
            ("subcritical_spindle_p3.5_r1.005", "(3.5, 7), r = 1.005"), ("subcritical_spindle_p3.0", "(3, 6), r = 1.02"),
            ("subcritical_spindle_p3.0_r1.01", "(3, 6), r = 1.01")]
    sers = [series(RUNS / r) for r, _ in runs]
    ts = times_union(sers, hold=8, maxn=80)
    cells = []
    for (run, lab), ser in zip(runs, sers):
        cells.append(cell(lambda t, ser=ser: [{"obj": str(at(ser, t)[1])}],
                          lambda t, ser=ser, lab=lab: f"Spindle torus, {lab}   step {at(ser, t)[0]}",
                          ts, az=-35, el=30, track=False))
    video("subcritical_spindles", cells, 3, 640, 540, fps=8, ps=17, gif_width=960,
          title="Subcritical exponents on a slightly too fat spindle torus (R = 1)")


# ------------------------------------------------------------------------------------------- eversion
def do_eversion():
    d = RUNS / "eversion"
    fs = sorted(d.glob("evert_*.obj"))
    n = len(fs)
    half = n // 2
    ts = list(range(n)) + [n - 1] * 10
    def lab(j):
        k = min(j, n - 1)
        if k < half:
            return f"Punctured torus, volume penalty: forward flow   frame {k + 1} of {half}"
        return f"Mirrored and rotated copy of the flow, played back (constructed, not a flow result)   frame {k + 1 - half} of {n - half}"
    c = cell(lambda j: [{"obj": str(fs[min(j, n - 1)])}], lab, ts, az=-35, el=25)
    video("eversion", [c], 1, 1280, 720, fps=20, ps=22, gif_width=800,
          title="NOT REPRODUCED: the flow did not complete the eversion (the end state is not inside out)",
          title_color="#b3261e")


ALL = {"dodeca": do_dodeca, "dodeca_perturbed": do_dodeca_perturbed, "genus": do_genus, "genus_plane": do_genus_plane,
       "envelopes": do_envelopes, "shrinkwrap": do_shrinkwrap, "h_pins": do_h_pins, "slideoff": do_slideoff,
       "cobblestones": do_cobblestones, "klein": do_klein, "spindles": do_spindles, "eversion": do_eversion}


def main():
    sel = sys.argv[1:] or list(ALL)
    MEDIA.mkdir(parents=True, exist_ok=True)
    for n in sel:  # one video at a time: its cells already run JOBS blender processes in parallel
        try:
            ALL[n]()
        except Exception as e:
            print("FAILED", n, repr(e), flush=True)


if __name__ == "__main__":
    main()
