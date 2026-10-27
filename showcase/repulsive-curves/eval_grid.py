# /// script
# requires-python = ">=3.10"
# dependencies = ["numpy", "matplotlib"]
# ///
"""Grid video of the Sec. 7 comparison (SPEC.md 14.2, Q6): one column per method, one knot, the
curves at equal wall-clock times, from the `.rcf` frame files of `examples/eval_knots.rs`.

  uv run showcase/repulsive-curves/eval_grid.py --frames-dir runs_tmp/rc-eval/runB/frames \
      --knot knot0102 --out showcase/media/rc_eval_knot0102.mp4 [--fps 30] [--seconds 60]
      [--width 2400 --height 420] [--methods hs-pg,l2-pg,...] [--still-at 0.5]

Every column is the last record with wall time <= t (the file format is in the header of
eval_knots.rs), so all columns show the same wall-clock instant; a run that has ended shows its
final curve with "done". The video lasts `--seconds` and covers the slowest run, unless
`--wall-max S` caps the wall time shown. A fixed camera (the union box of all records of the
knot) turning slowly; tubes are drawn as depth-shaded line segments (matplotlib Agg, no Blender),
which renders thousands of frames in minutes. `render.py` (Blender Cycles tubes) is for single
curves and per-run stills; it takes OBJ frames and one scene kind, which does not fit nine
panels x hundreds of frames, so this script is separate. `--still-at F` writes one PNG of the grid
at fraction F of the video and exits.
"""
import argparse
import glob
import os
import struct
import subprocess

import numpy as np

ORDER = ["hs-pg", "l2-pg", "h1-pg", "h2-pg", "hs-armijo",
         "soft-hs-ncg", "soft-hs-lbfgs", "soft-h1-ncg", "soft-h1-lbfgs"]
TITLE = {"hs-pg": "Hs grad.", "l2-pg": "L2 grad.", "h1-pg": "H1 grad.", "h2-pg": "H2 grad.",
         "hs-armijo": "Hs Armijo", "soft-hs-ncg": "Hs NCG", "soft-hs-lbfgs": "Hs L-BFGS",
         "soft-h1-ncg": "H1 NCG", "soft-h1-lbfgs": "H1 L-BFGS"}


def read_rcf(path):
    b = open(path, "rb").read()
    if b[:4] != b"RCF1":
        raise ValueError(f"{path}: not an RCF1 file")
    n, m = struct.unpack_from("<II", b, 4)
    off = 20
    edges = np.frombuffer(b, dtype="<u4", count=2 * m, offset=off).reshape(m, 2)
    off += 8 * m
    rec = 12 + 12 * n
    k = (len(b) - off) // rec
    T, I, F, P = np.empty(k), np.empty(k, int), np.empty(k, int), np.empty((k, n, 3), np.float32)
    for i in range(k):
        t, it, fl = struct.unpack_from("<fII", b, off + i * rec)
        T[i], I[i], F[i] = t, it, fl
        P[i] = np.frombuffer(b, dtype="<f4", count=3 * n, offset=off + i * rec + 12).reshape(n, 3)
    return dict(edges=edges, t=T, it=I, final=F, pos=P)


def rotation(az, el):
    ca, sa, ce, se = np.cos(az), np.sin(az), np.cos(el), np.sin(el)
    Rz = np.array([[ca, -sa, 0], [sa, ca, 0], [0, 0, 1]])
    Rx = np.array([[1, 0, 0], [0, ce, -se], [0, se, ce]])
    return Rx @ Rz


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--frames-dir", required=True)
    ap.add_argument("--knot", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--fps", type=int, default=30)
    ap.add_argument("--seconds", type=float, default=60.0)
    ap.add_argument("--wall-max", type=float)
    ap.add_argument("--width", type=int, default=2700)
    ap.add_argument("--height", type=int, default=420)
    ap.add_argument("--methods")
    ap.add_argument("--results-dir", help="results/ of the run: adds the energy to the labels")
    ap.add_argument("--still-at", type=float)
    a = ap.parse_args()
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.collections import LineCollection

    want = a.methods.split(",") if a.methods else ORDER
    runs = {}
    for m in want:
        f = os.path.join(a.frames_dir, f"{a.knot}.{m}.rcf")
        if os.path.exists(f):
            runs[m] = read_rcf(f)
    import json
    energy = {}
    if a.results_dir:
        for m in list(runs):
            f = os.path.join(a.results_dir, f"{a.knot}.{m}.jsonl")
            if os.path.exists(f):
                h = json.loads(open(f).readline())["history"]
                energy[m] = {row[0]: row[3] for row in h}
                runs[m]["status"] = json.loads(open(f).readline())["end_reason"]
    if not runs:
        raise SystemExit(f"no .rcf files for {a.knot} in {a.frames_dir}")
    methods = [m for m in want if m in runs]
    wall = a.wall_max or max(r["t"][-1] for r in runs.values())
    allp = np.concatenate([r["pos"].reshape(-1, 3) for r in runs.values()])
    ctr = (allp.min(0) + allp.max(0)) / 2
    rad = float(np.linalg.norm(allp.max(0) - allp.min(0))) / 2

    ncol = len(methods)
    dpi = 100
    fig, axs = plt.subplots(1, ncol, figsize=(a.width / dpi, a.height / dpi), dpi=dpi, squeeze=False)
    fig.patch.set_facecolor("white")
    nframes = int(a.seconds * a.fps)
    arts = []
    for ax, m in zip(axs[0], methods):
        ax.set_xlim(-rad * 0.75, rad * 0.75)
        ax.set_ylim(-rad * 0.75, rad * 0.75)
        ax.set_aspect("equal")
        ax.axis("off")
        lc = LineCollection([], linewidths=1.2)
        ax.add_collection(lc)
        title = ax.text(0.5, 1.0, TITLE[m], transform=ax.transAxes, ha="center", va="top", fontsize=10)
        info = ax.text(0.5, 0.0, "", transform=ax.transAxes, ha="center", va="bottom", fontsize=8,
                       family="monospace")
        arts.append((lc, info))
    fig.subplots_adjust(left=0, right=1, bottom=0, top=1, wspace=0)

    def draw(fi):
        t = wall * fi / max(nframes - 1, 1)
        az = 0.4 + 2 * np.pi * 0.5 * fi / max(nframes, 1)
        R = rotation(az, 0.9)
        for m, (lc, info) in zip(methods, arts):
            r = runs[m]
            k = int(np.searchsorted(r["t"], t, side="right")) - 1
            k = max(k, 0)
            P = (r["pos"][k] - ctr) @ R.T
            E = r["edges"]
            seg = np.stack([P[E[:, 0], :2], P[E[:, 1], :2]], 1)
            z = (P[E[:, 0], 2] + P[E[:, 1], 2]) / 2
            order = np.argsort(z)
            zz = (z[order] / rad + 1) / 2
            shade = 0.15 + 0.7 * (1 - np.clip(zz, 0, 1))
            lc.set_segments(seg[order])
            lc.set_color(np.stack([shade * 0.3, shade * 0.4, np.full_like(shade, 0.5) + 0.5 * (1 - shade) * 0, np.ones_like(shade)], 1))
            done = k == len(r["t"]) - 1 and r["final"][k] == 1 and t >= r["t"][-1]
            e = energy.get(m, {}).get(int(r["it"][k]))
            es = f"\nE {e:9.2f}" if e is not None else ""
            tail = ("  done: " + r.get("status", "")) if done else ""
            info.set_text(f"it {r['it'][k]:5d}  {min(t, r['t'][-1]):6.1f}s{es}{tail}")
        fig.canvas.draw()
        return np.asarray(fig.canvas.buffer_rgba())[:, :, :3]

    if a.still_at is not None:
        img = draw(int(a.still_at * (nframes - 1)))
        plt.imsave(a.out, img)
        print("wrote", a.out)
        return
    os.makedirs(os.path.dirname(os.path.abspath(a.out)), exist_ok=True)
    h, w = draw(0).shape[:2]
    ff = subprocess.Popen(["ffmpeg", "-y", "-loglevel", "error", "-f", "rawvideo", "-pix_fmt", "rgb24",
                           "-s", f"{w}x{h}", "-r", str(a.fps), "-i", "-", "-c:v", "libx264",
                           "-pix_fmt", "yuv420p", "-crf", "20", a.out], stdin=subprocess.PIPE)
    for fi in range(nframes):
        ff.stdin.write(np.ascontiguousarray(draw(fi)).tobytes())
    ff.stdin.close()
    ff.wait()
    print(f"wrote {a.out}: {nframes} frames, wall time shown {wall:.1f} s, {len(methods)} columns")


main()
