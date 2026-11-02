#!/usr/bin/env python3
"""Vertex statistics of OBJ files: count, mean, bounding box, smallest distance to the z axis and
the number of vertices closer to it than R (the pole of the HandcuffsPole scene has radius 0.075).

    python3 -I obj_stats.py [--radius R] FILE.obj ...
"""
import math
import sys


def stats(path, radius):
    pts = []
    with open(path) as f:
        for line in f:
            if line.startswith("v "):
                pts.append(tuple(float(t) for t in line.split()[1:4]))
    n = len(pts)
    mean = [sum(p[k] for p in pts) / n for k in range(3)]
    lo = [min(p[k] for p in pts) for k in range(3)]
    hi = [max(p[k] for p in pts) for k in range(3)]
    rad = [math.hypot(p[0], p[1]) for p in pts]
    close = sum(1 for r in rad if r < radius)
    print(
        f"{path}: n={n} mean=({mean[0]:.3f},{mean[1]:.3f},{mean[2]:.3f}) "
        f"box=[{lo[0]:.2f},{hi[0]:.2f}]x[{lo[1]:.2f},{hi[1]:.2f}]x[{lo[2]:.2f},{hi[2]:.2f}] "
        f"min_dist_to_z_axis={min(rad):.4f} vertices_within_{radius}={close}"
    )


if __name__ == "__main__":
    args = sys.argv[1:]
    radius = 0.1
    if args and args[0] == "--radius":
        radius = float(args[1])
        args = args[2:]
    for a in args:
        stats(a, radius)
