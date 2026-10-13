#!/usr/bin/env python3
"""Small generated meshes for the intrinsic-triangulation reference runs (standard library only).

usage: uv run gen_meshes.py OUTDIR
Writes OUTDIR/{icosphere,icosphere_jitter,thin_cylinder,long_torus,sheared_grid}.obj. Faces are counter-clockwise seen from outside
(outward normals for the closed ones), vertices are written with 17 significant digits.
"""
import math
import sys
from pathlib import Path


def write_obj(path, verts, faces):
    with open(path, "w") as f:
        for v in verts:
            f.write("v %.17g %.17g %.17g\n" % tuple(v))
        for t in faces:
            f.write("f %d %d %d\n" % tuple(i + 1 for i in t))


def icosphere(levels=1):
    # Coarse icosphere: icosahedron subdivided `levels` times (1 -> 42 vertices, 80 faces), projected to the unit sphere.
    p = (1 + 5 ** 0.5) / 2
    v = [(-1, p, 0), (1, p, 0), (-1, -p, 0), (1, -p, 0), (0, -1, p), (0, 1, p), (0, -1, -p), (0, 1, -p),
         (p, 0, -1), (p, 0, 1), (-p, 0, -1), (-p, 0, 1)]
    v = [tuple(c / math.sqrt(sum(x * x for x in q)) for c in q) for q in v]
    f = [(0, 11, 5), (0, 5, 1), (0, 1, 7), (0, 7, 10), (0, 10, 11), (1, 5, 9), (5, 11, 4), (11, 10, 2), (10, 7, 6),
         (7, 1, 8), (3, 9, 4), (3, 4, 2), (3, 2, 6), (3, 6, 8), (3, 8, 9), (4, 9, 5), (2, 4, 11), (6, 2, 10),
         (8, 6, 7), (9, 8, 1)]
    for _ in range(levels):
        cache, nf = {}, []

        def mid(a, b):
            k = (min(a, b), max(a, b))
            if k not in cache:
                m = [(v[a][i] + v[b][i]) / 2 for i in range(3)]
                n = math.sqrt(sum(x * x for x in m))
                v.append(tuple(x / n for x in m))
                cache[k] = len(v) - 1
            return cache[k]

        for a, b, c in f:
            ab, bc, ca = mid(a, b), mid(b, c), mid(c, a)
            nf += [(a, ab, ca), (b, bc, ab), (c, ca, bc), (ab, bc, ca)]
        f = nf
    return v, f


def jittered_icosphere(levels=1, amp=0.3, seed=12345):
    # Same sphere with vertices moved radially and sideways by a deterministic LCG, so that it is not Delaunay.
    v, f = icosphere(levels)
    state = [seed]

    def rnd():
        state[0] = (state[0] * 6364136223846793005 + 1442695040888963407) % (1 << 64)
        return (state[0] >> 11) / float(1 << 53) * 2 - 1

    v = [tuple(c * (1 + amp * rnd()) + 0.5 * amp * rnd() * 0.5 for c in p) for p in v]
    return v, f


def thin_cylinder(n_around=4, n_along=12, radius=0.1, length=10.0, twist=0.37):
    # Open tube (two boundary loops), long and thin: edges along the axis are ~ 0.83, around the axis ~ 0.14.
    # Ring i is rotated by twist * i steps: without the twist all quads are rectangles, i.e. cocircular ties, and
    # Delaunay flipping would do nothing.
    v = [(radius * math.cos(2 * math.pi * (j + twist * i) / n_around),
          radius * math.sin(2 * math.pi * (j + twist * i) / n_around),
          length * i / n_along) for i in range(n_along + 1) for j in range(n_around)]
    f = []
    for i in range(n_along):
        for j in range(n_around):
            a, b = i * n_around + j, i * n_around + (j + 1) % n_around
            c, d = a + n_around, b + n_around
            f += [(a, b, d), (a, d, c)]
    return v, f


def long_torus(n_major=8, n_minor=3, R=4.0, r=0.2, twist=0.37):
    # Coarse torus with a thin tube: n_minor vertices round the tube, n_major round the big circle, each ring
    # rotated by twist * i steps (breaks the cocircular ties). Delaunay flipping produces self-edges and double edges.
    v = []
    for i in range(n_major):
        u = 2 * math.pi * i / n_major
        for j in range(n_minor):
            w = 2 * math.pi * (j + twist * i) / n_minor
            rr = R + r * math.cos(w)
            v.append((rr * math.cos(u), rr * math.sin(u), r * math.sin(w)))
    f = []
    for i in range(n_major):
        for j in range(n_minor):
            a = i * n_minor + j
            b = i * n_minor + (j + 1) % n_minor
            c = ((i + 1) % n_major) * n_minor + j
            d = ((i + 1) % n_major) * n_minor + (j + 1) % n_minor
            f += [(a, c, d), (a, d, b)]
    return v, f


def sheared_grid(n=10, shear=6.0):
    # Flat (n x n)-cell grid on the unit square, sheared x' = x + shear * y: triangles become very obtuse
    # (the shear is 6 grid cells wide over the whole height), so Delaunay flipping flips a lot.
    v = [((i / n) + shear * (j / n), j / n, 0.0) for j in range(n + 1) for i in range(n + 1)]
    f = []
    for j in range(n):
        for i in range(n):
            a = j * (n + 1) + i
            b, c, d = a + 1, a + n + 1, a + n + 2
            f += [(a, b, d), (a, d, c)]
    return v, f


if __name__ == "__main__":
    out = Path(sys.argv[1])
    out.mkdir(parents=True, exist_ok=True)
    for name, (v, f) in {"icosphere": icosphere(1), "icosphere_jitter": jittered_icosphere(), "thin_cylinder": thin_cylinder(), "long_torus": long_torus(),
                         "sheared_grid": sheared_grid()}.items():
        write_obj(out / (name + ".obj"), v, f)
        print(name, len(v), "vertices", len(f), "faces")
