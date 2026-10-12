# Dumps Blender's own triangulation of random polygons, for comparison with
# ddg-core's `triangulate_polygons` (see scripts/blender-tessellation-diff.sh).
#
#   blender -b --factory-startup --python scripts/blender_tessellation_dump.py -- OUT.txt [SEED] [SCALE]
#
# SCALE multiplies the number of polygons per group (default 1.0). The script only creates
# mesh data blocks in memory; it reads and writes nothing but OUT.txt.
#
# Output format (text, whitespace separated; one "group" per polygon family):
#
#   # comment lines
#   group <name> <expect>      expect is "exact" (non-degenerate input, bit-exact agreement is
#                              the goal) or "degenerate" (collinear / zero-area input)
#   nv <n>                     then n lines "x y z": positions as stored by Blender (float32)
#   np <p>                     then p lines "k v0 .. v(k-1)": polygon vertex indices, loop order
#   nt <t>                     then t lines "l0 l1 l2": Blender's loop_triangles[i].loops
#   end

import math
import random
import sys

import bpy


def rotation(rng):
    """A random rotation matrix (rows), from a random unit quaternion."""
    while True:
        q = [rng.gauss(0.0, 1.0) for _ in range(4)]
        norm = math.sqrt(sum(c * c for c in q))
        if norm > 1e-6:
            break
    w, x, y, z = (c / norm for c in q)
    return (
        (1 - 2 * (y * y + z * z), 2 * (x * y - z * w), 2 * (x * z + y * w)),
        (2 * (x * y + z * w), 1 - 2 * (x * x + z * z), 2 * (y * z - x * w)),
        (2 * (x * z - y * w), 2 * (y * z + x * w), 1 - 2 * (x * x + y * y)),
    )


def place(rng, points2d, noise=0.0, rotate=True):
    """Embeds 2D points in 3D: random rotation and offset, optional out-of-plane noise."""
    rot = rotation(rng) if rotate else ((1, 0, 0), (0, 1, 0), (0, 0, 1))
    offset = [rng.uniform(-2.0, 2.0) for _ in range(3)] if rotate else (0.0, 0.0, 0.0)
    out = []
    for (x, y) in points2d:
        p = (x, y, rng.uniform(-noise, noise) if noise else 0.0)
        out.append(tuple(sum(rot[r][c] * p[c] for c in range(3)) + offset[r] for r in range(3)))
    return out


def maybe_reverse(rng, pts):
    return list(reversed(pts)) if rng.random() < 0.5 else pts


def rotate_start(rng, pts):
    k = rng.randrange(len(pts))
    return pts[k:] + pts[:k]


def spread_angles(rng, n):
    """n increasing angles around the circle, one in each of n equal sectors, so that no two are
    closer than a tenth of a sector (and for n >= 4 no gap reaches pi)."""
    return [2.0 * math.pi * (i + rng.uniform(0.05, 0.95)) / n for i in range(n)]


def convex_ngon(rng, n):
    """Points on a circle at increasing angles: strictly convex (n >= 4)."""
    r = rng.uniform(0.5, 2.0)
    return [(r * math.cos(a), r * math.sin(a)) for a in spread_angles(rng, n)]


def star_ngon(rng, n, r_min=0.3):
    """Star-shaped (generally concave) simple polygon: increasing angles, random radii."""
    return [(r * math.cos(a), r * math.sin(a))
            for a, r in ((a, rng.uniform(r_min, 1.0)) for a in spread_angles(rng, n))]


def dart(rng):
    """A concave quad (arrow head)."""
    tip = (rng.uniform(1.0, 2.0), 0.0)
    notch = (rng.uniform(0.1, 0.8), rng.uniform(-0.1, 0.1))
    a = (rng.uniform(-0.5, 0.0), rng.uniform(0.5, 1.5))
    b = (rng.uniform(-0.5, 0.0), -rng.uniform(0.5, 1.5))
    return [tip, a, notch, b]


def rectilinear(rng):
    """A rectilinear "histogram" polygon on the integer grid: many exactly collinear corners
    and exact ties in every orientation test."""
    columns = rng.randrange(2, 7)
    heights = [rng.randrange(1, 5) for _ in range(columns)]
    pts = [(0, 0), (columns, 0)]
    for i in reversed(range(columns)):
        pts.append((i + 1, heights[i]))
        pts.append((i, heights[i]))
    # Drop consecutive duplicates (equal neighbouring heights give repeated corners).
    out = []
    for p in pts:
        if not out or out[-1] != p:
            out.append(p)
    if out[0] == out[-1]:
        out.pop()
    return [(float(x), float(y)) for (x, y) in out]


def comb(rng):
    """A comb: teeth of random width and height on a common base. Many reflex corners, and
    ears that contain other corners."""
    teeth = rng.randrange(2, 9)
    upper = []
    x = 0.0
    for t in range(teeth):
        width = rng.uniform(0.2, 1.0)
        height = rng.uniform(1.0, 3.0)
        lean = rng.uniform(-0.08, 0.08)
        upper.append((x + rng.uniform(0.0, 0.05), 0.4 + rng.uniform(0.0, 0.3)))
        upper.append((x + lean, height))
        upper.append((x + width + lean, height + rng.uniform(-0.3, 0.3)))
        upper.append((x + width - rng.uniform(0.0, 0.05), 0.4 + rng.uniform(0.0, 0.3)))
        x += width + rng.uniform(0.2, 1.0)
    return [(-0.2, 0.0), (x, 0.0)] + list(reversed(upper))


def spiral(rng):
    """A spiral band: out along one arm, back along a parallel one."""
    turns = rng.uniform(0.75, 3.0)
    steps = rng.randrange(6, 40)
    pitch = rng.uniform(0.5, 1.0)
    width = rng.uniform(0.2, 0.6) * pitch
    outer, inner = [], []
    for i in range(steps + 1):
        theta = 2.0 * math.pi * turns * i / steps
        r = 1.0 + pitch * theta / (2.0 * math.pi) + rng.uniform(-0.02, 0.02)
        outer.append((r * math.cos(theta), r * math.sin(theta)))
        inner.append(((r - width) * math.cos(theta), (r - width) * math.sin(theta)))
    return outer + list(reversed(inner))


def subdivide_edges(rng, pts):
    """Inserts points on the edges of a polygon: collinear corners."""
    out = []
    n = len(pts)
    for i in range(n):
        a, b = pts[i], pts[(i + 1) % n]
        out.append(a)
        for k in range(rng.randrange(0, 3)):
            t = (k + 1) / 3.0
            out.append((a[0] + t * (b[0] - a[0]), a[1] + t * (b[1] - a[1])))
    return out


def make_groups(rng, scale):
    def count(n):
        return max(1, int(n * scale))

    groups = []

    def group(name, expect, n, make):
        polys = []
        for _ in range(count(n)):
            polys.append(make())
        groups.append((name, expect, polys))

    # ---- quads ----
    group("quad_random_3d", "exact", 4000, lambda: [tuple(rng.uniform(-1, 1) for _ in range(3)) for _ in range(4)])
    group("quad_convex_planar", "exact", 3000,
          lambda: place(rng, rotate_start(rng, maybe_reverse(rng, convex_ngon(rng, 4)))))
    group("quad_concave_planar", "exact", 3000,
          lambda: place(rng, rotate_start(rng, maybe_reverse(rng, dart(rng)))))
    group("quad_concave_nonplanar", "exact", 3000,
          lambda: place(rng, rotate_start(rng, maybe_reverse(rng, dart(rng))), noise=0.2))
    group("quad_convex_nonplanar", "exact", 3000,
          lambda: place(rng, rotate_start(rng, maybe_reverse(rng, convex_ngon(rng, 4))), noise=0.5))
    group("quad_axis_aligned", "exact", 500,
          lambda: place(rng, rotate_start(rng, maybe_reverse(rng, [(0.0, 0.0), (float(rng.randrange(1, 4)), 0.0),
                                                                     (float(rng.randrange(1, 4)), float(rng.randrange(1, 4))),
                                                                     (0.0, float(rng.randrange(1, 4)))])), rotate=False))

    def degenerate_quad():
        kind = rng.randrange(4)
        if kind == 0:  # three collinear corners on the integer grid
            pts = [(0.0, 0.0), (1.0, 0.0), (2.0, 0.0), (float(rng.randrange(0, 3)), float(rng.randrange(1, 3)))]
        elif kind == 1:  # two coincident corners
            pts = [(0.0, 0.0), (1.0, 0.0), (1.0, 0.0), (0.0, 1.0)]
        elif kind == 2:  # all collinear
            pts = [(float(i), 0.0) for i in rng.sample(range(8), 4)]
        else:  # all coincident
            pts = [(1.0, 1.0)] * 4
        return place(rng, rotate_start(rng, pts), rotate=rng.random() < 0.5)

    group("quad_degenerate", "degenerate", 1000, degenerate_quad)

    # ---- n-gons ----
    group("ngon_convex_planar", "exact", 3000,
          lambda: place(rng, maybe_reverse(rng, convex_ngon(rng, rng.randrange(5, 17)))))
    group("ngon_star_planar", "exact", 6000,
          lambda: place(rng, maybe_reverse(rng, star_ngon(rng, rng.randrange(5, 25)))))
    group("ngon_star_nonplanar", "exact", 4000,
          lambda: place(rng, maybe_reverse(rng, star_ngon(rng, rng.randrange(5, 25))), noise=0.1))
    group("ngon_star_xy_plane", "exact", 2000,
          lambda: place(rng, maybe_reverse(rng, star_ngon(rng, rng.randrange(5, 25))), rotate=False))
    group("ngon_star_large", "exact", 150,
          lambda: place(rng, maybe_reverse(rng, star_ngon(rng, rng.randrange(50, 400), r_min=0.6))))
    group("ngon_convex_large", "exact", 100,
          lambda: place(rng, maybe_reverse(rng, [(math.cos(2 * math.pi * i / n), math.sin(2 * math.pi * i / n))
                                                  for n in [rng.randrange(32, 257)] for i in range(n)])))
    group("ngon_random_3d", "exact", 1000,
          lambda: [tuple(rng.uniform(-1, 1) for _ in range(3)) for _ in range(rng.randrange(5, 10))])

    group("ngon_comb", "exact", 3000, lambda: place(rng, rotate_start(rng, maybe_reverse(rng, comb(rng)))))
    group("ngon_comb_nonplanar", "exact", 1500,
          lambda: place(rng, rotate_start(rng, maybe_reverse(rng, comb(rng))), noise=0.05))
    group("ngon_spiral", "exact", 3000, lambda: place(rng, rotate_start(rng, maybe_reverse(rng, spiral(rng)))))
    # Self-intersecting planar polygons: not valid faces, but Blender still has to answer.
    group("ngon_random_planar", "exact", 2000,
          lambda: place(rng, [(rng.uniform(-1, 1), rng.uniform(-1, 1)) for _ in range(rng.randrange(5, 12))]))

    def quantized_star():
        # A small polygon far from the origin: single precision leaves only a few bits per
        # coordinate, so corners collapse onto a coarse grid (collinear and duplicate corners).
        pts = place(rng, star_ngon(rng, rng.randrange(5, 16)))
        shift = [rng.uniform(-1.0, 1.0) * 4096.0 for _ in range(3)]
        return [tuple(shift[k] + 0.002 * p[k] for k in range(3)) for p in pts]

    group("ngon_star_quantized", "degenerate", 2000, quantized_star)

    # Collinear corners: exact ties, where bit-exact agreement is desirable but not required.
    group("ngon_rectilinear_xy", "degenerate", 2000,
          lambda: place(rng, rotate_start(rng, maybe_reverse(rng, rectilinear(rng))), rotate=False))
    group("ngon_rectilinear_rotated", "degenerate", 2000,
          lambda: place(rng, rotate_start(rng, maybe_reverse(rng, rectilinear(rng)))))
    group("ngon_subdivided_edges", "degenerate", 2000,
          lambda: place(rng, maybe_reverse(rng, subdivide_edges(rng, star_ngon(rng, rng.randrange(3, 10))))))

    def degenerate_ngon():
        n = rng.randrange(5, 10)
        kind = rng.randrange(3)
        if kind == 0:  # all collinear
            pts = [(float(i), 0.0) for i in rng.sample(range(20), n)]
        elif kind == 1:  # repeated corners
            base = star_ngon(rng, n)
            pts = []
            for p in base:
                pts.append(p)
                if rng.random() < 0.3:
                    pts.append(p)
        else:  # all coincident
            pts = [(0.5, 0.5)] * n
        return place(rng, pts, rotate=rng.random() < 0.5)

    group("ngon_degenerate", "degenerate", 1000, degenerate_ngon)
    return groups


def dump_group(out, name, expect, polys, normals_first):
    verts = []
    faces = []
    for poly in polys:
        start = len(verts)
        verts.extend(poly)
        faces.append(tuple(range(start, start + len(poly))))

    mesh = bpy.data.meshes.new("ddg_tessellation_" + name)
    mesh.from_pydata(verts, [], faces)
    if normals_first:
        # Fill the face-normal cache so that Blender tessellates with precomputed normals,
        # the code path used once normals have been requested.
        _ = [p.normal[:] for p in mesh.polygons]
    mesh.calc_loop_triangles()

    out.write("group %s %s\n" % (name, expect))
    out.write("nv %d\n" % len(mesh.vertices))
    for v in mesh.vertices:
        out.write("%r %r %r\n" % (v.co[0], v.co[1], v.co[2]))
    out.write("np %d\n" % len(mesh.polygons))
    for p in mesh.polygons:
        loops = range(p.loop_start, p.loop_start + p.loop_total)
        out.write("%d %s\n" % (p.loop_total, " ".join(str(mesh.loops[i].vertex_index) for i in loops)))
    out.write("nt %d\n" % len(mesh.loop_triangles))
    for t in mesh.loop_triangles:
        out.write("%d %d %d\n" % (t.loops[0], t.loops[1], t.loops[2]))
    out.write("end\n")
    bpy.data.meshes.remove(mesh)


def main():
    argv = sys.argv[sys.argv.index("--") + 1:] if "--" in sys.argv else []
    if not argv:
        raise SystemExit("usage: blender -b --factory-startup --python blender_tessellation_dump.py -- OUT.txt [SEED] [SCALE]")
    path = argv[0]
    seed = int(argv[1]) if len(argv) > 1 else 1
    scale = float(argv[2]) if len(argv) > 2 else 1.0
    rng = random.Random(seed)
    groups = make_groups(rng, scale)
    total = 0
    with open(path, "w") as out:
        out.write("# ddg-core tessellation dump, seed %d, scale %g\n" % (seed, scale))
        out.write("# blender %s\n" % bpy.app.version_string)
        for i, (name, expect, polys) in enumerate(groups):
            # Alternate between the two code paths (with and without cached face normals).
            dump_group(out, name, expect, polys, normals_first=(i % 2 == 1))
            total += len(polys)
        # The same n-gon families once more with the other normal-cache state.
        for i, (name, expect, polys) in enumerate(groups):
            if name.startswith("ngon_"):
                dump_group(out, name + "_alt", expect, polys, normals_first=(i % 2 == 0))
                total += len(polys)
    print("wrote %d polygons in %s" % (total, path))


main()
