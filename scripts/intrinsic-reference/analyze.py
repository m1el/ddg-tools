#!/usr/bin/env python3
"""Derive statistics from the demo's faceInds.dmat / faceLengths.dmat (standard library only).

usage: analyze.py OPDIR      OPDIR holds faceInds.dmat, faceLengths.dmat, stats.txt; ../input.obj is the input mesh.
       analyze.py --table ROOT   print a markdown table of all OPDIR/analysis.json under ROOT
Writes OPDIR/analysis.json.
"""
import json
import math
import sys
from collections import Counter, defaultdict
from pathlib import Path


def read_dmat(p):
    rows = []
    with open(p) as f:
        head = f.readline().split()
        assert head[:2] == ["#", "dense"], head
        for line in f:
            if line.strip():
                rows.append(line.split())
    assert len(rows) == int(head[2])
    return rows


def angle(a, b, c):
    """Angle opposite side a in a triangle with sides a, b, c (half-angle formula, accurate for tiny angles)."""
    s = (a + b + c) / 2
    num, den = (s - b) * (s - c), s * (s - a)
    if num <= 0 or den <= 0:
        return 0.0 if num <= 0 and den > 0 else math.pi
    return 2 * math.atan(math.sqrt(num / den))


def read_obj(p):
    v, f = [], []
    for line in open(p):
        t = line.split()
        if not t:
            continue
        if t[0] == "v":
            v.append(tuple(float(x) for x in t[1:4]))
        elif t[0] == "f":
            idx = [int(x.split("/")[0]) for x in t[1:]]
            idx = [i - 1 if i > 0 else len(v) + i for i in idx]
            for k in range(1, len(idx) - 1):
                f.append((idx[0], idx[k], idx[k + 1]))
    return v, f


def sub(a, b):
    return (a[0] - b[0], a[1] - b[1], a[2] - b[2])


def input_geometry(p):
    v, f = read_obj(p)
    angsum, area = defaultdict(float), 0.0
    for t in f:
        for k in range(3):
            o, a, b = v[t[k]], v[t[(k + 1) % 3]], v[t[(k + 2) % 3]]
            x, y = sub(a, o), sub(b, o)
            cr = (x[1] * y[2] - x[2] * y[1], x[2] * y[0] - x[0] * y[2], x[0] * y[1] - x[1] * y[0])
            n = math.sqrt(sum(c * c for c in cr))
            angsum[t[k]] += math.atan2(n, sum(x[i] * y[i] for i in range(3)))
            if k == 0:
                area += n / 2
    return len(v), len(f), area, angsum


def analyze(opdir):
    opdir = Path(opdir)
    inds = [[int(x) for x in r] for r in read_dmat(opdir / "faceInds.dmat")]
    lens = [[float(x) for x in r] for r in read_dmat(opdir / "faceLengths.dmat")]
    nF = len(inds)
    nV = max(max(r) for r in inds) + 1
    nV_in, nF_in, area_in, angsum_in = input_geometry(opdir.parent / "input.obj")

    # corner angle at face vertex k is opposite the edge k+1 -> k+2, i.e. lens[(k+1)%3]
    amin, amax, area, bad_tri = math.inf, 0.0, 0.0, 0
    angsum = defaultdict(float)
    hist = [0] * 18  # 10-degree bins of all corner angles
    for t, l in zip(inds, lens):
        a_, b_, c_ = l
        if not (a_ < b_ + c_ and b_ < a_ + c_ and c_ < a_ + b_):
            bad_tri += 1
        s = (a_ + b_ + c_) / 2
        area += math.sqrt(max(s * (s - a_) * (s - b_) * (s - c_), 0.0))
        for k in range(3):
            ang = angle(l[(k + 1) % 3], l[k], l[(k + 2) % 3])
            amin, amax = min(amin, ang), max(amax, ang)
            angsum[t[k]] += ang
            hist[min(int(math.degrees(ang) // 10), 17)] += 1

    # connectivity: directed halfedge counts
    he = Counter()
    for t in inds:
        for k in range(3):
            he[(t[k], t[(k + 1) % 3])] += 1
    self_edges = sum(c // 2 for (i, j), c in he.items() if i == j)
    faces_repeated_vertex = sum(1 for t in inds if len(set(t)) < 3)
    parallel = {}  # unordered pair -> number of distinct edges between them
    for (i, j), c in he.items():
        if i < j:
            parallel[(i, j)] = max(c, he.get((j, i), 0))
        elif (j, i) not in he:
            parallel[(j, i)] = c
    double_edges = sum(1 for m in parallel.values() if m > 1)
    # boundary halfedges: those i->j for which fewer j->i exist than i->j
    boundary_he = 0
    for (i, j), c in he.items():
        if i != j and he.get((j, i), 0) < c:
            boundary_he += c - he.get((j, i), 0)
    nE = (3 * nF + boundary_he) // 2

    # Delaunay check on the edges that can be matched unambiguously (directed pair occurs once each way)
    where = {}
    for fi, t in enumerate(inds):
        for k in range(3):
            where[(t[k], t[(k + 1) % 3], fi)] = k
    first = {}
    for fi, t in enumerate(inds):
        for k in range(3):
            first.setdefault((t[k], t[(k + 1) % 3]), []).append((fi, k))
    worst = -math.inf
    n_checked = 0
    for (i, j), lst in first.items():
        if i >= j or len(lst) != 1 or len(first.get((j, i), [])) != 1:
            continue
        (f1, k1), (f2, k2) = lst[0], first[(j, i)][0]
        # opposite corner of the halfedge k in its face is at vertex k+2, angle opposite edge k
        l1, l2 = lens[f1], lens[f2]
        a1 = angle(l1[k1], l1[(k1 + 1) % 3], l1[(k1 + 2) % 3])
        a2 = angle(l2[k2], l2[(k2 + 1) % 3], l2[(k2 + 2) % 3])
        worst = max(worst, a1 + a2 - math.pi)
        n_checked += 1

    # angle sums at the input vertices
    dev = max((abs(angsum[v] - angsum_in[v]) for v in range(nV_in)), default=0.0)
    # angle sums at inserted vertices (2 pi interior, pi on a boundary edge)
    ins = [angsum[v] for v in range(nV_in, nV)]
    stats = {}
    for line in open(opdir / "stats.txt"):
        k, v = line.split()
        stats[k] = float(v)
    out = {
        "nV_input": nV_in, "nF_input": nF_in, "nV": nV, "nF": nF, "nE": nE, "boundary_halfedges": boundary_he,
        "euler_characteristic": nV - nE + nF,
        "flips": int(stats["flips"]), "face_insertions": int(stats["face_insertions"]),
        "edge_splits": int(stats["edge_splits"]), "inserted_vertices_net": nV - nV_in,
        "is_delaunay_demo": int(stats["is_delaunay"]), "op_seconds_in_process": stats["op_seconds"],
        "min_angle_deg": math.degrees(amin), "max_angle_deg": math.degrees(amax),
        "min_angle_deg_demo": stats["min_angle_deg_gc"],
        "self_edges": self_edges, "faces_with_repeated_vertex": faces_repeated_vertex,
        "double_edge_pairs": double_edges, "max_edges_between_a_pair": max(parallel.values()),
        "triangle_inequality_violations": bad_tri,
        "area_input": area_in, "area_intrinsic": area, "area_rel_diff": (area - area_in) / area_in,
        "max_angle_sum_deviation_input_vertices": dev,
        "inserted_angle_sum_min": min(ins) if ins else None, "inserted_angle_sum_max": max(ins) if ins else None,
        "delaunay_excess_max_rad": worst if n_checked else None, "delaunay_edges_checked": n_checked,
        "angle_hist_10deg": hist,
    }
    assert nF == int(stats["nF"]) and nV == int(stats["nV"]), "helper and demo disagree on the size"
    (opdir / "analysis.json").write_text(json.dumps(out, indent=1) + "\n")
    return out


def table(root):
    rows = []
    for p in sorted(Path(root).glob("*/*/analysis.json")):
        a = json.loads(p.read_text())
        t = p.parent / "time.txt"
        wall = t.read_text().split()[0] if t.exists() else "?"
        rows.append((p.parent.parent.name, p.parent.name, a, wall))
    print("| input | op | nV in | nF in | nV out | nF out | flips | inserted | min angle | max angle | self-edges | "
          "double-edge pairs | Delaunay | wall s |")
    print("|" + "---|" * 14)
    for name, op, a, wall in rows:
        print("| %s | %s | %d | %d | %d | %d | %d | %d | %.3f | %.2f | %d | %d | %s | %s |" % (
            name, op, a["nV_input"], a["nF_input"], a["nV"], a["nF"], a["flips"], a["inserted_vertices_net"],
            a["min_angle_deg"], a["max_angle_deg"], a["self_edges"], a["double_edge_pairs"],
            "yes" if a["is_delaunay_demo"] else "NO", wall))


if __name__ == "__main__":
    if sys.argv[1] == "--table":
        table(sys.argv[2])
    else:
        r = analyze(sys.argv[1])
        print(json.dumps({k: v for k, v in r.items() if k != "angle_hist_10deg"}))
