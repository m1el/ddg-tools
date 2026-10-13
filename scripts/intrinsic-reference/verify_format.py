#!/usr/bin/env python3
"""Check the claims of data/intrinsic-reference/FORMAT.md against the files (standard library only).

usage: verify_format.py ROOT      ROOT = data/intrinsic-reference; every ROOT/<input>/<op>/ is checked.
For each run it checks:
  * headers of the five files and their sizes;
  * vertexPositions rows 0..nV_in-1 equal the input vertices; faceInds are in range, every vertex is used;
  * laplace.spmat equals the cotan Laplacian assembled here from faceInds/faceLengths, L_ij = -1/2 sum cot(opposite
    angles), L_ii = -sum_{j != i} L_ij (a self-edge contributes to the diagonal only), and is symmetric;
  * interpolate.spmat: rows sum to 1, columns < nV_in, identity on the input vertices, and
    sum_j W_ij * input_position_j reproduces vertexPositions row i (on a flat or curved input: it does for the
    inserted vertices only up to the position written by the demo, which is the point on the input face);
  * the length of a halfedge i->j equals the length of its twin j->i (where the pair is unambiguous).
Prints one line per run with the largest deviations.
"""
import math
import sys
from collections import defaultdict
from pathlib import Path


def dense(p):
    with open(p) as f:
        h = f.readline().split()
        assert h[:2] == ["#", "dense"], (p, h)
        rows = [line.split() for line in f if line.strip()]
    assert len(rows) == int(h[2]) and all(len(r) == int(h[3]) for r in rows), p
    return rows, int(h[2]), int(h[3])


def sparse(p):
    ent = []
    with open(p) as f:
        h = f.readline().split()
        assert h[:2] == ["#", "sparse"], (p, h)
        for line in f:
            if line.strip():
                i, j, v = line.split()
                ent.append((int(i) - 1, int(j) - 1, float(v)))
    return ent, int(h[2]), int(h[3])


def cot(a, b, c):
    """cotangent of the angle opposite side a (sides a, b, c): (b^2 + c^2 - a^2) / (4 * area)."""
    s = (a + b + c) / 2
    area = math.sqrt(max(s * (s - a) * (s - b) * (s - c), 0.0))
    return (b * b + c * c - a * a) / (4 * area)


def check(d):
    inp = [tuple(float(x) for x in line.split()[1:4]) for line in open(d.parent / "input.obj") if line.startswith("v ")]
    fi, nF, _ = dense(d / "faceInds.dmat")
    fl, nF2, _ = dense(d / "faceLengths.dmat")
    vp, nV, c3 = dense(d / "vertexPositions.dmat")
    assert nF == nF2 and c3 == 3
    inds = [[int(x) for x in r] for r in fi]
    lens = [[float(x) for x in r] for r in fl]
    pos = [[float(x) for x in r] for r in vp]
    nin = len(inp)
    assert min(min(r) for r in inds) == 0 and max(max(r) for r in inds) == nV - 1
    assert len({v for r in inds for v in r}) == nV
    e_pos = max(abs(pos[i][k] - inp[i][k]) for i in range(nin) for k in range(3))

    L = defaultdict(float)
    for t, l in zip(inds, lens):
        for k in range(3):
            # corner at t[k] is opposite the edge t[k+1] -> t[k+2], whose length is l[k+1]
            c = cot(l[(k + 1) % 3], l[k], l[(k + 2) % 3])
            i, j = t[(k + 1) % 3], t[(k + 2) % 3]
            L[(i, j)] -= c / 2
            L[(j, i)] -= c / 2
            L[(i, i)] += c / 2
            L[(j, j)] += c / 2
    # a self-edge i->i adds -c/2 twice to (i,i) in the loop above and +c/2 twice: net zero. Correct (L_ii = -sum_{j != i}).
    ent, r, c = sparse(d / "laplace.spmat")
    assert r == c == nV
    got = defaultdict(float)
    for i, j, v in ent:
        got[(i, j)] += v
    e_lap = max(abs(got.get(k, 0.0) - L.get(k, 0.0)) for k in set(got) | set(L))
    e_sym = max(abs(got.get((j, i), 0.0) - v) for (i, j), v in got.items())

    ent, r, c = sparse(d / "interpolate.spmat")
    assert r == c == nV
    rows = defaultdict(dict)
    for i, j, v in ent:
        rows[i][j] = v
    e_sum = max(abs(sum(rows[i].values()) - 1) for i in range(nV))
    assert all(j < nin for i in rows for j in rows[i])
    assert all(rows[i] == {i: 1.0} for i in range(nin)), "not identity on input vertices"
    e_ip = 0.0
    for i in range(nin, nV):
        q = [sum(w * inp[j][k] for j, w in rows[i].items()) for k in range(3)]
        e_ip = max(e_ip, max(abs(q[k] - pos[i][k]) for k in range(3)))
    maxnnz = max(len(rows[i]) for i in range(nV))

    # twin lengths
    he = defaultdict(list)
    for t, l in zip(inds, lens):
        for k in range(3):
            he[(t[k], t[(k + 1) % 3])].append(l[k])
    e_tw = 0.0
    for (i, j), ls in he.items():
        if i != j and len(ls) == 1 and len(he.get((j, i), [])) == 1:
            e_tw = max(e_tw, abs(ls[0] - he[(j, i)][0]))
    print("%-34s nV_in %6d nV %6d  |pos-input| %.1e  |L-L_own| %.1e  asym %.1e  |W 1-1| %.1e  |W X - pos| %.1e  "
          "max nnz/row %d  twin-len %.1e" % (d.parent.name + "/" + d.name, nin, nV, e_pos, e_lap, e_sym, e_sum, e_ip,
                                             maxnnz, e_tw))


if __name__ == "__main__":
    for d in sorted(Path(sys.argv[1]).glob("*/*/")):
        if (d / "faceInds.dmat").exists():
            check(d)
