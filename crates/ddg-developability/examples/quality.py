# /// script
# dependencies = ["numpy", "scipy", "matplotlib"]
# ///
"""Quality measures of a developability result against its input (QUALITY.md).

    uv run examples/quality.py INPUT.obj RESULT.obj [--png OUT.png] [--view az,el] [--tag NAME]

Prints one `QUALITY key=value ...` line:

* `e30`: interior edges with a dihedral angle (angle between the face normals) above 30 degrees;
  `folds`: above 170 degrees (a face folded back onto its neighbour); `e30_frac` of all edges.
* seams: the edges above `--seam-deg` (default 20) as a graph; `seam_comp` connected components,
  `seam_mean_len` their mean length (sum of edge lengths) in units of the bounding-box diagonal,
  `seam_short` components of at most 2 edges (short creases), `seam_long_frac` the fraction of the
  seam length in components longer than 5 % of the diagonal (clean curves).
* `d_mean`, `d_haus`: mean and max distance result -> input and the symmetric max (input <-> result),
  in units of the input's bounding-box diagonal, from area-uniform samples (approximate).
* `vol_ratio`: enclosed volume of the result / input (closed meshes; squashing shows here).
* `min_angle` (degrees), `faces_lt1`, `faces_lt5`: faces with a corner angle below 1 / 5 degrees.
* `sum_abs_defect`: sum of |angle defect| over interior vertices.

`--png` renders input and result side by side (flat shading, orthographic, painter's order).
"""

from __future__ import annotations

import argparse
import sys

import numpy as np
from scipy.spatial import cKDTree


def read_obj(path):
    V, F = [], []
    with open(path) as fh:
        for line in fh:
            if line.startswith("v "):
                V.append([float(x) for x in line.split()[1:4]])
            elif line.startswith("f "):
                idx = [int(t.split("/")[0]) - 1 for t in line.split()[1:]]
                for k in range(1, len(idx) - 1):
                    F.append([idx[0], idx[k], idx[k + 1]])
    return np.array(V, dtype=float), np.array(F, dtype=np.int64)


def face_normals(P, F):
    n = np.cross(P[F[:, 1]] - P[F[:, 0]], P[F[:, 2]] - P[F[:, 0]])
    a = np.linalg.norm(n, axis=1)
    return n / np.maximum(a, 1e-300)[:, None], 0.5 * a


def edges_of(F):
    """Unique edges (m, 2), and for each the (up to) two faces (m, 2), -1 when missing."""
    e = np.concatenate([F[:, [0, 1]], F[:, [1, 2]], F[:, [2, 0]]])
    fid = np.tile(np.arange(len(F)), 3)
    k = np.sort(e, axis=1)
    order = np.lexsort((k[:, 1], k[:, 0]))
    k, fid = k[order], fid[order]
    new = np.ones(len(k), dtype=bool)
    new[1:] = np.any(k[1:] != k[:-1], axis=1)
    starts = np.flatnonzero(new)
    E = k[starts]
    ef = -np.ones((len(E), 2), dtype=np.int64)
    ef[:, 0] = fid[starts]
    second = np.flatnonzero(~new)
    grp = np.cumsum(new) - 1
    ef[grp[second], 1] = fid[second]
    return E, ef


def corner_angles(P, F):
    out = np.zeros(F.shape)
    for i in range(3):
        a = P[F[:, (i + 1) % 3]] - P[F[:, i]]
        b = P[F[:, (i + 2) % 3]] - P[F[:, i]]
        out[:, i] = np.arctan2(np.linalg.norm(np.cross(a, b), axis=1), np.einsum("ij,ij->i", a, b))
    return out


def sample(P, F, n, rng):
    _, A = face_normals(P, F)
    f = rng.choice(len(F), size=n, p=A / A.sum())
    r1, r2 = rng.random(n), rng.random(n)
    s = np.sqrt(r1)
    a, b, c = P[F[f, 0]], P[F[f, 1]], P[F[f, 2]]
    return (1 - s)[:, None] * a + (s * (1 - r2))[:, None] * b + (s * r2)[:, None] * c


def volume(P, F):
    return float(np.einsum("ij,ij->i", P[F[:, 0]], np.cross(P[F[:, 1]], P[F[:, 2]])).sum() / 6.0)


def components(n, E):
    parent = np.arange(n)

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    for a, b in E:
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[ra] = rb
    return np.array([find(a) for a, _ in E]) if len(E) else np.zeros(0, dtype=int)


def measures(P0, F0, P, F, seam_deg=20.0, samples=60000, seed=1):
    rng = np.random.default_rng(seed)
    diag = float(np.linalg.norm(P0.max(0) - P0.min(0)))
    N, _ = face_normals(P, F)
    E, ef = edges_of(F)
    inner = ef[:, 1] >= 0
    dih = np.zeros(len(E))
    dih[inner] = np.degrees(np.arccos(np.clip(np.einsum("ij,ij->i", N[ef[inner, 0]], N[ef[inner, 1]]), -1, 1)))
    L = np.linalg.norm(P[E[:, 0]] - P[E[:, 1]], axis=1)
    m = {}
    m["faces"] = len(F)
    m["edges"] = len(E)
    m["e30"] = int((dih > 30).sum())
    m["e30_frac"] = round(float((dih > 30).mean()), 4)
    m["folds"] = int((dih > 170).sum())
    sel = dih > seam_deg
    SE, SL = E[sel], L[sel]
    comp = components(len(P), SE)
    if len(SE):
        _, inv, cnt = np.unique(comp, return_inverse=True, return_counts=True)
        clen = np.bincount(inv, weights=SL) / diag
        m["seam_comp"] = len(cnt)
        m["seam_mean_len"] = round(float(clen.mean()), 4)
        m["seam_short"] = int((cnt <= 2).sum())
        m["seam_long_frac"] = round(float(clen[clen > 0.05].sum() / clen.sum()), 3)
        m["seam_total_len"] = round(float(clen.sum()), 3)
    else:
        m.update(seam_comp=0, seam_mean_len=0.0, seam_short=0, seam_long_frac=0.0, seam_total_len=0.0)
    S = sample(P, F, samples, rng)
    S0 = sample(P0, F0, samples, rng)
    d_out = cKDTree(S0).query(S)[0] / diag
    d_in = cKDTree(S).query(S0)[0] / diag
    m["d_mean"] = round(float(d_out.mean()), 5)
    m["d_haus"] = round(float(max(d_out.max(), d_in.max())), 4)
    v0 = volume(P0, F0)
    m["vol_ratio"] = round(volume(P, F) / v0, 4) if abs(v0) > 0 and inner.all() else float("nan")
    ang = np.degrees(corner_angles(P, F))
    amin = ang.min(axis=1)
    m["min_angle"] = round(float(amin.min()), 5)
    m["faces_lt1"] = int((amin < 1).sum())
    m["faces_lt5"] = int((amin < 5).sum())
    tot = np.zeros(len(P))
    for i in range(3):
        np.add.at(tot, F[:, i], np.radians(ang[:, i]))
    bnd = np.zeros(len(P), dtype=bool)
    bnd[E[~inner].ravel()] = True
    used = np.zeros(len(P), dtype=bool)
    used[F.ravel()] = True
    m["sum_abs_defect"] = round(float(np.abs(2 * np.pi - tot)[used & ~bnd].sum()), 3)
    return m, dih, E


def render(ax, P, F, view, up_axis="y", title="", crease_deg=None):
    from matplotlib.collections import LineCollection, PolyCollection

    if up_axis == "y":  # rotate y-up models (Crane's repository) to z-up
        P = P[:, [2, 0, 1]]

    az, el = np.radians(view[0]), np.radians(view[1])
    fwd = np.array([np.cos(el) * np.cos(az), np.cos(el) * np.sin(az), np.sin(el)])  # towards the viewer
    up0 = np.array([0.0, 0.0, 1.0]) if abs(fwd[2]) < 0.95 else np.array([0.0, 1.0, 0.0])
    right = np.cross(up0, fwd)
    right /= np.linalg.norm(right)
    up = np.cross(fwd, right)
    X = np.stack([P @ right, P @ up], axis=1)
    depth = P @ fwd
    N, _ = face_normals(P, F)
    light = fwd + 0.6 * up + 0.3 * right
    light /= np.linalg.norm(light)
    facing = N @ fwd
    shade = np.abs(N @ light)
    shade = 0.15 + 0.85 * shade
    fd = depth[F].mean(axis=1)
    order = np.argsort(fd)
    cols = np.stack([shade * 0.92, shade * 0.93, shade * 0.97], axis=1)
    cols[facing < 0] *= np.array([1.0, 0.75, 0.75])  # back faces tinted (folds)
    pc = PolyCollection(X[F[order]], facecolors=np.clip(cols[order], 0, 1), edgecolors="none", antialiased=False)
    ax.add_collection(pc)
    if crease_deg is not None:
        E, ef = edges_of(F)
        inner = ef[:, 1] >= 0
        c = np.ones(len(E))
        c[inner] = np.einsum("ij,ij->i", N[ef[inner, 0]], N[ef[inner, 1]])
        sel = c < np.cos(np.radians(crease_deg))
        # front-most: drawn on top (no occlusion test; a sketch of where the creases are)
        vis = sel & ((facing[ef[:, 0]] > 0) | (inner & (facing[np.maximum(ef[:, 1], 0)] > 0)))
        ax.add_collection(LineCollection(X[E[vis]], colors="crimson", linewidths=0.5))
    ax.set_xlim(X[:, 0].min(), X[:, 0].max())
    ax.set_ylim(X[:, 1].min(), X[:, 1].max())
    ax.set_aspect("equal")
    ax.axis("off")
    ax.set_title(title, fontsize=9)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("input")
    ap.add_argument("result")
    ap.add_argument("--png")
    ap.add_argument("--view", default="-50,20", help="azimuth,elevation in degrees")
    ap.add_argument("--tag", default="")
    ap.add_argument("--up", default="y", choices=["y", "z"])
    ap.add_argument("--seam-deg", type=float, default=20.0)
    ap.add_argument("--crease-deg", type=float, help="draw the result's edges above this dihedral angle")
    ap.add_argument("--reference", help="mesh the distances are measured to (default: INPUT)")
    a = ap.parse_args()
    P0, F0 = read_obj(a.input)
    P, F = read_obj(a.result)
    Pr, Fr = read_obj(a.reference) if a.reference else (P0, F0)
    m, dih, E = measures(Pr, Fr, P, F, seam_deg=a.seam_deg)
    if a.png:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        view = [float(x) for x in a.view.split(",")]
        fig, axs = plt.subplots(1, 2, figsize=(12, 6), dpi=110)
        render(axs[0], P0, F0, view, a.up, title=f"input {len(F0)} faces")
        render(axs[1], P, F, view, a.up, crease_deg=a.crease_deg, title=f"{a.tag}  {len(F)} faces  e30 {m['e30']}  folds {m['folds']}  "
                                          f"seams {m['seam_comp']} (short {m['seam_short']})  vol {m['vol_ratio']}")
        fig.tight_layout()
        fig.savefig(a.png)
    print("QUALITY " + " ".join(f"{k}={v}" for k, v in ([("tag", a.tag)] + list(m.items()))))


if __name__ == "__main__":
    sys.exit(main())
