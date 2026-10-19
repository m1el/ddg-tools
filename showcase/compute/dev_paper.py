"""The developability pictures made the way the paper runs the method (crates/ddg-developability/QUALITY.md):
a uniform coarse input, E_P on it, E_lambda after each 4-1 subdivision. The flow is the library's
`develop` example (the Python module has no coarse-energy option); this script measures a result and
turns it into the arrays the Blender scenes read.

  dev_paper.py analyze TAG INPUT.obj RESULT.obj REFERENCE.obj [LOG]   -> work/dp_TAG.npz, .json
  dev_paper.py states  TAG                                            -> work/dp_snap_TAG.npz, .json
        (the OBJ files of the snapshot runs of dev_paper_runs.sh, in work/dp/snap_TAG_*.obj)

The seam curves are the paper's picture of where the surface is not smooth: the interior edges with a
dihedral angle above SEAM_DEG, grouped into connected components; components shorter than
MIN_LEN_FRAC of the bounding-box diagonal are dropped, the others are chained and smoothed a little.
"""
from __future__ import annotations

import importlib.util
import json
import re
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import ROOT, WORK, chain_edges  # noqa: E402

SEAM_DEG = 20.0
MIN_LEN_FRAC = 0.05

_spec = importlib.util.spec_from_file_location(
    "quality", ROOT / "crates/ddg-developability/examples/quality.py")
Q = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(Q)


def edge_dihedral(P, F):
    """Unique edges, the dihedral angle in degrees on each (0 on boundary edges)."""
    N, _ = Q.face_normals(P, F)
    E, ef = Q.edges_of(F)
    inner = ef[:, 1] >= 0
    dih = np.zeros(len(E))
    dih[inner] = np.degrees(np.arccos(np.clip(np.einsum("ij,ij->i", N[ef[inner, 0]], N[ef[inner, 1]]), -1, 1)))
    return E, dih


def corner_dihedral(P, F):
    """Per face corner i the dihedral angle (radians) at the edge opposite corner i."""
    N, _ = Q.face_normals(P, F)
    E, ef = Q.edges_of(F)
    out = np.zeros((len(F), 3))
    key = {(a, b): k for k, (a, b) in enumerate(E)}
    for f in range(len(F)):
        for i in range(3):
            a, b = F[f, (i + 1) % 3], F[f, (i + 2) % 3]
            k = key[(min(a, b), max(a, b))]
            if ef[k, 1] >= 0:
                out[f, i] = np.arccos(np.clip(N[ef[k, 0]] @ N[ef[k, 1]], -1, 1))
    return out


def angle_defect(P, F):
    """2 pi minus the corner angle sum at interior vertices, 0 at boundary vertices."""
    ang = Q.corner_angles(P, F)
    tot = np.zeros(len(P))
    for i in range(3):
        np.add.at(tot, F[:, i], ang[:, i])
    E, ef = Q.edges_of(F)
    bnd = np.zeros(len(P), dtype=bool)
    bnd[E[ef[:, 1] < 0].ravel()] = True
    d = 2 * np.pi - tot
    d[bnd] = 0.0
    return d


def smooth_chain(pts, passes=4):
    pts = pts.copy()
    closed = np.allclose(pts[0], pts[-1])
    for _ in range(passes):
        if closed:
            core = pts[:-1]
            new = 0.25 * np.roll(core, 1, 0) + 0.5 * core + 0.25 * np.roll(core, -1, 0)
            pts = np.concatenate([new, new[:1]])
        else:
            new = pts.copy()
            new[1:-1] = 0.25 * pts[:-2] + 0.5 * pts[1:-1] + 0.25 * pts[2:]
            pts = new
    return pts


def seam_curves(P, F, seam_deg=SEAM_DEG, min_len_frac=MIN_LEN_FRAC, diag=None):
    """Smoothed polylines of the long connected components of the edges above `seam_deg`; also
    the number of components and of the long ones."""
    diag = float(np.linalg.norm(P.max(0) - P.min(0))) if diag is None else diag
    E, dih = edge_dihedral(P, F)
    SE = E[dih > seam_deg]
    if not len(SE):
        return [], 0, 0
    comp = Q.components(len(P), SE)
    L = np.linalg.norm(P[SE[:, 0]] - P[SE[:, 1]], axis=1)
    ids, inv = np.unique(comp, return_inverse=True)
    clen = np.bincount(inv, weights=L) / diag
    keep = clen >= min_len_frac
    sel = keep[inv]
    curves = []
    if sel.any():
        for ch in chain_edges(SE[sel]):
            if len(ch) >= 2:
                curves.append(smooth_chain(P[np.asarray(ch)]))
    return curves, len(ids), int(keep.sum())


def pack_curves(curves):
    if not curves:
        return np.zeros((0, 3)), np.zeros(0, dtype=np.int32)
    return np.concatenate(curves), np.array([len(c) for c in curves], dtype=np.int32)


def parse_log(path):
    rounds, res = [], {}
    if not path or not Path(path).exists():
        return rounds, res
    for line in Path(path).read_text().splitlines():
        m = re.match(r"round (\d+): (\d+) v / (\d+) f, E (\S+) -> (\S+), .*?(\d+) iterations .*?, ([\d.]+)s$", line)
        if m:
            rounds.append(dict(round=int(m[1]), vertices=int(m[2]), faces=int(m[3]), e_start=float(m[4]),
                               e_end=float(m[5]), iterations=int(m[6]), seconds=float(m[7])))
        if line.startswith("RESULT"):
            res = dict(kv.split("=", 1) for kv in line.split()[1:])
    return rounds, res


def analyze(tag, inp, out, ref, log=None):
    P0, F0 = Q.read_obj(inp)
    P, F = Q.read_obj(out)
    Pr, Fr = Q.read_obj(ref)
    m, dih, E = Q.measures(Pr, Fr, P, F, seam_deg=SEAM_DEG)
    curves, n_comp, n_long = seam_curves(P, F)
    sc, ss = pack_curves(curves)
    rounds, res = parse_log(log)
    arrays = dict(in_P=P0, in_F=F0, out_P=P, out_F=F, out_defect=angle_defect(P, F),
                  out_dih=corner_dihedral(P, F), in_dih=corner_dihedral(P0, F0),
                  seam_pts=sc, seam_sizes=ss)
    np.savez_compressed(WORK / f"dp_{tag}.npz", **arrays)
    info = dict(tag=tag, input=str(inp), faces_in=len(F0), faces_out=len(F), rounds=rounds,
                result={k: res.get(k) for k in ("e_initial", "e_final", "seconds", "status", "iterations",
                                                "flips", "collapses", "bad_faces", "min_angle")},
                seam_curves=len(curves), seam_components=n_comp, seam_long_components=n_long,
                seam_deg=SEAM_DEG, min_len_frac=MIN_LEN_FRAC, **m)
    (WORK / f"dp_{tag}.json").write_text(json.dumps(info, indent=1))
    print(json.dumps(info))


def states(tag):
    """Snapshot runs: files work/dp/snap_<tag>_<k>.obj with <k>.json (round, iteration, energy)."""
    d = WORK / "dp"
    metas = sorted(d.glob(f"snap_{tag}_*.json"), key=lambda p: int(p.stem.split("_")[-1]))
    out, meta = {}, []
    for k, mp in enumerate(metas):
        mt = json.loads(mp.read_text())
        P, F = Q.read_obj(d / f"snap_{tag}_{mt['k']}.obj")
        curves, _, _ = seam_curves(P, F)
        sc, ss = pack_curves(curves)
        out[f"s{k}_P"], out[f"s{k}_F"] = P, F
        out[f"s{k}_seam_pts"], out[f"s{k}_seam_sizes"] = sc, ss
        meta.append(dict(mt, n_vertices=len(P), n_triangles=len(F), n_curves=len(curves)))
    out["n_states"] = len(meta)
    np.savez_compressed(WORK / f"dp_snap_{tag}.npz", **out)
    (WORK / f"dp_snap_{tag}.json").write_text(json.dumps(dict(tag=tag, states=meta), indent=1))
    print(tag, len(meta), "states")


if __name__ == "__main__":
    a = sys.argv[1:]
    if a[0] == "analyze":
        analyze(a[1], a[2], a[3], a[4], a[5] if len(a) > 5 else None)
    elif a[0] == "states":
        states(a[1])
