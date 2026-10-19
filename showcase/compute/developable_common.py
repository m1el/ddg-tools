"""Helpers of the developability scenes: triangle meshes from results, per-face dihedral angles,
angle defect, and the library calls with the arrays they take. Numpy only."""

from __future__ import annotations

import numpy as np

import ddg_tools as d


def mesh_args(P, F):
    """The four arrays the library takes for a triangle mesh P (n, 3), F (m, 3)."""
    P = np.ascontiguousarray(P, dtype=np.float64)
    F = np.asarray(F, dtype=np.int32).reshape(-1, 3)
    cv = np.ascontiguousarray(F.ravel(), dtype=np.int32)
    st = np.arange(0, len(cv), 3, dtype=np.int32)
    tr = np.arange(len(cv), dtype=np.int32).reshape(-1, 3)
    return P, cv, st, tr


def result_mesh(r):
    P = np.asarray(r.positions, dtype=np.float64).reshape(-1, 3)
    F = np.asarray(r.faces, dtype=np.int32).reshape(-1, 3)
    return P, F


def vertex_energy(P, F, energy="covariance", **kw):
    r = d.developability(*mesh_args(P, F), energy=energy, **kw)
    return np.asarray(r.vertex_energy, dtype=np.float64), r.energy


def face_normals(P, F):
    n = np.cross(P[F[:, 1]] - P[F[:, 0]], P[F[:, 2]] - P[F[:, 0]])
    return n / np.maximum(np.linalg.norm(n, axis=1, keepdims=True), 1e-300)


def dihedral(P, F):
    """Per face corner i the dihedral angle (radians, 0 = flat) at the edge opposite corner i;
    0 at boundary edges. Also the edge list (m, 2) and the angle per edge."""
    nf = len(F)
    N = face_normals(P, F)
    key = {}
    ang_f = np.zeros((nf, 3))
    for f in range(nf):
        for i in range(3):
            a, b = F[f, (i + 1) % 3], F[f, (i + 2) % 3]
            key.setdefault((min(a, b), max(a, b)), []).append((f, i))
    edges, angs = [], []
    for (a, b), lst in key.items():
        if len(lst) == 2:
            (f0, i0), (f1, i1) = lst
            c = np.clip(np.dot(N[f0], N[f1]), -1, 1)
            ang = float(np.arccos(c))
            ang_f[f0, i0] = ang_f[f1, i1] = ang
            edges.append((a, b))
            angs.append(ang)
    return ang_f, np.array(edges, dtype=np.int32), np.array(angs)


def angle_defect(P, F):
    """2 pi - sum of corner angles at interior vertices, 0 at boundary vertices."""
    n = len(P)
    tot = np.zeros(n)
    for i in range(3):
        a = P[F[:, (i + 1) % 3]] - P[F[:, i]]
        b = P[F[:, (i + 2) % 3]] - P[F[:, i]]
        c = np.einsum("ij,ij->i", a, b) / np.maximum(
            np.linalg.norm(a, axis=1) * np.linalg.norm(b, axis=1), 1e-300)
        np.add.at(tot, F[:, i], np.arccos(np.clip(c, -1, 1)))
    cnt = {}
    for f in F:
        for i in range(3):
            a, b = f[(i + 1) % 3], f[(i + 2) % 3]
            k = (min(a, b), max(a, b))
            cnt[k] = cnt.get(k, 0) + 1
    boundary = np.zeros(n, dtype=bool)
    for (a, b), c in cnt.items():
        if c == 1:
            boundary[a] = boundary[b] = True
    out = 2 * np.pi - tot
    out[boundary] = 0.0
    return out, boundary


def seams(P, F, **kw):
    s = d.developable_seams(*mesh_args(P, F), **kw)
    return (np.asarray(s.seam_edges, dtype=np.int32).reshape(-1, 2),
            np.asarray(s.terminals, dtype=np.int32), s)


def scalar_report(rep):
    return {k: v for k, v in rep.items() if isinstance(v, (int, float, str, bool, type(None)))}


def pack_state(prefix, P, F, energy_name="covariance"):
    """Arrays describing one mesh state: positions, faces, vertex energy, angle defect,
    dihedral per face corner and per edge."""
    e, total = vertex_energy(P, F, energy_name)
    defect, boundary = angle_defect(P, F)
    ang_f, edges, ang_e = dihedral(P, F)
    return {
        f"{prefix}_P": P, f"{prefix}_F": F, f"{prefix}_energy": e, f"{prefix}_defect": defect,
        f"{prefix}_boundary": boundary, f"{prefix}_dih": ang_f, f"{prefix}_edges": edges,
        f"{prefix}_edge_dih": ang_e,
    }, total
