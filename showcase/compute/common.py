"""Shared helpers of the compute scripts: OBJ input, polygon -> corner arrays, output paths."""

from __future__ import annotations

import os
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
SHOW = ROOT / "showcase"
WORK = SHOW / "work"
DATA = ROOT / "data"
DATASET = DATA / "optimal-cones" / "optimal-cones-data" / "data"
MODELS = DATA / "models"

WORK.mkdir(parents=True, exist_ok=True)


def model(name: str) -> Path:
    """Path of a mesh by short name: 'spot', 'bob', 'blub', or a dataset mesh like 'bunny'."""
    table = {
        "spot": MODELS / "spot/spot/spot_triangulated.obj",
        "spot_quad": MODELS / "spot/spot/spot_quadrangulated.obj",
        "spot_control": MODELS / "spot/spot/spot_control_mesh.obj",
        "bob": MODELS / "bob/bob/bob_tri.obj",
        "bob_quad": MODELS / "bob/bob/bob_quad.obj",
        "bob_control": MODELS / "bob/bob/bob_controlmesh.obj",
        "blub": MODELS / "blub/blub/blub_triangulated.obj",
        "blub_quad": MODELS / "blub/blub/blub_quadrangulated.obj",
        "blub_control": MODELS / "blub/blub/blub_control_mesh.obj",
    }
    return table.get(name, DATASET / f"{name}.obj")


def load_obj(path: Path | str):
    """Positions (n, 3) and polygons (list of vertex index lists). Texture coordinates are
    returned too when present: (uv (m, 2), polygon uv indices) or None."""
    pos, uv, polys, polys_uv = [], [], [], []
    with open(path) as f:
        for line in f:
            if line.startswith("v "):
                pos.append([float(x) for x in line.split()[1:4]])
            elif line.startswith("vt "):
                uv.append([float(x) for x in line.split()[1:3]])
            elif line.startswith("f "):
                p, q = [], []
                for tok in line.split()[1:]:
                    parts = tok.split("/")
                    p.append(int(parts[0]) - 1)
                    if len(parts) > 1 and parts[1]:
                        q.append(int(parts[1]) - 1)
                polys.append(p)
                polys_uv.append(q)
    pos = np.array(pos, dtype=np.float64)
    tex = None
    if uv and all(len(q) == len(p) for p, q in zip(polys, polys_uv)):
        tex = (np.array(uv, dtype=np.float64), polys_uv)
    return pos, polys, tex


def corners(polys):
    """The arrays `unwrap`/`CutSession` take: corner_vertices, polygon_starts, and a fan
    triangulation in corner indices."""
    cv, starts, tris = [], [], []
    for p in polys:
        s = len(cv)
        starts.append(s)
        cv.extend(p)
        for k in range(1, len(p) - 1):
            tris.append([s, s + k, s + k + 1])
    return (
        np.array(cv, dtype=np.int32),
        np.array(starts, dtype=np.int32),
        np.array(tris, dtype=np.int32),
    )


def triangles(polys):
    """Vertex-index triangles (fan triangulation of every polygon)."""
    t = []
    for p in polys:
        for k in range(1, len(p) - 1):
            t.append([p[0], p[k], p[k + 1]])
    return np.array(t, dtype=np.int32)


def flat_polys(polys):
    """Polygons as (flat vertex list, sizes) for npz storage."""
    return np.array([v for p in polys for v in p], dtype=np.int32), np.array(
        [len(p) for p in polys], dtype=np.int32
    )


def save(name: str, **arrays):
    out = WORK / f"{name}.npz"
    np.savez_compressed(out, **{k: np.asarray(v) for k, v in arrays.items()})
    print("wrote", out, file=sys.stderr)
    return out


def edges_of(polys):
    """Unique undirected polygon edges (m, 2)."""
    s = set()
    for p in polys:
        for a, b in zip(p, p[1:] + p[:1]):
            s.add((min(a, b), max(a, b)))
    return np.array(sorted(s), dtype=np.int32)


def chain_edges(edges):
    """Join undirected edges (m, 2) into maximal vertex chains (lists of vertex ids)."""
    from collections import defaultdict

    adj = defaultdict(list)
    for i, (a, b) in enumerate(edges):
        adj[a].append((b, i))
        adj[b].append((a, i))
    used = np.zeros(len(edges), dtype=bool)
    chains = []
    # Start at vertices whose degree is not 2 first, then close the remaining loops.
    starts = [v for v in adj if len(adj[v]) != 2] + list(adj)
    for v0 in starts:
        for w, i in adj[v0]:
            if used[i]:
                continue
            used[i] = True
            chain = [v0, w]
            cur = w
            while len(adj[cur]) == 2:
                nxt = [(x, j) for x, j in adj[cur] if not used[j]]
                if not nxt:
                    break
                x, j = nxt[0]
                used[j] = True
                chain.append(x)
                cur = x
            chains.append(chain)
    return chains


def env_int(name, default):
    return int(os.environ.get(name, default))
