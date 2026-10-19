"""Snapshots of the developability flow: positions, faces, vertex energy, angle defect and
dihedral angles at chosen iterations of every round. The flow is deterministic and independent of
the thread count, so a snapshot is the flow of the one-call run stopped at (round, iteration)
through the progress callback; the snapshots are computed in parallel, one process each.

developable_snap.py <tag> <mesh> [sub=2] [iters=500] [per_round=16] [energy=covariance] [procs=48]
"""
import json
import sys
import time
from multiprocessing import Pool

import ddg_tools as d
import numpy as np
from common import MODELS, WORK, corners, load_obj, model, save, triangles
from developable_common import pack_state, result_mesh

a = sys.argv[1:]
tag, name = a[0], a[1]
kv = dict(x.split("=") for x in a[2:])
sub = int(kv.get("sub", 2))
iters = int(kv.get("iters", 500))
per_round = int(kv.get("per_round", 16))
energy = kv.get("energy", "covariance")
procs = int(kv.get("procs", 48))

path = MODELS / "crumpleddevelopable/CrumpledDevelopable.obj" if name == "crumpled" else model(name)
pos, polys, _ = load_obj(path)
cv, st, tr = corners(polys)
F0 = triangles(polys)

targets = []
for rd in range(sub + 1):
    its = np.unique(np.round(np.logspace(0, np.log10(iters), per_round)).astype(int))
    targets += [(rd, int(i)) for i in its]


def run(target):
    rd, it = target

    def prog(r, i, e, g):
        return not (r == rd and i >= it)

    t0 = time.time()
    r = d.make_developable(pos, cv, st, tr, energy=energy, subdivisions=sub, iterations=iters,
                           return_faces=True, threads=1, progress=prog, progress_interval=0)
    P, F = result_mesh(r)
    arrays, total = pack_state("s", P, F, energy)
    arrays = {k[2:]: v for k, v in arrays.items()}
    return target, arrays, total, dict(cancelled=r.report["cancelled"], seconds=time.time() - t0)


if __name__ == "__main__":
    t0 = time.time()
    with Pool(procs) as pool:
        res = pool.map(run, targets, chunksize=1)
    out = {}
    first, e0 = pack_state("s", pos, F0, energy)
    first = {k[2:]: v for k, v in first.items()}
    states = [((-1, 0), first, e0)] + res
    n_it = 0
    meta = []
    for k, (tg, arrs, tot) in enumerate((s[0], s[1], s[2]) for s in states):
        for key, v in arrs.items():
            out[f"s{k}_{key}"] = v
        meta.append(dict(round=tg[0], iteration=tg[1], energy=tot,
                         cumulative=(0 if tg[0] < 0 else tg[0] * iters + tg[1]),
                         n_vertices=len(arrs["P"]), n_triangles=len(arrs["F"])))
    out["n_states"] = len(states)
    save(f"developable_snap_{tag}", **out)
    info = dict(tag=tag, mesh=name, energy_name=energy, subdivisions=sub, iterations_per_round=iters,
                states=meta, seconds=time.time() - t0)
    (WORK / f"developable_snap_{tag}.json").write_text(json.dumps(info, indent=1))
    print(json.dumps({k: v for k, v in info.items() if k != "states"}))
    print([ (m["round"], m["iteration"], round(m["energy"], 4), m["n_vertices"]) for m in meta])
