"""Developability flow (Stein, Grinspun, Crane 2018) through ddg_tools.make_developable: the
input and the result with per-vertex energy, angle defect, dihedral angles and the cut.

developable.py <tag> <mesh> [key=value ...]
keys: energy (covariance), sub (2), iters (500), seams_eps (1e-3), threads (1)
mesh: a showcase name (spot_control, bob_control, ...) or crumpled.
"""
import json
import sys
import time

import ddg_tools as d
import numpy as np
from common import MODELS, WORK, corners, load_obj, model, save, triangles
from developable_common import pack_state, result_mesh, scalar_report, seams

a = sys.argv[1:]
tag, name = a[0], a[1]
kv = dict(x.split("=") for x in a[2:])
energy = kv.get("energy", "covariance")
sub = int(kv.get("sub", 2))
iters = int(kv.get("iters", 500))
eps = float(kv.get("seams_eps", 1e-3))
threads = int(kv.get("threads", 1))

path = MODELS / "crumpleddevelopable/CrumpledDevelopable.obj" if name == "crumpled" else model(name)
pos, polys, _ = load_obj(path)
cv, st, tr = corners(polys)
F0 = triangles(polys)

t0 = time.time()
r = d.make_developable(pos, cv, st, tr, energy=energy, subdivisions=sub, iterations=iters,
                       return_faces=True, threads=threads)
secs = time.time() - t0
P, F = result_mesh(r)
rep = r.report
out = dict(mesh_name=name)
st0, e0 = pack_state("in", pos, F0, energy)
st1, e1 = pack_state("out", P, F, energy)
out.update(st0)
out.update(st1)
se, term, sres = seams(P, F, energy=energy, epsilon=eps)
out.update(seam_edges=se, terminals=term)
save(f"developable_{tag}", **out)
info = dict(tag=tag, mesh=name, energy=energy, subdivisions=sub, iterations_per_round=iters,
            n_vertices_in=len(pos), n_triangles_in=len(F0), n_vertices_out=len(P),
            n_triangles_out=len(F), energy_in=e0, energy_out=e1,
            energy_initial=rep["energy_initial"], energy_final=rep["energy_final"],
            iterations=rep["iterations"], seconds=secs, seams_epsilon=eps,
            n_seam_edges=len(se), n_terminals=len(term), bad_faces=rep.get("bad_faces"),
            stopped=rep.get("stopped"), min_angle=rep.get("min_angle"),
            sum_abs_defect_in=float(np.abs(st0[f"in_defect"]).sum()),
            sum_abs_defect_out=float(np.abs(st1["out_defect"]).sum()),
            max_dihedral_deg_out=float(np.degrees(st1["out_edge_dih"].max())),
            edges_above_30deg_out=int((st1["out_edge_dih"] > np.radians(30)).sum()),
            n_edges_out=len(st1["out_edge_dih"]))
(WORK / f"developable_{tag}.json").write_text(json.dumps(info, indent=1, default=str))
print(json.dumps(info, default=str))
