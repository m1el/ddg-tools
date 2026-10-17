"""Variational surface cutting (Sharp, Crane 2018) through ddg_tools.CutSession: record the
state, the cut polylines and the energy at every step, then insert the cut and unwrap.

cuts.py <tag> <mesh> circle <vertex> <radius> <steps> [wavelength] [energy] [penalty-camera]
cuts.py <tag> <mesh> clusters <n_regions> <iterations> <steps> [wavelength] [energy]

penalty-camera: "x,y,z,strength" (mesh coordinates): the visibility penalty of paper Sec. 7.1
from that camera position, scaled by strength.
"""
import json
import sys
import time

import ddg_tools as d
import numpy as np
from common import WORK, corners, flat_polys, load_obj, model, save

a = sys.argv[1:]
tag, name, start = a[0], a[1], a[2]
pos, polys, _ = load_obj(model(name))
cv, st, tr = corners(polys)

t0 = time.time()
if start == "circle":
    s = d.CutSession(pos, cv, st, tr, circle_vertex=int(a[3]), circle_radius=float(a[4]))
else:
    s = d.CutSession(pos, cv, st, tr, n_regions=int(a[3]), cluster_iterations=int(a[4]))
t_setup = time.time() - t0
steps = int(a[5])
wl = float(a[6]) if len(a) > 6 else 0.25
energy = a[7] if len(a) > 7 else "hencky"
controls = dict(energy=energy)
if len(a) > 8 and a[8] != "-":
    cx, cy, cz, strength = (float(x) for x in a[8].split(","))
    eta = np.asarray(d.visibility_penalty(pos, cv, st, tr, [cx, cy, cz]))
    controls["penalty"] = np.ascontiguousarray(eta * strength)
omega = s.roughening_rate(**controls)
alpha = d.length_weight_for_ripples(omega, wl)
controls["length_weight"] = alpha

N = s.n_regions


def record():
    p, starts, closed, regs = s.preview()
    return (np.asarray(p, dtype=np.float64).reshape(-1, 3), np.asarray(starts, dtype=np.int64),
            np.asarray(closed, dtype=bool), np.asarray(s.state, dtype=np.float64).copy())


recs = [record()]
energies, lengths, taus = [], [], []
t1 = time.time()
for i in range(steps):
    r = s.step(**controls)
    energies.append(r["energy_total"])
    lengths.append(r["cut_length"])
    taus.append(r["tau"])
    recs.append(record())
t_flow = time.time() - t1

# Ragged storage: points of all records, polyline starts per record.
pts = np.concatenate([x[0] for x in recs])
pt_off = np.cumsum([0] + [len(x[0]) for x in recs])
st_all = np.concatenate([x[1] for x in recs])
st_off = np.cumsum([0] + [len(x[1]) for x in recs])
closed_all = np.concatenate([x[2] for x in recs])
cl_off = np.cumsum([0] + [len(x[2]) for x in recs])
states = np.stack([x[3] for x in recs]).astype(np.float32)  # (n_rec, N, n_corners)

out = dict(positions=pos, corner_vertices=cv, polygon_starts=st, triangles=tr,
           pts=pts, pt_off=pt_off, starts=st_all, st_off=st_off, closed=closed_all,
           cl_off=cl_off, states=states, energy=np.array(energies), cut_length=np.array(lengths),
           tau=np.array(taus))
pv, ps = flat_polys(polys)
out.update(poly_verts=pv, poly_sizes=ps)

t2 = time.time()
try:
    u = s.unwrap()
    ins = u.inserted
    out.update(
        u_positions=np.asarray(ins.positions).reshape(-1, 3),
        u_corner_vertices=np.asarray(ins.corner_vertices),
        u_polygon_starts=np.asarray(ins.polygon_starts),
        u_triangles=np.asarray(ins.triangles).reshape(-1, 3),
        u_uvs=np.asarray(u.uvs).reshape(-1, 2),
        u_charts=np.asarray(u.polygon_charts),
        u_seams=np.asarray(u.seams).reshape(-1, 2),
        u_polygon_region=np.asarray(ins.polygon_region),
        u_polygon_parent=np.asarray(ins.polygon_parent),
    )
    # Region of every new polygon: where the insertion does not say, the region of its parent
    # polygon (argmin of the mean region functions over the parent's corners).
    reg = np.asarray(ins.polygon_region).astype(np.int64)
    par = np.asarray(ins.polygon_parent).astype(np.int64)
    fin = np.asarray(s.state, dtype=np.float64)
    sums = np.zeros((len(st), N))
    owner = np.repeat(np.arange(len(st)), np.diff(np.concatenate([st, [len(cv)]])))
    for k in range(N):
        sums[:, k] = np.bincount(owner, weights=fin[k], minlength=len(st))
    parent_region = np.argmin(sums, axis=1)
    out["u_region"] = np.where(reg >= 0, reg, parent_region[par])
    rep = u.report
    unwrap_info = {k: rep[k] for k in rep if isinstance(rep[k], (int, float, str, bool))}
except Exception as e:  # report, keep the flow data
    unwrap_info = {"error": str(e)}
t_unwrap = time.time() - t2

save(f"cuts_{tag}", **out)
info = dict(tag=tag, mesh=name, start=start, n_regions=N, steps=steps, wavelength=wl,
            energy=energy, omega=omega, alpha=alpha, energy_first=energies[0],
            energy_last=energies[-1], cut_length_first=lengths[0], cut_length_last=lengths[-1],
            seconds_setup=t_setup, seconds_flow=t_flow, seconds_unwrap=t_unwrap,
            n_vertices=len(pos), n_polygons=len(polys), unwrap=unwrap_info)
(WORK / f"cuts_{tag}.json").write_text(json.dumps(info, indent=1, default=str))
print(json.dumps({k: v for k, v in info.items() if k != "unwrap"}, default=str))
print({k: unwrap_info.get(k) for k in list(unwrap_info)[:40]})
