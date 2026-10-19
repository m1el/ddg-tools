"""Cut and flatten (paper Sec. 4.4, Figs. 19/20) of a result of developable.py: seams and the
paper's own flattening (developable_unwrap) for several thresholds epsilon.

developable_cut.py <tag> [eps1,eps2,...]
"""
import json
import sys
import time

import ddg_tools as d
import numpy as np
from common import WORK
from developable_common import mesh_args, scalar_report

tag = sys.argv[1]
epss = [float(x) for x in (sys.argv[2] if len(sys.argv) > 2 else "1e-1,1e-2,1e-3").split(",")]
D = np.load(WORK / f"developable_{tag}.npz")
P, F = D["out_P"], D["out_F"]
args = mesh_args(P, F)
out, infos = {}, []
for k, eps in enumerate(epss):
    t0 = time.time()
    u = d.developable_unwrap(*args, epsilon=eps)
    secs = time.time() - t0
    rep = u.report
    out[f"e{k}_seam_edges"] = np.asarray(u.seam_edges, dtype=np.int32).reshape(-1, 2)
    out[f"e{k}_terminals"] = np.asarray(u.terminals, dtype=np.int32)
    out[f"e{k}_uv"] = np.asarray(u.corner_uvs, dtype=np.float64).reshape(-1, 2)
    out[f"e{k}_chart"] = np.asarray(u.polygon_chart, dtype=np.int32)
    out[f"e{k}_logscale"] = np.asarray(u.corner_log_scale, dtype=np.float64)
    out[f"e{k}_defect"] = np.asarray(u.vertex_defect, dtype=np.float64)
    ch = rep["charts"]
    infos.append(dict(epsilon=eps, seconds=secs, n_seam_edges=len(out[f"e{k}_seam_edges"]),
                      n_terminals=len(out[f"e{k}_terminals"]), n_charts=rep["n_charts"],
                      flipped=rep["flipped_triangles"], mean_abs_log_scale=rep["mean_abs_log_scale"],
                      overlap_fraction_max=max(c["overlapped_fraction"] for c in ch),
                      overlap_fraction_mean=float(np.mean([c["overlapped_fraction"] for c in ch])),
                      qc_mean=float(np.mean([c["quasi_conformal_mean"] for c in ch])),
                      qc_max=max(c["quasi_conformal_max"] for c in ch),
                      is_injective=rep["is_injective"], n_triangles=len(F)))
    print(infos[-1], flush=True)
np.savez_compressed(WORK / f"developable_cut_{tag}.npz", **out)
(WORK / f"developable_cut_{tag}.json").write_text(json.dumps(infos, indent=1, default=str))
