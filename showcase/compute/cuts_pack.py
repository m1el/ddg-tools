"""Pack the charts of a cut unwrap (cuts.py) into the unit square: work/cuts_<tag>_packed.npy."""
import sys
import ddg_tools as d
import numpy as np
from common import WORK

tag = sys.argv[1]
D = np.load(WORK / f"cuts_{tag}.npz")
uv, rep = d.pack_charts(np.ascontiguousarray(D["u_uvs"].astype(np.float32)), D["u_polygon_starts"],
                        D["u_charts"], margin=0.01)
uv = np.asarray(uv, dtype=np.float64).reshape(-1, 2)
np.save(WORK / f"cuts_{tag}_packed.npy", uv)
print(tag, "charts", len(set(D["u_charts"].tolist())), "uv range", uv.min(0), uv.max(0), rep)
