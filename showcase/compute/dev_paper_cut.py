"""Cut and flatten the result of a dev_paper.py run (paper Sec. 4.4 with the library's automatic chart
splitting): the seams of `developable_seams` (epsilon 1e-3, the paper's value), then the general
unwrap with these edges as seams, once without and once with `split_overlaps`.

  dev_paper_cut.py TAG [epsilon]     -> work/dp_cut_<TAG>_nosplit.npz/json and work/dp_cut_<TAG>_split.npz/json
Arrays: P, F, face_chart, uv (3 F, 2) per triangle corner, seam_pts/seam_sizes (the cut of the seams),
split_pts/split_sizes (the cuts added by the splitting), face_overlap (triangle covered more than once).
"""
import json
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import ddg_tools as d  # noqa: E402
from common import WORK, chain_edges  # noqa: E402
from developable_common import mesh_args  # noqa: E402


def chains(P, E):
    if not len(E):
        return np.zeros((0, 3)), np.zeros(0, dtype=np.int32)
    cs = [P[np.asarray(c)] for c in chain_edges(np.asarray(E)) if len(c) >= 2]
    return np.concatenate(cs), np.array([len(c) for c in cs], dtype=np.int32)


def overlap_faces(uv, res=4096):
    """Per triangle: does one of its sample points lie where another triangle also is (coverage grid)."""
    T = uv.reshape(-1, 3, 2)
    lo, hi = T.reshape(-1, 2).min(0), T.reshape(-1, 2).max(0)
    s = (res - 2) / float((hi - lo).max())
    X = (T - lo) * s + 0.5
    X = X.mean(1, keepdims=True) + 0.6 * (X - X.mean(1, keepdims=True))   # shrunk: neighbours do not touch
    cnt = np.zeros((res, res), dtype=np.int16)
    inside = []
    for x in X:
        a, b, c = x
        x0, x1 = int(np.floor(x[:, 0].min())), int(np.ceil(x[:, 0].max()))
        y0, y1 = int(np.floor(x[:, 1].min())), int(np.ceil(x[:, 1].max()))
        gx, gy = np.meshgrid(np.arange(x0, x1 + 1), np.arange(y0, y1 + 1))
        p = np.stack([gx.ravel() + 0.5, gy.ravel() + 0.5], 1)
        den = (b[1] - c[1]) * (a[0] - c[0]) + (c[0] - b[0]) * (a[1] - c[1])
        if abs(den) < 1e-12:
            inside.append((None, None))
            continue
        l1 = ((b[1] - c[1]) * (p[:, 0] - c[0]) + (c[0] - b[0]) * (p[:, 1] - c[1])) / den
        l2 = ((c[1] - a[1]) * (p[:, 0] - c[0]) + (a[0] - c[0]) * (p[:, 1] - c[1])) / den
        m = (l1 >= 0) & (l2 >= 0) & (l1 + l2 <= 1)
        ix, iy = gx.ravel()[m], gy.ravel()[m]
        cnt[iy, ix] += 1
        inside.append((ix, iy))
    flag = np.zeros(len(X), dtype=np.int8)
    for f, (ix, iy) in enumerate(inside):
        if ix is not None and len(ix) and (cnt[iy, ix] > 1).mean() > 0.3:
            flag[f] = 1
    return flag


def shelf_pack(uv, face_chart, gap=0.02):
    """Charts come out of the unwrap in their own frames (overlapping each other): translate them into
    rows, largest first, one common scale."""
    T = uv.reshape(-1, 3, 2)
    ids = np.unique(face_chart)
    boxes = {c: (T[face_chart == c].reshape(-1, 2).min(0), T[face_chart == c].reshape(-1, 2).max(0)) for c in ids}
    size = {c: boxes[c][1] - boxes[c][0] for c in ids}
    order = sorted(ids, key=lambda c: -size[c][1])
    area = sum(size[c][0] * size[c][1] for c in ids)
    row_w = max(np.sqrt(area * 1.6), max(size[c][0] for c in ids))
    x = y = rowh = 0.0
    off = {}
    for c in order:
        w, h = size[c]
        if x + w > row_w and x > 0:
            x, y, rowh = 0.0, y + rowh + gap * row_w, 0.0
        off[c] = np.array([x, y]) - boxes[c][0]
        x += w + gap * row_w
        rowh = max(rowh, h)
    out = T.copy()
    for c in ids:
        out[face_chart == c] += off[c]
    return out.reshape(-1, 2)


def main(tag, eps=1e-3):
    D = np.load(WORK / f"dp_{tag}.npz")
    P, F = D["out_P"], D["out_F"]
    a = mesh_args(P, F)
    s = d.developable_seams(*a, epsilon=eps)
    se = np.ascontiguousarray(np.asarray(s.seam_edges, dtype=np.int32).reshape(-1, 2))
    seam_set = {(min(x, y), max(x, y)) for x, y in se}
    for split in (False, True):
        t0 = time.time()
        u = d.unwrap(*a, seam_edges=se.ravel(), split_overlaps=split)
        secs = time.time() - t0
        r = u.report
        ce = np.asarray(u.cut_edges, dtype=np.int32).reshape(-1, 2)
        extra = np.array([e for e in ce if (min(e), max(e)) not in seam_set], dtype=np.int32).reshape(-1, 2)
        uv = np.asarray(u.uvs, dtype=np.float64).reshape(-1, 2)
        fch = np.asarray(u.polygon_charts, dtype=np.int32)
        uv = shelf_pack(uv, fch)
        flag = overlap_faces(uv)
        sp, ss = chains(P, se)
        xp, xs = chains(P, extra)
        name = f"{tag}_{'split' if split else 'nosplit'}"
        np.savez_compressed(WORK / f"dp_cut_{name}.npz", P=P, F=F, face_chart=np.asarray(u.polygon_charts, dtype=np.int32),
                            uv=uv, seam_pts=sp, seam_sizes=ss, split_pts=xp, split_sizes=xs, face_overlap=flag)
        info = dict(tag=name, epsilon=eps, seconds=secs, n_triangles=len(F), n_seam_edges=len(se), n_terminals=len(s.terminals),
                    n_extra_cut_edges=len(extra), n_charts=r["n_charts"], overlapped_fraction=r["overlapped_fraction"],
                    max_chart_overlapped_fraction=r["max_chart_overlapped_fraction"], flipped=r["n_flipped_faces"],
                    area_distortion=r["area_distortion"], qc_max=r["max_quasi_conformal_error"],
                    overlap_faces=int(flag.sum()),
                    **{k: r[k] for k in r if "split" in k})
        (WORK / f"dp_cut_{name}.json").write_text(json.dumps(info, indent=1))
        print(json.dumps(info), flush=True)


if __name__ == "__main__":
    main(sys.argv[1], float(sys.argv[2]) if len(sys.argv) > 2 else 1e-3)
