# /// script
# dependencies = ["pillow", "numpy", "scikit-image", "scipy"]
# ///
"""Trace the Table 1 diagrams of Ishii-Kishimoto-Moriuchi-Suzuki 2012 (a 500 dpi render of page 3, not in the repo)
into spatial trivalent graphs: skeleton -> strand pieces -> gap pairing (under-crossings) -> edges -> lifted 3D polylines.
Usage: uv run ishii_extract.py PAGE.png OUTDIR [NAME...]   (writes OUTDIR/<name>.json and overlay PNGs)."""
import sys, json, math
import numpy as np
from PIL import Image, ImageDraw
from skimage.morphology import skeletonize
from scipy import ndimage as ndi

Image.MAX_IMAGE_PIXELS = None
NAMES = ['41','51','52','53','54','61','62','63','64','65','66','67','68','69','610','611','612','613','614','615','616']
XS = [480, 652, 822, 995]; YS = [500, 655, 820, 975, 1145, 1320]
S_PX = 90.0   # pixels per unit
H = 0.22      # lift of over (+) / under (-) strands

def load(page, idx):
    im = Image.open(page).convert('L'); s = im.size[0] / 1413
    r, c = divmod(idx, 4); cx, cy = XS[c] * s, YS[r] * s; h = 95 * s
    a = np.array(im.crop((int(cx - h * 1.05), int(cy - h), int(cx + h * 1.05), int(cy + h)))) < 128
    rows = np.where(a.any(1))[0]
    blocks = []; start = rows[0]; prev = rows[0]
    for y in rows[1:]:
        if y - prev > 15: blocks.append((start, prev)); start = y
        prev = y
    blocks.append((start, prev))
    b = max(blocks, key=lambda t: t[1] - t[0])
    a[:b[0]] = False; a[b[1] + 1:] = False
    return a

NB = [(-1,-1),(-1,0),(-1,1),(0,-1),(0,1),(1,-1),(1,0),(1,1)]

def analyse(a):
    a = ndi.binary_closing(np.pad(a, 5), np.ones((3, 3)), iterations=2)[5:-5, 5:-5]   # 1 px render seams
    sk = skeletonize(a)
    pts = set(map(tuple, np.argwhere(sk)))
    def nbrs(p): return [(p[0]+d[0], p[1]+d[1]) for d in NB if (p[0]+d[0], p[1]+d[1]) in pts]
    for _ in range(3):   # prune short spurs
        deg = {p: len(nbrs(p)) for p in pts}
        for p in [p for p, d in deg.items() if d == 1]:
            path = [p]; cur = p; prv = None
            while True:
                nx = [q for q in nbrs(cur) if q != prv and q not in path]
                if len(nx) != 1: break
                prv = cur; cur = nx[0]; path.append(cur)
                if len(path) > 10: break
            if len(path) <= 10 and len(nbrs(cur)) >= 3:
                for q in path[:-1]: pts.discard(q)
    deg = {p: len(nbrs(p)) for p in pts}
    mask = np.zeros(a.shape, bool)
    for p, d in deg.items():
        if d >= 3: mask[p] = True
    jl, nj = ndi.label(ndi.binary_dilation(mask, np.ones((3, 3)), iterations=2), np.ones((3, 3)))
    jid = {p: int(jl[p]) for p in pts if jl[p] > 0}
    rest = pts - set(jid)
    seen = set(); chains = []
    for p in rest:
        if p in seen: continue
        comp = []; st = [p]; seen.add(p)
        while st:
            q = st.pop(); comp.append(q)
            for d in NB:
                r2 = (q[0]+d[0], q[1]+d[1])
                if r2 in rest and r2 not in seen: seen.add(r2); st.append(r2)
        cs = set(comp)
        ends = [q for q in comp if sum(((q[0]+d[0], q[1]+d[1]) in cs) for d in NB) <= 1]
        if len(comp) < 12: continue
        if not ends: ends = [comp[0]]
        path = [ends[0]]; ps = {ends[0]}
        while True:
            cur = path[-1]
            nx = [(cur[0]+d[0], cur[1]+d[1]) for d in NB if (cur[0]+d[0], cur[1]+d[1]) in cs and (cur[0]+d[0], cur[1]+d[1]) not in ps]
            if not nx: break
            nx.sort(key=lambda t: abs(t[0]-cur[0]) + abs(t[1]-cur[1]))
            path.append(nx[0]); ps.add(nx[0])
        chains.append(path)
    cent = {}
    for j in range(1, nj + 1):
        ps = [p for p, k in jid.items() if k == j]
        if ps: cent[j] = np.mean(ps, axis=0)
    cends = []
    for ci, ch in enumerate(chains):
        for w, p in ((0, ch[0]), (1, ch[-1])):
            jn = None
            for d in [(i, j) for i in range(-3, 4) for j in range(-3, 4)]:
                q = (p[0]+d[0], p[1]+d[1])
                if q in jid: jn = jid[q]; break
            cends.append((ci, w, jn))
    return chains, cent, cends

def tangent(ch, w, k=12):
    if w == 0: p, q = np.array(ch[0], float), np.array(ch[min(k, len(ch)-1)], float)
    else: p, q = np.array(ch[-1], float), np.array(ch[max(len(ch)-1-k, 0)], float)
    t = p - q; return t / (np.linalg.norm(t) + 1e-9)

def build(a, name):
    chains, cent, cends = analyse(a)
    free = [(ci, w) for ci, w, jn in cends if jn is None]
    arms = [(ci, w, jn) for ci, w, jn in cends if jn is not None]
    allpx = np.array([p for ch in chains for p in ch], float)
    def endpt(ci, w): return np.array(chains[ci][0 if w == 0 else -1], float)
    cand = []
    for i in range(len(free)):
        for j in range(i + 1, len(free)):
            (c1, w1), (c2, w2) = free[i], free[j]
            p1, p2 = endpt(c1, w1), endpt(c2, w2)
            d = np.linalg.norm(p2 - p1)
            if d > 90: continue
            u = (p2 - p1) / max(d, 1e-9)
            cos1 = float(tangent(chains[c1], w1) @ u); cos2 = float(tangent(chains[c2], w2) @ -u)
            if cos1 < 0.6 or cos2 < 0.6: continue
            if d < 8:   # rendering seam in a line, not an under-crossing
                if cos1 > 0.9 and cos2 > 0.9: cand.append((d * 0.01, i, j, True))
                continue
            mid = (p1 + p2) / 2
            dist = np.linalg.norm(allpx - mid, axis=1)
            near = (dist < 14) & (np.linalg.norm(allpx - p1, axis=1) > 10) & (np.linalg.norm(allpx - p2, axis=1) > 10)
            if not near.any(): continue
            cand.append((d * (3 - cos1 - cos2), i, j, False))
    cand.sort(); used = set(); link = {}; seam = set()
    for cost, i, j, sm in cand:
        if i in used or j in used: continue
        used.add(i); used.add(j); link[free[i]] = free[j]; link[free[j]] = free[i]
        if sm: seam.add(free[i]); seam.add(free[j])
    unpaired = [f for f in free if f not in link]
    edges = []
    armset = {(ci, w) for ci, w, jn in arms}; jof = {(ci, w): jn for ci, w, jn in arms}
    done = set()
    for ci, w, jn in arms:
        if (ci, w) in done: continue
        poly = [(tuple(cent[jn]), False)]
        cur_c, cur_w = ci, w
        while True:
            done.add((cur_c, cur_w))
            ch = chains[cur_c]; seq = ch if cur_w == 0 else ch[::-1]
            poly += [(tuple(map(float, p)), False) for p in seq]
            out_w = 1 - cur_w
            if (cur_c, out_w) in armset:
                done.add((cur_c, out_w)); jb = jof[(cur_c, out_w)]
                poly.append((tuple(cent[jb]), False)); break
            nxt = link.get((cur_c, out_w))
            if nxt is None: raise RuntimeError(f'{name}: unpaired free end of chain {cur_c}')
            poly[-1] = (poly[-1][0], (cur_c, out_w) not in seam)
            cur_c, cur_w = nxt
        edges.append((jn, jb, poly))
    return edges, cent, chains, unpaired, (len(link) - len(seam)) // 2

def resample(poly, ds=4.0, smooth=6):
    pts = np.array([p for p, g in poly]); gap = np.array([g for p, g in poly])
    seg = np.linalg.norm(np.diff(pts, axis=0), axis=1); s = np.r_[0, np.cumsum(seg)]
    n = max(int(s[-1] / ds), 4); t = np.linspace(0, s[-1], n + 1)
    q = np.c_[np.interp(t, s, pts[:, 0]), np.interp(t, s, pts[:, 1])]
    gs = np.zeros(len(t) - 1, bool)
    for i in range(len(pts) - 1):
        if gap[i]:
            gs |= (t[:-1] >= s[i] - 1e-9) & (t[1:] <= s[i + 1] + 1e-9)
    for _ in range(smooth):
        q[1:-1] = 0.25 * q[:-2] + 0.5 * q[1:-1] + 0.25 * q[2:]
    return q, gs

def seg_inter(p, d, q, e):
    den = d[0]*e[1] - d[1]*e[0]
    if abs(den) < 1e-12: return None
    w = q - p
    u = (w[0]*e[1] - w[1]*e[0]) / den; v = (w[0]*d[1] - w[1]*d[0]) / den
    if 0 <= u < 1 and 0 <= v < 1: return u, v
    return None

def crossings(curves, gaps):
    res = []
    for a in range(len(curves)):
        for b in range(a, len(curves)):
            A, B = curves[a], curves[b]
            for i in range(len(A) - 1):
                for j in range(len(B) - 1):
                    if a == b and j <= i + 1: continue
                    r = seg_inter(A[i], A[i+1]-A[i], B[j], B[j+1]-B[j])
                    if r is None: continue
                    pt = A[i] + r[0] * (A[i+1] - A[i])
                    if min(np.linalg.norm(pt - A[0]), np.linalg.norm(pt - A[-1])) < 8 and min(np.linalg.norm(pt - B[0]), np.linalg.norm(pt - B[-1])) < 8: continue
                    ga, gb = gaps[a][i], gaps[b][j]
                    under = 'a' if ga and not gb else 'b' if gb and not ga else None
                    res.append((a, i + r[0], b, j + r[1], under))
    return res

def lift(curves, cr):
    nodes = [[(0.0, 0.0)] for _ in curves]
    for a, sa, b, sb, under in cr:
        nodes[a].append((sa, -H if under == 'a' else H)); nodes[b].append((sb, -H if under == 'b' else H))
    out = []
    for k, c in enumerate(curves):
        nodes[k].append((len(c) - 1.0, 0.0)); nl = sorted(nodes[k]); z = np.zeros(len(c))
        for idx in range(len(c)):
            for (s0, z0), (s1, z1) in zip(nl[:-1], nl[1:]):
                if s0 <= idx <= s1:
                    t = 0 if s1 == s0 else (idx - s0) / (s1 - s0)
                    z[idx] = z0 + (z1 - z0) * (1 - math.cos(math.pi * t)) / 2; break
        out.append(z)
    return out

def process(page, idx, outdir):
    name = NAMES[idx]; a = load(page, idx)
    edges, cent, chains, unpaired, ngaps = build(a, name)
    cx, cy = a.shape[1] / 2, a.shape[0] / 2
    curves, gaps, ends = [], [], []
    for ja, jb, poly in edges:
        q, g = resample(poly)
        curves.append(q); gaps.append(g); ends.append((ja, jb))
    cr = crossings(curves, gaps)
    info = dict(name=name, gap_links=ngaps, unpaired=len(unpaired), crossings=len(cr), undecided=sum(c[4] is None for c in cr),
                n_edges=len(edges), junctions=sorted({j for e in ends for j in e}))
    zs = lift(curves, cr)
    jids = sorted({j for e in ends for j in e}); vid = {j: i for i, j in enumerate(jids)}
    out = []
    for k, (c, z) in enumerate(zip(curves, zs)):
        pts = [[(x - cx) / S_PX, -(y - cy) / S_PX, float(zz)] for (y, x), zz in zip(c, z)]
        out.append(dict(from_=vid[ends[k][0]], to=vid[ends[k][1]], pts=pts))
    seqs = [[] for _ in curves]
    for ca, sa, cb, sb, under in cr:
        seqs[ca].append((sa, 'U' if under == 'a' else 'O', cb)); seqs[cb].append((sb, 'U' if under == 'b' else 'O', ca))
    info['sequences'] = [[(o, e) for s, o, e in sorted(x)] for x in seqs]
    info['edges'] = out
    json.dump(info, open(f'{outdir}/{name}.json', 'w'))
    img = Image.fromarray((~a * 255).astype(np.uint8)).convert('RGB'); d = ImageDraw.Draw(img)
    cols = [(220, 0, 0), (0, 140, 0), (0, 0, 220), (200, 120, 0)]
    for k, (c, g) in enumerate(zip(curves, gaps)):
        for i in range(len(c) - 1):
            d.line([(c[i][1], c[i][0]), (c[i+1][1], c[i+1][0])], fill=cols[k % 4], width=2 if not g[i] else 6)
    img.save(f'{outdir}/{name}_ov.png')
    return info

def emit(outdir, path, decimate=2):
    """Write the Rust data file `ishii_data.rs` from OUTDIR/<name>.json (run after the extraction)."""
    L = ['// Generated by tools/ishii_extract.py from the diagrams of Table 1 (Ishii et al. 2012); do not edit.',
         '/// `(name, edges)`; an edge is `(from, to, crossing sequence, points)`; the sequence lists `O`ver / `U`nder and the',
         '/// index of the other edge for each crossing along the edge, from `from` to `to`.',
         'pub(super) static DATA: &[(&str, &[(u8, u8, &str, &[[f64; 3]])])] = &[']
    for nm in NAMES:
        d = json.load(open(f'{outdir}/{nm}.json'))
        L.append(f'    ("{nm}", &[')
        for e, seq in zip(d['edges'], d['sequences']):
            pts = e['pts']; keep = pts[::decimate]
            if keep[-1] != pts[-1]: keep.append(pts[-1])
            s = ' '.join(f'{o}{k}' for o, k in seq)
            body = ', '.join(f'[{p[0]:.4f}, {p[1]:.4f}, {p[2]:.4f}]' for p in keep)
            L.append(f'        ({e["from_"]}, {e["to"]}, "{s}", &[{body}]),')
        L.append('    ]),')
    L.append('];')
    open(path, 'w').write('\n'.join(L) + '\n')


if __name__ == '__main__':
    page, outdir = sys.argv[1], sys.argv[2]
    if page == "--emit":
        emit(outdir, sys.argv[3]); sys.exit(0)
    sel = sys.argv[3:] or NAMES
    for nm in sel:
        try:
            i = process(page, NAMES.index(nm), outdir)
            print(nm, {k: v for k, v in i.items() if k not in ('edges',)})
        except Exception as e:
            print(nm, 'ERROR', repr(e))
