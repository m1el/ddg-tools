"""Generate the inputs of the saddle demo once (numpy); both the authors' program and our crate read them.

usage: make_inputs.py <out_root>

<out_root>/grid_axis64.npy, grid_axis128.npy       linspace(-1, 1, n); grids are meshgrid(ax, ax, ax, 'ij') C order
<out_root>/verify/                                  points.npy normals.npy coeffs_exact.npy (f32) + expected.json
<out_root>/cases/<case>/
    points.npy normals.npy            (n,3) f64, exact unit outward normals
    coeffs_exact.npy                  (n,6) f32, exact Monge coefficients (a00,a01,a10,a11,a02,a20)
    kappa_true.npy                    (n,2) f64, convex-positive (kmin, kmax), |kmin| <= |kmax|
    gauss_true.npy                    (n,) f64
    truth_sdf64.npy, gauss64.npy      (64,64,64) true SDF and K at the nearest surface point
    patch64.npy                       (64,64,64) bool: nearest surface point inside the sampled region (open patch only)
    truth_sdf128.npy                  (128,128,128) f32 true SDF
    surf_samples.npy, surf_gauss.npy  dense area-uniform surface samples (100k) and their K (for level-set distances)
    surf_patch.npy                    bool: dense sample inside the sampled region
    mesh_vertices.npy mesh_faces.npy  fine closed mesh of the true surface (renders only)
    meta.json
"""
import json, os, sys, time
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))  # -I drops the script dir
from shapes import Torus, Spool, exact_coefficients, kappa_from_coefficients, local_frame

OUT = sys.argv[1]
os.makedirs(f'{OUT}/cases', exist_ok=True)
AX64, AX128 = np.linspace(-1, 1, 64), np.linspace(-1, 1, 128)
np.save(f'{OUT}/grid_axis64.npy', AX64)
np.save(f'{OUT}/grid_axis128.npy', AX128)
grid = lambda ax: np.stack(np.meshgrid(ax, ax, ax, indexing='ij'), -1).reshape(-1, 3)
G64, G128 = grid(AX64), grid(AX128)

TOR = Torus(0.54, 0.36)  # R/r = 1.5 (the requested fat torus R .30 r .20, scaled x1.8 to fill the box)
SPOOL = Spool()


def orient(V, F):
    vol = np.einsum('ij,ij->i', V[F[:, 0]], np.cross(V[F[:, 1]], V[F[:, 2]])).sum() / 6
    return (F if vol > 0 else F[:, [0, 2, 1]]), abs(vol)


def kappa_sorted(k1, k2):
    sw = np.abs(k1) > np.abs(k2)
    return np.where(sw, k2, k1), np.where(sw, k1, k2)


def write_case(name, shape, kind, p, n, e1, k1, e2, k2, patch_fn, meta):
    d = f'{OUT}/cases/{name}'
    os.makedirs(d, exist_ok=True)
    c = exact_coefficients(n, e1, k1, e2, k2).astype(np.float32)
    kmin, kmax = kappa_sorted(k1, k2)
    np.save(f'{d}/points.npy', p)
    np.save(f'{d}/normals.npy', n)
    np.save(f'{d}/coeffs_exact.npy', c)
    np.save(f'{d}/kappa_true.npy', np.stack([kmin, kmax], -1))
    np.save(f'{d}/gauss_true.npy', k1 * k2)
    # self-check: curvatures from the stored coefficients
    q1, q2 = kappa_from_coefficients(c)
    meta['coeff_kappa_maxrelerr'] = float(np.max(np.abs(np.stack([q1 - kmin, q2 - kmax])) / np.abs(kmax)))
    for nm, (_, arr) in shape_cache[kind].items():
        np.save(f"{d}/{nm}.npy", np.ascontiguousarray(arr))
    if patch_fn is not None:
        np.save(f'{d}/patch64.npy', patch_fn(G64).reshape(64, 64, 64))
        np.save(f'{d}/surf_patch.npy', patch_fn(shape_cache[kind]['surf_samples'][1]))
    meta.update(kind=kind, n=len(p), saddle_point_fraction=float(np.mean(k1 * k2 < 0)))
    json.dump(meta, open(f'{d}/meta.json', 'w'), indent=1)
    print(name, json.dumps(meta), flush=True)


# ------------------------------------------------------------------ shape-level arrays (shared by cases)
shape_cache = {}
t0 = time.time()
for kind, shape in [('torus', TOR), ('spool', SPOOL)]:
    sd64 = shape.sdf(G64).reshape(64, 64, 64)
    K64, _ = shape.gauss_at_nearest(G64)
    sd128 = shape.sdf(G128).reshape(128, 128, 128).astype(np.float32)
    rng = np.random.default_rng(99)
    if kind == 'torus':
        th, ph = shape.sample(100000, rng)
        ps, _, _, km, _, kp = shape.point(th, ph)
    else:
        u, ph = shape.sample(100000, rng)
        ps, _, _, km, _, kp = shape.point(u, ph)
    V, F = shape.mesh()
    F, vol = orient(V, F)
    shape_cache[kind] = dict(truth_sdf64=(0, sd64), gauss64=(0, K64.reshape(64, 64, 64)), truth_sdf128=(0, sd128),
                             surf_samples=(0, ps), surf_gauss=(0, km * kp), mesh_vertices=(0, V), mesh_faces=(0, F))
    print(f'{kind}: arrays in {time.time() - t0:.1f}s; mesh volume {vol:.5f}; dense-sample saddle fraction {np.mean(km * kp < 0):.4f}', flush=True)
    if kind == 'torus':
        print(f'  torus analytic saddle area fraction {shape.saddle_area_fraction():.4f}')
    else:
        print(f'  spool area {shape.area:.4f}; check sdf on dense samples max |phi| {np.abs(shape.sdf(ps[:2000])).max():.2e}')

# ------------------------------------------------------------------ cases
for N in (512, 2048):
    rng = np.random.default_rng(1000 + N)
    # 1a. full torus, area uniform
    th, ph = TOR.sample(N, rng)
    p, n, em, km, ep, kp = TOR.point(th, ph)
    write_case(f'torus-uniform-{N}', TOR, 'torus', p, n, em, km, ep, kp, None, dict(R=TOR.R, r=TOR.r, sampling='area-uniform'))
    # 1b. full torus, 80 % of the points on the inner (saddle) half
    ni = int(round(0.8 * N))
    thi, phi_ = TOR.sample(ni, rng, (-1.0, 0.0))
    tho, pho = TOR.sample(N - ni, rng, (0.0, 1.0))
    th, ph = np.concatenate([thi, tho]), np.concatenate([phi_, pho])
    p, n, em, km, ep, kp = TOR.point(th, ph)
    write_case(f'torus-inner80-{N}', TOR, 'torus', p, n, em, km, ep, kp, None, dict(R=TOR.R, r=TOR.r, sampling='80% inner half, 20% outer half, area-uniform within each'))
    # 1c. inner half only (open patch: "tracing the inside of a torus")
    th, ph = TOR.sample(N, rng, (-1.0, 0.0))
    p, n, em, km, ep, kp = TOR.point(th, ph)
    patch = lambda x: np.cos(TOR.nearest_param(x)[0]) < -0.2
    write_case(f'torus-inneronly-{N}', TOR, 'torus', p, n, em, km, ep, kp, patch, dict(R=TOR.R, r=TOR.r, sampling='inner half only (cos theta < 0), open patch; evaluated where cos theta* < -0.2'))
    # 2. spool
    u, ph = SPOOL.sample(N, rng)
    p, n, em, km, ep, kp = SPOOL.point(u, ph)
    write_case(f'spool-{N}', SPOOL, 'spool', p, n, em, km, ep, kp, None, dict(h=SPOOL.h, a=SPOOL.a, b=SPOOL.b, sampling='area-uniform'))

# ------------------------------------------------------------------ verification points (closed forms)
vd = f'{OUT}/verify'
os.makedirs(vd, exist_ok=True)
P, Nn, E1, K1, E2, K2, exp = [], [], [], [], [], [], []
# sphere r = 0.5 at the origin, normals including |n_x| > 0.9 (the other frame branch)
dirs = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1], [0.95, 0.3, 0.1], [-0.3, 0.2, 0.93], [0.5, -0.5, 0.7], [-0.99, 0.1, 0.0], [0.2, 0.9, -0.4]], float)
dirs /= np.linalg.norm(dirs, axis=1, keepdims=True)
for dvec in dirs:
    tmp = np.array([0.3, 0.7, 0.2]) - dvec * np.dot([0.3, 0.7, 0.2], dvec)
    e1 = tmp / np.linalg.norm(tmp)
    e2 = np.cross(dvec, e1)
    P.append(0.5 * dvec); Nn.append(dvec); E1.append(e1); K1.append(2.0); E2.append(e2); K2.append(2.0)
    exp.append(dict(kind='sphere', center=[0, 0, 0], R=0.0, r=0.5, kmin=2.0, kmax=2.0))
# The torus of the verification is rotated by a generic rotation Q: with the axis along z, the local
# frame's s_hat is exactly the meridian direction at the equators, a11 = 0 and v2 = 0, and the
# reference eigenvector pair (both programs' EigenRule::Reference) degenerates (NOTES-WP2 item 1).
ax_ = np.array([1.0, 2.0, 3.0]) / np.sqrt(14.0)
Kx = np.array([[0, -ax_[2], ax_[1]], [ax_[2], 0, -ax_[0]], [-ax_[1], ax_[0], 0]])
Q = np.eye(3) + np.sin(0.7) * Kx + (1 - np.cos(0.7)) * Kx @ Kx
rot = lambda v: v @ Q.T
# torus outer equator (convex): the osculating torus is the torus itself
for ph in np.linspace(0, 2 * np.pi, 6, endpoint=False) + 0.1:
    p, n, em, km, ep, kp = TOR.point(np.array([0.0]), np.array([ph]))
    P.append(rot(p[0])); Nn.append(rot(n[0])); E1.append(rot(em[0])); K1.append(km[0]); E2.append(rot(ep[0])); K2.append(kp[0])
    exp.append(dict(kind='torus-outer-equator', center=[0, 0, 0], axis=list(Q[:, 2]), R=TOR.R, r=TOR.r, kmin=1 / (TOR.R + TOR.r), kmax=1 / TOR.r))
# torus inner equator (saddle): the osculating torus has tube radius R - r about the parallel circle and core radius
# r + (R - r) = R, centred at o - r n* with n* = -e_rho (pointing at the axis): centre R e_rho. The paper's rule mirrors
# the centre through o: (R - 2r) e_rho.
for ph in np.linspace(0, 2 * np.pi, 6, endpoint=False) + 0.3:
    p, n, em, km, ep, kp = TOR.point(np.array([np.pi]), np.array([ph]))
    erho = np.array([np.cos(ph), np.sin(ph), 0])
    P.append(rot(p[0])); Nn.append(rot(n[0])); E1.append(rot(em[0])); K1.append(km[0]); E2.append(rot(ep[0])); K2.append(kp[0])
    exp.append(dict(kind='torus-inner-equator (saddle)', center_geometric=list(rot(TOR.R * erho)), center_reference=list(rot((TOR.R - 2 * TOR.r) * erho)),
                    R=TOR.R, r=TOR.R - TOR.r, kmin=1 / TOR.r, kmax=-1 / (TOR.R - TOR.r)))
# the same outer-equator points WITHOUT the rotation (frame-aligned, a11 = 0): documents the eigenvector degeneracy
for ph in np.linspace(0, 2 * np.pi, 3, endpoint=False) + 0.1:
    p, n, em, km, ep, kp = TOR.point(np.array([0.0]), np.array([ph]))
    P.append(p[0]); Nn.append(n[0]); E1.append(em[0]); K1.append(km[0]); E2.append(ep[0]); K2.append(kp[0])
    exp.append(dict(kind='torus-outer-equator, axis z (a11 = 0: degenerate eigenvector pair)', center=[0, 0, 0], axis=[0, 0, 1], R=TOR.R, r=TOR.r, kmin=1 / (TOR.R + TOR.r), kmax=1 / TOR.r))
P, Nn, E1, E2 = map(np.array, (P, Nn, E1, E2))
K1, K2 = np.array(K1), np.array(K2)
np.save(f'{vd}/points.npy', P); np.save(f'{vd}/normals.npy', Nn)
np.save(f'{vd}/coeffs_exact.npy', exact_coefficients(Nn, E1, K1, E2, K2).astype(np.float32))
kmin, kmax = kappa_sorted(K1, K2)
np.save(f'{vd}/kappa_true.npy', np.stack([kmin, kmax], -1))
json.dump(exp, open(f'{vd}/expected.json', 'w'), indent=1)
print('verify points', len(P), 'total', f'{time.time() - t0:.1f}s')
