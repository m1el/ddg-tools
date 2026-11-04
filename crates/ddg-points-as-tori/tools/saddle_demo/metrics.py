"""Neutral metrics of the saddle demo: one numpy script reads every method's grids (the authors' program
in cases/<case>/theirs, our crate in cases/<case>/ours) and the shared truth.

usage: metrics.py <root> [case ...]     -> cases/<case>/metrics.json, <root>/summary.json
"""
import json, os, sys, time
import numpy as np
from scipy.spatial import cKDTree
from skimage.measure import marching_cubes

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))  # -I drops the script dir
from shapes import Torus, Spool, kappa_from_coefficients

ROOT = sys.argv[1]
cases = sys.argv[2:] or sorted(os.listdir(f'{ROOT}/cases'))
SHAPES = {'torus': Torus(0.54, 0.36), 'spool': Spool()}
BAND, NEAR_SIGN = 0.05, (0.005, 0.05)
THEIRS = ['N', 'N-all', 'Nfix', 'Eref', 'Eref-all', 'Efix', 'Efix-all', 'gwn']
OURS = ['Nref', 'Nref-all', 'Ngeo', 'Eref', 'Eref-all', 'Egeo', 'Egeo-all', 'gwn']


def eik(phi, h, mask):
	g = np.stack(np.gradient(phi, h), -1)
	e = np.abs(np.linalg.norm(g, axis=-1) - 1)
	inner = np.zeros_like(mask)
	inner[1:-1, 1:-1, 1:-1] = True
	m = mask & inner
	return float(e[m].mean()) if m.any() else float('nan')


def grid_metrics(phi, truth, K, mask, h):
	"""MAE over `mask`, band MAE, sign errors, eikonality; also split by region (K < 0 saddle, K > 0 convex)."""
	e = np.abs(phi - truth)
	band = mask & (np.abs(truth) < BAND)
	far = mask & (np.abs(truth) > BAND)
	near = mask & (np.abs(truth) > NEAR_SIGN[0]) & (np.abs(truth) < NEAR_SIGN[1])
	wrong = (phi < 0) != (truth < 0)
	out = dict(mae=float(e[mask].mean()), max=float(e[mask].max()), band_mae=float(e[band].mean()),
	           sign_err=float(wrong[far].mean()), band_sign_err=float(wrong[near].mean()),
	           eik=eik(phi, h, mask), eik_band=eik(phi, h, band), n=int(mask.sum()), band_n=int(band.sum()))
	mean = lambda v: float(v.mean()) if v.size else float('nan')
	vmax = lambda v: float(v.max()) if v.size else float('nan')
	for reg, rm in (('saddle', K < 0), ('convex', K > 0)):
		m, b, nr = mask & rm, band & rm, near & rm
		out[reg] = dict(mae=mean(e[m]), band_mae=mean(e[b]), band_max=vmax(e[b]),
		                band_sign_err=mean(wrong[nr]), eik_band=eik(phi, h, b), band_n=int(b.sum()))
	return out


def closest_on_triangles(p, a, b, c):
	"""Vectorized closest point on triangles (Ericson, Real-Time Collision Detection 5.1.5)."""
	ab, ac, ap = b - a, c - a, p - a
	d1, d2 = np.einsum('ij,ij->i', ab, ap), np.einsum('ij,ij->i', ac, ap)
	bp = p - b
	d3, d4 = np.einsum('ij,ij->i', ab, bp), np.einsum('ij,ij->i', ac, bp)
	cp = p - c
	d5, d6 = np.einsum('ij,ij->i', ab, cp), np.einsum('ij,ij->i', ac, cp)
	va, vb, vc = d3 * d6 - d5 * d4, d5 * d2 - d1 * d6, d1 * d4 - d3 * d2
	den = va + vb + vc
	den = np.where(np.abs(den) > 1e-300, den, 1e-300)
	v, w = vb / den, vc / den
	res = a + ab * v[:, None] + ac * w[:, None]
	sel = lambda m, q: np.where(m[:, None], q, res)
	# edge regions (later assignments take precedence in the order of the reference)
	with np.errstate(divide='ignore', invalid='ignore'):
		m = (va <= 0) & (d4 - d3 >= 0) & (d5 - d6 >= 0)
		res = sel(m, b + (c - b) * np.nan_to_num((d4 - d3) / ((d4 - d3) + (d5 - d6)))[:, None])
		m = (vb <= 0) & (d2 >= 0) & (d6 <= 0)
		res = sel(m, a + ac * np.nan_to_num(d2 / (d2 - d6))[:, None])
		m = (vc <= 0) & (d1 >= 0) & (d3 <= 0)
		res = sel(m, a + ab * np.nan_to_num(d1 / (d1 - d3))[:, None])
	res = sel((d6 >= 0) & (d5 <= d6), c)
	res = sel((d3 >= 0) & (d4 <= d3), b)
	res = sel((d1 <= 0) & (d2 <= 0), a)
	return res


def mesh_distance(q, V, F, k=24):
	cen = V[F].mean(1)
	_, idx = cKDTree(cen).query(q, k=k, workers=8)
	best = np.full(len(q), np.inf)
	for j in range(k):
		f = F[idx[:, j]]
		cp = closest_on_triangles(q, V[f[:, 0]], V[f[:, 1]], V[f[:, 2]])
		best = np.minimum(best, np.linalg.norm(cp - q, axis=1))
	return best


def stats(d):
	if len(d) == 0:
		return dict(n=0)
	return dict(n=int(len(d)), mean=float(d.mean()), p95=float(np.quantile(d, 0.95)), max=float(d.max()))


def level_set(phi128, ax, shape, surf, surf_K, surf_mask, vert_ok):
	"""Zero level set by marching cubes on the 128^3 grid; distances to the true surface split by region."""
	h = ax[1] - ax[0]
	if not (phi128.min() < 0 < phi128.max()):
		return dict(empty=True)
	V, F, _, _ = marching_cubes(phi128.astype(np.float64), 0.0, spacing=(h, h, h))
	V = V + ax[0]
	d_rt = np.abs(shape.sdf(V))  # recon -> truth (exact: the truth is an exact SDF)
	Kv, _ = shape.gauss_at_nearest(V)
	ok = vert_ok(V)
	d_tr = mesh_distance(surf[surf_mask], V, F)  # truth -> recon (exact point-triangle distance)
	Kt = surf_K[surf_mask]
	out = dict(vertices=int(len(V)), faces=int(len(F)))
	for reg, rv, rt in (('all', np.ones(len(V), bool), np.ones(len(Kt), bool)), ('saddle', Kv < 0, Kt < 0), ('convex', Kv > 0, Kt > 0)):
		a, b = stats(d_rt[ok & rv]), stats(d_tr[rt])
		out[reg] = dict(recon_to_truth=a, truth_to_recon=b, hausdorff=float(max(a.get('max', 0), b.get('max', 0))),
		                mean_sym=float(0.5 * (a.get('mean', 0) + b.get('mean', 0))))
	return out


def torus_sdf_grad(x, c, axis, R, r):
	"""The torus SDF of both programs, sign(minor)(|(|d x u| - R, d.u)| - |minor|), and its gradient by central differences."""
	def f(y):
		d = y - c
		v = np.linalg.norm(np.cross(d, axis), axis=-1)
		return np.sign(r) * (np.hypot(v - R, np.sum(d * axis, -1)) - np.abs(r))
	g = np.stack([(f(x + e) - f(x - e)) / 2e-6 for e in np.eye(3) * 1e-6], -1)
	return f(x), g


def axis_of(t):
	if 'axes' in t and t['axes'].shape[1] == 2:
		a = t['axes']
		return np.stack([np.cos(a[:, 0]) * np.sin(a[:, 1]), np.sin(a[:, 0]) * np.sin(a[:, 1]), np.cos(a[:, 1])], -1)
	return t['axes']


def torus_stats(t, P, N, Ktrue):
	ax = axis_of(t)
	g, gr = torus_sdf_grad(P, t['centers'], ax, t['major_radii'], t['minor_radii'])
	dn = np.sum(gr * N, -1)
	R, r = t['major_radii'], np.abs(t['minor_radii'])
	out = {}
	for reg, m in (('all', np.ones(len(P), bool)), ('saddle', Ktrue < 0), ('convex', Ktrue > 0)):
		if not m.any():
			continue
		out[reg] = dict(inside_out=float(np.mean(dn[m] < 0)), contact_max=float(np.abs(g[m]).max()), spindle=float(np.mean(R[m] < r[m])),
		                R_over_r_median=float(np.median(R[m] / r[m])), r_median=float(np.median(r[m])))
	return out


def kappa_stats(c, sigma, kt):
	kmin, kmax = kappa_from_coefficients(c)
	ts = kt[:, 0] * kt[:, 1] < 0
	an = lambda a, b: np.abs(a) / np.maximum(np.abs(b), 1e-300)
	out = dict(pred_saddle_frac=float(np.mean(kmin * kmax < 0)), true_saddle_frac=float(np.mean(ts)),
	           pred_saddle_at_true_saddle=float(np.mean(kmin[ts] * kmax[ts] < 0)) if ts.any() else None,
	           pred_anisotropy_median=float(np.median(an(kmin, kmax))), true_anisotropy_median=float(np.median(an(kt[:, 0], kt[:, 1]))))
	if ts.any():
		out.update(true_kmin_sigma_med_saddle=float(np.median(kt[ts, 0] * sigma[ts])), true_kmax_sigma_med_saddle=float(np.median(kt[ts, 1] * sigma[ts])),
		           pred_kmin_sigma_med_saddle=float(np.median(kmin[ts] * sigma[ts])), pred_kmax_sigma_med_saddle=float(np.median(kmax[ts] * sigma[ts])),
		           pred_anisotropy_median_saddle=float(np.median(an(kmin[ts], kmax[ts]))), true_anisotropy_median_saddle=float(np.median(an(kt[ts, 0], kt[ts, 1]))),
		           pred_H_sigma_med_saddle=float(np.median(0.5 * (kmin[ts] + kmax[ts]) * sigma[ts])), true_H_sigma_med_saddle=float(np.median(0.5 * (kt[ts, 0] + kt[ts, 1]) * sigma[ts])))
	return out


summary = {}
ax64, ax128 = np.load(f'{ROOT}/grid_axis64.npy'), np.load(f'{ROOT}/grid_axis128.npy')
for case in cases:
	t0 = time.time()
	d = f'{ROOT}/cases/{case}'
	meta = json.load(open(f'{d}/meta.json'))
	shape = SHAPES[meta['kind']]
	L = lambda f: np.load(f'{d}/{f}')
	truth, K64 = L('truth_sdf64.npy'), L('gauss64.npy')
	P, N, kt, Kp = L('points.npy'), L('normals.npy'), L('kappa_true.npy'), L('gauss_true.npy')
	surf, surf_K = L('surf_samples.npy'), L('surf_gauss.npy')
	open_patch = os.path.exists(f'{d}/patch64.npy')
	tree = cKDTree(P)
	if open_patch:
		# evaluate only near the sampled patch: nearest surface point inside the patch, |phi_true| < 0.1
		mask = L('patch64.npy') & (np.abs(truth) < 0.1)
		surf_mask = L('surf_patch.npy')
		vert_ok = lambda V: (np.cos(shape.nearest_param(V)[0]) < -0.2) & (tree.query(V, workers=8)[0] < 0.1)
	else:
		mask = np.ones_like(truth, bool)
		surf_mask = np.ones(len(surf), bool)
		vert_ok = lambda V: np.ones(len(V), bool)
	h64 = ax64[1] - ax64[0]
	res = dict(case=case, meta=meta, open_patch=open_patch, grid={}, level_set={}, crosscheck={})
	grids = {}
	for src, names in (('theirs', THEIRS), ('ours', OURS)):
		for m in names:
			f = f'{d}/{src}/grid64_{m}.npy'
			if os.path.exists(f):
				grids[(src, m)] = np.load(f)
				res['grid'][f'{src}:{m}'] = grid_metrics(grids[(src, m)], truth, K64, mask, h64)
			f128 = f'{d}/{src}/grid128_{m}.npy'
			if os.path.exists(f128):
				res['level_set'][f'{src}:{m}'] = level_set(np.load(f128), ax128, shape, surf, surf_K, surf_mask, vert_ok)
	# cross-check: our crate vs their program on the same inputs
	for a, b in (('N', 'Nref'), ('N-all', 'Nref-all'), ('Eref', 'Eref'), ('Eref-all', 'Eref-all'), ('Efix', 'Egeo'), ('Efix-all', 'Egeo-all'), ('gwn', 'gwn')):
		if ('theirs', a) in grids and ('ours', b) in grids:
			res['crosscheck'][f'theirs:{a} vs ours:{b}'] = float(np.abs(grids[('theirs', a)] - grids[('ours', b)]).max())
	for a, b in (('N', 'Nref'), ('Eref', 'Eref'), ('Efix', 'Egeo')):
		fa, fb = f'{d}/theirs/grid128_{a}.npy', f'{d}/ours/grid128_{b}.npy'
		if os.path.exists(fa) and os.path.exists(fb):
			res['crosscheck'][f'128: theirs:{a} vs ours:{b}'] = float(np.abs(np.load(fa).astype(float) - np.load(fb).astype(float)).max())
	nc_t, nc_o = L('theirs/net_coeffs.npy'), L('ours/net_coeffs.npy')
	res['crosscheck']['net_coeffs theirs vs ours'] = float(np.abs(nc_t - nc_o).max())
	# curvature statistics: network (their coefficients) and exact
	sig = L('theirs/net_sigma.npy').astype(float)
	res['kappa'] = dict(network=kappa_stats(nc_t.astype(float), sig, kt), exact=kappa_stats(L('coeffs_exact.npy').astype(float), sig, kt))
	dE = -L('ours/diag_E.npy')  # fit convention convex negative -> convex positive
	rel = np.max(np.abs(dE - kt), 1) / np.abs(kt[:, 1])
	res['kappa']['exact_fit_diagnostics_vs_truth'] = dict(max_rel=float(rel.max()), n_over_1em3=int((rel > 1e-3).sum()))
	# the tori themselves (authors' program): osculation side, spindle fraction
	res['tori'] = {m: torus_stats(dict(np.load(f'{d}/theirs/tori_{m}.npz')), P, N, Kp) for m in ('N', 'Eref', 'Efix')}
	json.dump(res, open(f'{d}/metrics.json', 'w'), indent=1)
	summary[case] = res
	print(f'{case}: {time.time() - t0:.1f}s', flush=True)
json.dump(summary, open(f'{ROOT}/summary.json', 'w'), indent=1)
