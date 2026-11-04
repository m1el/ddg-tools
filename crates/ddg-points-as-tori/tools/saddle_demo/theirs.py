"""The authors' program (points-as-tori, MIT, Nicole Feng) on the saddle-demo inputs. Nothing of theirs is modified.

usage: theirs.py <root> [case ...]        (run with /var/data/discrete-geometry/PointsAsTori-build/.venv/bin/python)

Per case, into <root>/cases/<case>/theirs/:
  (N)     their network + their fit + their SDF          = PointsAsTori(P, N)
  (E-ref) the EXACT fundamental forms (coeffs_exact.npy) through THEIR fit_tori_from_forms + their SDF
  (E-fix) (E-ref) with the torus centre moved to the osculating side at saddle points (one numpy line, below)
  (N-fix) (N) with the same one-line fix
  gwn     their winding number (their CGAL point areas) for the sign, nearest-sample distance for the magnitude
Grids: grid64_<m>.npy (f64; 32-NN "accelerate" evaluation, their default) and grid64_<m>-all.npy (all points),
grid128_<m>.npy (f32, 32-NN) for level-set extraction. Tori: tori_<m>.npz.
"""
import json, os, sys, time
import numpy as np
import jax, jax.numpy as jnp
from flax import nnx
from scipy.spatial import cKDTree

import pointsastori as pat
from pointsastori import infer, shape_3d
from pointsastori.network import FundamentalFormPredictor

ROOT = sys.argv[1]
GWN_ONLY = '--gwn-only' in sys.argv  # recompute only the GWN baseline
cases = [c for c in sys.argv[2:] if not c.startswith('--')] or sorted(os.listdir(f'{ROOT}/cases'))
grid = lambda ax: np.stack(np.meshgrid(ax, ax, ax, indexing='ij'), -1).reshape(-1, 3)
G64, G128 = grid(np.load(f'{ROOT}/grid_axis64.npy')), grid(np.load(f'{ROOT}/grid_axis128.npy'))
SL = np.linspace(-1, 1, 384)
GS = np.stack([np.repeat(SL, 384), np.zeros(384 * 384), np.tile(SL, 384)], -1)  # the plane y = 0, [ix, iz]
model, K_NB = FundamentalFormPredictor.load_saved_model(infer.DEFAULT_MODEL_PATH)
nnx.eval_mode(model)


@nnx.jit
def sigma_of(model, positions, normals):
	_, sigma = model.input_encoder(positions, normals)
	return sigma


def fix_saddle_centres(P, N, coeffs, centers):
	"""The one-line sign fix (crate SignRule::Geometric, Lean pat-fit-contact): at a saddle (det II < 0, i.e.
	4 a20 a02 - a11^2 < 0) the paper's centre o - s n*/|kmin| is the mirror image of the osculating one
	o + n*/kmin through the contact point o = p + a00 n; reflect it back. Radii, axis and sign are unchanged."""
	a00, a11, a02, a20 = coeffs[:, 0:1], coeffs[:, 3], coeffs[:, 4], coeffs[:, 5]
	saddle = 4 * a20 * a02 - a11 * a11 < 0
	return np.where(saddle[:, None], 2 * (P + a00 * N) - centers, centers), int(saddle.sum())


def tori_of(obj):
	c, a, M, m = obj._tdf.get_tori()
	return dict(centers=c, axes=a, major_radii=M, minor_radii=m)


def evaluate(obj, out, tag, fine=True, all_points=True):
	t = time.time()
	np.save(f'{out}/grid64_{tag}.npy', obj.signed_distance(G64, accelerate=True).reshape(64, 64, 64))
	if all_points:
		np.save(f'{out}/grid64_{tag}-all.npy', obj.signed_distance(G64, accelerate=False).reshape(64, 64, 64))
	np.save(f'{out}/slice_{tag}.npy', obj.signed_distance(GS, accelerate=True).reshape(384, 384))
	if fine:
		np.save(f'{out}/grid128_{tag}.npy', obj.signed_distance(G128, accelerate=True).reshape(128, 128, 128).astype(np.float32))
	return time.time() - t


def verify():
	d = f'{ROOT}/verify'
	P, N = np.load(f'{d}/points.npy'), np.load(f'{d}/normals.npy')
	C = np.load(f'{d}/coeffs_exact.npy').astype(np.float64)
	exp = json.load(open(f'{d}/expected.json'))
	c, a, M, m = shape_3d.fit_tori_from_forms(P, N, C)
	cf, _ = fix_saddle_centres(P, N, C, c)
	axis = np.stack([np.cos(a[:, 0]) * np.sin(a[:, 1]), np.sin(a[:, 0]) * np.sin(a[:, 1]), np.cos(a[:, 1])], -1)
	rows = []
	for i, e in enumerate(exp):
		r = dict(kind=e['kind'], R=float(M[i]), r=float(m[i]), center=c[i].tolist(), center_fixed=cf[i].tolist(), axis=axis[i].tolist())
		r['err_R'] = abs(M[i] - e['R']); r['err_r'] = abs(abs(m[i]) - e['r'])
		if 'center' in e:
			r['err_center'] = float(np.linalg.norm(c[i] - e['center']))
		else:
			r['err_center_reference'] = float(np.linalg.norm(c[i] - e['center_reference']))
			r['err_center_fixed_vs_geometric'] = float(np.linalg.norm(cf[i] - e['center_geometric']))
		if 'axis' in e:
			r['err_axis'] = float(min(np.linalg.norm(axis[i] - e['axis']), np.linalg.norm(axis[i] + np.array(e['axis']))))
		rows.append(r)
	json.dump(rows, open(f'{d}/verify_theirs.json', 'w'), indent=1)
	for kind in dict.fromkeys(r['kind'] for r in rows):
		rs = [r for r in rows if r['kind'] == kind]
		worst = {k: float(max(r[k] for r in rs)) for k in rs[0] if k.startswith('err_')}
		print(f'verify (theirs) {kind}:', json.dumps(worst), flush=True)


def gwn(objN, P, out, log):
	"""Their winding number (their CGAL point areas) for the sign, nearest-sample distance for the magnitude.
	Their winding_number sums a_i n_i.(q - p_i)/(4 pi |q - p_i|^3), i.e. MINUS the usual generalized winding
	number (-1 inside a closed surface, 0 outside), so inside is w < -1/2."""
	areas = objN._tdf.get_point_areas()
	log['area_total'] = float(np.sum(areas))
	tree = cKDTree(P)
	for G, tag, dt in ((G64, 'grid64_gwn', np.float64), (G128, 'grid128_gwn', np.float32), (GS, 'slice_gwn', np.float64)):
		w = objN._tdf.winding_number(G)
		dist, _ = tree.query(G, workers=8)
		shape = (384, 384) if tag.startswith('slice') else (round(len(G) ** (1 / 3)),) * 3
		np.save(f'{out}/{tag}.npy', np.where(w < -0.5, -dist, dist).reshape(shape).astype(dt))
	log['gwn_minus_w_inside_median'] = float(np.median(w[w < -0.5])) if (w < -0.5).any() else None


if GWN_ONLY:
	for case in cases:
		d = f'{ROOT}/cases/{case}'
		P, N = np.load(f'{d}/points.npy'), np.load(f'{d}/normals.npy')
		log = json.load(open(f'{d}/theirs/log.json'))
		gwn(pat.PointsAsTori(P, N, tori=dict(np.load(f'{d}/theirs/tori_N.npz'))), P, f'{d}/theirs', log)
		json.dump(log, open(f'{d}/theirs/log.json', 'w'), indent=1)
		print(case, 'gwn', log['area_total'], flush=True)
	sys.exit(0)

verify()
for case in cases:
	t0 = time.time()
	d = f'{ROOT}/cases/{case}'
	out = f'{d}/theirs'
	os.makedirs(out, exist_ok=True)
	P, N = np.load(f'{d}/points.npy'), np.load(f'{d}/normals.npy')
	C = np.load(f'{d}/coeffs_exact.npy').astype(np.float64)
	log = dict(case=case, n=len(P), jax_backend=jax.default_backend())
	# (N): their public API end to end
	objN = pat.PointsAsTori(P, N)
	idx = shape_3d.get_neighbors(P, K_NB, np.array([]))
	net = np.array(model(P[idx], N[idx]))
	np.save(f'{out}/net_coeffs.npy', net)
	np.save(f'{out}/net_sigma.npy', np.array(sigma_of(model, P[idx], N[idx])))
	tN = tori_of(objN)
	c2, a2, M2, m2 = shape_3d.fit_tori_from_forms(P, N, net)
	log['refit_vs_api_maxdiff'] = float(max(np.abs(c2 - tN['centers']).max(), np.abs(M2 - tN['major_radii']).max(), np.abs(m2 - tN['minor_radii']).max()))
	np.savez(f'{out}/tori_N.npz', **tN)
	log['sec_N'] = evaluate(objN, out, 'N')
	# (E-ref): exact forms through their fit
	c, a, M, m = shape_3d.fit_tori_from_forms(P, N, C)
	tE = dict(centers=c, axes=a, major_radii=M, minor_radii=m)
	np.savez(f'{out}/tori_Eref.npz', **tE)
	log['sec_Eref'] = evaluate(pat.PointsAsTori(P, N, tori=tE), out, 'Eref')
	# (E-fix), (N-fix): their tori with the one-line centre fix
	cf, ns = fix_saddle_centres(P, N, C, c)
	tF = dict(tE, centers=cf)
	np.savez(f'{out}/tori_Efix.npz', **tF)
	log['saddles_exact'] = ns
	log['sec_Efix'] = evaluate(pat.PointsAsTori(P, N, tori=tF), out, 'Efix')
	cnf, nsn = fix_saddle_centres(P, N, net.astype(np.float64), tN['centers'])
	log['saddles_network'] = nsn
	tNF = dict(tN, centers=cnf)
	np.savez(f'{out}/tori_Nfix.npz', **tNF)
	log['sec_Nfix'] = evaluate(pat.PointsAsTori(P, N, tori=tNF), out, 'Nfix', fine=nsn > 0, all_points=nsn > 0)
	gwn(objN, P, out, log)
	log['seconds'] = time.time() - t0
	json.dump(log, open(f'{out}/log.json', 'w'), indent=1)
	print(json.dumps(log), flush=True)
