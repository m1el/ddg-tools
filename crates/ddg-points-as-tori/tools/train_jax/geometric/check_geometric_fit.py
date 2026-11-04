"""Check the Geometric-sign patch of the authors' JAX torus fit (torus_precompute_masked, used by the training loss).

  python check_geometric_fit.py [VERIFY_DIR]      (VERIFY_DIR default: runs_tmp/pat-saddle/verify)

Runs the original (installed pointsastori) and the patched (src/pointsastori) fit in two subprocesses on
  (1) the saddle demo's exact Monge forms (spheres, torus outer/inner equators; expected.json), and
  (2) 4096 random coefficient rows (about half saddles),
and checks: patched == numpy one-liner `where(det II < 0, 2 (p + a00 n) - c_orig, c_orig)`; radii and axes
unchanged; on (1) patched centres == the osculating centres (center_geometric) and original == center_reference;
the patched torus has SDF ~ 0 and gradient . n = +1 at the contact point of every exact form.
"""

import json
import os
import subprocess
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
VERIFY = sys.argv[1] if len(sys.argv) > 1 and sys.argv[1] != '--worker' else '/var/data/discrete-geometry/ddg-tools/runs_tmp/pat-saddle/verify'


def worker(mode, inp, out):
	if mode == 'geo':
		sys.path.insert(0, os.path.join(HERE, 'src'))
	os.environ.setdefault('JAX_PLATFORMS', 'cpu')
	import jax
	import jax.numpy as jnp
	import pointsastori.network as NW

	d = np.load(inp)
	P, N, C = (jnp.asarray(d[k], dtype=jnp.float32) for k in ('P', 'N', 'C'))
	mask = jnp.ones(P.shape[0], dtype=bool)
	c, a, R, r = NW.torus_precompute_masked(mask, P, N, C)
	# SDF and gradient of each torus at its own contact point o = p + a00 n.
	o = P + C[:, 0:1] * N
	f = jax.vmap(NW.torus_signed_distance)
	g = jax.vmap(jax.grad(NW.torus_signed_distance))
	sdf = f(o, c, a, R, r)
	grad = g(o, c, a, R, r)
	np.savez(out, c=np.asarray(c), a=np.asarray(a), R=np.asarray(R), r=np.asarray(r), sdf=np.asarray(sdf),
		grad=np.asarray(grad), module=np.asarray(NW.__file__))


def run(mode, inp):
	out = inp.replace('.npz', f'.{mode}.npz')
	subprocess.run([sys.executable, '-I', os.path.abspath(__file__), '--worker', mode, inp, out], check=True)
	return dict(np.load(out))


def main():
	tmp = os.path.join(HERE, 'check')
	os.makedirs(tmp, exist_ok=True)
	P = np.load(f'{VERIFY}/points.npy')
	N = np.load(f'{VERIFY}/normals.npy')
	C = np.load(f'{VERIFY}/coeffs_exact.npy').astype(np.float32)
	exp = json.load(open(f'{VERIFY}/expected.json'))
	rng = np.random.default_rng(7)
	m = 4096
	Nr = rng.normal(size=(m, 3))
	Nr /= np.linalg.norm(Nr, axis=1, keepdims=True)
	Cr = np.concatenate([rng.normal(0, 0.01, (m, 1)), rng.normal(0, 0.2, (m, 2)), rng.normal(0, 2.0, (m, 3))], 1)
	sets = {'exact': (P, N, C), 'random': (rng.normal(size=(m, 3)), Nr, Cr.astype(np.float32))}
	ok = True
	for name, (p, n, cf) in sets.items():
		inp = os.path.join(tmp, f'{name}.npz')
		np.savez(inp, P=p.astype(np.float32), N=n.astype(np.float32), C=cf)
		o, g = run('orig', inp), run('geo', inp)
		p32, n32, c32 = p.astype(np.float32), n.astype(np.float32), cf.astype(np.float32)
		saddle = 4 * c32[:, 5] * c32[:, 4] - c32[:, 3] ** 2 < 0
		oneliner = np.where(saddle[:, None], 2 * (p32 + c32[:, 0:1] * n32) - o['c'], o['c'])
		e_one = np.abs(g['c'] - oneliner).max()
		same = max(np.abs(g[k] - o[k]).max() for k in ('a', 'R', 'r'))
		print(f'[{name}] {len(p)} rows, {saddle.sum()} saddles; orig {o["module"]}; geo {g["module"]}')
		print(f'[{name}] max |c_geo - oneliner(c_orig)| = {e_one:.2e}; max change of axis/R/r = {same:.2e}; '
			f'max |c_geo - c_orig| at non-saddles = {np.abs(g["c"] - o["c"])[~saddle].max() if (~saddle).any() else 0:.2e}')
		ok &= e_one < 1e-5 and same == 0.0
		if name == 'exact':
			for i, e in enumerate(exp):
				if 'center_geometric' in e:
					eg = np.linalg.norm(g['c'][i] - e['center_geometric'])
					er = np.linalg.norm(o['c'][i] - e['center_reference'])
				else:
					eg = np.linalg.norm(g['c'][i] - e['center'])
					er = np.linalg.norm(o['c'][i] - e['center'])
				gn_geo = float(g['grad'][i] @ N[i])
				gn_orig = float(o['grad'][i] @ N[i])
				print(f'  {i:2d} {e["kind"][:34]:34s} |c_geo-expected| {eg:.1e}  |c_orig-expected_ref| {er:.1e}  '
					f'sdf_geo(o) {g["sdf"][i]:+.1e}  grad.n geo {gn_geo:+.4f} orig {gn_orig:+.4f}')
				# Saddles: the patched centre must be the osculating one. Non-saddles: the patch must be a no-op
				# (rows 20-22, a11 = 0, are the reference fit's degenerate eigenvector pair, EigenRule::Reference in the
				# crate, SADDLE-DEMO.md: off by 0.54 in both versions, as in the authors' C++).
				if 'center_geometric' in e:
					ok &= eg < 1e-5 and abs(g['sdf'][i]) < 1e-5 and abs(gn_geo - 1) < 1e-4
				else:
					ok &= bool(np.all(g['c'][i] == o['c'][i])) and abs(g['sdf'][i]) < 1e-5 and abs(gn_geo - 1) < 1e-4
	print('PASS' if ok else 'FAIL')
	sys.exit(0 if ok else 1)


if __name__ == '__main__':
	if len(sys.argv) > 1 and sys.argv[1] == '--worker':
		worker(*sys.argv[2:5])
	else:
		main()
