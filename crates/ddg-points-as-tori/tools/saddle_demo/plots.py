"""Matplotlib figures of the saddle demo, made directly from the arrays of the authors' program
(cases/<case>/theirs) and the shared truth. Neutral: no crate code involved.

usage: plots.py <root> <out_dir>
"""
import json, os, sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))  # -I drops the script dir
from shapes import Torus, Spool

ROOT, OUT = sys.argv[1], sys.argv[2]
os.makedirs(OUT, exist_ok=True)
SHAPES = {'torus': Torus(0.54, 0.36), 'spool': Spool()}
# reference categorical palette (dataviz skill), fixed order; truth in ink
INK, C_N, C_EREF, C_EGEO, C_GWN = '#0b0b0b', '#2a78d6', '#eb6834', '#1baf7a', '#6250d6'
plt.rcParams.update({'font.size': 9, 'axes.edgecolor': '#52514e', 'axes.labelcolor': '#52514e', 'xtick.color': '#52514e', 'ytick.color': '#52514e'})
SL = np.linspace(-1, 1, 384)
GS = np.stack([np.repeat(SL, 384), np.zeros(384 * 384), np.tile(SL, 384)], -1)
LABEL = {'truth': 'truth (exact SDF)', 'N': '(N) network (theirs)', 'Eref': '(E-ref) exact forms,\npaper rule (theirs)',
         'Efix': '(E-fix) exact forms + one-line\ncentre fix (their code)', 'gwn': 'GWN sign + distance (theirs)'}


def slices(case, methods=('N', 'Eref', 'Efix')):
	d = f'{ROOT}/cases/{case}'
	meta = json.load(open(f'{d}/meta.json'))
	shape = SHAPES[meta['kind']]
	truth = shape.sdf(GS).reshape(384, 384)
	P = np.load(f'{d}/points.npy')
	near = np.abs(P[:, 1]) < 0.03
	fields = [('truth', truth)] + [(m, np.load(f'{d}/theirs/slice_{m}.npy')) for m in methods]
	fig, axs = plt.subplots(2, len(fields), figsize=(10.24, 5.4), dpi=100, constrained_layout=True)
	ext = (-1, 1, -1, 1)
	norm = TwoSlopeNorm(vmin=-0.4, vcenter=0, vmax=0.4)
	for j, (m, phi) in enumerate(fields):
		ax = axs[0, j]
		im = ax.imshow(phi.T, origin='lower', extent=ext, cmap='RdBu_r', norm=norm)
		ax.contour(SL, SL, phi.T, levels=np.arange(-0.4, 0.41, 0.05), colors='#52514e', linewidths=0.3, alpha=0.5)
		ax.contour(SL, SL, truth.T, levels=[0], colors=INK, linewidths=0.9, linestyles='--')
		ax.contour(SL, SL, phi.T, levels=[0], colors=INK, linewidths=1.6)
		ax.plot(P[near, 0], P[near, 2], '.', ms=2.5, color='#fcfcfb', mec=INK, mew=0.3)
		ax.set_title(LABEL[m], fontsize=8.5)
		ax.set_xticks([]); ax.set_yticks([])
		ax2 = axs[1, j]
		if m == 'truth':
			ax2.axis('off')
			ax2.text(0.02, 0.95, f'{case}\nplane y = 0 (through the axis)\n\ntop: phi, isolines every 0.05;\nsolid: phi = 0 of the method,\ndashed: true surface;\ndots: samples with |y| < 0.03\n\nbottom: |phi - phi_true|', va='top', fontsize=8.5, transform=ax2.transAxes)
			continue
		e = np.abs(phi - truth)
		im2 = ax2.imshow(e.T, origin='lower', extent=ext, cmap='Greys', vmin=0, vmax=0.1)
		ax2.contour(SL, SL, truth.T, levels=[0], colors=C_EREF, linewidths=0.7, linestyles='--')
		ax2.set_xticks([]); ax2.set_yticks([])
		band = np.abs(truth) < 0.05
		lab = 'band MAE'
		if os.path.exists(f'{d}/patch64.npy'):  # open patch: only where the nearest surface point is in the sampled region
			band &= (np.cos(shape.nearest_param(GS)[0]) < -0.2).reshape(384, 384)
			lab = 'band MAE (sampled patch)'
		ax2.set_title(f'{lab} {e[band].mean():.4f}', fontsize=8.5)
	fig.colorbar(im, ax=axs[0, :], shrink=0.8, label='phi')
	fig.colorbar(im2, ax=axs[1, 1:], shrink=0.8, label='|error| (clipped at 0.1)')
	fig.savefig(f'{OUT}/slice_{case}.png')
	plt.close(fig)


def torus_sdf(x, c, axis, R, r):
	d = x - c
	v = np.linalg.norm(np.cross(d, axis), axis=-1)
	return np.sign(r) * (np.hypot(v - R, d @ axis) - abs(r))


def axis_vec(a):
	return np.array([np.cos(a[0]) * np.sin(a[1]), np.sin(a[0]) * np.sin(a[1]), np.cos(a[1])])


def principal_frame(shape, kind, p):
	if kind == 'torus':
		th, ph = shape.nearest_param(p[None])
		_, n, em, km, ep, kp = shape.point(th, ph)
	else:
		_, u = shape.sdf(p[None], return_u=True)
		ph = np.arctan2(p[None, 1], p[None, 0])
		_, n, em, km, ep, kp = shape.point(u, ph)
	return n[0], (em[0], km[0]), (ep[0], kp[0])


def tori_plot(picks, fname, W=0.42):
	fig, axs = plt.subplots(len(picks), 2, figsize=(7.6, 3.6 * len(picks)), dpi=100, constrained_layout=True)
	a = np.linspace(-W, W, 401)
	A, B = np.meshgrid(a, a, indexing='ij')
	for row, (case, i, title) in enumerate(picks):
		d = f'{ROOT}/cases/{case}'
		kind = json.load(open(f'{d}/meta.json'))['kind']
		shape = SHAPES[kind]
		P, N = np.load(f'{d}/points.npy'), np.load(f'{d}/normals.npy')
		p = P[i]
		n, (em, km), (ep, kp) = principal_frame(shape, kind, p)
		tor = {m: dict(np.load(f'{d}/theirs/tori_{m}.npz')) for m in ('N', 'Eref', 'Efix')}
		for col, (e, k, dname) in enumerate(((em, km, 'meridian'), (ep, kp, 'parallel'))):
			ax = axs[row, col]
			X = p + A[..., None] * e + B[..., None] * n
			ts = shape.sdf(X.reshape(-1, 3)).reshape(A.shape)
			ax.contourf(A, B, ts, levels=[-10, 0], colors=['#e8e7e2'])
			ax.contour(A, B, ts, levels=[0], colors=INK, linewidths=2.0)
			handles = []
			for m, c, ls in (('N', C_N, '-'), ('Eref', C_EREF, '--'), ('Efix', C_EGEO, '-')):
				t = tor[m]
				g = torus_sdf(X.reshape(-1, 3), t['centers'][i], axis_vec(t['axes_angles'][i] if 'axes_angles' in t else t['axes'][i]), t['major_radii'][i], t['minor_radii'][i]).reshape(A.shape)
				cs = ax.contour(A, B, g, levels=[0], colors=c, linewidths=1.4, linestyles=ls)
				# shade the torus's own inside (g < 0) faintly so an inside-out torus is visible
				ax.contourf(A, B, g, levels=[-10, 0], colors=[c], alpha=0.10)
				handles.append(plt.Line2D([], [], color=c, ls=ls, lw=1.4, label={'N': '(N) network torus', 'Eref': '(E-ref) paper rule', 'Efix': '(E-fix) one-line fix'}[m] + f'  R {t["major_radii"][i]:.3f}, r {t["minor_radii"][i]:+.3f}'))
			ax.annotate('', xy=(0, 0.12), xytext=(0, 0), arrowprops=dict(arrowstyle='->', color=INK, lw=1))
			ax.plot([0], [0], 'o', color=INK, ms=4)
			ax.set_aspect('equal'); ax.set_xlim(-W, W); ax.set_ylim(-W, W)
			ax.set_xlabel(f'{dname} direction (true kappa {k:+.2f})'); ax.set_ylabel('outward normal')
			ax.set_title(f'{title}: {dname} plane', fontsize=9)
			if col == 1:
				handles.insert(0, plt.Line2D([], [], color=INK, lw=2, label='true surface (grey: inside)'))
				ax.legend(handles=handles, fontsize=7, loc='lower left', framealpha=0.9)
	fig.savefig(f'{OUT}/{fname}')
	plt.close(fig)


def summary_chart(fname):
	S = json.load(open(f'{ROOT}/summary.json'))
	cases = list(S)
	series = [('theirs:N', '(N) network', C_N, 'o'), ('theirs:Eref', '(E-ref) exact forms, paper rule', C_EREF, 's'),
	          ('theirs:Efix', '(E-fix) exact forms + one-line fix', C_EGEO, 'D'), ('theirs:gwn', 'GWN sign + nearest distance', C_GWN, '^')]
	fig, axs = plt.subplots(1, 2, figsize=(10.24, 4.6), dpi=100, layout='constrained', sharey=True)
	for ax, reg in zip(axs, ('saddle', 'convex')):
		for j, (k, lab, c, mk) in enumerate(series):
			y = [S[cs]['grid'][k][reg]['band_mae'] for cs in cases]
			x = np.arange(len(cases)) + (j - 1.5) * 0.12
			ax.plot(x, y, mk, color=c, ms=7, mec='#fcfcfb', mew=1.0, label=lab, ls='none')
		ax.set_yscale('log')
		ax.set_xticks(np.arange(len(cases)), cases, rotation=35, ha='right')
		ax.grid(axis='y', color='#e8e7e2', lw=0.6)
		ax.set_axisbelow(True)
		for sp in ('top', 'right'):
			ax.spines[sp].set_visible(False)
		ax.set_title(f'narrow-band MAE (|phi_true| < 0.05), {reg} region', fontsize=9.5)
	axs[0].set_ylabel('MAE (log scale)')
	h, l = axs[0].get_legend_handles_labels()
	fig.legend(h, l, loc='outside lower center', ncol=4, fontsize=8.5, frameon=False)
	fig.suptitle("Authors' program on the same inputs: field error near the surface, split by the sign of K", fontsize=10)
	fig.savefig(f'{OUT}/{fname}')
	plt.close(fig)


def pick(case, rule):
	d = f'{ROOT}/cases/{case}'
	P = np.load(f'{d}/points.npy'); K = np.load(f'{d}/gauss_true.npy')
	return int(rule(P, K))


if __name__ == '__main__':
	for c in sorted(os.listdir(f'{ROOT}/cases')):
		if os.path.exists(f'{ROOT}/cases/{c}/theirs/slice_N.npy'):
			slices(c)
	t = 'torus-uniform-2048'
	tori_plot([
		(t, pick(t, lambda P, K: np.argmin(K)), 'torus, inner equator'),
		(t, pick(t, lambda P, K: np.argmin(np.abs(np.arctan2(P[:, 2], np.hypot(P[:, 0], P[:, 1]) - 0.54) - 2.2))), 'torus, theta ~ 126 deg'),
		('spool-2048', pick('spool-2048', lambda P, K: np.argmin(np.abs(P[:, 2]))), 'spool, waist'),
	], 'tori_at_saddles.png')
	summary_chart('band_mae_by_region.png')
	print('figures in', OUT)
