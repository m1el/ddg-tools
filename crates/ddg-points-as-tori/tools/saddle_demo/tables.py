"""Markdown tables from summary.json (saddle demo).   usage: tables.py <root> > tables.md"""
import json, sys

S = json.load(open(f'{sys.argv[1]}/summary.json'))
f = lambda x, k=4: '—' if x is None or x != x else f'{x:.{k}f}'
pct = lambda x: '—' if x is None or x != x else f'{100 * x:.1f} %'
METHODS = [('theirs:N', '(N) network, theirs'), ('theirs:N-all', '(N) all points, theirs'), ('theirs:Eref', '(E-ref) exact, paper rule, theirs'),
           ('theirs:Efix', '(E-fix) exact + one-line fix, their code'), ('ours:Egeo', '(E-geo) exact, Geometric, ours'), ('theirs:gwn', 'GWN sign + distance, theirs')]

print('## Field error by region (64³ grid; band |φ_true| < 0.05; region = sign of K at the nearest true surface point)\n')
for case, r in S.items():
	print(f'### {case} ({r["meta"]["sampling"]}; {pct(r["meta"]["saddle_point_fraction"])} of the points are saddles)\n')
	print('| method | grid MAE | band MAE saddle | band MAE convex | band max saddle | near sign err saddle | near sign err convex | far sign err | eik band |')
	print('|---|---|---|---|---|---|---|---|---|')
	for k, lab in METHODS:
		g = r['grid'].get(k)
		if g is None:
			continue
		print(f'| {lab} | {f(g["mae"])} | {f(g["saddle"]["band_mae"])} | {f(g["convex"]["band_mae"])} | {f(g["saddle"]["band_max"], 3)} | {pct(g["saddle"]["band_sign_err"])} | {pct(g["convex"]["band_sign_err"])} | {pct(g["sign_err"])} | {f(g["eik_band"], 3)} |')
	print()
	print('| zero level set (MC 128³) | mean sym. dist. saddle | convex | Hausdorff saddle | convex | recon→truth max (all) |')
	print('|---|---|---|---|---|---|')
	for k, lab in METHODS:
		ls = r['level_set'].get(k)
		if ls is None or ls.get('empty'):
			continue
		print(f'| {lab} | {f(ls["saddle"]["mean_sym"])} | {f(ls["convex"]["mean_sym"])} | {f(ls["saddle"]["hausdorff"], 3)} | {f(ls["convex"]["hausdorff"], 3)} | {f(ls["all"]["recon_to_truth"].get("max"), 3)} |')
	print()

print('## Curvatures at the input points (authors\' network coefficients vs truth; convex positive, × σ)\n')
print('| case | true saddle frac | predicted saddle frac | true (κmin σ, κmax σ) at saddles | predicted | true H σ | predicted H σ | anisotropy true / pred (saddles) |')
print('|---|---|---|---|---|---|---|---|')
for case, r in S.items():
	k = r['kappa']['network']
	print(f'| {case} | {pct(k["true_saddle_frac"])} | {pct(k["pred_saddle_frac"])} | ({f(k.get("true_kmin_sigma_med_saddle"), 2)}, {f(k.get("true_kmax_sigma_med_saddle"), 2)}) | ({f(k.get("pred_kmin_sigma_med_saddle"), 2)}, {f(k.get("pred_kmax_sigma_med_saddle"), 2)}) | {f(k.get("true_H_sigma_med_saddle"), 2)} | {f(k.get("pred_H_sigma_med_saddle"), 2)} | {f(k.get("true_anisotropy_median_saddle"), 2)} / {f(k.get("pred_anisotropy_median_saddle"), 2)} |')
print()
print('## The tori (authors\' program): inside-out at the contact point (∇g·n < 0), spindle fraction, median R/|r|\n')
print('| case | N inside-out (saddle pts) | N spindle | N median R/r | E-ref inside-out (saddle pts) | E-fix inside-out (saddle pts) | E spindle (saddle pts) |')
print('|---|---|---|---|---|---|---|')
for case, r in S.items():
	t = r['tori']
	print(f'| {case} | {pct(t["N"]["saddle"]["inside_out"])} | {pct(t["N"]["all"]["spindle"])} | {f(t["N"]["all"]["R_over_r_median"], 2)} | {pct(t["Eref"]["saddle"]["inside_out"])} | {pct(t["Efix"]["saddle"]["inside_out"])} | {pct(t["Eref"]["saddle"]["spindle"])} |')
print()
print('## Cross-check: our crate vs the authors\' program on the same inputs (max |Δφ| over the grid)\n')
keys = list(next(iter(S.values()))['crosscheck'].keys())
print('| case | ' + ' | '.join(keys) + ' |')
print('|---|' + '---|' * len(keys))
for case, r in S.items():
	print(f'| {case} | ' + ' | '.join(f'{r["crosscheck"].get(k, float("nan")):.1e}' for k in keys) + ' |')
print()
print('exact-coefficient fit_diagnostics vs true κ (max relative error):', ', '.join(f'{c}: {r["kappa"]["exact_fit_diagnostics_vs_truth"]["max_rel"]:.1e}' for c, r in S.items()))
