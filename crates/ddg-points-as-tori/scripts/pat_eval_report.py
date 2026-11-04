# /// script
# requires-python = ">=3.11"
# dependencies = ["numpy", "matplotlib", "pillow"]
# ///
"""Aggregate WP14's per-case results into the tables of EVALUATION.md / COMPARISON.md and draw the
figures (SPEC.md §16).

    uv run --no-project --with numpy --with matplotlib --with pillow python -I \
        crates/ddg-points-as-tori/scripts/pat_eval_report.py runs_tmp/pat-wp14/cases runs_tmp/pat-wp14/report

Reads `<cases>/<dataset>/<case>/{result,continuity,orient}.json` and `<cases>/density*/density.json`;
writes `<out>/tables.md` (all tables) and PNGs to `<out>/fig/`.
"""

import glob
import json
import math
import os
import sys

import numpy as np

METHODS = [
    ('pat', 'PAT (32-NN, default)'),
    ('pat-all', 'PAT (all points)'),
    ('pat-paper', 'PAT (paper λ, radius)'),
    ('pat-taper', 'PAT (tapered 32-NN, ours)'),
    ('pat-geo', 'PAT Geometric (32-NN)'),
    ('pat-geo-all', 'PAT Geometric (all)'),
    ('authors', "authors' program (accel)"),
    ('authors-all', "authors' program (all)"),
    ('shc', 'SHC'),
    ('sspd', 'SSPD'),
    ('gwn-nearest', 'GWN-sign + distance (ours)'),
    ('rwn-nearest', 'RWN(0.02)-sign + distance'),
    ('lse-gwn', 'LSE + GWN-sign'),
    ('spsr', 'SPSR mesh SDF'),
]
DATASETS = [('abc', 'ABC (100)'), ('blobs', 'blobs (100)'), ('thingi', 'Thingi10K (45)'), ('analytic', 'analytic (70)'), ('ref', "reference data/ (3)")]


def load(cases, ds, fname='result.json'):
    out = []
    for p in sorted(glob.glob(os.path.join(cases, ds, '*', fname))):
        try:
            out.append(json.load(open(p)))
        except Exception as e:  # a broken file is reported, not fatal
            print(f'skip {p}: {e}', file=sys.stderr)
    return out


def m(r, meth, key='mae'):
    x = r['methods'].get(meth, {}).get('m', {}).get(key)
    return float('nan') if x is None else float(x)


def fmt(x, d=4):
    if x is None or (isinstance(x, float) and not math.isfinite(x)):
        return '–'
    return f'{x:.{d}f}'


def table(head, rows):
    s = '| ' + ' | '.join(head) + ' |\n|' + '|'.join('---' for _ in head) + '|\n'
    for r in rows:
        s += '| ' + ' | '.join(str(c) for c in r) + ' |\n'
    return s


def headline(res, label, methods=METHODS):
    rows = []
    for key, name in methods:
        rs = [r for r in res if key in r['methods'] and 'm' in r['methods'][key]]
        if not rs:
            continue
        mae = np.array([m(r, key) for r in rs])
        band = np.array([m(r, key, 'band_mae') for r in rs])
        eik = np.array([m(r, key, 'eik') for r in rs])
        se = np.array([m(r, key, 'sign_err') for r in rs])
        cat = int((se > 0.01).sum())
        big = int((mae > 0.05).sum())
        sem = mae.std(ddof=1) / math.sqrt(len(mae)) if len(mae) > 1 else float('nan')
        rows.append([name, len(rs), f'{np.nanmean(mae):.4f} ± {sem:.4f}', fmt(np.nanmedian(mae)), fmt(np.nanmean(band)), fmt(np.nanmean(eik), 3), cat, big])
    return f'**{label}**\n\n' + table(['method', 'shapes', 'grid MAE (mean ± s.e.)', 'median', 'band MAE', 'eikonality', 'wrong sign > 1 %', 'MAE > 0.05'], rows)


def pairwise(res, a, b):
    xs = [(m(r, a), m(r, b)) for r in res if a in r['methods'] and b in r['methods']]
    if not xs:
        return None
    w = sum(1 for x, y in xs if x < y)
    return w, len(xs)


def main():
    cases, out = sys.argv[1], sys.argv[2]
    os.makedirs(os.path.join(out, 'fig'), exist_ok=True)
    T = []
    data = {ds: load(cases, ds) for ds, _ in DATASETS}
    # exact duplicates of training shapes among the validation shapes (scripts/pat_eval_data.py)
    ov = os.path.join(cases, 'src', 'overlap.json')
    if os.path.isfile(ov):
        per = json.load(open(ov)).get('per_shape', {})
        for r in data['abc']:
            if per.get(r['name'], {}).get('exact'):
                r['training'] = 'duplicate-of-training'
    # WP10's retrained weights, when the case stage ran with --retrained label=path
    extra = sorted({k for rs in data.values() for r in rs for k in r['methods'] if k.startswith('retrained')})
    METHODS.extend((k, f'PAT {k} weights') for k in extra if k not in dict(METHODS))
    # ---------------------------------------------------------------- headline
    T.append('## Headline accuracy (512 points, 64³ grid over [−1, 1]³)\n')
    for ds, label in DATASETS:
        if data[ds]:
            T.append(headline(data[ds], label))
    allres = sum((data[d] for d in ('abc', 'blobs', 'thingi', 'analytic')), [])
    T.append(headline(allres, f'all four datasets ({len(allres)} shapes)'))
    # PAT vs GWN-sign + distance, per shape
    rows = []
    for ds, label in DATASETS[:4]:
        for a, b in [('pat', 'gwn-nearest'), ('pat', 'lse-gwn'), ('pat', 'spsr'), ('pat', 'shc'), ('pat-all', 'pat')]:
            pw = pairwise(data[ds], a, b)
            if pw:
                rows.append([label, f'{a} < {b}', f'{pw[0]}/{pw[1]}'])
    T.append('**Per-shape wins (grid MAE)**\n\n' + table(['dataset', 'comparison', 'shapes where the first is lower'], rows))
    # ABC by training status
    for st in ['unseen', 'duplicate-of-training', 'seen']:
        rs = [r for r in data['abc'] if r.get('training') == st]
        if st == 'unseen':
            rs = rs
        if rs:
            T.append(headline(rs, f'ABC, {st} ({len(rs)})', [x for x in METHODS if x[0] in ('pat', 'pat-all', 'gwn-nearest', 'lse-gwn', 'shc', 'spsr', 'authors')]))
    # SHC degenerate counts, saddles predicted
    for ds, label in DATASETS:
        if data[ds]:
            z = sum(r['extra'].get('shc_zero_w', 0) for r in data[ds])
            s = sum(r['extra'].get('saddles', 0) for r in data[ds])
            sh = sum(1 for r in data[ds] if r['extra'].get('saddles', 0) > 0)
            T.append(f'- {label}: SHC grid points with w = 0 (set to 0): {z}; saddle tori predicted: {s} on {sh} shapes\n')
    # ---------------------------------------------------------------- sign rule
    T.append('\n## Sign rule (Reference vs Geometric)\n')
    rows = []
    for ds in ['saddle', 'saddle2048', 'analyticsad', 'analyticsad2048']:
        for r in load(cases, ds):
            rows.append([ds, r['name'], r['n'], r['extra'].get('saddles', 0), fmt(m(r, 'pat'), 5), fmt(m(r, 'pat-geo'), 5), fmt(m(r, 'pat', 'band_mae'), 5), fmt(m(r, 'pat-geo', 'band_mae'), 5),
                         fmt(m(r, 'pat-all'), 5), fmt(m(r, 'pat-geo-all'), 5), fmt(m(r, 'gwn-nearest'), 5)])
    T.append(table(['set', 'shape', 'n', 'saddle tori', 'MAE Ref', 'MAE Geo', 'band Ref', 'band Geo', 'MAE Ref all', 'MAE Geo all', 'GWN-sign+dist'], rows))
    rows = []
    for ds, label in DATASETS + [('saddle', 'saddle set'), ('saddle2048', 'saddle set 2048'), ('analyticsad', 'App. E saddle 13'), ('analyticsad2048', 'App. E saddle 13, 2048')]:
        rs = data.get(ds) or load(cases, ds)
        if not rs:
            continue
        d = np.array([m(r, 'pat-geo') - m(r, 'pat') for r in rs])
        da = np.array([m(r, 'pat-geo-all') - m(r, 'pat-all') for r in rs])
        rows.append([label, len(rs), sum(r['extra'].get('saddles', 0) for r in rs), int((d < 0).sum()), int((d > 0).sum()), f'{d.mean():+.2e}', f'{d.min():+.2e} / {d.max():+.2e}', f'{da.mean():+.2e}'])
    T.append(table(['set', 'shapes', 'saddle tori', 'Geo better', 'Ref better', 'mean ΔMAE (Geo − Ref, 32-NN)', 'min / max Δ', 'mean Δ (all)'], rows))
    # ---------------------------------------------------------------- §16.3 sweep
    T.append('\n## §16.3 choices (first 50 ABC + 50 blobs)\n')
    sub = data['abc'][:50] + data['blobs'][:50]
    sw = [('pat', 'C 64, k 32 (default)'), ('pat-c32', 'C 32'), ('pat-c128', 'C 128'), ('pat-k16', 'k 16'), ('pat-k64', 'k 64'), ('pat-paper', 'paper λ = 10³/D, radius 2C/λ'),
          ('pat-all', 'all, C 64'), ('pat-all-c32', 'all, C 32'), ('pat-all-c128', 'all, C 128'), ('pat-taper', 'tapered 32-NN')]
    T.append(headline([r for r in sub if 'pat-c32' in r['methods']], 'sweep', sw))
    # ---------------------------------------------------------------- black box
    T.append('\n## Black box: ours vs the authors\' program\n')
    rows = []
    keys = ['coeff_rel', 'tori_center_rel', 'tori_radius_rel', 'tori_axis_sin', 'sdf_accel', 'sdf_all']
    for ds in ['abc', 'blobs', 'thingi', 'analytic', 'ref', 'saddle', 'analyticsad', 'saddle2048', 'analyticsad2048', 'reffull']:
        rs = [r for r in (data.get(ds) or load(cases, ds)) if 'black_box' in r]
        if not rs:
            continue
        row = [ds, len(rs)]
        for k in keys:
            v = [r['black_box'][k]['max'] for r in rs if k in r['black_box'] and r['black_box'][k]['max'] is not None]
            med = [r['black_box'][k]['median'] for r in rs if k in r['black_box'] and r['black_box'][k]['median'] is not None]
            row.append(f'{max(v):.1e} ({np.median(med):.0e})' if v else '–')
        row.append(sum(r['black_box'].get('knn_rows_differ', 0) for r in rs))
        row.append(sum(r['black_box'].get('tori_sign_flips', 0) for r in rs))
        row.append(sum(r['black_box'].get('tori_flat_skipped', 0) for r in rs))
        rows.append(row)
    T.append(table(['set', 'clouds', 'coeff |Δ|/(1+|a|) max (median)', 'center', 'radii', 'axis sin', 'SDF 32-NN |Δ|', 'SDF all |Δ|', 'kNN rows differing', 'minor-sign flips', 'flat tori skipped'], rows))
    # worst SDF cases
    worst = []
    for ds in ['abc', 'blobs', 'thingi', 'analytic', 'ref', 'saddle', 'analyticsad', 'saddle2048', 'analyticsad2048', 'reffull']:
        for r in (data.get(ds) or load(cases, ds)):
            bb = r.get('black_box', {})
            if 'sdf_accel' in bb:
                worst.append((bb['sdf_accel']['max'] or 0, ds, r['name'], bb['sdf_accel'].get('p999'), bb.get('sdf_all', {}).get('max'), bb.get('coeff_rel', {}).get('max')))
    worst.sort(reverse=True)
    T.append('**Largest SDF differences**\n\n' + table(['set', 'cloud', 'SDF 32-NN max', 'p99.9', 'SDF all max', 'coeff max'], [[w[1], w[2], f'{w[0]:.2e}', f'{w[3]:.1e}', f'{w[4]:.1e}' if w[4] else '–', f'{w[5]:.1e}' if w[5] else '–'] for w in worst[:8]]))
    # timings (512 points)
    tt = []
    for ds in ['abc', 'blobs', 'thingi', 'analytic']:
        for r in data[ds]:
            at = r.get('black_box', {}).get('authors_timing')
            if at:
                tt.append((r['pat_build_s'], r['methods']['pat']['seconds'], r['methods']['pat-all']['seconds'], at['constructor_warm_s'], at['grid64_accel_s'], at['grid64_noaccel_s']))
    if tt:
        a = np.array(tt)
        T.append('**Timings, 512 points, CPU only, shared machine (median over clouds; ours 4 threads, theirs OMP 4 + XLA 4 intra-op)**\n\n' + table(
            ['stage', 'ours (s)', "authors' (s)", 'ratio'],
            [['build (kNN + network + fit)', f'{np.median(a[:, 0]):.3f}', f'{np.median(a[:, 3]):.3f}', f'{np.median(a[:, 3]) / np.median(a[:, 0]):.1f}×'],
             ['64³ grid, 32-NN', f'{np.median(a[:, 1]):.3f}', f'{np.median(a[:, 4]):.3f}', f'{np.median(a[:, 4]) / np.median(a[:, 1]):.1f}×'],
             ['64³ grid, all points', f'{np.median(a[:, 2]):.3f}', f'{np.median(a[:, 5]):.3f}', f'{np.median(a[:, 5]) / np.median(a[:, 2]):.1f}×']]))
    # ---------------------------------------------------------------- continuity
    T.append('\n## Continuity per evaluation rule\n')
    cont = []
    for p in sorted(glob.glob(os.path.join(cases, '*', '*', 'continuity.json'))):
        try:
            cont.append(json.load(open(p)))
        except Exception:
            pass
    rules = ['knn32', 'paper', 'all', 'taper']
    for res in ['64', '128', '256']:
        rows = []
        # only clouds with a truth and all four rules at this resolution, so the rows compare like with like
        cl = [c for c in cont if res in c['res'] and 'truth_mc_components' in c['res'][res] and all(r in c['res'][res] for r in rules)]
        for rule in rules:
            vals = [c['res'][res][rule] for c in cl]
            if not vals:
                continue
            mx = [v['jumps']['max'] for v in vals]
            p999 = [v['jumps']['p999'] for v in vals]
            o5 = [v['jumps']['over5'] / v['jumps']['edges'] * 1e3 for v in vals]
            comps = [v['mc_components'] for v in vals]
            det = [v['mc_detached'] for v in vals]
            rows.append([rule, len(vals), fmt(np.median(mx), 1), fmt(max(mx), 1), fmt(np.median(p999), 2), fmt(np.mean(o5), 2),
                         fmt(np.mean([v['m']['mae'] for v in vals]), 4), fmt(np.mean([v['m']['band_mae'] for v in vals]), 4),
                         fmt(np.mean([v['m']['eik'] for v in vals]), 3), sum(comps), sum(c['res'][res]['truth_mc_components'] for c in cl), sum(det)])
        T.append(f'**{res}³**\n\n' + table(['rule', 'clouds', 'max |Δφ|/h (median)', 'max (worst)', 'p99.9 (median)', 'edges > 5 (‰, mean)', 'MAE', 'band MAE', 'eikonality', 'MC components (sum)', 'truth components (sum)', 'detached (sum)'], rows))
    # per cloud components at 256 vs truth
    rows = []
    for c in cont:
        top = max(c['res'].keys(), key=int)
        r = c['res'][top]
        row = [c.get('dataset'), c.get('name'), c['n'], top, r.get('truth_mc_components', c.get('truth_components', '–'))]
        for rule in rules:
            v = r.get(rule)
            has_truth = 'truth_mc_components' in r or 'truth_components' in c
            row.append((f"{v['mc_components']} ({v['mc_detached']})" if has_truth else f"{v['mc_components']}") if v else '–')
        for rule in rules:
            v = r.get(rule)
            row.append(fmt(v['jumps']['max'], 1) if v else '–')
        rows.append(row)
    T.append('**Per cloud, finest grid: zero-level-set components (detached: no vertex within 2h of the truth surface) and max slope**\n\n' + table(
        ['set', 'cloud', 'n', 'res', 'truth comps'] + [f'comps {x}' for x in rules] + [f'slope {x}' for x in rules], rows))
    # mechanism
    rows = []
    for c in cont:
        top = max(c['res'].keys(), key=int)
        for rule in ['knn32', 'paper']:
            mm = c['res'][top].get(rule, {}).get('mechanism')
            if mm:
                rows.append([c.get('name'), rule, mm['edges'], fmt(mm['continuous_frac'], 2), fmt(mm['set_change_frac'], 2), fmt(mm['far_field_frac'], 2), fmt(mm['rho_min_over_sigma_median'], 2), fmt(mm['truth_dist_median'], 3), fmt(mm['w_entering_median'], 3)])
    T.append('**Mechanism at the largest jumps (finest grid, top edges by slope)**\n\n' + table(
        ['cloud', 'rule', 'edges', 'continuous (no jump after bisection)', 'neighbor set changes', 'far field (ρmin/σ > 1)', 'ρmin/σ median', '|φ_true| median', 'weight share entering'], rows))
    # ---------------------------------------------------------------- orientation
    ors = sum((load(cases, d, 'orient.json') for d in ['orient', 'abc2048', 'blobs2048', 'thingi2048']), [])
    if ors:
        rows = [[os.path.basename(o['case']), o['n'], o['flipped'], fmt(100 * o['it10']['residual'], 2), fmt(100 * o['it30']['residual'], 2), o['it10']['converged'], o['it30']['converged'], len(o['it30']['flips'])] for o in ors]
        r10 = np.array([o['it10']['residual'] for o in ors])
        r30 = np.array([o['it30']['residual'] for o in ors])
        rows.append(['**mean**', '', '', fmt(100 * r10.mean(), 2), fmt(100 * r30.mean(), 2), '', '', ''])
        rows.append(['**median**', '', '', fmt(100 * np.median(r10), 2), fmt(100 * np.median(r30), 2), '', '', ''])
        T.append('\n## Orientation (30 % flipped)\n\n' + table(['cloud', 'n', 'flipped', 'wrong after 10 (%)', 'wrong after 30 (%)', 'converged 10', 'converged 30', 'iterations used (30)'], rows))
    # ---------------------------------------------------------------- density
    for p in sorted(glob.glob(os.path.join(cases, '..', 'density*', 'density.json'))):
        d = json.load(open(p))
        rows = [[r['n'], fmt(r['m']['mae'], 5), fmt(r['m']['band_mae'], 5), fmt(r['interior_mae'], 5), fmt(r['exterior_mae'], 5), fmt(r['m']['eik'], 3), fmt(r['m']['sign_err'] * 100, 2),
                 fmt(r['major_over_sigma_median'], 1), fmt(r['minor_over_sigma_median'], 1), fmt(r['minor_over_sigma_p10'], 2), fmt(100 * r['negative_minor_frac'], 1), r.get('mc_components', '–'), fmt(r['pat_build_s'], 1)] for r in d['rows']]
        T.append(f"\n## Density: {os.path.basename(os.path.dirname(p))} ({d['mesh']}, {d['grid']}³)\n\n" + table(
            ['points', 'MAE', 'band', 'interior', 'exterior', 'eikonality', 'wrong sign %', 'R/σ median', '|r|/σ median', '|r|/σ p10', 'r < 0 %', 'MC comps', 'build s'], rows))
    open(os.path.join(out, 'tables.md'), 'w').write('\n'.join(T))
    print('\n'.join(T))
    figures(cases, out, data, cont)
    # 8-bit palette PNGs (small enough to commit)
    from PIL import Image
    for p in glob.glob(os.path.join(out, 'fig', '*.png')):
        Image.open(p).convert('RGB').quantize(colors=96, method=Image.Quantize.MEDIANCUT).save(p, optimize=True)
    return 0


def figures(cases, out, data, cont):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    blue, orange, aqua, violet, gray = '#2a78d6', '#eb6834', '#1baf7a', '#6250d6', '#8a8984'
    fam = {'pat': blue, 'pat-all': blue, 'pat-paper': blue, 'pat-taper': blue, 'authors': aqua, 'authors-all': aqua, 'shc': orange, 'sspd': orange,
           'gwn-nearest': gray, 'rwn-nearest': gray, 'lse-gwn': gray, 'spsr': violet}
    keys = ['pat', 'pat-all', 'pat-paper', 'authors', 'shc', 'sspd', 'gwn-nearest', 'lse-gwn', 'spsr']
    names = dict(METHODS)
    # Fig. 18 analogue: per-dataset mean grid MAE (top) and band MAE (bottom), small multiples
    fig, axs = plt.subplots(2, 4, figsize=(14, 6.4), sharey=True)
    for row, (key_m, xl) in enumerate([('mae', 'grid MAE (64³)'), ('band_mae', 'band MAE (|φ_true| < 0.05)')]):
        for ax, (ds, label) in zip(axs[row], DATASETS[:4]):
            rs = data[ds]
            vals = [np.nanmean([m(r, k, key_m) for r in rs if k in r['methods']]) if rs else np.nan for k in keys]
            ax.barh(range(len(keys)), vals, color=[fam[k] for k in keys], height=0.7)
            for i, v in enumerate(vals):
                if math.isfinite(v):
                    ax.text(v, i, f' {v:.3f}', va='center', fontsize=7, color='#52514e')
            if row == 0:
                ax.set_title(label, fontsize=10)
            ax.set_yticks(range(len(keys)))
            ax.set_yticklabels([names[k] for k in keys], fontsize=8)
            ax.set_xlabel(xl, fontsize=8)
            ax.spines[['top', 'right']].set_visible(False)
            ax.tick_params(axis='x', labelsize=7)
            ax.set_xlim(0, max(v for v in vals if math.isfinite(v)) * 1.35)
    axs[0, 0].invert_yaxis()
    fig.suptitle('512 points; mean over shapes (blue: PAT, aqua: authors\' program, orange: paper baselines, gray: winding-number baselines, violet: SPSR)', fontsize=9)
    fig.tight_layout()
    fig.savefig(os.path.join(out, 'fig', 'accuracy.png'), dpi=100)
    plt.close(fig)
    # Fig. 17 analogue: z-slices of truth and methods for a few shapes
    picks = []
    for ds in ['abc', 'blobs', 'thingi', 'analytic']:
        rs = data[ds]
        if rs:
            r = rs[len(rs) // 3]
            picks.append((ds, r))
    show = ['truth', 'pat', 'pat-all', 'shc', 'sspd', 'gwn-nearest', 'spsr']
    if picks:
        fig, axs = plt.subplots(len(picks), len(show), figsize=(2.1 * len(show), 2.2 * len(picks)))
        for i, (ds, r) in enumerate(picks):
            sl = np.load(os.path.join(cases, ds, r['name'], 'slices.npy'))
            nm = r['slice_names']
            for j, k in enumerate(show):
                ax = axs[i, j]
                ax.set_xticks([])
                ax.set_yticks([])
                if k not in nm:
                    ax.axis('off')
                    continue
                z = sl[nm.index(k), 2].T
                ax.imshow(z, origin='lower', cmap='RdBu_r', vmin=-0.5, vmax=0.5, extent=(-1, 1, -1, 1))
                ax.contour(np.linspace(-1, 1, 64), np.linspace(-1, 1, 64), z, levels=[0], colors='k', linewidths=0.8)
                if 'truth' in nm and k != 'truth':
                    ax.contour(np.linspace(-1, 1, 64), np.linspace(-1, 1, 64), sl[nm.index('truth'), 2].T, levels=[0], colors='#eda100', linewidths=0.6, linestyles='dashed')
                if i == 0:
                    ax.set_title(names.get(k, k), fontsize=7)
                if j == 0:
                    ax.set_ylabel(f"{ds}\n{r['name']}", fontsize=7)
        fig.suptitle('z ≈ 0 slices (red outside, blue inside, ±0.5; black: zero level, dashed: truth)', fontsize=8)
        fig.tight_layout()
        fig.savefig(os.path.join(out, 'fig', 'slices.png'), dpi=60)
        plt.close(fig)
    # continuity: max slope vs resolution per rule (median over clouds)
    if cont:
        fig, ax = plt.subplots(figsize=(5, 3.4))
        for rule, col in [('knn32', blue), ('paper', orange), ('all', aqua), ('taper', violet)]:
            xs, ys = [], []
            for res in ['64', '128', '256']:
                v = [c['res'][res][rule]['jumps']['max'] for c in cont if res in c['res'] and rule in c['res'][res]]
                if v:
                    xs.append(int(res))
                    ys.append(np.median(v))
            ax.plot(xs, ys, marker='o', color=col, lw=2, ms=6, label=rule)
            if ys:
                ax.annotate(rule, (xs[-1], ys[-1]), textcoords='offset points', xytext=(6, 0), fontsize=8, color='#52514e', va='center')
        ax.set_xscale('log', base=2)
        ax.set_yscale('log')
        ax.set_xlabel('grid points per axis')
        ax.set_ylabel('max |Δφ|/h (median over clouds)')
        ax.legend(fontsize=8, frameon=False)
        ax.spines[['top', 'right']].set_visible(False)
        fig.tight_layout()
        fig.savefig(os.path.join(out, 'fig', 'continuity.png'), dpi=110)
        plt.close(fig)
    # density slices
    for p in sorted(glob.glob(os.path.join(cases, '..', 'density*', 'density.json'))):
        d = json.load(open(p))
        sl = np.load(os.path.join(os.path.dirname(p), 'slices.npy'))
        k = sl.shape[0]
        fig, axs = plt.subplots(1, k, figsize=(2.6 * k, 2.8))
        ax1 = np.linspace(-1, 1, sl.shape[1])
        for i in range(k):
            ax = axs[i]
            ax.imshow(sl[i], origin='lower', cmap='RdBu_r', vmin=-0.4, vmax=0.4, extent=(-1, 1, -1, 1))
            ax.contour(ax1, ax1, sl[i], levels=[0], colors='k', linewidths=0.7)
            ax.set_xticks([])
            ax.set_yticks([])
            ax.set_title('truth' if i == 0 else f"{d['rows'][i - 1]['n']:,} points\nMAE {d['rows'][i - 1]['m']['mae']:.4f}", fontsize=8)
        fig.tight_layout()
        fig.savefig(os.path.join(out, 'fig', f"{os.path.basename(os.path.dirname(p))}.png"), dpi=70)
        plt.close(fig)


if __name__ == '__main__':
    sys.exit(main())
