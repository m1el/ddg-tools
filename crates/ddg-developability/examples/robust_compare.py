# /// script
# dependencies = []
# ///
"""Before/after table of the robustness study (NOTES-flow.md): `robust_compare.py OLD NEW`, where OLD
is the directory of the cut study (flow/ and cut/ with the covariance / 1e-3 / 0.1 cut) and NEW that of
`robust_study.sh`. Prints a markdown table."""
import glob, os, sys, re

def res(fn):
    d = {}
    if not os.path.exists(fn):
        return d
    for l in open(fn):
        if l.startswith('RESULT '):
            for kv in l.split()[1:]:
                if '=' in kv:
                    k, v = kv.split('=', 1); d[k] = v
        if l.startswith('error:') or l.startswith('cut failed'):
            d['err'] = l.strip()[:60]
    return d

def cell(d, cutd):
    if 'e_final' not in d:
        return ['ERR ' + d.get('err', '?')[:40], '', '', '', '']
    st = d['status'].replace('linesearch_failed', 'lsf').replace('max_iter', 'max')
    cut = 'cut fails' if cutd.get('cut') != 'ok' else cutd.get('flipped_uv', '?')
    return [d['e_final'], d['iterations'], d['seconds'], d.get('faces_out', ''), cut + ' / ' + st]

old, new = sys.argv[1:3]
print('| flow | E_final old | E_final new | iterations old | new | s old | s new | faces old | new | flipped uv / status old | flipped uv / status new | new min angle | new bad faces |')
print('|---|---|---|---|---|---|---|---|---|---|---|---|---|')
for f in sorted(glob.glob(os.path.join(new, 'flow', '*.log'))):
    n = os.path.basename(f)[:-4]
    o = cell(res(os.path.join(old, 'flow', n + '.log')), res(os.path.join(old, 'cut', n + '_covariance_1e-3_0.1.log')))
    w = cell(res(f), res(os.path.join(new, 'cut', n + '.log')))
    nd = res(f)
    print('| ' + ' | '.join([n, o[0], w[0], o[1], w[1], o[2], w[2], o[3], w[3], o[4], w[4], nd.get('min_angle', ''), nd.get('bad_faces', '')]) + ' |')
