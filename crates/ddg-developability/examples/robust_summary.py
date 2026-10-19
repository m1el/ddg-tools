# /// script
# dependencies = ["numpy"]
# ///
"""Summary of flow runs (SPEC §12.4, SPEC D31/D32).

Usage: uv run robust_summary.py DIR [DIR ...]
Reads every DIR/*.log written by `develop` (RESULT line) and the matching .obj (if any); prints a
markdown table: run, status, E_final, seconds, inadmissible trials, smallest corner angle (rad),
faces below 1 degree, smallest fold gap (pi minus the angle between normals of adjacent faces).
"""
import glob, os, sys
from collections import defaultdict
import numpy as np


def read_obj(fn):
    V, F = [], []
    for l in open(fn):
        t = l.split()
        if not t:
            continue
        if t[0] == 'v':
            V.append([float(x) for x in t[1:4]])
        elif t[0] == 'f':
            F.append([int(x.split('/')[0]) - 1 for x in t[1:4]])
    return np.array(V), np.array(F)


def stats(fn):
    V, F = read_obj(fn)
    P = V[F]

    def ang(a, b, c):
        u, v = b - a, c - a
        return np.arctan2(np.linalg.norm(np.cross(u, v), axis=1), (u * v).sum(1))

    A = np.stack([ang(P[:, 0], P[:, 1], P[:, 2]), ang(P[:, 1], P[:, 2], P[:, 0]), ang(P[:, 2], P[:, 0], P[:, 1])], 1)
    m = A.min(1)
    n = np.cross(P[:, 1] - P[:, 0], P[:, 2] - P[:, 0])
    nl = np.linalg.norm(n, axis=1)
    N = n / np.maximum(nl, 1e-300)[:, None]
    ef = defaultdict(list)
    for f, t in enumerate(F):
        for k in range(3):
            ef[tuple(sorted((t[k], t[(k + 1) % 3])))].append(f)
    gap = np.pi
    for fs in ef.values():
        if len(fs) == 2:
            c = float(N[fs[0]] @ N[fs[1]])
            s = float(np.linalg.norm(np.cross(N[fs[0]], N[fs[1]])))
            gap = min(gap, np.pi - np.arctan2(s, c))
    return len(F), m.min(), int((m < np.pi / 180).sum()), gap


def main():
    print("| run | faces | status | E_final | s | inadmissible | min angle | faces < 1 deg | min fold gap |")
    print("|---|---|---|---|---|---|---|---|---|")
    for d in sys.argv[1:]:
        for log in sorted(glob.glob(os.path.join(d, '*.log'))):
            name = os.path.basename(log)[:-4]
            res = {}
            err = None
            for l in open(log):
                if l.startswith('RESULT '):
                    for kv in l.split()[1:]:
                        if '=' in kv:
                            k, v = kv.split('=', 1)
                            res[k] = v
                if l.startswith('error:') or l.startswith('cut failed'):
                    err = l.strip()[:90]
            obj = log[:-4] + '.obj'
            st = stats(obj) if os.path.exists(obj) else None
            if 'error' in res and 'status' not in res:
                print(f"| {name} | - | ERROR {err} | | | | | | |")
                continue
            row = [name, str(st[0]) if st else res.get('faces_out', '-'), res.get('status', '-'), res.get('e_final', '-'),
                   res.get('seconds', '-'), res.get('inadmissible_trials', '-')]
            row += [f"{st[1]:.2e}", str(st[2]), f"{st[3]:.2e}"] if st else ['-', '-', '-']
            print('| ' + ' | '.join(row) + ' |')


main()
