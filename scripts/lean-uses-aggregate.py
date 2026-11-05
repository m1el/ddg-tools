#!/usr/bin/env python3
"""Aggregate the TSV written by scripts/lean-uses.lean.

Input lines: module <TAB> declaration <TAB> used-module <TAB> used-constant.

usage: lean-uses-aggregate.py [uses.tsv] [--by module|part|depth:N] [--parts-map FILE]
                              [--keep-artifacts] [--min N] [--from PFX] [--to PFX] [--decls]

  --by module     one row per (module, used-module) pair (default)
  --by depth:N    group modules by their first N name components (DdgProofs.Cuts.Flow -> depth 2:
                  DdgProofs.Cuts)
  --by part       group by the longest matching prefix in --parts-map FILE, a TSV
                  `module-prefix <TAB> part` (lines starting with # ignored)
  --keep-artifacts  keep .eq_N / .match_N / _proof_N / splitter declarations (dropped by default,
                  on both the using and the used side)
  --min N         print only pairs with at least N uses
  --from/--to PFX only rows whose using / used module starts with PFX
  --decls         with --by module: also list the distinct (declaration -> used constant) pairs
Output: count <TAB> from <TAB> to, sorted by count descending; a summary goes to stderr.
Default input: runs_tmp/uses.tsv relative to the repository root.
"""
import argparse, os, re, sys
from collections import Counter, defaultdict

ARTIFACT = re.compile(r'(\.eq_\d+|\.eq_def|\.match_\d+|_proof_\d+|\.splitter|\.sunfold|\._unary|\.proof_\d+)(\.|$)|\.splitter')

def is_artifact(name):
    return ARTIFACT.search(name) is not None

def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    ap.add_argument('tsv', nargs='?', default=os.path.join(root, 'runs_tmp', 'uses.tsv'))
    ap.add_argument('--by', default='module')
    ap.add_argument('--parts-map')
    ap.add_argument('--keep-artifacts', action='store_true')
    ap.add_argument('--min', type=int, default=1)
    ap.add_argument('--from', dest='frm', default='')
    ap.add_argument('--to', default='')
    ap.add_argument('--decls', action='store_true')
    a = ap.parse_args()

    pmap = []
    if a.by == 'part':
        if not a.parts_map:
            ap.error('--by part needs --parts-map')
        for line in open(a.parts_map):
            line = line.rstrip('\n')
            if line and not line.startswith('#'):
                pre, part = line.split('\t')[:2]
                pmap.append((pre, part))
        pmap.sort(key=lambda x: -len(x[0]))

    def group(m):
        if a.by == 'module':
            return m
        if a.by.startswith('depth:'):
            n = int(a.by.split(':')[1])
            return '.'.join(m.split('.')[:n])
        if a.by == 'part':
            for pre, part in pmap:
                if m == pre or m.startswith(pre + '.'):
                    return part
            return '(unmapped)'
        ap.error('bad --by')

    cnt = Counter(); decls = defaultdict(set)
    total = dropped = 0
    for line in open(a.tsv):
        p = line.rstrip('\n').split('\t')
        if len(p) != 4:
            continue
        m, d, m2, u = p
        total += 1
        if not a.keep_artifacts and (is_artifact(d) or is_artifact(u)):
            dropped += 1
            continue
        if a.frm and not m.startswith(a.frm): continue
        if a.to and not m2.startswith(a.to): continue
        g1, g2 = group(m), group(m2)
        if g1 == g2:
            continue
        cnt[(g1, g2)] += 1
        if a.decls:
            decls[(g1, g2)].add((d, u))
    for (g1, g2), c in sorted(cnt.items(), key=lambda kv: (-kv[1], kv[0])):
        if c < a.min: continue
        print(f'{c}\t{g1}\t{g2}')
        if a.decls:
            for d, u in sorted(decls[(g1, g2)]):
                print(f'\t\t{d} -> {u}')
    print(f'lean-uses-aggregate: {total} uses read, {dropped} artifact uses dropped, '
          f'{sum(cnt.values())} cross-group uses in {len(cnt)} pairs', file=sys.stderr)

main()
