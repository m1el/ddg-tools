#!/usr/bin/env python3
"""Check the part DAG of the Lean library `DDG` (LIBRARY-PLAN.md, section 4).

usage: check-lean-dag.py [--repo DIR] [--parts FILE] [--self-test]

Reads lean/ddg/parts.toml:

    # one table per part; the part's modules are DDG.<name>[.*], files lean/ddg/DDG/<name>/**.lean
    # (and lean/ddg/DDG/<name>.lean); `dir` overrides <name> when the directory differs and may be a
    # nested path such as "Delaunay/Planar" (the longest matching dir wins).
    # `deps` lists the DIRECT dependencies; the allowed set is their transitive closure.
    [parts.Geometry]
    deps = []

    [parts.Mesh]
    deps = ["Geometry"]

    [parts.Laplace]
    dir = "Laplace"          # optional
    deps = ["Mesh"]

Fails (exit 1) if
  * a file of part P imports DDG.Q (Q a different part) with Q not in the transitive closure of
    P's deps;
  * a file under lean/ddg (outside .lake) imports DdgProofs or DdgOracle (module prefix `DdgProofs`,
    `DdgOracle`);
  * a file under lean/ddg/DDG/ lies in no declared part, or imports a DDG.Q whose Q is not declared
    (top-level DDG/*.lean umbrella files are exempt: they may import every part);
  * parts.toml has a dependency on an undeclared part or a cycle.
If lean/ddg does not exist, prints a note and exits 0. `--self-test` builds synthetic trees in a
temporary directory and checks that good trees pass and bad ones fail.
"""
import argparse, os, re, sys, tempfile, shutil
try:
    import tomllib
except ImportError:  # python < 3.11
    tomllib = None

FORBIDDEN = ('DdgProofs', 'DdgOracle')
IMPORT_RE = re.compile(r'^\s*(?:(?:public|meta|private)\s+)*import\s+(.+?)\s*$')


def strip_comments(text):
    """Remove block comments (nested) and line comments; keeps line structure."""
    out, i, depth, n = [], 0, 0, len(text)
    while i < n:
        if text.startswith('/-', i):
            depth += 1; i += 2
        elif depth and text.startswith('-/', i):
            depth -= 1; i += 2
        elif depth:
            out.append('\n' if text[i] == '\n' else ' '); i += 1
        elif text.startswith('--', i):
            while i < n and text[i] != '\n': i += 1
        else:
            out.append(text[i]); i += 1
    return ''.join(out)


def imports_of(path):
    """Module names imported by a .lean file (the header only: imports precede all commands)."""
    text = strip_comments(open(path, encoding='utf-8').read())
    mods = []
    for line in text.splitlines():
        if not line.strip():
            continue
        m = IMPORT_RE.match(line)
        if m:
            mods.extend(m.group(1).split())
        elif line.strip() == 'prelude' or line.strip().startswith('module'):
            continue
        else:
            break
    return mods


def load_parts(path):
    if tomllib is None:
        sys.exit('check-lean-dag: needs Python >= 3.11 (tomllib)')
    with open(path, 'rb') as f:
        data = tomllib.load(f)
    parts = {}
    for name, cfg in data.get('parts', {}).items():
        parts[name] = {'dir': cfg.get('dir', name), 'deps': list(cfg.get('deps', []))}
    return parts


def closure(parts):
    errors = []
    for p, c in parts.items():
        for d in c['deps']:
            if d not in parts:
                errors.append(f'parts.toml: {p} depends on undeclared part {d}')
    clo, state = {}, {}

    def visit(p, stack):
        if p in clo: return clo[p]
        if state.get(p) == 1:
            errors.append('parts.toml: dependency cycle: ' + ' -> '.join(stack + [p]))
            return set()
        state[p] = 1
        s = set()
        for d in parts[p]['deps']:
            if d in parts:
                s.add(d); s |= visit(d, stack + [p])
        state[p] = 2
        clo[p] = s
        return s
    for p in parts: visit(p, [])
    return clo, errors


def match_part(segs, parts):
    """Part whose `dir` (a path like `Delaunay/Planar`, or one name) is the longest prefix of segs
    (the segments after `DDG`), or None."""
    best, blen = None, 0
    for p, c in parts.items():
        d = c['dir'].split('/')
        if segs[:len(d)] == d and len(d) > blen:
            best, blen = p, len(d)
    return best


def module_part(mod, parts):
    """The part of module `DDG.X.Y...`, or None."""
    segs = mod.split('.')
    if len(segs) < 2 or segs[0] != 'DDG':
        return None
    return match_part(segs[1:], parts)


def check(ddg_dir, parts_file):
    errors = []
    if not os.path.isfile(parts_file):
        return [f'{parts_file} is missing']
    parts = load_parts(parts_file)
    clo, errors = closure(parts)
    n_files = 0
    for dp, dns, fns in os.walk(ddg_dir):
        dns[:] = [d for d in dns if d not in ('.lake', '.git')]
        for fn in fns:
            if not fn.endswith('.lean'): continue
            path = os.path.join(dp, fn)
            rel = os.path.relpath(path, ddg_dir)
            n_files += 1
            mods = imports_of(path)
            for m in mods:
                if m.split('.')[0] in FORBIDDEN:
                    errors.append(f'{rel}: imports {m} (the library must not depend on the app or the oracle)')
            segs = rel[:-5].split(os.sep)
            if segs[0] != 'DDG' or len(segs) == 1:
                continue  # umbrella / other roots: only the forbidden-import rule applies
            me = match_part(segs[1:], parts)
            if me is None:
                errors.append(f'{rel}: not in any declared part')
                continue
            for m in mods:
                if m.split('.')[0] != 'DDG':
                    continue
                if len(m.split('.')) < 2:
                    continue  # `import DDG` (umbrella) inside a part is itself wrong
                q = module_part(m, parts)
                if q is None:
                    errors.append(f'{rel}: imports {m}, which is in no declared part')
                elif q != me and q not in clo[me]:
                    errors.append(f'{rel}: part {me} imports {m} (part {q}), outside its dependency closure '
                                  f'{sorted(clo[me])}')
            if 'DDG' in mods:
                errors.append(f'{rel}: imports the umbrella DDG')
    return errors, n_files


def run(repo, parts_override=None):
    ddg = os.path.join(repo, 'lean', 'ddg')
    if not os.path.isdir(ddg):
        print(f'check-lean-dag: {ddg} does not exist yet (before S4); nothing to check')
        return 0
    r = check(ddg, parts_override or os.path.join(ddg, 'parts.toml'))
    if isinstance(r, list):
        errs, n = r, 0
    else:
        errs, n = r
    if errs:
        for e in sorted(set(errs)): print('FAIL: ' + e, file=sys.stderr)
        print(f'check-lean-dag: {len(set(errs))} violation(s)', file=sys.stderr)
        return 1
    print(f'check-lean-dag: ok ({n} files)')
    return 0


def self_test():
    parts_toml = '[parts.A]\ndeps = []\n[parts.B]\ndeps = ["A"]\n[parts.C]\ndeps = ["B"]\n[parts.D]\ndeps = ["A"]\n'
    cases = {
        'good': ({'DDG/A/X.lean': 'import Mathlib\n', 'DDG/B/Y.lean': 'import DDG.A.X\n',
                  'DDG/C/Z.lean': '/- import DdgProofs.Foo -/\nimport DDG.A.X -- c\nimport DDG.B.Y\n',
                  'DDG.lean': 'import DDG.A.X\nimport DDG.D.W\n', 'DDG/D/W.lean': 'import DDG.A.X\n'}, 0),
        'upward': ({'DDG/A/X.lean': 'import DDG.B.Y\n', 'DDG/B/Y.lean': 'import Mathlib\n'}, 1),
        'sideways (D not in closure of C... via B)': ({'DDG/D/W.lean': 'import Mathlib\n',
                  'DDG/C/Z.lean': 'import DDG.D.W\n'}, 1),
        'app import': ({'DDG/A/X.lean': 'import DdgProofs.Foo\n'}, 1),
        'oracle import': ({'DDG/B/Y.lean': 'import Mathlib\nimport DdgOracle.Scalar\n'}, 1),
        'unplaced file': ({'DDG/Q/Y.lean': 'import Mathlib\n'}, 1),
        'umbrella imports app': ({'DDG.lean': 'import DdgProofs\n'}, 1),
    }
    bad = 0
    for name, (files, want) in cases.items():
        d = tempfile.mkdtemp(prefix='dagtest-')
        try:
            ddg = os.path.join(d, 'lean', 'ddg')
            for rel, content in files.items():
                p = os.path.join(ddg, rel); os.makedirs(os.path.dirname(p), exist_ok=True)
                open(p, 'w').write(content)
            open(os.path.join(ddg, 'parts.toml'), 'w').write(parts_toml)
            saved = sys.stderr; sys.stderr = open(os.devnull, 'w')
            try: got = run(d)
            finally: sys.stderr = saved
            ok = got == want
            print(f'self-test {name}: expected exit {want}, got {got}: {"ok" if ok else "WRONG"}')
            bad += not ok
        finally:
            shutil.rmtree(d)
    d = tempfile.mkdtemp(prefix='dagtest-')
    try:
        print('self-test no lean/ddg: ', end=''); got = run(d); bad += got != 0
    finally:
        shutil.rmtree(d)
    # cycle detection
    d = tempfile.mkdtemp(prefix='dagtest-')
    try:
        ddg = os.path.join(d, 'lean', 'ddg'); os.makedirs(ddg)
        open(os.path.join(ddg, 'parts.toml'), 'w').write('[parts.A]\ndeps=["B"]\n[parts.B]\ndeps=["A"]\n')
        saved = sys.stderr; sys.stderr = open(os.devnull, 'w')
        try: got = run(d)
        finally: sys.stderr = saved
        print(f'self-test cycle: expected exit 1, got {got}: {"ok" if got == 1 else "WRONG"}'); bad += got != 1
    finally:
        shutil.rmtree(d)
    print('self-test: ' + ('PASSED' if not bad else f'{bad} FAILED'))
    return 1 if bad else 0


if __name__ == '__main__':
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--repo', default=os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    ap.add_argument('--parts')
    ap.add_argument('--self-test', action='store_true')
    a = ap.parse_args()
    sys.exit(self_test() if a.self_test else run(a.repo, a.parts))
