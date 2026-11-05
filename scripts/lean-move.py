#!/usr/bin/env python3
"""Move Lean modules between packages/directories and fix every reference.

usage: lean-move.py MOVES.tsv [--renames RENAMES.tsv] [--repo DIR] [--dry-run] [--diff]
                    [--names-in EXT,EXT,...] [--allow-dirty]

MOVES.tsv    one move per line:  old-module <TAB> new-module   (e.g. DdgProofs.Cuts.Weighted <TAB>
             DDG.Mesh.Weighted). '#' starts a comment. The package of a module is chosen by its root:
               DdgProofs -> lean/proofs, DDG -> lean/ddg, DdgOracle -> lean/oracle
             (other roots: add to PACKAGES below). File path = <package>/<Mod/As/Path>.lean.
RENAMES.tsv  optional, one per line:  old-full-name <TAB> new-full-name  (declarations or
             namespaces). Replaced as whole identifiers: not preceded by an identifier char or '.',
             and followed by anything except an identifier char; a trailing `.suffix` is allowed, so
             a namespace rename also rewrites `Old.NS.foo` and a declaration rename also rewrites
             `Old.decl.eq_1` / `Old.decl.symm`. Names written relative to an open namespace are
             NOT found (the build reports them). All renames apply in a single pass (no chaining).

What it does, in order
  1. `git mv` of each module's .lean file (destination directories are created).
  2. In every tracked/untracked-unignored .lean file: `import` lines rewritten from the move map;
     name renames applied (also in --names-in files, default md,rs,sh,py).
  3. In .lean/.md/.rs/.sh/.py/.toml files: path mentions `DdgProofs/Foo/Bar.lean` (optionally
     preceded by `lean/proofs/`) are rewritten to the new path; in .md/.rs/.sh/.py/.toml files the
     module name `DdgProofs.Foo.Bar` (when not followed by `.ident`) is rewritten too.
Not handled (do by hand): lakefile `globs`/`roots`, lake manifests, mentions of a module's directory
without a file name, relative names in docs.
--dry-run prints the planned moves and per-file replacement counts (--diff: unified diffs) and
changes nothing. The tree must be clean unless --allow-dirty (so a bad move is one `git revert`).
"""
import argparse, difflib, os, re, subprocess, sys

PACKAGES = {'DdgProofs': 'lean/proofs', 'DDG': 'lean/ddg', 'DdgOracle': 'lean/oracle'}
IDCH = r"[\w'!?À-￿]"   # identifier characters (rough; includes unicode letters/subscripts)


def sh(repo, *cmd):
    return subprocess.run(cmd, cwd=repo, check=True, capture_output=True, text=True).stdout


def read_tsv(path):
    rows = []
    for ln, line in enumerate(open(path, encoding='utf-8'), 1):
        line = line.split('#', 1)[0].rstrip('\n').rstrip()
        if not line.strip():
            continue
        parts = line.split('\t')
        if len(parts) < 2:
            sys.exit(f'{path}:{ln}: expected two tab-separated columns')
        rows.append((parts[0].strip(), parts[1].strip()))
    return rows


def mod_path(mod):
    root = mod.split('.')[0]
    if root not in PACKAGES:
        sys.exit(f'unknown module root {root} in {mod}')
    return os.path.join(PACKAGES[root], *mod.split('.')) + '.lean'


def alternation(keys):
    return '|'.join(re.escape(k) for k in sorted(keys, key=len, reverse=True))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('moves')
    ap.add_argument('--renames')
    ap.add_argument('--repo', default=os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    ap.add_argument('--dry-run', action='store_true')
    ap.add_argument('--diff', action='store_true')
    ap.add_argument('--names-in', default='md,rs,sh,py')
    ap.add_argument('--allow-dirty', action='store_true')
    a = ap.parse_args()
    repo = os.path.abspath(a.repo)

    moves = dict(read_tsv(a.moves))
    renames = dict(read_tsv(a.renames)) if a.renames else {}
    if len(moves) != len(read_tsv(a.moves)):
        sys.exit('duplicate old module in moves')
    if not a.dry_run and not a.allow_dirty and sh(repo, 'git', 'status', '--porcelain').strip():
        sys.exit('working tree is not clean (commit first, or pass --allow-dirty)')

    # --- plan the file moves
    file_moves = {}
    for old, new in moves.items():
        src, dst = mod_path(old), mod_path(new)
        if not os.path.isfile(os.path.join(repo, src)):
            sys.exit(f'missing source file {src}')
        if os.path.exists(os.path.join(repo, dst)):
            sys.exit(f'destination exists: {dst}')
        file_moves[src] = dst
    if len(set(file_moves.values())) != len(file_moves):
        sys.exit('two modules move to the same destination')
    print(f'{len(file_moves)} module move(s), {len(renames)} name rename(s)')
    for s, d in file_moves.items():
        print(f'  git mv {s} {d}')

    # --- the file list after the moves
    listed = sh(repo, 'git', 'ls-files', '-co', '--exclude-standard', '-z').split('\0')
    files = []
    for f in listed:
        if not f or '/.lake/' in '/' + f or not os.path.isfile(os.path.join(repo, f)):
            continue
        files.append(file_moves.get(f, f))
    if not a.dry_run:
        for s, d in file_moves.items():
            os.makedirs(os.path.dirname(os.path.join(repo, d)), exist_ok=True)
            sh(repo, 'git', 'mv', s, d)

    # --- regexes (single pass, longest first)
    name_re = None
    if renames:
        name_re = re.compile(r'(?<!' + IDCH + r')(?<!\.)(' + alternation(renames) + r')(?!' + IDCH + r')')
    path_map = {}
    for old, new in moves.items():
        op, np_ = mod_path(old), mod_path(new)
        # relative-to-package form and package-prefixed form
        path_map[op] = np_
        path_map[op.split('/', 2)[2]] = np_.split('/', 2)[2]
    path_re = re.compile(r'(?<![\w/.-])(' + alternation(path_map) + r')') if path_map else None
    mod_re = (re.compile(r'(?<!' + IDCH + r')(?<!\.)(' + alternation(moves) + r')(?!' + IDCH + r')(?!\.[A-Za-z_])')
              if moves else None)
    imp_re = re.compile(r'^(\s*(?:(?:public|meta|private)\s+)*import\s+)(.*)$')
    names_ext = {'.' + e for e in a.names_in.split(',') if e}

    total = 0
    for rel in sorted(files):
        ext = os.path.splitext(rel)[1]
        if rel == 'scripts/lean-move.py':
            continue
        if ext not in ('.lean', '.md', '.rs', '.sh', '.py', '.toml'):
            continue
        full = os.path.join(repo, rel if not a.dry_run else next((s for s, d in file_moves.items() if d == rel), rel))
        try:
            text = open(full, encoding='utf-8').read()
        except (UnicodeDecodeError, FileNotFoundError):
            continue
        new, n = text, 0

        def sub(rx, fn, s):
            nonlocal n
            out, k = rx.subn(fn, s)
            n += k
            return out
        if ext == '.lean':
            lines = new.split('\n')
            in_header = True
            for i, line in enumerate(lines):
                m = imp_re.match(line)
                if m:
                    toks = m.group(2).split(' ')
                    t2 = [moves.get(t, t) for t in toks]
                    if t2 != toks:
                        n += 1
                        lines[i] = m.group(1) + ' '.join(t2)
            new = '\n'.join(lines)
        if name_re and (ext == '.lean' or ext in names_ext):
            new = sub(name_re, lambda m: renames[m.group(1)], new)
        if path_re:
            def pr(m, _pm=path_map):
                return _pm[m.group(1)]
            # `lean/proofs/` prefix before a relative path: map to the new package prefix
            def pfx(m):
                full_old = m.group(1) + m.group(2)
                return path_map[full_old] if full_old in path_map else m.group(0)
            new = sub(re.compile(r'(?<![\w/.-])(lean/(?:proofs|ddg|oracle)/)(' + alternation(
                [k.split('/', 2)[2] for k in path_map if k.startswith('lean/')]) + r')'), pfx, new)
            new = sub(path_re, pr, new)
        if mod_re and ext in ('.md', '.rs', '.sh', '.py', '.toml'):
            new = sub(mod_re, lambda m: moves[m.group(1)], new)
        if new != text:
            total += n
            print(f'  {rel}: {n} replacement(s)')
            if a.diff:
                sys.stdout.writelines(difflib.unified_diff(text.splitlines(True), new.splitlines(True),
                                                           'a/' + rel, 'b/' + rel, n=0))
            if not a.dry_run:
                with open(full, 'w', encoding='utf-8') as fh:
                    fh.write(new)
    print(('would make' if a.dry_run else 'made') + f' {total} replacement(s)')


main()
