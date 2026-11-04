# /// script
# requires-python = ">=3.11"
# dependencies = ["numpy"]
# ///
"""WP14 data extraction for the Points as Tori evaluation (SPEC.md §16.2).

Extracts the ABC meshes of the evaluation from the authors' training data and determines the
train/test overlap that can be determined from the shipped files.

* "unseen" ABC shapes: the 10 + 10 shapes of `ABC0_validation_0.npz` and `ABC17_validation_3.npz`
  (the authors' test splits of chunks 2/17/55/58 are not shipped).
* "seen" ABC shapes: a seeded random pick from `nice/ABC*_train_mode0_*.npz` (used in training
  phases 1-4) to fill the set to 100; reported separately in EVALUATION.md.
* Overlap: every validation shape is fingerprinted (exact: sorted vertex coordinates rounded to
  1e-6; coarse: face count, area and volume rounded to 3 digits) and compared with every shape in
  `nice/` and `simulated/`.

Writes `<out>/abc/<name>/{vertices,faces}.npy` (float64 / int32, padding removed, unreferenced
vertices dropped), `<out>/abc/manifest.tsv` and `<out>/overlap.json`.

    uv run --no-project --with numpy python -I crates/ddg-points-as-tori/scripts/pat_eval_data.py \
        --train-root /var/data/discrete-geometry/PointsAsTori/PAT_training_data/PAT_training_data \
        --out runs_tmp/pat-wp14/src --seen 140
"""

import argparse
import glob
import hashlib
import json
import os
import sys

import numpy as np


def compact(v, f, fm):
    f = f[fm.astype(bool)] if fm is not None else f
    f = f[(f >= 0).all(1)]
    used = np.unique(f.reshape(-1))
    remap = -np.ones(len(v), dtype=np.int64)
    remap[used] = np.arange(len(used))
    return v[used].astype(np.float64), remap[f].astype(np.int32)


def fingerprints(v, f):
    if len(f) == 0:
        return None, None
    key = np.round(v, 6)
    key = key[np.lexsort(key.T[::-1])]
    exact = hashlib.sha1(key.tobytes()).hexdigest()
    a, b, c = v[f[:, 0]], v[f[:, 1]], v[f[:, 2]]
    cr = np.cross(b - a, c - a)
    area = 0.5 * np.linalg.norm(cr, axis=1).sum()
    vol = np.einsum('ij,ij->i', a, np.cross(b, c)).sum() / 6.0
    coarse = (len(f), float(f"{area:.3g}"), float(f"{vol:.3g}"))
    return exact, coarse


def shapes(path):
    z = np.load(path)
    V, F, FM = z['mesh_vertices'], z['mesh_faces'], z['face_masks']
    for s in range(len(V)):
        yield s, compact(V[s], F[s], FM[s])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--train-root', required=True)
    ap.add_argument('--out', required=True)
    ap.add_argument('--seen', type=int, default=140, help='candidates drawn from the training files (validated later in Rust)')
    ap.add_argument('--seed', type=int, default=14)
    args = ap.parse_args()
    root = args.train_root
    out_abc = os.path.join(args.out, 'abc')
    os.makedirs(out_abc, exist_ok=True)
    manifest = []
    val = {}
    for fname, tag in [('ABC0_validation_0.npz', 'abc0val'), ('ABC17_validation_3.npz', 'abc17val')]:
        for s, (v, f) in shapes(os.path.join(root, fname)):
            name = f'{tag}_{s}'
            d = os.path.join(out_abc, name)
            os.makedirs(d, exist_ok=True)
            np.save(os.path.join(d, 'vertices.npy'), v)
            np.save(os.path.join(d, 'faces.npy'), f)
            val[name] = fingerprints(v, f)
            manifest.append((name, fname, s, 'unseen', len(v), len(f)))
    # fingerprints of every training shape
    train_files = sorted(glob.glob(os.path.join(root, 'nice', 'ABC*_train_mode0_*.npz'))) + sorted(
        glob.glob(os.path.join(root, 'nice', 'blobs_train_mode0_*.npz'))) + sorted(glob.glob(os.path.join(root, 'simulated', '*.npz')))
    exact_idx, coarse_idx = {}, {}
    abc_nice = []
    for path in train_files:
        rel = os.path.relpath(path, root)
        for s, (v, f) in shapes(path):
            e, c = fingerprints(v, f)
            if e is None:
                continue
            exact_idx.setdefault(e, []).append(f'{rel}:{s}')
            coarse_idx.setdefault(c, []).append(f'{rel}:{s}')
            if rel.startswith('nice/ABC'):
                abc_nice.append((rel, s, len(f)))
        print(f'fingerprinted {rel}', flush=True)
    overlap = {}
    for name, (e, c) in val.items():
        overlap[name] = {'exact': exact_idx.get(e, []), 'coarse': coarse_idx.get(c, []) if c else []}
    n_exact = sum(1 for o in overlap.values() if o['exact'])
    n_coarse = sum(1 for o in overlap.values() if o['coarse'])
    # seen candidates: seeded pick among non-degenerate nice ABC shapes (>= 100 faces)
    rng = np.random.default_rng(args.seed)
    cands = [x for x in abc_nice if x[2] >= 100]
    pick = rng.choice(len(cands), size=min(args.seen, len(cands)), replace=False)
    by_file = {}
    for i in sorted(pick):
        rel, s, _ = cands[i]
        by_file.setdefault(rel, []).append(s)
    for rel, idx in sorted(by_file.items()):
        z = np.load(os.path.join(root, rel))
        V, F, FM = z['mesh_vertices'], z['mesh_faces'], z['face_masks']
        stem = os.path.basename(rel).replace('_train_mode0_', 't').replace('.npz', '')
        for s in idx:
            v, f = compact(V[s], F[s], FM[s])
            name = f'{stem}_{s}'.lower()
            d = os.path.join(out_abc, name)
            os.makedirs(d, exist_ok=True)
            np.save(os.path.join(d, 'vertices.npy'), v)
            np.save(os.path.join(d, 'faces.npy'), f)
            manifest.append((name, rel, s, 'seen', len(v), len(f)))
    with open(os.path.join(out_abc, 'manifest.tsv'), 'w') as fh:
        fh.write('name\tsource\tindex\ttraining\tvertices\tfaces\n')
        for row in manifest:
            fh.write('\t'.join(str(x) for x in row) + '\n')
    summary = {
        'training_files_fingerprinted': [os.path.relpath(p, root) for p in train_files],
        'training_shapes': sum(len(x) for x in exact_idx.values()),
        'nice_abc_candidates': len(cands),
        'validation_shapes': len(val),
        'validation_with_exact_match_in_training': n_exact,
        'validation_with_coarse_match_in_training': n_coarse,
        'per_shape': overlap,
        'seen_picked': int(len(pick)),
        'seed': args.seed,
    }
    with open(os.path.join(args.out, 'overlap.json'), 'w') as fh:
        json.dump(summary, fh, indent=1)
    print(json.dumps({k: v for k, v in summary.items() if k != 'per_shape' and k != 'training_files_fingerprinted'}))
    return 0


if __name__ == '__main__':
    sys.exit(main())
