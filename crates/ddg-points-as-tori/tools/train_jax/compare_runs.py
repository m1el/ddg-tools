"""Compare two runs of the JAX trainer (resume test): [PAT] progress lines, their epoch/validation lines, and
the resume checkpoints (npz, no pickle).  python -I compare_runs.py LOG_A CKPT_A LOG_B CKPT_B"""

import re
import sys

import numpy as np


def lines(path):
	steps, epochs = {}, []
	for ln in open(path):
		m = re.search(r'\[PAT\] \S+ phase (\S+) .* step (\d+)/\d+ .* loss (\S+) d (\S+) e (\S+) ', ln)
		if m:
			steps[(m.group(1), int(m.group(2)))] = tuple(m.group(i) for i in (3, 4, 5))
		elif re.match(r'\s*(total loss:|validation:)', ln) or re.match(r'\s*Epoch \d+, shape', ln):
			epochs.append(ln.strip())
	return steps, epochs


la, ca, lb, cb = sys.argv[1:5]
sa, ea = lines(la)
sb, eb = lines(lb)
common = sorted(set(sa) & set(sb))
diff = [k for k in common if sa[k] != sb[k]]
print(f'progress lines: {len(sa)} vs {len(sb)}, common {len(common)}, differing {len(diff)}')
for k in diff[:5]:
	print('  ', k, sa[k], sb[k])
# Epoch/validation lines of B (resumed) must be a suffix of A's lines (the resumed run prints only the rest).
suffix = ea[-len(eb):] if eb else []
print(f'epoch/validation lines: A {len(ea)}, B {len(eb)}; B == tail of A: {suffix == eb}')
if suffix != eb:
	for x, y in zip(suffix, eb):
		if x != y:
			print('   A:', x, '\n   B:', y)
A, B = np.load(ca), np.load(cb)
keys = sorted(set(A.files) | set(B.files))
missing = [k for k in keys if k not in A.files or k not in B.files]
bad = [(k, float(np.max(np.abs(A[k].astype(np.float64) - B[k].astype(np.float64))))) for k in keys
	if k not in missing and (A[k].shape != B[k].shape or not np.array_equal(A[k], B[k]))]
print(f'checkpoint arrays: {len(keys)}, missing in one: {len(missing)}, not bitwise equal: {len(bad)}')
for k, v in sorted(bad, key=lambda x: -x[1])[:8]:
	print(f'   {k}: max |diff| {v:.3e}')
