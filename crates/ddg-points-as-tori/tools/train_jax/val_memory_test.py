"""GPU-memory test: one phase-0 train step (batch 1024) followed by the authors' epoch-end validation code path
(model forward on a validation shape + validation_loss at batch 1024), in one process, as at an epoch end.
Fresh weights (memory only). Run with the same XLA_PYTHON_CLIENT_MEM_FRACTION as the run; cwd = a work dir."""

import sys
import time

sys.path.insert(0, '/var/data/discrete-geometry/PointsAsTori-build/points-as-tori/training')
sys.path.insert(0, '/var/data/discrete-geometry/ddg-tools/runs_tmp/pat-train-jax/src')
import train as T  # noqa: E402
import jax  # noqa: E402
import jax.numpy as jnp  # noqa: E402
from flax import nnx  # noqa: E402

model, optimizer, step = T.FundamentalFormPredictorTrainer(embed_dim=128, n_freq=0, n_layers=8, n_heads=8, mlp_dim=512,
	rngs=nnx.Rngs(42), learning_rate=1e-5, warmup_steps=100, total_steps=1e8)
P, N, NB, M, V, F, FM = T.load_mesh_training_data('datasets/nice/blobs_train_mode0_0.npz')
key = jax.random.PRNGKey(0)
pos, nor, nb = jnp.array(P[0]), jnp.array(N[0]), jnp.array(NB[0])
q, d = T.generate_shape_features(pos, nb, jnp.array(V[0]), jnp.array(F[0]), jnp.array(FM[0]), key, n_queries=120,
	narrow_band_offset=0.2, bbox_extension=1.0)
mask = jnp.array(M[0])
for b in range(2):
	s = b * 1024
	loss, _ = step(model, optimizer, mask, jax.lax.dynamic_slice_in_dim(mask, s, 1024), pos, nor, nb,
		jax.lax.dynamic_slice_in_dim(q, s, 1024), jax.lax.dynamic_slice_in_dim(d, s, 1024))
	print('train step', b, float(loss), flush=True)
vp, vn, vnb, vm, vv, vf, vfm = T.load_mesh_training_data('datasets/ABC17_validation_3.npz')
nnx.eval_mode(model)
t = time.time()
for i in range(len(vp)):
	vq, vd = T.generate_shape_features(vp[i], vnb[i], vv[i], vf[i], vfm[i], jax.random.fold_in(key, 10000 + i), n_queries=120)
	w = model(vp[i][vnb[i]], vn[i][vnb[i]])
	for b in range(2):
		s = b * 1024
		bm = jax.lax.dynamic_slice_in_dim(vm[i], s, 1024)
		if jnp.sum(bm) > 0:
			out = T.validation_loss(vm[i], bm, vp[i], vn[i], w, jax.lax.dynamic_slice_in_dim(vq, s, 1024),
				jax.lax.dynamic_slice_in_dim(vd, s, 1024))
	print('val shape', i, [float(x) for x in out], flush=True)
print(f'validation {time.time() - t:.1f} s; peak {jax.devices()[0].memory_stats()["peak_bytes_in_use"] / 2**30:.2f} GiB')
print('OK')
