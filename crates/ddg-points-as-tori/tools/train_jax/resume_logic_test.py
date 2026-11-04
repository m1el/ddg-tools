"""Mid-epoch resume test across shapes (the phase-0..4 path: many shapes per epoch, cursor inside the shuffled
permutation). Calls the authors' TrainFundamentalFormPredictor (patched copy) directly on one blob file with
batch 256 (8 steps per shape; the phase-0 batch 1024 does not fit with deterministic ops on 16 GB), lr 1e-5,
warmup 100, from fresh weights, and stops after PAT_MAX_STEPS. cwd = a work dir with datasets/ and models/.
Run with XLA_FLAGS=--xla_gpu_deterministic_ops=true."""

import sys

sys.path.insert(0, '/var/data/discrete-geometry/PointsAsTori-build/points-as-tori/training')
sys.path.insert(0, '/var/data/discrete-geometry/ddg-tools/runs_tmp/pat-train-jax/src')
import train as T  # noqa: E402

try:
	T.TrainFundamentalFormPredictor(['datasets/nice/blobs_train_mode0_0.npz'], num_epochs=2, learning_rate=1e-5,
		batch_size=256, warmup_steps=1e2, n_epochs_to_save=1, n_shapes_to_save=1, load_model_path=None,
		final_save_filepath='resume_logic')
except T.PatStopPhase as e:
	print('[PAT] stop:', e)
