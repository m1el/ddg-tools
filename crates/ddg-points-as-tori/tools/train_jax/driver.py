"""Driver for the authors' Points-as-Tori trainer (Feng, Gkioulekas, Crane, TOG 2026; code MIT, (c) 2026 Nicole Feng).

Runs `train()` from the authors' `training/train.py` verbatim (their curriculum, hyperparameters, seeds, f32),
using a copy of train.py patched only for periodic checkpoints, resume and progress lines (see PATCHES.diff).
`mesh_sampling_jax.py`, `mesh_data.py` and the `pointsastori` package are imported unmodified from their clone.

Their code uses paths relative to the working directory (`datasets/...`, `models/`); this driver chdirs into a
work directory that holds a `datasets` symlink to the released training data and the `models/` output folder.

Environment:
  PAT_WORK         work directory (required)
  PAT_SRC          directory with the patched train.py (default: <this dir>/src)
  PAT_PHASES       comma list of phase tags to run; the driver stops before the first phase not listed
                   (tags: FundamentalFormPredictor_pretrained, _phase0 ... _phase4; default: all)
  PAT_RESUME_DIR, PAT_CKPT_SECONDS, PAT_LOG_EVERY, PAT_MAX_STEPS: see the patch in src/train.py
"""

import hashlib
import os
import sys
import time

AUTHORS_TRAINING = '/var/data/discrete-geometry/PointsAsTori-build/points-as-tori/training'
HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.environ.get('PAT_SRC', os.path.join(HERE, 'src'))
WORK = os.environ['PAT_WORK']

# Patched train.py first, then the authors' training dir (mesh_data, mesh_sampling_jax).
sys.path.insert(0, AUTHORS_TRAINING)
sys.path.insert(0, SRC)
os.makedirs(os.path.join(WORK, 'models'), exist_ok=True)
os.makedirs(os.path.join(WORK, 'plots'), exist_ok=True)
os.chdir(WORK)
assert os.path.isdir('datasets/nice'), f'{WORK}/datasets must link to PAT_training_data/PAT_training_data'

import jax  # noqa: E402

# Persistent compilation cache (shared by the run, its restarts, tests and the Geometric probe): compiled
# executables (and XLA's per-fusion autotune results) are reused instead of re-running ptxas. No effect on the math.
CACHE = os.environ.get('PAT_JAX_CACHE', '/var/data/discrete-geometry/ddg-tools/runs_tmp/pat-train-jax/jax-cache')
if CACHE:
	jax.config.update('jax_compilation_cache_dir', CACHE)
	jax.config.update('jax_persistent_cache_min_compile_time_secs', 0.0)
	jax.config.update('jax_persistent_cache_min_entry_size_bytes', -1)
if os.environ.get('PAT_LOG_COMPILES'):
	jax.config.update('jax_log_compiles', True)  # one WARNING line per compilation, with the function name

# Compile counter: backend compiles (including persistent-cache retrievals) and cache hits; shown in the progress lines.
COMPILES = {'n': 0, 's': 0.0, 'hits': 0}


def _on_duration(event, duration, **kw):
	if event == '/jax/core/compile/backend_compile_duration':
		COMPILES['n'] += 1
		COMPILES['s'] += duration
		if os.environ.get('PAT_LOG_COMPILES'):
			print(f'[PAT] compile #{COMPILES["n"]} {kw.get("fun_name", "?")} {duration:.2f} s', flush=True)


def _on_event(event, **kw):
	if event == '/jax/compilation_cache/cache_hits':
		COMPILES['hits'] += 1


jax.monitoring.register_event_duration_secs_listener(_on_duration)
jax.monitoring.register_event_listener(_on_event)

import train as T  # noqa: E402

T.PAT_COMPILES = COMPILES
import flax  # noqa: E402
import optax  # noqa: E402
import numpy as np  # noqa: E402

CURRICULUM = [  # steps per phase (NOTES-REFERENCE.md §9.3), for the global step / ETA in the progress lines
	('FundamentalFormPredictor_pretrained', 4096),
	('FundamentalFormPredictor_phase0', 96000),
	('FundamentalFormPredictor_phase1', 64000),
	('FundamentalFormPredictor_phase2', 119344),
	('FundamentalFormPredictor_phase3', 124590),
	('FundamentalFormPredictor_phase4', 89828),
]
T.PAT_CURRICULUM = CURRICULUM
PHASES = os.environ.get('PAT_PHASES', '')
PHASES = [p for p in PHASES.split(',') if p] or [t for t, _ in CURRICULUM]


class StopCurriculum(Exception):
	pass


def sha(path):
	with open(path, 'rb') as f:
		return hashlib.sha256(f.read()).hexdigest()[:16]


print(
	f'[PAT] driver start {time.strftime("%Y-%m-%dT%H:%M:%S")} pid {os.getpid()} work {WORK}\n'
	f'[PAT] train.py {T.__file__} sha {sha(T.__file__)}; authors train.py sha '
	f'{sha(os.path.join(AUTHORS_TRAINING, "train.py"))}; mesh_sampling_jax {sys.modules["mesh_sampling_jax"].__file__}\n'
	f'[PAT] pointsastori.network {sys.modules["pointsastori.network"].__file__} sha '
	f'{sha(sys.modules["pointsastori.network"].__file__)}\n'
	f'[PAT] jax {jax.__version__} flax {flax.__version__} optax {optax.__version__} numpy {np.__version__} '
	f'x64 {jax.config.jax_enable_x64} devices {jax.devices()} phases {PHASES}',
	flush=True,
)

_orig = T.TrainFundamentalFormPredictor


def phase(*args, final_save_filepath=None, **kwargs):
	tag = final_save_filepath
	done = os.path.join(T.PAT_RESUME_DIR, tag + '.done')
	if os.path.exists(done):
		print(f'[PAT] phase {tag} already done ({done}); skipping', flush=True)
		return
	if tag not in PHASES:
		raise StopCurriculum(f'phase {tag} not in PAT_PHASES')
	t0 = time.time()
	print(f'[PAT] phase {tag} start {time.strftime("%Y-%m-%dT%H:%M:%S")}', flush=True)
	_orig(*args, final_save_filepath=final_save_filepath, **kwargs)
	if os.path.exists('PlotTrainingHistory.png'):
		os.replace('PlotTrainingHistory.png', os.path.join('plots', tag + '.png'))
	os.makedirs(T.PAT_RESUME_DIR, exist_ok=True)
	with open(done, 'w') as f:
		f.write(f'{time.strftime("%Y-%m-%dT%H:%M:%S")} {time.time() - t0:.1f} s (this process)\n')
	print(f'[PAT] phase {tag} done after {(time.time() - t0) / 3600:.2f} h (this process)', flush=True)


T.TrainFundamentalFormPredictor = phase

try:
	T.train()
	print(f'[PAT] curriculum finished {time.strftime("%Y-%m-%dT%H:%M:%S")}', flush=True)
except StopCurriculum as e:
	print(f'[PAT] stop: {e}', flush=True)
except T.PatStopPhase as e:
	print(f'[PAT] stop: {e}', flush=True)
