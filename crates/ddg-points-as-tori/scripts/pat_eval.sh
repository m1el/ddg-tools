#!/usr/bin/env bash
# WP14: the whole Points as Tori evaluation (SPEC.md §16.2–16.4, §16.7), CPU only.
#
#   crates/ddg-points-as-tori/scripts/pat_eval.sh [stage...]
#
# Stages (default: all, in this order): data prepare spsr ref case continuity orient density report.
# Run from the repository root. Everything goes to runs_tmp/pat-wp14/ (one log per job). Jobs run
# under `nice -n 19`; JOBS × THREADS bounds the threads in use (default 12 × 4 = 48). Environment:
#   JOBS, THREADS            parallel jobs and threads per job
#   RETRAINED                optional extra `--retrained label=path[:geometric]` arguments of the
#                            case stage, e.g.
#                            RETRAINED="--retrained retrained=data/points-as-tori/retrained-reference.safetensors
#                                       --retrained retrained-geo=data/points-as-tori/retrained-geometric.safetensors:geometric"
#   PAT_TRAIN_ROOT           the authors' training data (PAT_training_data/PAT_training_data)
#   PAT_REF_VENV             the reference build's venv (runs_tmp/pat-ref/README.md)
set -euo pipefail
JOBS=${JOBS:-12}
THREADS=${THREADS:-4}
W=runs_tmp/pat-wp14
C=$W/cases
L=$W/logs
B=${CARGO_TARGET_DIR:-$PWD/target-pat14}/release/examples
PAT_TRAIN_ROOT=${PAT_TRAIN_ROOT:-/var/data/discrete-geometry/PointsAsTori/PAT_training_data/PAT_training_data}
PAT_REF_VENV=${PAT_REF_VENV:-/var/data/discrete-geometry/PointsAsTori-build/.venv}
S=crates/ddg-points-as-tori/scripts
mkdir -p "$C" "$L"
stages=("$@")
[ ${#stages[@]} -eq 0 ] && stages=(data prepare spsr ref case continuity orient density report)

par() { # par <list-file> <log-prefix> <command with {}>
	local list=$1 pre=$2
	shift 2
	mkdir -p "$L/$pre"
	xargs -P "$JOBS" -I{} sh -c "nice -n 19 $* > $L/$pre/\$(echo {} | tr / _).log 2>&1 || echo FAILED {}" < "$list"
}

for st in "${stages[@]}"; do
	case $st in
	data)
		nice -n 19 uv run --no-project --with numpy python -I $S/pat_eval_data.py --train-root "$PAT_TRAIN_ROOT" --out $W/src --seen 140 > $L/data.log 2>&1
		ln -sfn ../src $C/src
		;;
	prepare)
		cargo build --release -p ddg-points-as-tori --example eval --example density > $L/build.log 2>&1
		for d in abc blobs thingi analytic saddle ref reffull orient analyticsad; do nice -n 19 $B/eval prepare $d $C --threads $THREADS > $L/prep-$d.log 2>&1; done
		for d in saddle analyticsad abc blobs thingi; do nice -n 19 $B/eval prepare $d $C --points 2048 --threads $THREADS > $L/prep-$d-2048.log 2>&1; done
		ls -d $C/{abc,blobs,thingi,analytic,saddle,ref,analyticsad,reffull,saddle2048,analyticsad2048}/*/ | sed 's:/$::' > $W/all_cases.txt
		grep -v reffull $W/all_cases.txt > $W/ref_cases.txt
		;;
	spsr)
		par $W/all_cases.txt spsr "uv run $S/pat_spsr.py --cloud {} --out {}/spsr"
		;;
	ref)
		rm -f $W/refchunk_*
		split -n r/6 -d $W/ref_cases.txt $W/refchunk_
		for c in $W/refchunk_*; do
			JAX_PLATFORMS=cpu OMP_NUM_THREADS=4 XLA_FLAGS="--xla_cpu_multi_thread_eigen=false intra_op_parallelism_threads=4" \
				nice -n 19 uv run --no-project --python $PAT_REF_VENV/bin/python -- python -I $S/pat_eval_ref.py $c > $L/ref-$(basename $c).log 2>&1 &
		done
		wait
		;;
	case)
		# the §16.3 sweep on the first 50 ABC and 50 blob shapes
		(ls -d $C/abc/*/ | head -50; ls -d $C/blobs/*/ | head -50) | sed 's:/$::' > $W/sweep_cases.txt
		grep -vxF -f $W/sweep_cases.txt $W/all_cases.txt > $W/nosweep_cases.txt || true
		par $W/sweep_cases.txt case "$B/eval case {} --threads $THREADS --sweep ${RETRAINED:-}"
		par $W/nosweep_cases.txt case "$B/eval case {} --threads $THREADS ${RETRAINED:-}"
		;;
	continuity)
		(ls -d $C/reffull/*/; for d in abc blobs thingi analytic; do ls -d $C/$d/*/ | head -8; done) | sed 's:/$::' > $W/cont_cases.txt
		par $W/cont_cases.txt cont "$B/eval continuity {} --threads $THREADS"
		;;
	orient)
		(ls -d $C/orient/*/; for d in abc2048 blobs2048 thingi2048; do ls -d $C/$d/*/ | head -10; done) | sed 's:/$::' > $W/orient_cases.txt
		par $W/orient_cases.txt orient "$B/eval orient {} --threads $THREADS"
		;;
	density)
		nice -n 19 $B/density data/models/spot/spot/spot_triangulated.obj $W/density-fig7-spot --counts 1024,8192,65536 --grid 128 --mc 128 --threads 16 > $L/density-fig7.log 2>&1
		nice -n 19 $B/density data/thingi10k/converted/59340.obj $W/density-fig29-thingi59340 --counts 110000,880000,7300000 --grid 128 --mc 256 --threads 32 > $L/density-fig29.log 2>&1
		;;
	report)
		uv run --no-project --with numpy --with matplotlib --with pillow python -I $S/pat_eval_report.py $C $W/report > $L/report.log 2>&1
		;;
	*) echo "unknown stage $st" >&2; exit 2 ;;
	esac
done
