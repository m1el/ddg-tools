#!/usr/bin/env bash
# Repulsive Curves, Sec. 7 comparison (crates/ddg-repulsive-curves/SPEC.md 14.2): job runner.
#
#   scripts/rc-eval.sh jobs  [all|subset]  > jobs.txt   job list: "<knot.obj> <method>" per line
#   scripts/rc-eval.sh A JOBS.txt     phase A: parallel (xargs -P $RC_P), one thread per job, audit on,
#                                     final curves; one log + one result file per job
#   scripts/rc-eval.sh B JOBS.txt     phase B: sequential, timed alone, $RC_THREADS threads, no audit,
#                                     frames every $RC_FRAME_S seconds (0 = none)
#
# Environment: RC_OUT (output dir, default runs_tmp/rc-eval/run), RC_P (default 60), RC_MAX_ITERS,
# RC_MAX_SECONDS (caps; empty = none), RC_HW (Hencky weight, default 1e3), RC_THREADS (B, default 64;
# A is always 1), RC_FRAME_S (B, default 0), RC_BACKEND (default auto), RC_BIN.
# Finished jobs (result file present) are skipped, so a run can be resumed.
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"
DATA=data/repulsive-curves/KnotDatasets
BIN="${RC_BIN:-$ROOT/target-rc-wp12/release/examples/eval_knots}"
OUT="${RC_OUT:-runs_tmp/rc-eval/run}"
METHODS="hs-pg l2-pg h1-pg h2-pg hs-armijo soft-hs-ncg soft-hs-lbfgs soft-h1-ncg soft-h1-lbfgs"

run_job() { # knot method mode
  local knot="$1" method="$2" mode="$3"
  local base; base="$(basename "$knot" .obj).$method"
  local res="$OUT/results/$base.jsonl"
  [ -s "$res" ] && return 0
  local args=(--knot "$knot" --method "$method" --hencky-w "${RC_HW:-1e3}" --backend "${RC_BACKEND:-auto}"
              --out "$res.tmp" --tag "$mode")
  [ -n "${RC_MAX_ITERS:-}" ] && args+=(--max-iters "$RC_MAX_ITERS")
  [ -n "${RC_MAX_SECONDS:-}" ] && args+=(--max-seconds "$RC_MAX_SECONDS")
  if [ "$mode" = A ]; then
    args+=(--threads 1 --final-obj "$OUT/final/$base.obj")
  else
    args+=(--threads "${RC_THREADS:-64}" --no-audit --final-obj "$OUT/final/$base.obj")
    if [ "${RC_FRAME_S:-0}" != 0 ]; then
      args+=(--frames-every-seconds "$RC_FRAME_S" --frames-out "$OUT/frames/$base.rcf")
    fi
  fi
  if "$BIN" "${args[@]}" > "$OUT/logs/$base.log" 2>&1; then
    mv "$res.tmp" "$res"
  else
    echo "FAILED $base (see $OUT/logs/$base.log)" >&2
    return 0
  fi
}
export -f run_job
export BIN OUT

case "${1:-}" in
  jobs)
    which="${2:-all}"
    if [ "$which" = subset ]; then
      list="$(cat "${RC_SUBSET:-scripts/rc-eval-subset.txt}")"
    else
      list="$(ls $DATA/Knot128/perturbed/*.obj $DATA/Trefoil100/*.obj)"
    fi
    for k in $list; do for m in $METHODS; do echo "$k $m"; done; done
    ;;
  A)
    mkdir -p "$OUT"/{logs,results,final}
    # slowest first (larger knots first) is up to the job list; xargs keeps the order
    xargs -a "$2" -P "${RC_P:-60}" -L 1 bash -c 'run_job "$0" "$1" A'
    ;;
  B)
    mkdir -p "$OUT"/{logs,results,final,frames}
    while read -r knot method; do run_job "$knot" "$method" B; done < "$2"
    ;;
  *)
    sed -n 2,14p "$0"; exit 2 ;;
esac
