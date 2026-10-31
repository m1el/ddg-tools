#!/usr/bin/env bash
# WP10: the Thingi10k evaluation of integer coordinates and signposts (ddg-integer-coords SPEC 10.2).
#
#   scripts/ic-eval.sh convert            # raw meshes -> data/thingi10k/converted/ + converted.tsv
#   scripts/ic-eval.sh run [BACKENDS...]  # one process per (mesh, backend), in parallel
#   scripts/ic-eval.sh summary            # runs_tmp/wp10/summary.md from the job results
#   scripts/ic-eval.sh cmp [MAX_FACES]    # integer runs with the IDT edge lengths dumped, for the demo comparison
#
# BACKENDS: integer sp-ruppert sp-chew (default: all three); integer-ref = integers with the authors' Delaunay test (cotan weight < -1e-6). Env: IC_JOBS (parallel jobs in the
# first phase, default 56), IC_TIMEOUT (seconds per job, default 600), IC_IDS (file with ids,
# default: every kept mesh). Results: runs_tmp/wp10/jobs/<backend>/<id>.json and .log. A job
# killed by the time cap is recorded as {"status":"timeout"}, one that dies without a result
# (address-space limit reached, abort) as "memory" or "crash".
# Meshes are run in three phases by size so that memory stays bounded:
#   <= 100k faces: IC_JOBS jobs, 4 GB each; <= 1M faces: 12 jobs, 12 GB; larger: 3 jobs, 24 GB.
set -u
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
export CARGO_TARGET_DIR="${CARGO_TARGET_DIR:-$ROOT/target-wp10}"
BIN="$CARGO_TARGET_DIR/release/examples"
CONV="$ROOT/data/thingi10k/converted"
OUT="$ROOT/runs_tmp/wp10"
TIMEOUT="${IC_TIMEOUT:-600}"

job() { # job ID BACKEND MEM_KB
  local id="$1" backend="$2" mem="$3" d="$OUT/jobs/$2"
  mkdir -p "$d"
  [ -s "$d/$id.json" ] && return 0
  ( ulimit -v "$mem"; timeout --signal=KILL "$TIMEOUT" "$BIN/thingi_eval" "$CONV/$id.obj" --backend "$backend" --id "$id" \
      > "$d/$id.out" 2> "$d/$id.log"; echo "exit=$?" >> "$d/$id.log" )
  if grep -q '^exit=0' "$d/$id.log" && [ -s "$d/$id.out" ]; then
    mv "$d/$id.out" "$d/$id.json"
  else
    local code status
    code=$(grep '^exit=' "$d/$id.log" | tail -1 | cut -d= -f2)
    if [ "$code" = 137 ] || [ "$code" = 124 ]; then status=timeout
    elif grep -qiE 'memory allocation|allocat|capacity overflow' "$d/$id.log"; then status=memory
    else status=crash; fi
    echo "{\"id\":\"$id\",\"backend\":\"$backend\",\"status\":\"$status\",\"exit\":$code}" > "$d/$id.json"
    rm -f "$d/$id.out"
  fi
}

case "${1:-}" in
  --job) shift; job "$@"; exit 0 ;;
  convert)
    mkdir -p "$OUT"
    cargo build --release -p ddg-integer-coords --example thingi_convert > "$OUT/build-convert.log" 2>&1 || exit 1
    "$BIN/thingi_convert" > "$OUT/convert.log" 2>&1 ;;
  run)
    shift; backends=("$@"); [ ${#backends[@]} -eq 0 ] && backends=(integer sp-ruppert sp-chew)
    cargo build --release -p ddg-integer-coords --example thingi_eval > "$OUT/build-eval.log" 2>&1 || exit 1
    mkdir -p "$OUT/jobs"
    tsv="$ROOT/data/thingi10k/converted.tsv"
    for phase in 1 2 3; do
      case $phase in
        1) lo=0; hi=100000; par="${IC_JOBS:-56}"; mem=4000000 ;;
        2) lo=100000; hi=1000000; par=12; mem=12000000 ;;
        3) lo=1000000; hi=999999999; par=3; mem=24000000 ;;
      esac
      for b in "${backends[@]}"; do
        awk -F'\t' -v lo=$lo -v hi=$hi -v b="$b" -v m=$mem 'NR>1 && $3=="kept" && $8>lo && $8<=hi {print $1, b, m}' "$tsv" |
          { if [ -n "${IC_IDS:-}" ]; then grep -wFf "$IC_IDS"; else cat; fi; } 
      done > "$OUT/phase$phase.list"
      echo "phase $phase: $(wc -l < "$OUT/phase$phase.list") jobs, $par parallel"
      xargs -P "$par" -L 1 "$0" --job < "$OUT/phase$phase.list"
    done ;;
  cmp)
    # Integer runs with the IDT edge lengths dumped (IC_DUMP_IDT), for the comparison with the demo
    # (scripts/ic-demo-eval.sh): meshes up to $2 faces (default 20000), results in runs_tmp/wp10/cmp/.
    mkdir -p "$OUT/cmp"
    awk -F'\t' -v m="${2:-20000}" 'NR>1 && $3=="kept" && $8<=m {print $1}' "$ROOT/data/thingi10k/converted.tsv" |
      xargs -P "${IC_JOBS:-56}" -I{} sh -c 'IC_DUMP_IDT="$1/cmp/{}.idt" timeout 600 "$2/thingi_eval" "$3/{}.obj" --backend "${IC_CMP_BACKEND:-integer-ref}" --id {} > "$1/cmp/{}.json" 2>/dev/null' _ "$OUT" "$BIN" "$CONV" ;;
  summary)
    uv run --no-project python -I "$ROOT/scripts/ic_eval_summary.py" "$ROOT" ;;
  *) sed -n 2,17p "$0"; exit 2 ;;
esac
