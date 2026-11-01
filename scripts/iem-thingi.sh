#!/usr/bin/env bash
# WP9b: the Thingi10k benchmark of intrinsic simplification (ddg-intrinsic-simplification SPEC 10.2).
#
#   scripts/iem-thingi.sh run [MODES...]   # one process per (mesh, mode), in parallel
#   scripts/iem-thingi.sh summary          # $IEM_OUT/summary.md from the job results
#
# MODES: refine (pre-refinement, integer coordinates, 25 deg; default) noref (no pre-refinement),
# refine-sp (pre-refinement with signposts). Default: refine noref.
# Env: IEM_JOBS (parallel jobs, default and maximum 24), IEM_TIMEOUT (seconds per job, default 900),
# IEM_IDS (file with ids, default: every kept mesh), IEM_OUT (default runs_tmp/iem/thingi).
# Results: $IEM_OUT/jobs/<mode>/<id>.json and .log. A job killed by the time cap is recorded as
# {"status":"timeout"}, one that dies without a result as "memory" or "crash".
# Phases by size: <= 100k faces: IEM_JOBS jobs, 4 GB each; <= 1M faces: 12 jobs, 12 GB; larger: 3, 24 GB.
set -u
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
export CARGO_TARGET_DIR="${CARGO_TARGET_DIR:-$ROOT/target-iem-thingi}"
BIN="$CARGO_TARGET_DIR/release/examples"
CONV="$ROOT/data/thingi10k/converted"
OUT="${IEM_OUT:-$ROOT/runs_tmp/iem/thingi}"
TIMEOUT="${IEM_TIMEOUT:-900}"
JOBS="${IEM_JOBS:-24}"
[ "$JOBS" -gt 24 ] && JOBS=24

job() { # job ID MODE MEM_KB
  local id="$1" mode="$2" mem="$3" d="$OUT/jobs/$2" flags=""
  case "$mode" in noref) flags="--no-refine" ;; refine-sp) flags="--signposts" ;; esac
  mkdir -p "$d"
  [ -s "$d/$id.json" ] && return 0
  ( ulimit -v "$mem"; timeout --signal=KILL "$TIMEOUT" "$BIN/thingi_simplify" "$CONV/$id.obj" --id "$id" $flags \
      > "$d/$id.out" 2> "$d/$id.log"; echo "exit=$?" >> "$d/$id.log" )
  if grep -q '^exit=0' "$d/$id.log" && [ -s "$d/$id.out" ]; then
    mv "$d/$id.out" "$d/$id.json"
  else
    local code status
    code=$(grep '^exit=' "$d/$id.log" | tail -1 | cut -d= -f2)
    if [ "$code" = 137 ] || [ "$code" = 124 ]; then status=timeout
    elif grep -qiE 'memory allocation|allocat|capacity overflow' "$d/$id.log"; then status=memory
    else status=crash; fi
    echo "{\"id\":\"$id\",\"mode\":\"$mode\",\"status\":\"$status\",\"exit\":$code}" > "$d/$id.json"
    rm -f "$d/$id.out"
  fi
}

case "${1:-}" in
  --job) shift; job "$@"; exit 0 ;;
  run)
    shift; modes=("$@"); [ ${#modes[@]} -eq 0 ] && modes=(refine noref)
    mkdir -p "$OUT/jobs"
    cargo build --release -p ddg-intrinsic-simplification --example thingi_simplify > "$OUT/build.log" 2>&1 || exit 1
    tsv="$ROOT/data/thingi10k/converted.tsv"
    for phase in 1 2 3; do
      case $phase in
        1) lo=0; hi=100000; par="$JOBS"; mem=4000000 ;;
        2) lo=100000; hi=1000000; par=$(( JOBS < 12 ? JOBS : 12 )); mem=12000000 ;;
        3) lo=1000000; hi=999999999; par=$(( JOBS < 3 ? JOBS : 3 )); mem=24000000 ;;
      esac
      for b in "${modes[@]}"; do
        awk -F'\t' -v lo=$lo -v hi=$hi -v b="$b" -v m=$mem 'NR>1 && $3=="kept" && $8>lo && $8<=hi {print $1, b, m}' "$tsv" |
          { if [ -n "${IEM_IDS:-}" ]; then grep -wFf "$IEM_IDS"; else cat; fi; }
      done > "$OUT/phase$phase.list"
      echo "phase $phase: $(wc -l < "$OUT/phase$phase.list") jobs, $par parallel"
      xargs -P "$par" -L 1 "$0" --job < "$OUT/phase$phase.list"
    done ;;
  summary)
    IEM_OUT="$OUT" uv run --no-project python -I "$ROOT/scripts/iem_thingi_summary.py" "$ROOT" "$OUT" ;;
  *) sed -n 2,14p "$0"; exit 2 ;;
esac
