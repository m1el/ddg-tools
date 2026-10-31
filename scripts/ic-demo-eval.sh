#!/usr/bin/env bash
# WP10, SPEC 10.3: runs the authors' demo (built headless in runs_tmp/ic-demo-build/, see
# scripts/ic-demo-build.sh) on the converted Thingi10k meshes, as the paper's benchmark/run_benchmark.py
# does: `int_tri MESH --{flip,refine}Delaunay --noGUI --logStats --refineMaxInsertions=0`, once per
# backend (integer, signpost) and operation, one process per job, in parallel, 10-minute cap.
# For meshes up to IC_SMALL faces (default 20000) the flip run also writes faceLengths.dmat and the
# common subdivision, reduced here to sorted IDT edge lengths and V/E/F (the big files are deleted).
#   scripts/ic-demo-eval.sh [MAX_FACES]      (default 200000)     env: IC_JOBS (56), IC_TIMEOUT (600)
set -u
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
BIN="$ROOT/runs_tmp/ic-demo-build/build/bin/int_tri"
CONV="$ROOT/data/thingi10k/converted"
OUT="$ROOT/runs_tmp/wp10/demo"
TIMEOUT="${IC_TIMEOUT:-600}"
SMALL="${IC_SMALL:-20000}"

job() { # job ID BACKEND OP FACES
  local id="$1" b="$2" op="$3" nf="$4" d="$OUT/$2"
  mkdir -p "$d"
  local p="$d/${id}_${op}_"
  [ -e "${p}done" ] && return 0
  local extra=""
  if [ "$op" = flip ] && [ "$nf" -le "$SMALL" ]; then extra="--intrinsicFaces --commonSubdivision"; fi
  ( ulimit -v 12000000; timeout --signal=KILL "$TIMEOUT" "$BIN" "$CONV/$id.obj" --backend "$b" --${op}Delaunay --noGUI --logStats \
      --refineMaxInsertions=0 $extra --outputPrefix="$p" > "${p}out.log" 2>&1; echo "exit=$?" >> "${p}out.log" )
  if [ -e "${p}faceLengths.dmat" ]; then
    grep -v '^#' "${p}faceLengths.dmat" | tr ' ' '\n' | grep -v '^$' | sort -g > "${p}idt_lengths.txt"
    rm -f "${p}faceLengths.dmat" "${p}faceInds.dmat"
  fi
  if [ -e "${p}common_subdivision.obj" ]; then
    awk '$1=="v"{v++} $1=="f"{f++; for(i=2;i<=NF;i++){a=$i+0;b=$(i<NF?i+1:2)+0; k=(a<b)?a" "b:b" "a; e[k]=1}} END{n=0; for(k in e)n++; print v, n, f}' "${p}common_subdivision.obj" > "${p}vef.txt"
    rm -f "${p}common_subdivision.obj"
  fi
  touch "${p}done"
}

case "${1:-}" in
  --job) shift; job "$@"; exit 0 ;;
esac
MAXF="${1:-200000}"
tsv="$ROOT/data/thingi10k/converted.tsv"
mkdir -p "$OUT"
for op in flip refine; do for b in integer signpost; do
  awk -F'\t' -v m="$MAXF" -v b="$b" -v op="$op" 'NR>1 && $3=="kept" && $8<=m {print $1, b, op, $8}' "$tsv"
done; done > "$OUT/jobs.list"
echo "$(wc -l < "$OUT/jobs.list") jobs"
xargs -P "${IC_JOBS:-56}" -L 1 "$0" --job < "$OUT/jobs.list"
