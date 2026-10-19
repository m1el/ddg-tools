#!/bin/sh
# End-to-end study (SPEC §12.5): every OBJ under data/models x all six energies x iteration
# budgets, one process per job (xargs -P), one log per job in runs_tmp/dev-flow-study/.
# Usage (from the repository root):
#   cargo build --release -p ddg-developability --example develop
#   sh crates/ddg-developability/examples/study.sh [PARALLEL=60] [BUDGETS="100 400"]
set -e
BIN=${CARGO_TARGET_DIR:-target}/release/examples/develop
OUT=runs_tmp/dev-flow-study
PAR=${1:-60}
BUDGETS=${2:-"100 400"}
mkdir -p "$OUT"
rm -f "$OUT"/job-*.log
for f in $(find data/models -name '*.obj' -not -path '*__MACOSX*' | sort); do
  for e in combinatorial combinatorial-max covariance covariance-max intrinsic intrinsic-max; do
    for b in $BUDGETS; do
      echo "$f $e $b"
    done
  done
done | xargs -P "$PAR" -L 1 sh -c \
  'n=$(basename "$0" .obj); "'"$BIN"'" "$0" --energy "$1" --subdivisions 1 --iterations "$2" --threads 1 > "'"$OUT"'/job-$n-$1-$2.log" 2>&1 || true'
