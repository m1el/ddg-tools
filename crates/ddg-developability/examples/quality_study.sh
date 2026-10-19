#!/bin/sh
# Quality study (QUALITY.md): flows of a job list, each measured by examples/quality.py.
# Usage (repository root), after `cargo build --release -p ddg-developability --example develop`:
#   sh crates/ddg-developability/examples/quality_study.sh JOBS OUTDIR [PARALLEL=60]
# JOBS: one job per line, `tag input.obj reference.obj [develop options...]` (# comments allowed);
# the reference is the mesh distances are measured to (the dense model). One log per job:
# OUTDIR/<tag>.log with the develop RESULT line and the QUALITY line; OUTDIR/<tag>.{obj,png}.
# Aggregate: uv run crates/ddg-developability/examples/quality_summary.py OUTDIR
set -e
BIN=${CARGO_TARGET_DIR:-target}/release/examples/develop
JOBS=$1
OUT=$2
PAR=${3:-60}
Q=crates/ddg-developability/examples/quality.py
mkdir -p "$OUT"
# (trailing blanks are stripped: xargs -L would join such a line with the next)
grep -v '^\s*#' "$JOBS" | grep -v '^\s*$' | sed 's/[[:space:]]*$//' | xargs -P "$PAR" -L 1 sh -c '
  tag=$0; in=$1; ref=$2; shift 2
  "'"$BIN"'" "$in" "$@" --threads 1 --no-cut --out "'"$OUT"'/$tag.obj" > "'"$OUT"'/$tag.log" 2>&1 || true
  if [ -f "'"$OUT"'/$tag.obj" ]; then
    uv run --quiet "'"$Q"'" "$in" "'"$OUT"'/$tag.obj" --reference "$ref" --tag "$tag" --png "'"$OUT"'/$tag.png" >> "'"$OUT"'/$tag.log" 2>&1 || true
  fi'
