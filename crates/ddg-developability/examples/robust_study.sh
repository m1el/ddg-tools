#!/bin/sh
# Robustness study (SPEC §12, SPEC D31/D32): the flows of the cut study re-run with the admissibility floor, then the cut (covariance terminals, eps 1e-3, defect tolerance 0.1) on each result.
# One process and one log per job (xargs -P). Same flows as examples/cut_study.sh.
# Usage (repository root):
#   cargo build --release -p ddg-developability --example develop
#   sh crates/ddg-developability/examples/robust_study.sh [PARALLEL=56]
# Logs and the per-job RESULT lines: runs_tmp/dev-robust/study/{flow,cut}/ ; aggregate with
#   uv run crates/ddg-developability/examples/robust_summary.py runs_tmp/dev-robust/study/flow
set -e
BIN=${CARGO_TARGET_DIR:-target}/release/examples/develop
OUT=runs_tmp/dev-robust/study
PAR=${1:-56}
M=data/models
mkdir -p "$OUT/flow" "$OUT/cut"
rm -f "$OUT"/flow/* "$OUT"/cut/*
# mesh, subdivisions, budgets
cat > "$OUT/meshes.txt" <<LIST
$M/spot/spot/spot_control_mesh.obj 2 300 3000
$M/bob/bob/bob_controlmesh.obj 2 300 3000
$M/blub/blub/blub_control_mesh.obj 2 300 3000
$M/crumpleddevelopable/CrumpledDevelopable.obj 0 500 3000
$M/spot/spot/spot_triangulated.obj 1 300
$M/spot/spot/spot_quadrangulated.obj 1 300
$M/bob/bob/bob_tri.obj 1 300
$M/bob/bob/bob_quad.obj 1 300
$M/bob/bob/bob_isotropic.obj 1 300
$M/blub/blub/blub_triangulated.obj 1 300
$M/blub/blub/blub_quadrangulated.obj 1 300
LIST
# phase 1: flows
while read -r f sub budgets; do
  for e in covariance intrinsic combinatorial; do
    for b in $budgets; do echo "$f $sub $e $b"; done
  done
done < "$OUT/meshes.txt" | xargs -P "$PAR" -L 1 sh -c \
  'n=$(basename "$0" .obj)-$2-$3; "'"$BIN"'" "$0" --subdivisions "$1" --energy "$2" --iterations "$3" --threads 1 --no-cut --out "'"$OUT"'/flow/$n.obj" > "'"$OUT"'/flow/$n.log" 2>&1 || true'
# phase 2: one cut per flow (the variant of the cut study that is the default of D30)
for f in "$OUT"/flow/*.obj; do echo "$f"; done | xargs -P "$PAR" -L 1 sh -c \
  'n=$(basename $0 .obj); "'"$BIN"'" "$0" --cut-only --cut-energy covariance --epsilon 1e-3 --defect-tol 0.1 --threads 1 > "'"$OUT"'/cut/$n.log" 2>&1 || true'
