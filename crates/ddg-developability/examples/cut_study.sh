#!/bin/sh
# Cut study (SPEC §12.3, SPEC D30): the flow is run once per (mesh, flow energy,
# budget) and written to OBJ; every cut variant (terminal energy x epsilon x defect tolerance) is
# then a cheap `--cut-only` job on that OBJ. One process and one log per job (xargs -P).
# Usage (repository root):
#   cargo build --release -p ddg-developability --example develop
#   sh crates/ddg-developability/examples/cut_study.sh [PARALLEL=56]
# Logs and the per-job RESULT lines: runs_tmp/dev-cut-study/{flow,cut}/ ; aggregate with
#   uv run crates/ddg-developability/examples/cut_study.py runs_tmp/dev-cut-study
set -e
BIN=${CARGO_TARGET_DIR:-target}/release/examples/develop
OUT=${OUT:-runs_tmp/dev-cut-study}
# EXTRA: more develop options for every flow (e.g. EXTRA=--no-fallback)
EXTRA=${EXTRA:-}
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
done < "$OUT/meshes.txt" | EXTRA="$EXTRA" xargs -P "$PAR" -L 1 sh -c \
  'n=$(basename "$0" .obj)-$2-$3; "'"$BIN"'" "$0" --subdivisions "$1" --energy "$2" --iterations "$3" --threads 1 $EXTRA --no-cut --out "'"$OUT"'/flow/$n.obj" > "'"$OUT"'/flow/$n.log" 2>&1 || true'
# phase 2: cuts
for f in "$OUT"/flow/*.obj; do
  for ce in covariance intrinsic combinatorial; do
    for eps in 1e-2 1e-3 1e-4; do
      for tol in off 0.1 0.02; do
        echo "$f $ce $eps $tol"
      done
    done
  done
done | xargs -P "$PAR" -L 1 sh -c \
  'n=$(basename "$0" .obj)_$1_$2_$3; t=""; [ "$3" != off ] && t="--defect-tol $3"; "'"$BIN"'" "$0" --cut-only --cut-energy "$1" --epsilon "$2" $t --threads 1 > "'"$OUT"'/cut/$n.log" 2>&1 || true'
