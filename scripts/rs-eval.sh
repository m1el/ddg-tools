#!/usr/bin/env bash
# WP16B: method comparison (paper Fig. 5 and Fig. 12) for ddg-repulsive-surfaces.
#   scripts/rs-eval.sh gen      generate the meshes (runs_tmp/wp16b/meshes)
#   scripts/rs-eval.sh run      launch all runs in parallel (one thread each), wait, write logs/CSVs
# Env: BUDGET (Fig. 12 wall seconds, default 2400; tunnels 3600), F5_BUDGET (Fig. 5 wall cap, default 5400).
set -u
cd "$(dirname "$0")/.."
export CARGO_TARGET_DIR=${CARGO_TARGET_DIR:-$PWD/target-wp16b}
BIN=$CARGO_TARGET_DIR/release/examples/compare_methods
OUT=runs_tmp/wp16b
M=$OUT/meshes
mkdir -p "$M"
BUDGET=${BUDGET:-2400}
F5_BUDGET=${F5_BUDGET:-5400}
METHODS="hs l2 h1 h2 aqp h1-lbfgs bcqn"

if [ "${1:-}" = gen ]; then
  pids=()
  for i in 1 2 3 4 5; do "$BIN" gen fig12:$i $M/fig12_$i.obj > $OUT/gen_fig12_$i.log 2>&1 & pids+=($!); done
  for f in 11000 40000 135000; do "$BIN" gen genus:2:$f $M/genus2_$f.obj > $OUT/gen_genus2_$f.log 2>&1 & pids+=($!); done
  wait "${pids[@]}"
  cat $OUT/gen_*.log
  exit 0
fi

echo "load at start: $(uptime)" > $OUT/load.txt
pids=()
for i in 1 2 3 4 5; do
  b=$BUDGET; if [ $i -ge 4 ]; then b=$((BUDGET * 3 / 2)); fi
  for m in $METHODS; do
    "$BIN" $M/fig12_$i.obj $m 0 1 --budget $b --csv $OUT/fig12_${i}_$m.csv \
      --exact-time 60,300,1200,2400,3600 --exact-end > $OUT/fig12_${i}_$m.log 2>&1 &
    pids+=($!)
  done
done
for f in 11000 40000 135000; do
  for m in $METHODS; do
    ex=""; if [ $f -le 40000 ]; then ex="--exact-steps 100,300"; fi
    "$BIN" $M/genus2_$f.obj $m 300 1 --budget $F5_BUDGET --csv $OUT/fig5_${f}_$m.csv $ex > $OUT/fig5_${f}_$m.log 2>&1 &
    pids+=($!)
  done
done
echo "launched ${#pids[@]}"
wait "${pids[@]}"
echo "load at end: $(uptime)" >> $OUT/load.txt
echo done
