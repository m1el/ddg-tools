#!/usr/bin/env bash
# Fig. 12 analogues: sphere-traced renders with the reference weights, CPU, 8 threads.
# Usage: showcase/points-as-tori/render/render.sh   (from the repo root; builds the example first)
set -euo pipefail
cd "$(dirname "$0")/../../.."
D=${PAT_DATA:-/var/data/discrete-geometry/PointsAsTori/points-as-tori/data}
OUT=showcase/points-as-tori/render
LOG=runs_tmp/pat-wp11
mkdir -p "$LOG"
export CARGO_TARGET_DIR=${CARGO_TARGET_DIR:-$PWD/target-pat11}
nice -n 19 cargo build --release -p ddg-points-as-tori -j 8 --example render > "$LOG/render-build.log" 2>&1
R=$CARGO_TARGET_DIR/release/examples/render
r() { nice -n 19 "$R" "$@"; }
{
r $D/tooth.ply $OUT/tooth-level0.png res=1024
r $D/tooth.ply $OUT/tooth-offset-0.08.png res=1024 level=0.08
r $D/tooth.ply $OUT/tooth-offset-0.16.png res=1024 level=0.16
r $D/snail.ply $OUT/snail-level0.png res=1024 az=40 el=20
r $D/snail.ply $OUT/snail-offset-csg-0.05.png res=1024 az=40 el=20 offset=0.05
r $D/holey-bunny.ply $OUT/bunny-level0.png res=1024 az=-30 el=15
r $D/holey-bunny.ply $OUT/bunny-cut-plane.png res=1024 az=-30 el=15 cut=2 cut_offset=-0.1
r analytic:5 $OUT/analytic5-level0.png res=1024
r analytic:5 $OUT/boolean-subtract.png res=1024 op=subtract other=analytic:30 other_shift=0.35,0.1,0
r analytic:5 $OUT/boolean-union.png res=1024 op=union other=analytic:30 other_shift=0.45,0.1,0
r $D/tooth.ply $OUT/tooth-demo-constants.png res=1024 demo=1
} > "$LOG/showcase-render.log" 2>&1
cat "$LOG/showcase-render.log"
