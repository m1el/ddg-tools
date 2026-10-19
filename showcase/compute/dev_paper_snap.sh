#!/usr/bin/env bash
# States of the coarse-to-fine recipe for the animation: the flow is deterministic, so the state at
# (round r, iteration k) is the flow whose last round is cut at k iterations; one process per state.
#   compute/dev_paper_snap.sh TAG INPUT.obj
set -u
cd "$(dirname "$0")/../.."
BIN=${DEVELOP_BIN:-target-show-dev/release/examples/develop}
tag=$1; in=$2; D=showcase/work/dp; mkdir -p $D
rm -f $D/snap_${tag}_*
R0="combinatorial"; R1="covariance"
jobs=""; k=0
add() { jobs="$jobs$k $1 $2 $3\n"; k=$((k+1)); }
# k round iteration rounds-spec
for it in 1 3 8 20 50 120 250 500; do add 0 $it "combinatorial:$it"; done
for it in 5 20 60 150 500; do add 1 $it "combinatorial:500,+covariance:$it"; done
for it in 5 20 60 150 500; do add 2 $it "combinatorial:500,+covariance:500,+covariance:$it"; done
snap() { k=$1; r=$2; it=$3; spec=$4
  $BIN $in --rounds $spec --threads 1 --no-cut --out $D/snap_${tag}_$k.obj > $D/snap_${tag}_$k.log 2>&1
  e=$(grep -o "^RESULT.*" $D/snap_${tag}_$k.log | grep -o "e_final=[^ ]*" | cut -d= -f2)
  e0=$(grep -o "^RESULT.*" $D/snap_${tag}_$k.log | grep -o "e_initial=[^ ]*" | cut -d= -f2)
  echo "{\"k\": $k, \"round\": $r, \"iteration\": $it, \"energy\": $e, \"energy_initial\": $e0}" > $D/snap_${tag}_$k.json
}
export -f snap; export BIN in D tag
printf "$jobs" | xargs -P 20 -L 1 bash -c 'snap "$@"' _
uv run --no-project --with numpy --with scipy python showcase/compute/dev_paper.py states $tag
