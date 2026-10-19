#!/usr/bin/env bash
# The developability flows of the showcase, in parallel, one thread each (single-thread timings are
# only meaningful on a quiet machine: see `uptime` in work/dp/load.txt). Run from anywhere.
#   compute/dev_paper_runs.sh            flows + analysis
# Needs: cargo build --release -p ddg-developability --example develop (CARGO_TARGET_DIR below),
#        work/inputs/{spot,bob,blub}_iso1000.obj from
#        uv run crates/ddg-developability/examples/coarse_inputs.py showcase/work/inputs 1000
set -u
cd "$(dirname "$0")/../.."
BIN=${DEVELOP_BIN:-target-show-dev/release/examples/develop}
W=showcase/work; D=$W/dp; mkdir -p $D
M=data/models
# the library default: FlowOptions::coarse_to_fine(2, 500)
REC="combinatorial:500,+covariance:500,+covariance:500"
uptime > $D/load.txt
cat > $D/jobs.txt <<JOBS
spot_iso $W/inputs/spot_iso1000.obj $M/spot/spot/spot_triangulated.obj $REC
bob_iso $W/inputs/bob_iso1000.obj $M/bob/bob/bob_tri.obj $REC
blub_iso $W/inputs/blub_iso1000.obj $M/blub/blub/blub_triangulated.obj $REC
spot_cage $M/spot/spot/spot_control_mesh.obj $M/spot/spot/spot_triangulated.obj $REC
crumpled $M/crumpleddevelopable/CrumpledDevelopable.obj $M/crumpleddevelopable/CrumpledDevelopable.obj $REC
JOBS
for e in covariance combinatorial combinatorial-max covariance-max intrinsic intrinsic-max; do
  echo "spot_e_$e $W/inputs/spot_iso1000.obj $M/spot/spot/spot_triangulated.obj $e:500,+$e:500,+$e:500" >> $D/jobs.txt
done
run() { tag=$1; in=$2; ref=$3; rounds=$4
  $BIN $in --rounds $rounds --threads 1 --no-cut --out $D/$tag.obj --history $D/$tag.hist > $D/$tag.log 2>&1
  uv run --no-project --with numpy --with scipy python showcase/compute/dev_paper.py analyze $tag $in $D/$tag.obj $ref $D/$tag.log >> $D/$tag.log 2>&1
}
export -f run; export BIN D W
xargs -P 12 -L 1 bash -c 'run "$@"' _ < $D/jobs.txt
