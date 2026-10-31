#!/usr/bin/env bash
# WP10, SPEC 10.2: single-thread timings, run ALONE on a quiet machine, one job after the other, pinned
# to one core (CORE, default 7). For the 50-mesh subset (every 155th kept mesh by id) and the paper's
# meshes (44395 Fig. 15, 54674 Fig. 16, 719791 Fig. 19) it runs our integers, our signposts (Ruppert)
# and the authors' demo (integer and signpost backends, refinement), and records the load average
# before and after every job.   scripts/ic-timings.sh [subset|paper|all]
set -u
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
EX="${CARGO_TARGET_DIR:-$ROOT/target-wp10}/release/examples"
DEMO="$ROOT/runs_tmp/ic-demo-build/build/bin/int_tri"
CONV="$ROOT/data/thingi10k/converted"
OUT="$ROOT/runs_tmp/wp10/timings"
CORE="${CORE:-7}"
mkdir -p "$OUT"
awk -F'\t' 'NR>1 && $3=="kept"{print $1}' "$ROOT/data/thingi10k/converted.tsv" | sort -n | awk 'NR%155==1' > "$OUT/subset50.txt"
printf '44395\n54674\n719791\n' > "$OUT/paper.txt"
case "${1:-all}" in subset) lists=(subset50) ;; paper) lists=(paper) ;; *) lists=(paper subset50) ;; esac
for l in "${lists[@]}"; do
  while read -r id; do
    echo "$id load_before $(uptime | sed "s/.*load average: //")" >> "$OUT/load.log"
    for b in integer sp-ruppert; do
      taskset -c "$CORE" timeout --signal=KILL 1800 "$EX/thingi_eval" "$CONV/$id.obj" --backend "$b" --id "$id" > "$OUT/${id}_ours_$b.json" 2> "$OUT/${id}_ours_$b.err"
    done
    for b in integer signpost; do
      for op in flip refine; do
        p="$OUT/${id}_demo_${b}_${op}_"
        taskset -c "$CORE" timeout --signal=KILL 1800 "$DEMO" "$CONV/$id.obj" --backend "$b" --${op}Delaunay --noGUI --logStats --refineMaxInsertions=0 --outputPrefix="$p" > "${p}out.log" 2>&1
      done
    done
    echo "$id load_after $(uptime | sed "s/.*load average: //")" >> "$OUT/load.log"
  done < "$OUT/$l.txt"
done
