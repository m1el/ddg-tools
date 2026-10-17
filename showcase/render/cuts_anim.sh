#!/usr/bin/env bash
# Frames of both cut animations, one Blender at a time; SHOW_KEEP skips frames already there.
cd "$(dirname "$0")/.."
export SHOW_KEEP=1
anim() {  # log name args...
    local log=$1; shift
    for try in $(seq 1 40); do
        blender -b --python render/cuts.py -- "$@" > "$log" 2>&1 && ! grep -q "Traceback" "$log" && return 0
        echo "retry $try" >> "$log.retries"; sleep 10
    done
}
CUTS_NFLOW=192 anim ../runs_tmp/show_cuts_anim_flow.log anim spot2 cuts_flow 100 15 50
CUTS_NFLOW=230 anim ../runs_tmp/show_cuts_anim_atlas.log anim spot8 cuts_atlas -130 20 60
echo done
