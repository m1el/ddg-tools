#!/usr/bin/env bash
# Frames and videos of the two flow animations (one Blender each, retried on GPU out-of-memory).
cd "$(dirname "$0")/.."
export SHOW_KEEP=1
anim() {  # log frames-name args...
    local log=$1; shift
    for try in $(seq 1 40); do
        blender -b --python render/developable.py -- flow "$@" > "$log" 2>&1 && ! grep -q "Traceback" "$log" && return 0
        sleep 8
    done
}
anim work/dev_anim_spot.log spot_flow developable_flow_spot -50 16 40 4 &
anim work/dev_anim_crumpled.log crumpled developable_flow_crumpled -30 38 30 3 &
wait
python3 post.py video work/frames/developable_flow_spot media/developable_flow_spot work/labels_developable_flow_spot.json 30 720
python3 post.py video work/frames/developable_flow_crumpled media/developable_flow_crumpled work/labels_developable_flow_crumpled.json 30 720
