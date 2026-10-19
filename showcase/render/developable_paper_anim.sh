#!/usr/bin/env bash
# The coarse-to-fine animation: frames (retried on GPU out-of-memory), then mp4 + gif.
#   render/developable_paper_anim.sh  (needs work/dp_spot_iso.npz, work/dp_snap_spot.npz)
cd "$(dirname "$0")/.."
export SHOW_KEEP=1 DEV_FULL=spot_iso
for try in $(seq 1 40); do
    blender -b --python render/developable_paper.py -- flow spot developable_flow_spot -50 16 40 4 > work/dev2/anim.log 2>&1 && ! grep -q "Traceback" work/dev2/anim.log && break
    sleep 8
done
python3 post.py video work/frames/developable_flow_spot media/developable_flow_spot work/labels_developable_flow_spot.json 30 720
