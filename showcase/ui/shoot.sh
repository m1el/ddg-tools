#!/usr/bin/env bash
# Run one scene of showcase/ui/drive.py in a windowed Blender on a private Xvfb display.
#   showcase/ui/shoot.sh <scene> [display-number]
set -euo pipefail
source "$(dirname "${BASH_SOURCE[0]}")/../../scripts/blender-common.sh"
scene="$1"; disp="${2:-77}"
env_dir="$root/showcase/work/ui-env"
blender_scratch_env "$env_dir/user"
blank="$env_dir/blank.blend"
[ -f "$blank" ] || "$BLENDER" -b --factory-startup --python-expr \
    "import bpy; bpy.ops.wm.save_as_mainfile(filepath='$blank')" > /dev/null
Xvfb ":$disp" -screen 0 1920x1080x24 -nolisten tcp > "$env_dir/xvfb-$disp.log" 2>&1 &
xpid=$!
trap 'kill $xpid 2>/dev/null || true' EXIT
sleep 1
rm -f "$env_dir/shots-$disp/"*.png
shift 2 || true
DISPLAY=":$disp" timeout 600 "$BLENDER" --factory-startup --enable-event-simulate -p 0 0 1920 1080 "$blank" \
    --python "$root/showcase/ui/drive.py" -- "$scene" "$env_dir/shots-$disp" "$@"
