#!/usr/bin/env bash
# developable_blender.sh <log> <blender args after --> : run render/developable.py, retrying
# while the shared GPU is out of memory.
cd "$(dirname "$0")/.."
log="$1"; shift
for try in $(seq 1 30); do
    blender -b --python render/developable.py -- "$@" > "$log" 2>&1
    if ! grep -q "out of GPU memory\|Out of memory\|Traceback" "$log"; then exit 0; fi
    grep -qi "out of .*memory" "$log" || exit 1
    sleep $((5 + RANDOM % 10))
done
exit 1
