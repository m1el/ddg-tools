#!/usr/bin/env bash
# Run a cuts render, retrying while the shared GPU is out of memory. Args: log, then blender args.
log="$1"; shift
for try in $(seq 1 30); do
    blender -b --python render/cuts.py -- "$@" > "$log" 2>&1
    if ! grep -q "Out of memory\|Traceback" "$log"; then exit 0; fi
    grep -q "Out of memory" "$log" || exit 1
    sleep $((5 + RANDOM % 10))
done
exit 1
