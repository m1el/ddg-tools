#!/usr/bin/env bash
# retry_frames.sh <n-frames> <frames-dir> <blender args...>: rerun until every frame exists.
n=$1; dir=$2; shift 2
for i in $(seq 1 40); do
  have=$(ls $dir 2>/dev/null | grep -c '^[0-9]*\.png$')
  [ "$have" -ge "$n" ] && exit 0
  SHOW_KEEP=1 blender -b --python "$@"
  sleep 3
done
