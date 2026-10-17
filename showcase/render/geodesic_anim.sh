#!/usr/bin/env bash
# Expanding-front animation of the heat method: input mesh | intrinsic Delaunay refinement,
# rendered as two 640x800 halves, joined and labelled with the measured mean errors
# (compute/geodesic_error.py), then encoded by post.py.
#   render/geodesic_anim.sh <mesh> <out-basename>
set -euo pipefail
cd "$(dirname "$0")/.."
mesh=$1; out=$2
for t in input refined; do
  [ "${KEEP:-}" = 1 ] && continue
  rm -rf "work/frames/geodesic_${mesh}_$t"
  blender -b --python render/geodesic.py -- "$mesh" $t anim 640 800 > "../runs_tmp/show_geo_anim_$t.log" 2>&1 &
done
wait
for t in input refined; do
  n=$(ls "work/frames/geodesic_${mesh}_$t" | wc -l)
  [ "$n" = 150 ] || { echo "only $n frames for $t (see runs_tmp/show_geo_anim_$t.log)"; exit 1; }
done
read -r ein eref < <(python3 -c "import json; e=json.load(open('work/geodesic_error_$mesh.json')); print(f\"{100*e['input']:.0f} {100*e['refined']:.1f}\")")
j=work/frames/geodesic_joined; rm -rf $j; mkdir -p $j
lab="-font DejaVu-Sans -fill #1f2329 -pointsize 24 -gravity north"
ls "work/frames/geodesic_${mesh}_input" | xargs -P 32 -I{} convert \
  \( "work/frames/geodesic_${mesh}_input/{}" $lab -annotate +0+30 "input mesh  ·  mean error $ein %%" \) \
  \( "work/frames/geodesic_${mesh}_refined/{}" $lab -annotate +0+30 "intrinsic Delaunay refinement  ·  $eref %%" \) \
  +append "$j/{}"
python3 post.py video $j "$out" - 30 720
