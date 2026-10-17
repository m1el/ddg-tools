#!/usr/bin/env bash
# Side-by-side still of the heat method: input mesh | intrinsic Delaunay refinement, labelled with
# the mean relative error against the exact polyhedral distance (compute/geodesic_error.py).
#   render/geodesic_compose.sh <mesh> <input.png> <refined.png> <added-vertices> <out.png>
set -euo pipefail
cd "$(dirname "$0")/.."
mesh=$1; left=$2; right=$3; added=$4; out=$5
read -r ein eref < <(python3 -c "import json; e=json.load(open('work/geodesic_error_$mesh.json')); print(f\"{100*e['input']:.0f} {100*e['refined']:.1f}\")")
panel() { # in title subtitle out
  convert "$1" -background '#f3f4f6' -alpha remove -alpha off \
    -font DejaVu-Sans -fill '#1f2329' -gravity north -pointsize 40 -annotate +0+40 "$2" \
    -fill '#5a6270' -pointsize 28 -annotate +0+98 "$3" "$4"
}
panel "$left" "Heat method on the input mesh" "mean error $ein %%" work/geo_l.png
panel "$right" "… on its intrinsic Delaunay refinement" "mean error $eref %%  ·  +$added vertices" work/geo_r.png
convert work/geo_l.png work/geo_r.png +append "$out"
echo "wrote $out"
