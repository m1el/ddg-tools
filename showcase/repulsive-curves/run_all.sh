#!/usr/bin/env bash
# Runs the full-size scenes in parallel (one process per scene, threads split), 8 at a time.
# Frames and reports go to showcase/work/rc/<name>/, logs to runs_tmp/rc-wp13/<name>.log.
# Usage: run_all.sh [name ...]   (default: all)
set -u
here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
root="$(cd "$here/../.." && pwd)"
bin="${RC_SCENES_BIN:-$root/target-rc-wp13/release/examples/scenes}"
S="$root/data/repulsive-curves/scenes"
F="$root/data/repulsive-curves/FreedmanHeWang"
out="$root/showcase/work/rc"
logs="$root/runs_tmp/rc-wp13"
mkdir -p "$out" "$logs"
# name | scene | threads | frame-every | extra arguments
jobs=(
"freedman1200|$S/FreedmanUnknot/scene.txt|12|10|--curve $F/freedman-he-wang-1200.obj --backend mg --steps 1500"
"bunny|$S/BunnyCurvePacking/scene.txt|8|10|"
"capillary|$S/CapillaryNetwork/scene.txt|8|5|--backend mg"
"rings|$S/RingsPacking/scene.txt|8|10|--steps 300"
"cup|$S/SimpleCurvePacking/circle-scene.txt|6|10|"
"rect|$S/SimpleCurvePacking/rect-scene.txt|6|10|"
"torus|$S/ImplicitTorus/scene.txt|4|10|"
"torus-multi|$S/ImplicitTorus/multi-scene.txt|4|5|"
"doubletorus1|$S/ImplicitDoubleTorus/scene1.txt|2|2|"
"doubletorus2|$S/ImplicitDoubleTorus/scene2.txt|2|2|"
"sphere|$S/ImplicitSphere/scene.txt|4|10|"
"bluenoise|$S/BluenoiseCurves/scene.txt|4|5|"
"k5|$S/GraphDrawing/k5scene.txt|2|5|--steps 150"
"k33|$S/GraphDrawing/k33scene.txt|2|5|--steps 150"
"planar-graph|$S/PlanarGraphDrawing/scene.txt|2|2|"
"interp|$S/CurveInterpolation/scene.txt|2|5|--steps 150"
"tangent|$S/TangentInterpolation/scene.txt|1|2|--steps 40"
)
run() {
    IFS='|' read -r name scene th ev extra <<< "$1"
    start=$(date +%s)
    # shellcheck disable=SC2086
    "$bin" "$scene" --out "$out/$name" --every "$ev" --threads "$th" $extra > "$logs/$name.log" 2>&1
    echo "$name rc=$? wall=$(( $(date +%s) - start )) s" >> "$logs/run_all.summary"
}
export -f run
export bin out logs
sel=("$@")
for j in "${jobs[@]}"; do
    n="${j%%|*}"
    if [ ${#sel[@]} -eq 0 ] || [[ " ${sel[*]} " == *" $n "* ]]; then printf '%s\n' "$j"; fi
done | xargs -P 8 -d '\n' -I{} bash -c 'run "$1"' _ {}
