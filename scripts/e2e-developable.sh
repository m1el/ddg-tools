#!/usr/bin/env bash
# End-to-end measurement of the developability operators on real meshes, headless, in parallel.
#
#   scripts/e2e-developable.sh [--no-sweep] [OUTDIR]
#
# Needs the installed extension of scripts/blender-test.sh (runs_tmp/blender-user-test). For each
# mesh (spot, bob, blub: quad, triangle and control meshes; Blender primitives) and each number of
# subdivision rounds (2 = the defaults, 1) one Blender process runs Make Developable and then
# Developable Seams + UV (blender/tests/e2e_developable.py). Output: OUTDIR (default
# runs_tmp/e2e-developable)/<case>.jsonl and .log; the table: scripts/e2e_table.py OUTDIR.
set -euo pipefail
source "$(dirname "${BASH_SOURCE[0]}")/blender-common.sh"
sweep=(--sweep)
if [ "${1:-}" = "--no-sweep" ]; then sweep=(); shift; fi
out="${1:-$root/runs_tmp/e2e-developable}"
rm -rf "$out"; mkdir -p "$out"
blender_scratch_env "$root/runs_tmp/blender-user-test"
models="$root/data/models"
meshes=(
    "spot_quad:$models/spot/spot/spot_quadrangulated.obj"
    "spot_tri:$models/spot/spot/spot_triangulated.obj"
    "spot_control:$models/spot/spot/spot_control_mesh.obj"
    "bob_quad:$models/bob/bob/bob_quad.obj"
    "bob_tri:$models/bob/bob/bob_tri.obj"
    "bob_control:$models/bob/bob/bob_controlmesh.obj"
    "blub_quad:$models/blub/blub/blub_quadrangulated.obj"
    "blub_tri:$models/blub/blub/blub_triangulated.obj"
    "blub_control:$models/blub/blub/blub_control_mesh.obj"
)
jobs=()
for m in "${meshes[@]}"; do
    for s in 2 1 -1; do jobs+=("${m%%:*}|--mesh ${m#*:}|$s"); done
done
for p in cube sphere torus monkey cylinder plane; do
    for s in 2 1 -1; do jobs+=("prim_$p|--primitive $p|$s"); done
done
run_one() {
    IFS='|' read -r name args subs <<< "$1"
    tag="$name-s$subs"
    "$BLENDER" -b --factory-startup --python-exit-code 1 \
        --python "$root/blender/tests/e2e_developable.py" -- \
        --module bl_ext.user_default.ddg_uv_unwrap $args --label "$name" --subdivisions "$subs" \
        ${SWEEP:-} --out "$out/$tag.jsonl" > "$out/$tag.log" 2>&1 || echo "FAILED $tag"
}
export -f run_one
export RAYON_NUM_THREADS="${RAYON_NUM_THREADS:-2}"
export root BLENDER out SWEEP="${sweep[*]:-}"
printf '%s\n' "${jobs[@]}" | xargs -P "${E2E_JOBS:-24}" -I{} bash -c 'run_one "{}"'
echo "done: $out"
