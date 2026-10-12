#!/usr/bin/env bash
# Differential test of ddg-core's polygon triangulation against a real Blender.
#
#   scripts/blender-tessellation-diff.sh [SEED] [SCALE]
#
# 1. Blender (headless, factory settings, user resources redirected to a scratch directory so
#    that nothing touches the user's configuration) triangulates random quads and n-gons and
#    writes its `loop_triangles` to runs_tmp/blender-tessellation-<seed>.txt.
# 2. `cargo test -p ddg-core --test blender_tessellation` triangulates the same polygons and
#    prints the agreement per polygon family. The test fails if any non-degenerate polygon is
#    split differently.
#
# Environment: BLENDER (default: blender), CARGO_TARGET_DIR (default: target-core).
set -euo pipefail

seed="${1:-1}"
scale="${2:-1.0}"
root="$(cd "$(dirname "$0")/.." && pwd)"
blender="${BLENDER:-blender}"
dump="$root/runs_tmp/blender-tessellation-$seed.txt"
log="$root/runs_tmp/blender-tessellation-$seed.blender.log"
mkdir -p "$root/runs_tmp"

scratch="$(mktemp -d)"
trap 'rm -rf "$scratch"' EXIT
if ! BLENDER_USER_RESOURCES="$scratch" "$blender" -b --factory-startup \
    --python "$root/scripts/blender_tessellation_dump.py" -- "$dump" "$seed" "$scale" > "$log" 2>&1; then
    echo "Blender failed; see $log" >&2
    exit 1
fi
grep -E "^(Blender|wrote)" "$log" || true

export CARGO_TARGET_DIR="${CARGO_TARGET_DIR:-$root/target-core}"
cd "$root"
DDG_BLENDER_DUMP="$dump" cargo test --release -p ddg-core --test blender_tessellation \
    external_dump_agreement -- --nocapture
