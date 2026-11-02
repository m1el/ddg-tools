#!/usr/bin/env bash
# Copies the reference's scenes (repulsive-surfaces/scenes: scene files + OBJs, MIT, Copyright (c) 2019
# Nick Sharp) to data/repulsive-surfaces/scenes/ (gitignored). Idempotent. RS_REFERENCE overrides the
# location of the reference checkout.
set -euo pipefail
here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
root="$here/.."
ref="${RS_REFERENCE:-/var/data/discrete-geometry/RepulsiveSurfaces/repulsive-surfaces}"
[ -d "$ref/scenes" ] || ref="$root/../RepulsiveSurfaces/repulsive-surfaces"
src="$ref/scenes"
dst="$root/data/repulsive-surfaces/scenes"
[ -d "$src" ] || { echo "no scenes directory at $src (set RS_REFERENCE)" >&2; exit 1; }
mkdir -p "$dst"
for d in "$src"/*/; do
    n="$(basename "$d")"
    mkdir -p "$dst/$n"
    cp -u "$d"/*.obj "$d"/*.txt "$dst/$n/" 2>/dev/null || true
done
cp -u "$src/FORMAT.txt" "$dst/"
if [ -f "$ref/LICENSE" ]; then cp -u "$ref/LICENSE" "$dst/LICENSE-reference-MIT.txt"; fi
echo "scenes in $dst"
