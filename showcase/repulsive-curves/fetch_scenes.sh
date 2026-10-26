#!/usr/bin/env bash
# Copies the authors' scenes (repulsive-curves/scenes, MIT, Copyright (c) 2019 Christopher Yu)
# to data/repulsive-curves/scenes/ (gitignored). Their .gitignore'd media (flow/, frames/) is skipped.
set -euo pipefail
here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
root="$here/../.."
src="${RC_SCENES:-$root/../RepulsiveCurves/repulsive-curves/scenes}"
dst="$root/data/repulsive-curves/scenes"
mkdir -p "$dst"
for d in "$src"/*/; do
    n="$(basename "$d")"
    mkdir -p "$dst/$n"
    cp -u "$d"/*.obj "$d"/*.txt "$dst/$n/" 2>/dev/null || true
done
cp "$src/FORMATS.md" "$dst/"
echo "scenes in $dst"
