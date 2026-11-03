#!/usr/bin/env bash
# Build the Blender extension zip, with the native module bundled as a wheel.
#
#   scripts/blender-build-extension.sh
#
# 1. scripts/py-build-wheel.sh builds the wheel of `ddg_tools` (skipped with DDG_SKIP_WHEEL=1;
#    then the wheels already in dist/wheels are used).
# 2. scripts/blender_stage_extension.py stages blender/ddg_uv_unwrap with every
#    dist/wheels/ddg_tools-<version>-*.whl and a manifest that lists them.
# 3. `blender --command extension validate` and `... build` produce
#    dist/blender/ddg_uv_unwrap-<version>.zip. The last line printed is that path.
#
# With DDG_WHEEL_CUDA=1 the bundled wheel is the CUDA one (scripts/py-build-wheel.sh: cargo feature
# `cuda`, dist/wheels-cuda/) and the zip goes to dist/blender-cuda/, next to (not over) the CPU
# extension: the Repulsive Surfaces panel then offers Device CUDA / Auto when the machine has an
# NVIDIA GPU and driver. The wheel is not portable (glibc of this machine) unless
# DDG_WHEEL_PORTABLE=1 is set too.
#
# Blender runs with --factory-startup and a scratch user directory (runs_tmp/blender-user-build).
# Environment: BLENDER (default: blender), CARGO_TARGET_DIR (default: target-py).
set -euo pipefail
source "$(dirname "${BASH_SOURCE[0]}")/blender-common.sh"

version="$(ddg_version)"
out="$root/dist/blender"
wheel_dir="$root/dist/wheels"
if [ "${DDG_WHEEL_CUDA:-0}" = 1 ]; then
    out="$root/dist/blender-cuda"
    wheel_dir="$root/dist/wheels-cuda"
fi
stage="$out/stage/ddg_uv_unwrap"
log="$root/runs_tmp/blender-build-extension.log"
mkdir -p "$out" "$root/runs_tmp"

if [ "${DDG_SKIP_WHEEL:-0}" != 1 ]; then
    "$root/scripts/py-build-wheel.sh" > /dev/null
fi
shopt -s nullglob
wheels=("$wheel_dir"/ddg_tools-"$version"-*.whl)
shopt -u nullglob
if [ "${#wheels[@]}" -eq 0 ]; then
    echo "no wheel $wheel_dir/ddg_tools-$version-*.whl; run scripts/py-build-wheel.sh" >&2
    exit 1
fi

# Ask a wheel that runs on this machine whether it reports a stand-in pipeline (the current
# wheel never does: ddg_tools.IS_PLACEHOLDER is False). `-I` and a
# neutral working directory keep the repository out of the import path.
placeholder_flag=()
is_placeholder="$(cd "$root/runs_tmp" && uv run --quiet --no-project --python '>=3.11' --with "${wheels[0]}" \
    python -I -c 'import ddg_tools; print(int(ddg_tools.IS_PLACEHOLDER))')"
if [ "$is_placeholder" = 1 ]; then
    placeholder_flag=(--placeholder)
fi

uv run --quiet --no-project --python '>=3.11' "$root/scripts/blender_stage_extension.py" \
    --source "$root/blender/ddg_uv_unwrap" --out "$stage" --license "$root/LICENSE" \
    "${placeholder_flag[@]}" "${wheels[@]}"

blender_scratch_env "$root/runs_tmp/blender-user-build"
{
    "$BLENDER" --factory-startup --command extension validate "$stage"
    "$BLENDER" --factory-startup --command extension build \
        --source-dir "$stage" --output-dir "$out"
} > "$log" 2>&1 || { echo "building the extension failed; see $log" >&2; tail -n 30 "$log" >&2; exit 1; }

zip="$out/ddg_uv_unwrap-$version.zip"
if [ ! -f "$zip" ]; then
    echo "expected $zip; see $log" >&2
    exit 1
fi
echo "$zip"
