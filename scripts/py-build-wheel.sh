#!/usr/bin/env bash
# Build the wheel of the Python module `ddg_tools` (crates/ddg-py) into dist/wheels/.
#
#   scripts/py-build-wheel.sh [extra maturin arguments]
#
# One wheel serves every Python from 3.11 on (stable ABI, tag cp311-abi3). The last line printed
# is the path of the wheel.
#
# By default the wheel is linked against this machine's C library and tagged accordingly (for
# example manylinux_2_34: it needs glibc 2.34 or newer). With DDG_WHEEL_PORTABLE=1 the build
# uses zig as the linker (fetched by uv, as the `zig` extra of maturin) and targets glibc 2.28,
# which is what official Blender builds for Linux require; use that for a wheel to hand out.
#
# With DDG_WHEEL_CUDA=1 the wheel is built with the `cuda` feature (the Repulsive Surfaces flow
# can run its products on an NVIDIA GPU: device="cuda"/"auto"; nothing CUDA is needed to build,
# the wheel loads the driver and libnvrtc.so.12 at run time). It goes to dist/wheels-cuda/ (same
# file name as the CPU wheel, so the two never overwrite each other) and its log is
# runs_tmp/py-build-wheel-cuda.log.
#
# Environment: CARGO_TARGET_DIR (default: target-py), DDG_WHEEL_PORTABLE, DDG_WHEEL_CUDA. The
# build log is runs_tmp/py-build-wheel.log.
set -euo pipefail
source "$(dirname "${BASH_SOURCE[0]}")/blender-common.sh"

out="$root/dist/wheels"
log="$root/runs_tmp/py-build-wheel.log"
features=()
if [ "${DDG_WHEEL_CUDA:-0}" = 1 ]; then
    out="$root/dist/wheels-cuda"
    log="$root/runs_tmp/py-build-wheel-cuda.log"
    features=(--features python,cuda)
fi
mkdir -p "$out" "$root/runs_tmp"

maturin=(uvx maturin)
portable=()
if [ "${DDG_WHEEL_PORTABLE:-0}" = 1 ]; then
    maturin=(uvx --from 'maturin[zig]' maturin)
    portable=(--zig --compatibility manylinux_2_28)
fi
# One wheel per version and machine: a wheel left over from a build with other settings must
# not be bundled next to the new one.
rm -f "$out"/ddg_tools-"$(ddg_version)"-*_"$(uname -m)".whl

if ! "${maturin[@]}" build --release --manifest-path "$root/crates/ddg-py/Cargo.toml" \
        --out "$out" "${portable[@]}" "${features[@]}" "$@" > "$log" 2>&1; then
    echo "maturin failed; see $log" >&2
    tail -n 30 "$log" >&2
    exit 1
fi
wheel="$(sed -n 's/.*Built wheel.* to \(.*\.whl\)$/\1/p' "$log" | tail -n 1)"
if [ -z "$wheel" ] || [ ! -f "$wheel" ]; then
    echo "could not find the wheel maturin built; see $log" >&2
    exit 1
fi
case "$wheel" in
    /*) ;;
    *) wheel="$PWD/$wheel" ;;
esac
echo "$wheel"
