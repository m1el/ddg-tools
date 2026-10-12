#!/usr/bin/env bash
# Build the wheel of the Python module `ddg_tools` (crates/ddg-py) into dist/wheels/.
#
#   scripts/py-build-wheel.sh [extra maturin arguments]
#
# One wheel serves every Python from 3.11 on (stable ABI, tag cp311-abi3). The last line printed
# is the path of the wheel.
#
# Environment: CARGO_TARGET_DIR (default: target-py). The build log is
# runs_tmp/py-build-wheel.log.
set -euo pipefail
source "$(dirname "${BASH_SOURCE[0]}")/blender-common.sh"

out="$root/dist/wheels"
log="$root/runs_tmp/py-build-wheel.log"
mkdir -p "$out" "$root/runs_tmp"

if ! uvx maturin build --release --manifest-path "$root/crates/ddg-py/Cargo.toml" \
        --out "$out" "$@" > "$log" 2>&1; then
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
