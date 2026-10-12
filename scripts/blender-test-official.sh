#!/usr/bin/env bash
# Run the add-on tests in an official Blender build: Blender-as-a-Python-module from PyPI.
#
#   scripts/blender-test-official.sh [BPY_VERSION [PYTHON_VERSION]] [-k NAME]
#
# Defaults: bpy 5.0.1 on Python 3.11 (what official Blender 5.0 bundles), with numpy installed
# as in official builds. The extension zip in dist/blender is used as it is (build it with
# scripts/blender-build-extension.sh, preferably with DDG_WHEEL_PORTABLE=1). `uv` downloads the
# `bpy` wheel (about 400 MB) into its cache on first use.
#
# The zip is installed into a fresh scratch user directory (runs_tmp/blender-user-official) and
# blender/tests/test_addon.py runs unchanged. The real user configuration is never touched.
# Log: runs_tmp/blender-test-official-<version>.log.
set -euo pipefail
source "$(dirname "${BASH_SOURCE[0]}")/blender-common.sh"

bpy_version="5.0.1"
python_version="3.11"
positional=0
extra=()
while [ $# -gt 0 ]; do
    case "$1" in
        -k) extra=(-k "$2"); shift ;;
        -h|--help) sed -n '2,14p' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//'; exit 0 ;;
        *)
            if [ "$positional" = 0 ]; then bpy_version="$1"; else python_version="$1"; fi
            positional=$((positional + 1)) ;;
    esac
    shift
done

zip="$root/dist/blender/ddg_uv_unwrap-$(ddg_version).zip"
[ -f "$zip" ] || { echo "no extension zip at $zip; run scripts/blender-build-extension.sh" >&2; exit 1; }
scratch="$root/runs_tmp/blender-user-official"
rm -rf "$scratch"
blender_scratch_env "$scratch"

log="$root/runs_tmp/blender-test-official-$bpy_version.log"
status=0
(cd "$root/runs_tmp" && uv run --quiet --no-project --python "$python_version" \
    --with "bpy==$bpy_version" --with numpy \
    python -I "$root/blender/tests/run_in_bpy_module.py" "$zip" "${extra[@]}") > "$log" 2>&1 || status=$?
grep -E '^(Blender .* as a module|  ddg_tools|  numpy|PASS|FAIL|SKIP|RESULT|ALL TESTS PASSED|  failed)' "$log" || true
if [ "$status" != 0 ] || ! grep -q '^ALL TESTS PASSED$' "$log"; then
    echo "tests in bpy $bpy_version FAILED (exit status $status); see $log" >&2
    exit 1
fi
