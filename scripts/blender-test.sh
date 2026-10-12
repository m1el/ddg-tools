#!/usr/bin/env bash
# Build, package, install and test the Blender extension, headless.
#
#   scripts/blender-test.sh [--no-build] [--bench] [-k NAME]
#
# 1. scripts/blender-build-extension.sh builds the wheel and the extension zip
#    (--no-build: use the zip that is already in dist/blender).
# 2. The zip is installed into a fresh scratch Blender user directory
#    (runs_tmp/blender-user-test) with `blender --command extension install-file`.
# 3. blender/tests/test_addon.py runs in that Blender, against the installed extension and its
#    bundled wheel. `-k NAME` runs only the tests whose name contains NAME.
# 4. With --bench, blender/tests/bench_overhead.py measures the cost of the Blender side on
#    meshes of about 300k triangles.
#
# Every Blender call uses --factory-startup and the scratch user directory; the real user
# configuration is never read or written. Logs: runs_tmp/blender-test-*.log.
# Environment: BLENDER (default: blender), CARGO_TARGET_DIR (default: target-py).
set -euo pipefail
source "$(dirname "${BASH_SOURCE[0]}")/blender-common.sh"

build=1
bench=0
only=()
while [ $# -gt 0 ]; do
    case "$1" in
        --no-build) build=0 ;;
        --bench) bench=1 ;;
        -k) only=(-k "$2"); shift ;;
        -h|--help) sed -n '2,18p' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//'; exit 0 ;;
        *) echo "unknown argument: $1" >&2; exit 2 ;;
    esac
    shift
done

zip="$root/dist/blender/ddg_uv_unwrap-$(ddg_version).zip"
if [ "$build" = 1 ]; then
    zip="$("$root/scripts/blender-build-extension.sh" | tail -n 1)"
fi
[ -f "$zip" ] || { echo "no extension zip at $zip; run without --no-build" >&2; exit 1; }
echo "extension: $zip"
"$BLENDER" --version | head -n 1

scratch="$root/runs_tmp/blender-user-test"
rm -rf "$scratch"
blender_scratch_env "$scratch"
module="bl_ext.user_default.ddg_uv_unwrap"

log="$root/runs_tmp/blender-test-install.log"
if ! "$BLENDER" --factory-startup --command extension install-file -r user_default -e "$zip" \
        > "$log" 2>&1 || ! grep -q 'Installed "ddg_uv_unwrap"' "$log"; then
    echo "installing the extension failed; see $log" >&2
    tail -n 20 "$log" >&2
    exit 1
fi
echo "installed into $scratch"

log="$root/runs_tmp/blender-test.log"
status=0
"$BLENDER" -b --factory-startup --python-exit-code 1 \
    --python "$root/blender/tests/test_addon.py" -- --module "$module" --zip "$zip" "${only[@]}" \
    > "$log" 2>&1 || status=$?
grep -E '^(PASS|FAIL|SKIP|RESULT|ALL TESTS PASSED|  failed)' "$log" || true
if [ "$status" != 0 ] || ! grep -q '^ALL TESTS PASSED$' "$log"; then
    echo "Blender tests FAILED (exit status $status); see $log" >&2
    exit 1
fi

if [ "$bench" = 1 ]; then
    log="$root/runs_tmp/blender-test-bench.log"
    if ! "$BLENDER" -b --factory-startup --python-exit-code 1 \
            --python "$root/blender/tests/bench_overhead.py" -- --module "$module" > "$log" 2>&1; then
        echo "benchmark failed; see $log" >&2
        tail -n 30 "$log" >&2
        exit 1
    fi
    sed -n '/^BENCH BEGIN$/,/^BENCH END$/p' "$log" | sed '1d;$d'
fi
