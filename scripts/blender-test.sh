#!/usr/bin/env bash
# Build, package, install and test the Blender extension, headless.
#
#   scripts/blender-test.sh [--no-build] [--no-ui] [--bench] [-k NAME]
#
# 1. scripts/blender-build-extension.sh builds the wheel and the extension zip
#    (--no-build: use the zip that is already in dist/blender).
# 2. The zip is installed into a fresh scratch Blender user directory
#    (runs_tmp/blender-user-test) with `blender --command extension install-file`.
# 3. blender/tests/test_addon.py runs in that Blender, against the installed extension and its
#    bundled wheel. `-k NAME` runs only the tests whose name contains NAME.
# 3b. blender/tests/test_cut.py runs the same way (Start, Flow, Apply of the variational cut),
#    in parallel, once per variant (intrinsic, caller, refined, intrinsic and refined with
#    Optimal Delaunay on; logs
#    runs_tmp/blender-test-cut-<name>.log; -k applies to both files).
# 3c. blender/tests/test_developable.py runs the same way (Make Developable, Show Developability,
#    Developable Seams + UV; log runs_tmp/blender-test-developable.log), in parallel with 3b.
# 4. blender/tests/ui_smoke.py runs in a windowed Blender on a private virtual display
#    (xvfb-run): operator from the viewport, redo panel, menus, undo and redo. Skipped with a
#    notice when xvfb-run is not installed, or with --no-ui.
# 5. With --bench, blender/tests/bench_overhead.py measures the cost of the Blender side on
#    meshes of about 300k triangles.
#
# Every Blender call uses --factory-startup and the scratch user directory; the real user
# configuration is never read or written. Logs: runs_tmp/blender-test-*.log.
# Environment: BLENDER (default: blender), CARGO_TARGET_DIR (default: target-py).
set -euo pipefail
source "$(dirname "${BASH_SOURCE[0]}")/blender-common.sh"

build=1
bench=0
ui=1
only=()
while [ $# -gt 0 ]; do
    case "$1" in
        --no-build) build=0 ;;
        --bench) bench=1 ;;
        --no-ui) ui=0 ;;
        -k) only=(-k "$2"); shift ;;
        -h|--help) sed -n '2,20p' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//'; exit 0 ;;
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

# The cut tests run once per variant, the Blender processes in parallel; one log each: the
# triangulations intrinsic, caller, refined, and intrinsic and refined with Optimal Delaunay on
# by default ("intrinsic-odt", "refined-odt").
variants=(intrinsic caller refined intrinsic-odt refined-odt)
pids=()
for tri in "${variants[@]}"; do
    odt=()
    case "$tri" in *-odt) odt=(--optimal-delaunay) ;; esac
    "$BLENDER" -b --factory-startup --python-exit-code 1 \
        --python "$root/blender/tests/test_cut.py" -- --module "$module" \
        --triangulation "${tri%-odt}" "${odt[@]}" \
        "${only[@]}" > "$root/runs_tmp/blender-test-cut-$tri.log" 2>&1 &
    pids+=($!)
done
dev_log="$root/runs_tmp/blender-test-developable.log"
"$BLENDER" -b --factory-startup --python-exit-code 1 \
    --python "$root/blender/tests/test_developable.py" -- --module "$module" \
    "${only[@]}" > "$dev_log" 2>&1 &
dev_pid=$!
cut_failed=0
for i in "${!variants[@]}"; do
    tri="${variants[$i]}"
    status=0
    wait "${pids[$i]}" || status=$?
    log="$root/runs_tmp/blender-test-cut-$tri.log"
    echo "cut tests, triangulation $tri: $(grep -E '^RESULT' "$log")"
    grep -E '^(FAIL|  failed)' "$log" || true
    if [ "$status" != 0 ] || ! grep -q '^ALL TESTS PASSED$' "$log"; then
        echo "Blender cut tests FAILED for $tri (exit status $status); see $log" >&2
        cut_failed=1
    fi
done
dev_status=0
wait "$dev_pid" || dev_status=$?
echo "developable tests: $(grep -E '^RESULT' "$dev_log")"
grep -E '^(FAIL|  failed|SPOT)' "$dev_log" || true
if [ "$dev_status" != 0 ] || ! grep -q '^ALL TESTS PASSED$' "$dev_log"; then
    echo "Blender developable tests FAILED (exit status $dev_status); see $dev_log" >&2
    cut_failed=1
fi
[ "$cut_failed" = 0 ] || exit 1

if [ "$ui" = 1 ] && [ "${#only[@]}" -eq 0 ]; then
    if command -v xvfb-run > /dev/null; then
        log="$root/runs_tmp/blender-test-ui.log"
        status=0
        timeout 300 xvfb-run -a -s "-screen 0 1280x800x24" "$BLENDER" --factory-startup \
            --python "$root/blender/tests/ui_smoke.py" -- --module "$module" > "$log" 2>&1 || status=$?
        grep -E '^(UI: FAILED|UI RESULT|UI SMOKE TEST PASSED|  failed)' "$log" || true
        if [ "$status" != 0 ] || ! grep -q '^UI SMOKE TEST PASSED$' "$log" || grep -q 'Traceback' "$log"; then
            echo "UI smoke test FAILED (exit status $status); see $log" >&2
            exit 1
        fi
        echo "UI smoke test: $(grep -c '^UI: ok' "$log") checks passed"
    else
        echo "UI smoke test SKIPPED: xvfb-run is not installed"
    fi
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
    log="$root/runs_tmp/blender-test-bench-cut.log"
    if ! "$BLENDER" -b --factory-startup --python-exit-code 1 \
            --python "$root/blender/tests/bench_cut.py" -- --module "$module" > "$log" 2>&1; then
        echo "cut benchmark failed; see $log" >&2
        tail -n 30 "$log" >&2
        exit 1
    fi
    sed -n '/^BENCH BEGIN$/,/^BENCH END$/p' "$log" | sed '1d;$d'
fi
