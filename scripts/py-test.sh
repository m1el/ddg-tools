#!/usr/bin/env bash
# Test the wheel of `ddg_tools` under plain Python, without Blender.
#
#   scripts/py-test.sh [--no-build]
#
# Builds the wheel (scripts/py-build-wheel.sh; skipped with --no-build, then the newest
# dist/wheels/ddg_tools-<version>-*.whl is used) and runs
# crates/ddg-py/python-tests/test_ddg_tools.py and test_cut_session.py with `uv run` in throw-away environments:
#
#   Python 3.11 with numpy     what official Blender builds bundle
#   Python 3.14 without numpy  what this machine's distribution build of Blender runs
#   Python 3.12 without numpy  one in between, same binary
#
# The same wheel file is used for all three (stable ABI). Logs: runs_tmp/py-test-*.log.
# Environment: DDG_PY_TEST_MATRIX overrides the list, e.g. "3.13: 3.13:numpy".
set -euo pipefail
source "$(dirname "${BASH_SOURCE[0]}")/blender-common.sh"

if [ "${1:-}" = "--no-build" ]; then
    wheel="$(ls -t "$root"/dist/wheels/ddg_tools-"$(ddg_version)"-*.whl 2>/dev/null | head -n 1)"
    [ -n "$wheel" ] || { echo "no wheel in dist/wheels; run without --no-build" >&2; exit 1; }
else
    wheel="$("$root/scripts/py-build-wheel.sh" | tail -n 1)"
fi
echo "wheel: $wheel"

matrix="${DDG_PY_TEST_MATRIX:-3.11:numpy 3.14: 3.12:}"
status=0
for entry in $matrix; do
    python="${entry%%:*}"
    extra="${entry#*:}"
    with=(--with "$wheel")
    [ -z "$extra" ] || with+=(--with "$extra")
    log="$root/runs_tmp/py-test-$python${extra:+-$extra}.log"
    # --refresh-package: never test a stale copy of a wheel that was rebuilt under the same name.
    # -I and a neutral working directory keep the repository out of the import path.
    if (cd "$root/runs_tmp" && for test in test_ddg_tools test_cut_session; do
            DDG_WHEEL="$wheel" uv run --quiet --no-project --python "$python" \
                --refresh-package ddg-tools "${with[@]}" \
                python -I "$root/crates/ddg-py/python-tests/$test.py" || exit 1
        done) > "$log" 2>&1; then
        echo "PASS python $python ${extra:-(no numpy)}: $(grep -E '^Ran [0-9]+ tests' "$log" | sed 's/ in .*//' | tr '\n' ';') $(grep -cE '^OK' "$log") files OK"
    else
        echo "FAIL python $python ${extra:-(no numpy)}: see $log"
        status=1
    fi
done
exit $status
