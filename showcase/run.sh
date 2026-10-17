#!/usr/bin/env bash
# Run a compute script with the ddg_tools wheel from dist/wheels.
here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
wheel="$(ls "$here"/../dist/wheels/ddg_tools-*.whl | head -1)"
exec uv run --no-project --python 3.12 --with "$wheel" --with numpy --with scipy python "$@"
