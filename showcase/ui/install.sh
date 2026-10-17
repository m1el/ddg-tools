#!/usr/bin/env bash
# Stage blender/ddg_uv_unwrap with the newest wheel from dist/wheels into showcase/work/ui-env,
# build the extension zip there and install it into a scratch Blender user directory.
set -euo pipefail
source "$(dirname "${BASH_SOURCE[0]}")/../../scripts/blender-common.sh"
env_dir="$root/showcase/work/ui-env"
stage="$env_dir/stage/ddg_uv_unwrap"
rm -rf "$env_dir/stage" "$env_dir/user"
mkdir -p "$env_dir/stage" "$env_dir/zip"
wheel="$(ls -t "$root"/dist/wheels/ddg_tools-*.whl | head -1)"
uv run --quiet --no-project --python '>=3.11' "$root/scripts/blender_stage_extension.py" \
    --source "$root/blender/ddg_uv_unwrap" --out "$stage" --license "$root/LICENSE" "$wheel"
blender_scratch_env "$env_dir/user"
"$BLENDER" --factory-startup --command extension build --source-dir "$stage" --output-dir "$env_dir/zip"
"$BLENDER" --factory-startup --command extension install-file -r user_default -e "$env_dir"/zip/ddg_uv_unwrap-*.zip
