# Shared by the blender-*.sh and py-*.sh scripts. Source it; do not run it.
#
# Defines:
#   root                       the repository root
#   BLENDER                    the Blender executable (default: blender)
#   CARGO_TARGET_DIR           default: $root/target-py
#   blender_scratch_env DIR    point every Blender user directory and the temporary directory
#                              at DIR, so that nothing reads or writes the real user
#                              configuration (~/.config/blender) or session files in /tmp
#   ddg_version                version of the ddg-py crate (the wheel's version)

root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BLENDER="${BLENDER:-blender}"
export CARGO_TARGET_DIR="${CARGO_TARGET_DIR:-$root/target-py}"

blender_scratch_env() {
    local dir="$1"
    mkdir -p "$dir"
    export BLENDER_USER_RESOURCES="$dir"
    export BLENDER_USER_CONFIG="$dir/config"
    export BLENDER_USER_SCRIPTS="$dir/scripts"
    export BLENDER_USER_EXTENSIONS="$dir/extensions"
    export BLENDER_USER_DATAFILES="$dir/datafiles"
    # Blender's temporary directory: a windowed Blender writes its "recover last session" file
    # (quit.blend) there on exit, and the real one in /tmp belongs to the user.
    mkdir -p "$dir/tmp"
    export TMPDIR="$dir/tmp"
    export TMP="$dir/tmp"
}

ddg_version() {
    # The crate inherits its version from the workspace.
    sed -n '/^\[workspace.package\]/,/^\[/ s/^version *= *"\(.*\)"/\1/p' "$root/Cargo.toml"
}
