#!/usr/bin/env bash
# Crop the captures of shoot.sh into showcase/media/blender_* (stills, GIFs, MP4s).
# Expects: shots-77 (unwrap), shots-78 (cut, 2 regions), shots-80 (cut, 4 regions), shots-79 (sphere).
set -euo pipefail
here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
S="$here/../work/ui-env"
M="$here/../media"
mkdir -p "$M"
UV=1576x1055+0+0       # UV Editing workspace without properties editor and status bar
LAY=1576x977+0+0       # Layout workspace without properties editor and timeline
VIEW=1100x917+200+60   # the 3D viewport around the model
c() { convert "$1" -crop "$2" +repage "$3"; echo "wrote $3"; }
c "$S"/shots-77/unwrap_*_result.png  $UV  "$M/blender_unwrap.png"
c "$S"/shots-77/unwrap_*_stretch.png $UV  "$M/blender_unwrap_stretch.png"
c "$S"/shots-77/unwrap_*_menu.png    $UV  "$M/blender_unwrap_menu.png"
c "$S"/shots-77/unwrap_*_options.png $UV  "$M/blender_unwrap_options.png"
c "$S"/shots-78/cut_30_flow30.png    $LAY "$M/blender_cut.png"
c "$S"/shots-78/cut_*_options.png    $LAY "$M/blender_cut_options.png"
c "$S"/shots-78/cut_*_applied.png    $UV  "$M/blender_cut_applied.png"
c "$S"/shots-80/cut_30_flow30.png    $LAY "$M/blender_cut_atlas.png"
c "$S"/shots-80/cut_*_applied.png    $UV  "$M/blender_cut_atlas_applied.png"
c "$S"/shots-79/sphere_*_options.png $LAY "$M/blender_sphere.png"

anim() { # name, frames..., (forward, hold, back)
    local name="$1"; shift
    local tmp="$S/anim-$name"; rm -rf "$tmp"; mkdir -p "$tmp"
    local i=0 f
    local list=("$@")
    for f in "${list[@]}"; do convert "$f" -crop $VIEW +repage "$tmp/$(printf %04d $i).png"; i=$((i+1)); done
    local last="${list[-1]}"
    for _ in $(seq 12); do cp "$tmp/$(printf %04d $((${#list[@]}-1))).png" "$tmp/$(printf %04d $i).png"; i=$((i+1)); done
    for ((k=${#list[@]}-1; k>=0; k--)); do cp "$tmp/$(printf %04d $k).png" "$tmp/$(printf %04d $i).png"; i=$((i+1)); done
    for _ in $(seq 8); do cp "$tmp/0000.png" "$tmp/$(printf %04d $i).png"; i=$((i+1)); done
    ffmpeg -y -loglevel error -framerate "$FPS" -i "$tmp/%04d.png" -vf "scale=1100:-2" -c:v libx264 \
        -pix_fmt yuv420p -crf 20 -movflags +faststart "$M/$name.mp4"
    ffmpeg -y -loglevel error -framerate "$FPS" -i "$tmp/%04d.png" -filter_complex \
        "scale=720:-1:flags=lanczos,split[a][b];[a]palettegen=max_colors=128[p];[b][p]paletteuse=dither=sierra2_4a" \
        "$M/$name.gif"
    echo "wrote $M/$name.mp4 $M/$name.gif"
}
FPS=10 anim blender_cut_flow $(ls "$S"/shots-78/cut_*_flow*.png | sort)
FPS=10 anim blender_cut_atlas_flow $(ls "$S"/shots-80/cut_*_flow*.png | sort)
FPS=20 anim blender_sphere_morph $(ls "$S"/shots-79/sphere_*_morph*.png | sort)
