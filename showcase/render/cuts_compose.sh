#!/usr/bin/env bash
# Compose the cut stills in media/ from the renders in work/cuts_test/ (see render/cuts.py).
set -euo pipefail
cd "$(dirname "$0")/.."
T=work/cuts_test
M=media
BG="#f3f4f6"; INK="#1f2329"; F=DejaVu-Sans
flat() { convert "$1" -background "$BG" -alpha remove -alpha off "$2"; }
strip() {  # in out text pointsize: a label strip above the image
    convert "$1" -background "$BG" -gravity north -splice 0x$(( $4 * 2 )) -font $F -pointsize "$4" \
        -fill "$INK" -gravity northwest -annotate +28+$(( $4 / 2 )) "$3" "$2"
}
label() {  # in out text pointsize [gravity]
    convert "$1" -font $F -pointsize "$4" -fill "$INK" -gravity "${5:-northwest}" -annotate +28+22 "$3" "$2"
}
mkdir -p $M

# Atlas: converged network on spot, the same at the start, bob.
flat $T/atlas_spot8.png $M/cuts_atlas.png
flat $T/atlas_bob6.png $M/cuts_atlas_bob.png
flat $T/atlas0_spot8.png $T/a0.png; flat $T/atlas_spot8.png $T/a1.png
label $T/a0.png $T/a0l.png "start: normal clusters" 44
label $T/a1.png $T/a1l.png "after the cut flow" 44
convert $T/a0l.png $T/a1l.png +append -resize 1800x $M/cuts_atlas_start_end.png

# Unwrap: the cut surface with a checker and its UV layout.
flat $T/unwrap_spot8.png $T/u0.png; flat $T/layout_spot8.png $T/u1.png
label $T/u0.png $T/u0l.png "cut along the region boundaries" 40
label $T/u1.png $T/u1l.png "conformal UV layout" 40
convert $T/u0l.png $T/u1l.png +append -resize 1800x $M/cuts_unwrap.png

# Visibility: without / with the penalty, from the camera and from the other side.
for f in vis_spot2_100 vis_vis_10_100 vis_spot2_-80 vis_vis_10_-80; do flat $T/${f}.png $T/${f}_f.png; done
strip $T/vis_spot2_100_f.png $T/v00.png "no penalty: seen from the camera" 52
strip $T/vis_vis_10_100_f.png $T/v01.png "visibility penalty: seen from the camera" 52
strip $T/vis_spot2_-80_f.png $T/v10.png "no penalty: other side" 52
strip $T/vis_vis_10_-80_f.png $T/v11.png "visibility penalty: other side" 52
convert \( $T/v00.png $T/v01.png +append \) \( $T/v10.png $T/v11.png +append \) -append \
    -resize 1800x $M/cuts_visibility.png
ls -la $M/cuts_*
