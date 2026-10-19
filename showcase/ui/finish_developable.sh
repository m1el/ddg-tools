#!/usr/bin/env bash
# Crop the captures of `shoot.sh developable 81` into showcase/media/blender_developable_*.png.
set -euo pipefail
here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
S="$here/../work/ui-env/shots-81"
M="$here/../media"
LAY=1576x977+0+0       # Layout workspace without properties editor and timeline
UV=1576x1055+0+0       # UV Editing workspace without properties editor and status bar
c() { convert "$1" -crop "$2" +repage "$3"; echo "wrote $3"; }
c "$S"/developable_01_result.png         $LAY "$M/blender_developable_result.png"
c "$S"/developable_02_options.png        $LAY "$M/blender_developable_options.png"
c "$S"/developable_03_energy.png         $LAY "$M/blender_developable_energy.png"
c "$S"/developable_04_energy_options.png $LAY "$M/blender_developable_energy_options.png"
c "$S"/developable_05_seams.png          $UV  "$M/blender_developable_seams.png"
c "$S"/developable_06_seams_options.png  $UV  "$M/blender_developable_seams_options.png"
