#!/usr/bin/env bash
# Panels of the cut and flatten figures into work/dev2/c/ (compose with developable_paper_compose.py cuts).
cd "$(dirname "$0")/.."
T=work/dev2/c; mkdir -p $T
R() { while [ $(jobs -r | wc -l) -ge ${DEV_JOBS:-6} ]; do sleep 2; done; render/developable_paper.sh "$T/$1.log" "${@:2}" & }
for m in "spot_iso -50 16 130 16" "bob_iso -40 22 140 22"; do
  set -- $m
  R $1_front chart3d $1_split $T/$1_front.png $2 $3 1000 800
  R $1_back  chart3d $1_split $T/$1_back.png  $4 $5 1000 800
  R $1_lay0  layout  $1_nosplit $T/$1_lay0.png 1000
  R $1_lay1  layout  $1_split   $T/$1_lay1.png 1000
done
R crumpled_3d chart3d crumpled_split $T/crumpled_3d.png -30 38 1000 800
R crumpled_lay layout crumpled_split $T/crumpled_lay.png 1000
wait
