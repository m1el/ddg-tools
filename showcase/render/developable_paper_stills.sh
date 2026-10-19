#!/usr/bin/env bash
# Panels of the developability stills into work/dev2/p/ (compose with render/developable_paper_compose.py).
cd "$(dirname "$0")/.."
T=work/dev2/p; mkdir -p $T
R() { while [ $(jobs -r | wc -l) -ge ${DEV_JOBS:-6} ]; do sleep 2; done; render/developable_paper.sh "$T/$1.log" panel "${@:2}" & }
for m in "spot_iso -50 16" "bob_iso -40 22" "blub_iso -40 20" "spot_cage -50 16" "crumpled -30 38"; do
  set -- $m
  R $1_in  $1 in  $T/$1_in.png  $2 $3 plain 800 800
  R $1_out $1 out $T/$1_out.png $2 $3 plain 800 800
  R $1_def $1 out $T/$1_def.png $2 $3 defect 800 800
done
for e in covariance combinatorial combinatorial-max covariance-max intrinsic intrinsic-max; do
  R e_$e spot_e_$e out $T/e_$e.png -50 16 plain 700 700
done
wait
