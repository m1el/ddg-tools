#!/usr/bin/env bash
# Renders of the developability stills into work/dev/ (then developable_compose.sh), in parallel.
cd "$(dirname "$0")/.."
T=work/dev; mkdir -p $T
R() { while [ $(jobs -r | wc -l) -ge ${DEV_JOBS:-5} ]; do sleep 2; done; render/developable_blender.sh "$T/$1.log" "${@:2}" & }
R hero      pair f_spot_i6000 $T/hero.png -50 16 plain 2000 1000
R energy    pair f_spot_i6000 $T/energy.png -50 16 energy 2000 1000
R energy_bob pair f_bob_i4000 $T/energy_bob.png -40 24 energy 2000 1000
R hero_bob  pair f_bob_i4000 $T/hero_bob.png -40 24 plain 2000 1000
R hero_blub pair f_blub_i3000 $T/hero_blub.png -40 20 plain 2000 1000
R crumpled  pair hero_crumpled $T/crumpled.png -30 38 plain 2000 1000
R crumpled_e pair hero_crumpled $T/crumpled_energy.png -30 38 energy 2000 1000
DEV_COLOR=defect R defect cut3d f_spot_i6000 $T/defect.png -50 16 1 1800 1200
DEV_COLOR=defect R defect_back cut3d f_spot_i6000 $T/defect_back.png 130 16 1 1800 1200
for e in covariance combinatorial combinatorial_max covariance_max intrinsic intrinsic_max; do
    R cmp_$e single cmp_$e $T/cmp_$e.png -50 16 plain 1000 1000
done
R cmp_in single cmp_covariance $T/cmp_in.png -50 16 plain 1000 1000 in
R cut_bob0 cut3d f_bob_i4000 $T/cut_bob0.png -40 24 0 1400 1000
R cut_bob1 cut3d f_bob_i4000 $T/cut_bob1.png -40 24 1 1400 1000
R cut_bob2 cut3d f_bob_i4000 $T/cut_bob2.png -40 24 2 1400 1000
for k in 0 1 2; do R lay_bob$k layout f_bob_i4000 $T/lay_bob$k.png $k 1000; done
R cut_cr cut3d hero_crumpled $T/cut_cr.png -30 38 0 1400 1000
R lay_cr layout hero_crumpled $T/lay_cr.png 0 1000
wait
