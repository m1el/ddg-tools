#!/usr/bin/env bash
# Black-box comparison of ddg-repulsive-surfaces with the authors' program (WP18, SPEC 14.7).
#
#   scripts/rs-compare-reference.sh NAME SCENE STEPS [OBJ_EVERY] [OUTDIR]
#
# Runs, in parallel and with one log each: the reference (`rsurfaces_headless`, 1 thread, the
# `_nonorientable` binary when NAME starts with "kb"), ours with `reference_compat` and without CCD
# (`ours_NAME_compat`), and ours with the defaults (`ours_NAME_default`); for "kb*" also the
# variants without splits and collapses (`noedge`, `compat-noedge`). Output (default
# runs_tmp/wp18b/): ref_NAME/{steps.csv,step_*.obj}, ours_NAME_MODE.{csv,obj}, ours_NAME_MODE_NNNNN.obj,
# one .log per run. Afterwards `compare_reference stats` on the dumps and
# `scripts/rs_compare_tables.py OUTDIR NAME` give the tables of crates/ddg-repulsive-surfaces/COMPARISON.md.
# Timings taken this way are contended; for ms/step re-run each alone with `--threads 1`.
#
# Needs: scenes (scripts/rs-fetch-scenes.sh), the reference build in $RS_BUILD
# (default /var/data/discrete-geometry/RepulsiveSurfaces-build), and
#   CARGO_TARGET_DIR=$PWD/target-wp18b cargo build --release -p ddg-repulsive-surfaces --example compare_reference
set -euo pipefail
name="$1"; scene="$2"; steps="$3"; every="${4:-10}"; out="${5:-runs_tmp/wp18b}"
build="${RS_BUILD:-/var/data/discrete-geometry/RepulsiveSurfaces-build}"
tgt="${CARGO_TARGET_DIR:-target-wp18b}"
mkdir -p "$out"
ref="$build/build/bin/rsurfaces_headless"
modes="compat default"
case "$name" in kb*) ref="$build/build/bin/rsurfaces_headless_nonorientable"; modes="compat compat-noedge default noedge";; esac
"$ref" --scene "$scene" --iterations "$steps" --obj-every "$every" --out "$out/ref_$name" --threads 1 \
    > "$out/ref_$name.log" 2>&1 &
for m in $modes; do
    "$tgt/release/examples/compare_reference" run "$scene" --mode "$m" --steps "$steps" --obj-every "$every" \
        --out "$out/ours_${name}_$m" --threads 2 > "$out/ours_${name}_$m.log" 2>&1 &
done
wait
