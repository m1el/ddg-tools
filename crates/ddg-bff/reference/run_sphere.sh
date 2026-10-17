#!/usr/bin/env bash
# Runs the reference driver (see build.sh) in MODE sphere (`BFF::mapToSphere`) on the closed
# genus-0 meshes that are present (spot and blub from data/models, bunny, axes, chair and brain
# from the optimal-cones dataset, and two generated spheres), then compares this crate's
# `map_to_sphere` with the driver's output (`examples/reference_sphere_diff.rs`).
#
#   reference/run_sphere.sh WORK_DIR > comparison.txt
#
# Run from the repository root, after `reference/build.sh WORK_DIR`. The driver runs are
# independent and run in parallel; each is single-threaded.
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "$HERE/../../.." && pwd)"
W="$(cd "${1:?usage: run_sphere.sh WORK_DIR}" && pwd)"
BIN="$W/bin/bff_ref"
export OPENBLAS_NUM_THREADS=1
mkdir -p "$W/sphere"
cd "$ROOT"

cargo build --release -p ddg-bff --example reference_sphere_diff >&2
TARGET="${CARGO_TARGET_DIR:-$ROOT/target}"
DIFF="$TARGET/release/examples/reference_sphere_diff"
"$DIFF" write "$W/sphere" >&2

meshes=()
for f in data/models/spot/spot/spot_triangulated.obj data/models/blub/blub/blub_triangulated.obj \
	data/optimal-cones/optimal-cones-data/data/{bunny,axes,chair,brain}.obj \
	"$W"/sphere/generated-*.obj; do
	case "$f" in *-punctured.obj) continue ;; esac
	[ -f "$f" ] && meshes+=("$f")
done

pids=()
outs=()
for obj in "${meshes[@]}"; do
	name="$(basename "$obj" .obj)"
	"$BIN" sphere "$obj" "$W/sphere/$name.out.txt" > "$W/sphere/$name.log" 2>&1 &
	pids+=($!)
	outs+=("$W/sphere/$name.out.txt")
done
for p in "${pids[@]}"; do wait "$p" || echo "a driver run failed, see $W/sphere/*.log" >&2; done

echo "## map_to_sphere against BFF::mapToSphere"
"$DIFF" compare "${outs[@]}"

# The punctured meshes (written by the comparison) through the driver's own disk mode.
pids=()
douts=()
for obj in "$W"/sphere/*-punctured.obj; do
	name="$(basename "$obj" .obj)"
	"$BIN" disk "$obj" "$W/sphere/$name.out.txt" > "$W/sphere/$name.log" 2>&1 &
	pids+=($!)
	douts+=("$W/sphere/$name.out.txt")
done
for p in "${pids[@]}"; do wait "$p" || echo "a driver run failed, see $W/sphere/*.log" >&2; done
echo
echo "## The same meshes without the star of the pole, as disks (BFF::flattenToDisk)"
"$DIFF" disk "${douts[@]}"
