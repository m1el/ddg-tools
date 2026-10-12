#!/usr/bin/env bash
# Runs the reference driver (see build.sh) on
#   * the 23 meshes of the optimal-cones dataset with their published cones,
#   * the generated cases of `examples/reference_cases.rs`,
#   * five disks that ship with the reference implementation (automatic and unit-disk maps),
# then compares this crate with the driver's output (`examples/reference_diff.rs`).
#
#   reference/run.sh WORK_DIR > comparison.txt
#
# Run from the repository root, after `reference/build.sh WORK_DIR`. One thread throughout.
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "$HERE/../../.." && pwd)"
W="$(cd "${1:?usage: run.sh WORK_DIR}" && pwd)"
REF="${BFF_REF_SRC:-$ROOT/../ddg-course/boundary-first-flattening}"
DATA="$ROOT/data/optimal-cones/optimal-cones-data/data"
BIN="$W/bin/bff_ref"
export OPENBLAS_NUM_THREADS=1
mkdir -p "$W/dataset" "$W/cases" "$W/disks"
cd "$ROOT"

cargo build --release -p ddg-bff --example reference_cases --example reference_diff >&2
TARGET="${CARGO_TARGET_DIR:-$ROOT/target}"
"$TARGET/release/examples/reference_cases" "$W/cases" >&2

dataset=()
for obj in "$DATA"/*.obj; do
	name="$(basename "$obj" .obj)"
	"$BIN" cones "$obj" "$W/dataset/$name.out.txt" "$DATA/$name-cones.txt" > "$W/dataset/$name.log" 2>&1
	dataset+=("$W/dataset/$name.out.txt" "$DATA/$name-cones.txt" "$obj")
done
cases=()
for obj in "$W"/cases/*.obj; do
	name="$(basename "$obj" .obj)"
	"$BIN" cones "$obj" "$W/cases/$name.out.txt" "$W/cases/$name-cones.txt" > "$W/cases/$name.log" 2>&1
	cases+=("$W/cases/$name.out.txt" "$W/cases/$name-cones.txt" "-")
done
auto=()
disk=()
for name in cowhead face hemisphere maze spothead; do
	"$BIN" auto "$REF/input/$name.obj" "$W/disks/$name-auto.out.txt" > "$W/disks/$name-auto.log" 2>&1
	"$BIN" disk "$REF/input/$name.obj" "$W/disks/$name-disk.out.txt" > "$W/disks/$name-disk.log" 2>&1
	auto+=("$W/disks/$name-auto.out.txt")
	disk+=("$W/disks/$name-disk.out.txt")
done

DIFF="$TARGET/release/examples/reference_diff"
echo "## Dataset meshes, published cones (flatten_with_cones)"
"$DIFF" cones "${dataset[@]}"
echo
echo "## Generated cases (flatten_with_cones)"
"$DIFF" cones "${cases[@]}"
echo
echo "## Disks of the reference implementation, automatic flattening (u = 0 on the boundary)"
"$DIFF" auto "${auto[@]}"
echo
echo "## The same disks, flattened to the unit disk"
"$DIFF" disk "${disk[@]}"
