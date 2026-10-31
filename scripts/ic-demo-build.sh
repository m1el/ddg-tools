#!/usr/bin/env bash
# WP10, SPEC 13 Q2: builds the authors' demo (MIT, read-only reference tree) headless into
# runs_tmp/ic-demo-build/ (outside tracked files). Nothing is written into the reference tree.
# - geometry-central's Eigen: GC_EIGEN_LOCATION (Eigen 3.4.0, sha256
#   8586084f71f9bde545ee7fa6d00288b264a2b7ac3607b974e54d13e7162c1c72), default the scratch copy.
# - the old cmake_minimum_required of polyscope's glad needs CMAKE_POLICY_VERSION_MINIMUM=3.5.
# - geometry-central's own deps/happly is an empty submodule: the demo's deps/happly goes on the include path.
# - `setMarkedEdges` (declared, not defined at ab116f1) is not needed at link time by the demo.
set -eu
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
SRC=/var/data/discrete-geometry/IntegerCoordinates/intrinsic-triangulations-demo
EIGEN="${GC_EIGEN_LOCATION:?set GC_EIGEN_LOCATION to an Eigen 3.4.0 source tree}"
export CPLUS_INCLUDE_PATH="$SRC/deps/happly"
cmake -S "$SRC" -B "$ROOT/runs_tmp/ic-demo-build/build" -G Ninja -DGC_EIGEN_LOCATION="$EIGEN" \
  -DCMAKE_BUILD_TYPE=Release -DCMAKE_POLICY_VERSION_MINIMUM=3.5 > "$ROOT/runs_tmp/ic-demo-build/configure.log" 2>&1
ninja -C "$ROOT/runs_tmp/ic-demo-build/build" -j"${JOBS:-8}" int_tri > "$ROOT/runs_tmp/ic-demo-build/build.log" 2>&1
