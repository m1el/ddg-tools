#!/usr/bin/env bash
# Rebuild the authors' demo of "Navigating Intrinsic Triangulations" (Sharp, Soliman, Crane 2019; MIT) out of tree
# and regenerate the reference data for differential tests of crates/ddg-intrinsic.
#
#   scripts/intrinsic-reference.sh            build (if needed) + regenerate everything
#   NIT_SKIP_BUILD=1 scripts/intrinsic-reference.sh   only regenerate the data
#
# Output: data/intrinsic-reference/<input>/{input.obj, <op>/...}; see data/intrinsic-reference/FORMAT.md
# (FORMAT.md is hand-written and is not touched here). Logs: runs_tmp/nit-reference/*.log.
#
# Choices made, and why:
#  * The demo directory is read-only for us. Its only change is `git submodule update --init --recursive`
#    (geometry-central v1 and polyscope), needed by its build. All build products go to
#    runs_tmp/nit-reference/build.
#  * No root: CMake 3.31 comes from the micromamba prefix /var/data/ddg-quarantine/prefix/cm (CMake 4 of the system
#    rejects old bundled dependencies); -DCMAKE_POLICY_VERSION_MINIMUM=3.5 is passed as well. SuiteSparse is not
#    needed (geometry-central falls back to Eigen), nothing else had to be installed.
#  * The demo is run with --noGUI: it never opens a window, so neither xvfb nor a display is needed.
#  * The demo ignores --outputPrefix (parsed, never used): it writes faceInds.dmat etc. into the working directory
#    under fixed names. Hence every run is done in its own directory.
#  * The demo prints no counts. scripts/intrinsic-reference/stats.cpp repeats the same calls on the same library and
#    counts flips and insertions with geometry-central's callbacks; analyze.py derives angle statistics and the
#    self-edge / double-edge counts from the demo's own output files.
#  * Operations: --flipDelaunay (flip), --refineDelaunay with --refineAngle 25 (the demo's default; directory
#    refine25) and 30 (the paper's bound; refine30). Other refinement options (circumradius, insertion cap) are left
#    at the demo's defaults (no size bound; cap 10 * nVertices).
#  * All outputs are requested: --intrinsicFaces --vertexPositions --laplaceMat --interpolateMat.
#  * Inputs: three Crane-repository models, seven optimal-cones meshes (bump has boundary), five generated meshes
#    (scripts/intrinsic-reference/gen_meshes.py).
#  * Wall time is that of the whole demo process (reading the mesh, building the triangulation, the operation,
#    writing the files); op_seconds in stats.txt is the operation alone.
set -euo pipefail

ROOT=$(cd "$(dirname "$0")/.." && pwd)
DEMO=${NIT_DEMO:-/var/data/discrete-geometry/NavigatingIntrinsicTriangulations/navigating-intrinsic-triangulations-demo}
CMAKE_PREFIX=${NIT_CMAKE_PREFIX:-/var/data/ddg-quarantine/prefix/cm}
WORK=$ROOT/runs_tmp/nit-reference
BUILD=$WORK/build
OUT=$ROOT/data/intrinsic-reference
HELP=$ROOT/scripts/intrinsic-reference
OC=$ROOT/data/optimal-cones/optimal-cones-data/data
MODELS=$ROOT/data/models

mkdir -p "$WORK" "$OUT"
export PATH="$CMAKE_PREFIX/bin:$PATH"

if [ -z "${NIT_SKIP_BUILD:-}" ]; then
  echo "== submodules"
  (cd "$DEMO" && git submodule update --init --recursive) > "$WORK/submodule.log" 2>&1
  echo "== configure"
  command -v ninja >/dev/null || { echo "ninja is required"; exit 1; }
  GEN=(-G Ninja); BUILDCMD=(cmake --build "$BUILD" --target int_tri -j 32)
  cmake -S "$DEMO" -B "$BUILD" "${GEN[@]}" -DCMAKE_BUILD_TYPE=Release -DCMAKE_POLICY_VERSION_MINIMUM=3.5 \
    > "$WORK/cmake.log" 2>&1
  echo "== build (log: $WORK/build.log)"
  "${BUILDCMD[@]}" > "$WORK/build.log" 2>&1

  echo "== stats helper"
  # Reuse the compile flags and include paths the build used for the demo's main.cpp (needs ninja).
  CC=$(cd "$BUILD" && ninja -t commands int_tri | grep -- 'src/main.cpp' | grep -- ' -c ' | head -1)
  CC=$(echo "${CC%% -o *}" | sed -E 's/ -MD -MT [^ ]+ -MF [^ ]+//')   # keep the flags, drop "-o ..." and depfile options
  # shellcheck disable=SC2086
  (cd "$BUILD" && $CC -o "$WORK/stats.o" -c "$HELP/stats.cpp" && \
   c++ "$WORK/stats.o" -o "$BUILD/bin/nit_stats" lib/libgeometry-central.a -lm -ldl) > "$WORK/stats-build.log" 2>&1
fi

INT_TRI=$BUILD/bin/int_tri
STATS=$BUILD/bin/nit_stats
[ -x "$INT_TRI" ] && [ -x "$STATS" ] || { echo "demo not built (see $WORK/*.log)"; exit 1; }

echo "== generated meshes"
GEN_DIR=$WORK/generated
uv run "$HELP/gen_meshes.py" "$GEN_DIR" > "$WORK/gen.log" 2>&1

# name -> input path
declare -A INPUT=(
  [spot]=$MODELS/spot/spot/spot_triangulated.obj
  [bob]=$MODELS/bob/bob/bob_tri.obj
  [blub]=$MODELS/blub/blub/blub_triangulated.obj
  [bunny]=$OC/bunny.obj
  [camel]=$OC/camel.obj
  [camelbad]=$OC/camelbad.obj
  [camelbadnonuniform]=$OC/camelbadnonuniform.obj
  [camel_noise1]=$OC/camel_noise1.obj
  [camelshotnoise]=$OC/camelshotnoise.obj
  [bump]=$OC/bump.obj
  [icosphere]=$GEN_DIR/icosphere.obj
  [icosphere_jitter]=$GEN_DIR/icosphere_jitter.obj
  [thin_cylinder]=$GEN_DIR/thin_cylinder.obj
  [long_torus]=$GEN_DIR/long_torus.obj
  [sheared_grid]=$GEN_DIR/sheared_grid.obj
)
NAMES=(spot bob blub bunny camel camelbad camelbadnonuniform camel_noise1 camelshotnoise bump
       icosphere icosphere_jitter thin_cylinder long_torus sheared_grid)
# op directory -> demo arguments (also given to stats)
OPS=(flip refine25 refine30)

for name in "${NAMES[@]}"; do
  src=${INPUT[$name]}
  [ -f "$src" ] || { echo "missing input $src"; exit 1; }
  d=$OUT/$name
  mkdir -p "$d"
  cp "$src" "$d/input.obj"
  for op in "${OPS[@]}"; do
    case $op in
      flip)     DEMO_ARGS=(--flipDelaunay);                   STAT_ARGS=(flip) ;;
      refine25) DEMO_ARGS=(--refineDelaunay --refineAngle 25); STAT_ARGS=(refine 25) ;;
      refine30) DEMO_ARGS=(--refineDelaunay --refineAngle 30); STAT_ARGS=(refine 30) ;;
    esac
    od=$d/$op
    mkdir -p "$od"
    # The demo's own files are overwritten on every run; no other file in the directory is left stale.
    CMD=("$INT_TRI" ../input.obj --noGUI "${DEMO_ARGS[@]}" --intrinsicFaces --vertexPositions --laplaceMat
         --interpolateMat)
    echo "${CMD[*]}" | sed "s#$INT_TRI#int_tri#" > "$od/command.txt"
    echo "(cd $od && int_tri as above)  # run from this directory; input is ../input.obj" >> "$od/command.txt"
    t0=$(date +%s.%N)
    (cd "$od" && timeout 1800 "${CMD[@]}") > "$od/stdout.txt" 2>&1 || { echo "FAILED $name $op"; continue; }
    t1=$(date +%s.%N)
    printf '%.3f seconds wall\n' "$(echo "$t1 - $t0" | bc -l)" > "$od/time.txt"
    "$STATS" "$od/../input.obj" "${STAT_ARGS[@]}" > "$od/stats.txt" 2> "$od/stats.err"
    uv run "$HELP/analyze.py" "$od" > "$od/analysis.log" 2>&1 || echo "analysis failed $name $op"
    echo "done $name $op: $(cat "$od/time.txt")"
  done
done

uv run "$HELP/analyze.py" --table "$OUT" > "$OUT/SUMMARY.md"
echo "== summary in $OUT/SUMMARY.md"
