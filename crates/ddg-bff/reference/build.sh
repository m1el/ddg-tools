#!/usr/bin/env bash
# Builds the differential-test driver `bff_ref` around the reference C++ implementation of
# Boundary First Flattening (https://github.com/GeometryCollective/boundary-first-flattening,
# MIT, (c) 2017 Rohan Sawhney and Keenan Crane). The reference sources are compiled as they are;
# nothing of them is copied into this repository.
#
#   reference/build.sh WORK_DIR
#
# WORK_DIR is a scratch directory outside the repository; everything is written there:
#   WORK_DIR/ss-install   SuiteSparse 7.8.3 (CHOLMOD and what it needs), downloaded and built
#                         on first use
#   WORK_DIR/bin/bff_ref  the driver
#   WORK_DIR/logs         build logs
#
# Environment:
#   BFF_REF_SRC   checkout of the reference implementation
#                 (default: ../../../../ddg-course/boundary-first-flattening next to this repo)
#   CXX           C++ compiler (default g++)
#
# Needs cmake, ninja, curl and OpenBLAS (libopenblas-dev).
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
W="${1:?usage: build.sh WORK_DIR}"
mkdir -p "$W"
W="$(cd "$W" && pwd)"
REF="${BFF_REF_SRC:-$HERE/../../../../ddg-course/boundary-first-flattening}"
CXX="${CXX:-g++}"
SS="$W/ss-install"
mkdir -p "$W/dl" "$W/logs" "$W/ss-src" "$W/obj" "$W/bin"

if [ ! -f "$SS/lib/libcholmod.a" ]; then
	VER=7.8.3
	if [ ! -f "$W/dl/SuiteSparse-$VER.tar.gz" ]; then
		curl -L -sS -o "$W/dl/SuiteSparse-$VER.tar.gz" \
			"https://github.com/DrTimothyAldenDavis/SuiteSparse/archive/refs/tags/v$VER.tar.gz"
	fi
	[ -d "$W/ss-src/SuiteSparse-$VER" ] || tar -xzf "$W/dl/SuiteSparse-$VER.tar.gz" -C "$W/ss-src"
	cmake -S "$W/ss-src/SuiteSparse-$VER" -B "$W/ss-build" -G Ninja \
		-DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX="$SS" \
		-DSUITESPARSE_ENABLE_PROJECTS="suitesparse_config;amd;camd;colamd;ccolamd;cholmod" \
		-DSUITESPARSE_USE_FORTRAN=OFF -DSUITESPARSE_USE_CUDA=OFF -DSUITESPARSE_USE_OPENMP=OFF \
		-DBUILD_SHARED_LIBS=OFF -DBUILD_STATIC_LIBS=ON -DSUITESPARSE_DEMOS=OFF -DBUILD_TESTING=OFF \
		-DBLA_VENDOR=OpenBLAS -DCMAKE_POLICY_VERSION_MINIMUM=3.5 > "$W/logs/ss_configure.log" 2>&1
	cmake --build "$W/ss-build" > "$W/logs/ss_build.log" 2>&1
	cmake --install "$W/ss-build" > "$W/logs/ss_install.log" 2>&1
fi

# -ffp-contract=off: no fused multiply-add, whatever the target.
# -include cstdint: the reference headers use uint8_t without including <cstdint>, which
# current standard libraries no longer tolerate; force-including it avoids patching them.
CXXFLAGS="-std=c++11 -O2 -g -ffp-contract=off -include cstdint -I$REF/include -I$SS/include/suitesparse"
SRCS=(
	"$REF"/src/linear-algebra/*.cpp
	"$REF"/src/mesh/*.cpp
	"$REF"/src/project/{Bff,Cutter,Distortion,Generators,HoleFiller,ConePlacement}.cpp
	"$HERE/BinPackingStub.cpp"
	"$HERE/bff_ref.cpp"
)
pids=()
objs=()
for s in "${SRCS[@]}"; do
	name="$(basename "${s%.cpp}")"
	objs+=("$W/obj/$name.o")
	$CXX $CXXFLAGS -c "$s" -o "$W/obj/$name.o" > "$W/logs/cc_$name.log" 2>&1 &
	pids+=($!)
done
fail=0
for p in "${pids[@]}"; do wait "$p" || fail=1; done
if [ "$fail" != 0 ]; then
	echo "compile failed, see $W/logs/cc_*.log" >&2
	exit 1
fi
$CXX -g -o "$W/bin/bff_ref" "${objs[@]}" \
	"$SS/lib/libcholmod.a" "$SS/lib/libamd.a" "$SS/lib/libcamd.a" \
	"$SS/lib/libcolamd.a" "$SS/lib/libccolamd.a" "$SS/lib/libsuitesparseconfig.a" \
	-lopenblas -lm -lpthread > "$W/logs/link.log" 2>&1 || {
	echo "link failed, see $W/logs/link.log" >&2
	exit 1
}
echo "built $W/bin/bff_ref"
