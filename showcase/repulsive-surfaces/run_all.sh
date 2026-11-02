#!/usr/bin/env bash
# Runs every application of Sec. 9 of "Repulsive Surfaces" (Yu, Brakensiek, Schumacher, Crane 2021)
# in parallel (one job and one log per application).
#
#   showcase/repulsive-surfaces/run_all.sh [--smoke] [--threads T] [--jobs J] [--out DIR] [--frame-every K] [NAME...]
#
# Needs the reference scenes: scripts/rs-fetch-scenes.sh. Logs: runs_tmp/wp15/<mode>/<name>.log
# (RESULT lines: steps, time, energy, constraints, self-intersection). Frames: DIR/<name>/*.obj
# (default showcase/media/repulsive-surfaces/, untracked).
set -euo pipefail
here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
root="$here/../.."
mode=full; threads=4; jobs=14; out="$root/showcase/media/repulsive-surfaces"; every=10; names=()
while [ $# -gt 0 ]; do
  case "$1" in
    --smoke) mode=smoke ;;
    --threads) threads="$2"; shift ;;
    --jobs) jobs="$2"; shift ;;
    --out) out="$2"; shift ;;
    --frame-every) every="$2"; shift ;;
    *) names+=("$1") ;;
  esac
  shift
done
bin="${CARGO_TARGET_DIR:-$root/target}/release/examples/apps"
[ -x "$bin" ] || (cd "$root" && cargo build --release -p ddg-repulsive-surfaces --example apps)
log="$root/runs_tmp/wp15/$mode"
mkdir -p "$log"
[ ${#names[@]} -gt 0 ] || mapfile -t names < <("$bin" list)
extra=()
[ "$mode" = smoke ] && extra+=(--smoke) || extra+=(--out "$out" --frame-every "$every")
printf '%s\n' "${names[@]}" | xargs -P "$jobs" -I{} bash -c "'$bin' {} --threads $threads ${extra[*]} > '$log/{}.log' 2>&1 || true"
grep -h '^RESULT' "$log"/*.log
