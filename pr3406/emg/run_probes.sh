#!/bin/bash
# Claims M1, M2 and M4 of README_emg.md: build the four emg probes against
# each tree and write their outputs to out/emg/.  Run from pr3406/.
#
# Usage: TREES_DIR=<dir with one Stan Math tree per commit> \
#          emg/run_probes.sh <commit> [<commit> ...]
# Environment: as build.sh (TBB_LIB, OPT); JOBS = parallel builds (default 4).
set -eo pipefail
: "${TREES_DIR:?set TREES_DIR to the directory of the trees}"
JOBS="${JOBS:-4}"
mkdir -p out/emg

build_one() {
  local commit="$1" probe="$2"
  STAN_MATH="$TREES_DIR/$commit" ./build.sh "emg/$probe.cpp" \
    "out/emg/bin_${probe}_$commit" > "out/emg/build_${probe}_$commit.log" 2>&1 \
    || { echo "build failed: $probe $commit" >&2; return 1; }
}
export -f build_one
export TREES_DIR

for c in "$@"; do
  for p in probe_emg_phi probe_emg_dense check_ref_rows probe_deep_grid; do
    echo "$c $p"
  done
done | xargs -P "$JOBS" -n 2 bash -c 'build_one "$0" "$1"'

for c in "$@"; do
  for p in probe_emg_phi probe_emg_dense check_ref_rows probe_deep_grid; do
    timeout 300 "out/emg/bin_${p}_$c" > "out/emg/${p}_$c.txt"
  done
  echo "ran probes on $c"
done
