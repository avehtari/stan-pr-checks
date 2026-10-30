#!/bin/bash
# Claim M3 of README_emg.md: device code compiled on the host against the CPU
# kernel, per tree.  Run from pr3406/.
# Usage: TREES_DIR=<dir> emg/ocl/run_cmp.sh <commit> [<commit> ...]
set -eo pipefail
: "${TREES_DIR:?set TREES_DIR to the directory of the trees}"
for c in "$@"; do
  d="out/emg/ocl_$c"
  mkdir -p "$d"
  STAN_MATH="$TREES_DIR/$c" ./build.sh emg/ocl/gen_dev_src.cpp "$d/gen_dev_src"
  "$d/gen_dev_src" > "$d/dev_src.inc"
  STAN_MATH="$TREES_DIR/$c" EXTRA_CXXFLAGS="-I $d -Wno-unknown-pragmas" \
    ./build.sh emg/ocl/cmp_host.cpp "$d/cmp_host"
  for g in dense fine; do
    timeout 300 "$d/cmp_host" "$g" > "out/emg/ocl_cmp_${g}_$c.txt"
  done
  echo "compared on $c"
done
