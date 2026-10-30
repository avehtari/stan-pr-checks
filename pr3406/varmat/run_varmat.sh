#!/bin/bash
# Build and run each case of varmat/probe_varmat.cpp with AddressSanitizer
# on two trees; prints the compile result and the run output per case.
# Usage (from pr3406/): DEV=<tree> HEAD=<tree> TBB_LIB=<dir> varmat/run_varmat.sh
set -u
: "${DEV:?set DEV}" "${HEAD:?set HEAD}"
OUT=out/varmat
mkdir -p "$OUT"
ASAN="-fsanitize=address -fno-omit-frame-pointer -g"
for tree in dev head; do
  if [ "$tree" = dev ]; then SM=$DEV; else SM=$HEAD; fi
  for c in AS_TWICE AS_ARRAY REF_TYPE DANGLE; do
    bin=$OUT/probe_${c}_$tree
    log=$OUT/build_${c}_$tree.log
    echo "== $tree $c"
    if STAN_MATH=$SM OPT=-O1 EXTRA_CXXFLAGS="$ASAN -DCASE_$c" \
        timeout 600 ./build.sh varmat/probe_varmat.cpp "$bin" > "$log" 2>&1; then
      echo "compile: ok"
      ASAN_OPTIONS=detect_leaks=0 timeout 60 "$bin" > "$OUT/run_${c}_$tree.txt" 2>&1
      echo "exit: $?"
      grep -E "^(ERROR|SUMMARY): AddressSanitizer|expect|reference:" \
        "$OUT/run_${c}_$tree.txt" | head -8
    else
      echo "compile: FAILED"
      grep -m 3 -E "error:" "$log" | cut -c1-300
    fi
  done
done
