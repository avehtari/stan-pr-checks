#!/bin/bash
# Claim B1 of README_emg.md: build bitwise/compare.cpp against each commit and
# its parent (-O3 unless OPT is set) and diff the hex outputs.  Run from
# pr3406/.  Usage: TREES_DIR=<dir> bitwise/run_compare.sh
set -eo pipefail
: "${TREES_DIR:?set TREES_DIR to the directory of the trees}"
PAIRS="0d4ac36e34:f4cd10ee4a f4cd10ee4a:c30ee6f9ef 9b38ecbcb0:927099efb8"
mkdir -p out/bitwise
trees=$(for p in $PAIRS; do echo "${p%%:*}"; echo "${p#*:}"; done | sort -u)
for c in $trees; do
  STAN_MATH="$TREES_DIR/$c" ./build.sh bitwise/compare.cpp \
    "out/bitwise/compare_$c" > "out/bitwise/build_$c.log" 2>&1 &
done
wait
for c in $trees; do
  timeout 120 "out/bitwise/compare_$c" > "out/bitwise/out_$c.txt"
done
for p in $PAIRS; do
  a="out/bitwise/out_${p%%:*}.txt"
  b="out/bitwise/out_${p#*:}.txt"
  la=$(wc -l < "$a")
  lb=$(wc -l < "$b")
  nd=$(paste -d '\n' "$a" "$b" | paste - - | awk -F'\t' '$1 != $2' | wc -l)
  echo "${p%%:*} vs ${p#*:}: lines $la / $lb, differing lines $nd"
  diff "$a" "$b" | head -6 || true
done
