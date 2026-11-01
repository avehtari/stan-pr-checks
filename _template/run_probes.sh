#!/bin/bash
# Run every probe and scorer of this directory against every tree and
# write the outputs to out/.  A new run must give files identical to
# results/ (compare with cmp); copy them to results/ only after a
# deliberate rerun.  The cluster job logs in results/ come from
# tools/triton_fetch_logs.sh --scrub, not from this script.
#
# Usage: STAN_MATH_LIB=<clone>/lib TBB_LIB=<dir with libtbb.so> ./run_probes.sh
# Run ./make_trees.sh first.  TREES lists the trees (default: develop pr).
set -eo pipefail
cd "$(dirname "$0")"

export STAN_MATH_LIB TBB_LIB
OPT="${OPT:--O1}"
export OPT
TREES="${TREES:-develop pr}"
mkdir -p out
{
  "${CXX:-g++}" --version | head -1
  python3 --version
  python3 -c 'import mpmath; print("mpmath", mpmath.__version__)'
  echo "OPT=$OPT"
} > out/versions.txt

for t in $TREES; do
  STAN_MATH=out/trees/$t ./build.sh probes/probe.cpp out/probe_$t
  ./out/probe_$t > out/probe_$t.txt
  python3 refs/score.py < out/probe_$t.txt > out/score_$t.txt
done

echo "outputs in out/; compare with results/:"
for f in out/*.txt; do
  b=$(basename "$f")
  if [ -f "results/$b" ]; then
    cmp -s "$f" "results/$b" && echo "  same $b" || echo "  DIFF $b"
  else
    echo "  new  $b"
  fi
done
