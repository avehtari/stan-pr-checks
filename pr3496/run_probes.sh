#!/bin/bash
# Run the probes, references and checks of README sections 1 to 5. Writes
# the files of results/ (except the Triton logs) to out/.
#
# Usage: CLONE=<stan-dev/math clone> STAN_MATH_LIB=<clone>/lib \
#          TBB_LIB=<dir with libtbb.so> ./run_probes.sh
# Run ./make_trees.sh first.
set -eo pipefail
cd "$(dirname "$0")"

: "${CLONE:?set CLONE to a stan-dev/math clone}"
export STAN_MATH_LIB TBB_LIB
OPT="${OPT:--O1}"
export OPT
mkdir -p out
{
  echo "$("${CXX:-g++}" --version | head -1)"
  python3 --version
  python3 -c 'import mpmath; print("mpmath", mpmath.__version__)'
  echo "OPT=$OPT"
} > out/versions.txt

# Section 1 and 2: the edge probe
for t in develop pr; do
  STAN_MATH=out/trees/$t ./build.sh edges/probe_edges.cpp out/probe_edges_$t
  ./out/probe_edges_$t > out/edges_$t.txt
done
python3 refs/ref_claims.py > out/ref_claims.txt

# Section 3: UBSan
: > out/ubsan_nb_lcdf.txt
for t in develop pr; do
  STAN_MATH=out/trees/$t \
    EXTRA_CXXFLAGS="-fsanitize=undefined -fno-sanitize-recover=undefined" \
    ./build.sh edges/ubsan_nb_lcdf.cpp out/ubsan_nb_lcdf_$t
  echo "== $t" >> out/ubsan_nb_lcdf.txt
  # The develop binary stops at the first error, with exit status 1
  ./out/ubsan_nb_lcdf_$t 2>&1 | sed "s|$PWD/||" >> out/ubsan_nb_lcdf.txt || true
done

# Section 4: large counts (needs out/trees/pr3433)
for t in develop pr pr3433; do
  STAN_MATH=out/trees/$t ./build.sh large_n/probe_nb2_lcdf_large_n.cpp \
    out/probe_nb2_large_n_$t
  ./out/probe_nb2_large_n_$t > out/nb2_large_n_$t.txt
  python3 large_n/score_large_n.py < out/nb2_large_n_$t.txt \
    > out/nb2_large_n_score_$t.txt
done

# Section 5: overlap with the open PRs
./overlap/overlap.sh > out/overlap.txt

echo "outputs in out/; compare with results/"
