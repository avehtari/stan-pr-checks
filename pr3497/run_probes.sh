#!/bin/bash
# Run the local probes and reference generators of the README.  Writes the
# files of results/local/ and results/background/ to out/ (the Triton
# outputs in results/jobs/ and results/tests/ come from the job scripts).
#
# Usage: STAN_MATH_LIB=<clone>/lib TBB_LIB=<dir with libtbb.so> ./run_probes.sh
# Run first:
#   CLONE=<clone> ./make_trees.sh develop:dc66941969 branch:b40012df74 \
#       pr3433:709a1a4d76
set -eo pipefail
cd "$(dirname "$0")"

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

# The reference rows of the unit tests
python3 refs/gen_beta_reference.py > out/beta_reference.hpp \
  2> out/beta_reference_routes.txt
python3 refs/gen_caller_reference.py > out/caller_reference.txt

# Edge inputs and the beta_lcdf x derivative
for t in develop branch; do
  STAN_MATH=out/trees/$t ./build.sh probe_edges.cpp out/probe_edges_$t
  ./out/probe_edges_$t > out/edges_$t.txt
  STAN_MATH=out/trees/$t ./build.sh probe_caller.cpp out/probe_caller_$t
  ./out/probe_caller_$t > out/caller_$t.txt
  python3 refs/score_caller.py < out/caller_$t.txt > out/score_caller_$t.txt
done
python3 refs/compare_edges.py out/edges_develop.txt out/edges_branch.txt \
  > out/edges_compare.txt

# The table of the PR summary
for t in develop branch; do
  STAN_MATH=out/trees/$t ./build.sh probe_table.cpp out/probe_table_$t
  ./out/probe_table_$t > out/table_$t.txt
done
python3 refs/score_table.py out/table_develop.txt out/table_branch.txt \
  > out/score_table.txt

# Background: the probes of the issue draft (develop and #3433)
for t in develop pr3433; do
  STAN_MATH=out/trees/$t ./build.sh background/probe_beta_lcdf.cpp out/probe_$t
  ./out/probe_$t > out/beta_lcdf_$t.txt
  python3 background/score_beta_lcdf.py < out/beta_lcdf_$t.txt \
    > out/score_$t.txt
  STAN_MATH=out/trees/$t ./build.sh background/probe_beta_lccdf_dx.cpp \
    out/probe_lccdf_dx_$t
  ./out/probe_lccdf_dx_$t > out/beta_lccdf_dx_$t.txt
  python3 background/score_lccdf_dx.py < out/beta_lccdf_dx_$t.txt \
    > out/score_lccdf_dx_$t.txt
done
echo "done; compare out/ with results/local/ and results/background/"
