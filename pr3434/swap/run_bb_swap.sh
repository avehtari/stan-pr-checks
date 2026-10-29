#!/bin/bash
#SBATCH --job-name=bb_swap
#SBATCH --output=bb_swap_%j.log
#SBATCH --time=04:00:00
#SBATCH --mem=120G
#SBATCH --cpus-per-task=16

# beta_binomial_lpmf repair (stan-dev/math#3153), branch
# stable-beta-binomial-lpmf, tree <tree>.
#
#   stage patched_1  the new and the existing unit tests on the branch
#   stage develop    the develop beta_binomial_lpmf.hpp: the new tests must
#                    fail (STAN_MATH_DISTRIBUTIONS.md section 8)
#   stage patched_2  the branch header again: everything must pass
#   extra            digamma_diff prim and mix tests, the other existing
#                    beta_binomial unit tests, and the generated
#                    test/prob/beta_binomial tests (RUN_PROB=1)
#
# The develop header is deployed with this job as
# $HOME/bb_jobs/beta_binomial_lpmf.hpp.develop. Do not use `set -u` or
# `set -e`: /etc/profile reads unset variables, and nothing must abort.
set -o pipefail

source /etc/profile
module load triton/2025.1-gcc scicomp-r-env fftw/3.3.10

TREE="${TREE:?set TREE to the Stan Math tree}"
JOBDIR="$HOME/bb_jobs"
cd "$TREE" || exit 1

echo "== host: $(hostname)   date: $(date -Is)   g++: $(g++ --version | head -1)"
[ -f make/local ] && echo "== WARNING: make/local exists"
J="${SLURM_CPUS_PER_TASK:-4}"
J3="${J3:-6}"
RUN_PROB="${RUN_PROB:-1}"

HDR=stan/math/prim/prob/beta_binomial_lpmf.hpp
DEV_HDR="$JOBDIR/beta_binomial_lpmf.hpp.develop"
SWAP_TESTS=(
  test/unit/math/prim/prob/beta_binomial_test
  test/unit/math/rev/prob/beta_binomial_lpmf_test
)
EXTRA_TESTS=(
  test/unit/math/prim/fun/digamma_diff_test
  test/unit/math/mix/fun/digamma_diff_test
  test/unit/math/prim/prob/beta_binomial_cdf_log_test
  test/unit/math/prim/prob/beta_binomial_ccdf_log_test
  test/unit/math/prim/fun/lbeta_test
  test/unit/math/prim/fun/lgamma_stirling_diff_test
)

[ -f "$DEV_HDR" ] || { echo "SETUP FAILED: no $DEV_HDR"; exit 1; }
for t in "${SWAP_TESTS[@]}" "${EXTRA_TESTS[@]}"; do
  [ -f "$t.cpp" ] || { echo "SETUP FAILED: missing $t.cpp"; exit 1; }
done
mkdir -p swap_keep
cp -f "$HDR" swap_keep/beta_binomial_lpmf.hpp.patched || exit 1
echo "== header md5: branch $(md5sum < "$HDR" | cut -c1-32)" \
     " develop $(md5sum < "$DEV_HDR" | cut -c1-32)"

run_tests() {
  tag="$1"
  shift
  suites_failed=0
  for t in "$@"; do
    base=$(basename "$t")
    rm -f "$t" "$t.o" "$t.d"
    if ! make -j"$J" "$t" > "bb_build_${tag}_${base}.log" 2>&1; then
      echo "  BUILD FAILED  $base"
      tail -25 "bb_build_${tag}_${base}.log"
      suites_failed=$((suites_failed + 1))
      continue
    fi
    if ./"$t" > "bb_run_${tag}_${base}.log" 2>&1; then
      n=$(grep -c '^\[       OK \]' "bb_run_${tag}_${base}.log")
      echo "  PASS          $base  ($n cases)"
    else
      nfail=$(grep -c '^\[  FAILED  \].*(' "bb_run_${tag}_${base}.log")
      npass=$(grep -c '^\[       OK \]' "bb_run_${tag}_${base}.log")
      echo "  FAIL          $base  ($nfail of $((nfail + npass)) cases failed)"
      grep '^\[  FAILED  \]' "bb_run_${tag}_${base}.log" | head -12
      grep -B1 -A3 'Failure' "bb_run_${tag}_${base}.log" | head -40
      suites_failed=$((suites_failed + 1))
    fi
  done
  echo "== $tag: $suites_failed of $# suites failed"
}

echo ""
echo "######## stage patched_1 ########"
md5sum "$HDR"
run_tests patched_1 "${SWAP_TESTS[@]}"

echo ""
echo "######## stage develop ########"
cp -f "$DEV_HDR" "$HDR" || exit 1
md5sum "$HDR"
run_tests develop "${SWAP_TESTS[@]}"

echo ""
echo "######## stage patched_2 ########"
cp -f swap_keep/beta_binomial_lpmf.hpp.patched "$HDR" || exit 1
md5sum "$HDR"
run_tests patched_2 "${SWAP_TESTS[@]}"

echo ""
echo "######## extra unit tests ########"
run_tests extra "${EXTRA_TESTS[@]}"

if [ "$RUN_PROB" = "1" ]; then
  echo ""
  echo "######## test/prob/beta_binomial (generated), -j$J3 ########"
  echo "  start $(date -Is)"
  if python3 runTests.py -j"$J3" test/prob/beta_binomial \
       > bb_prob_beta_binomial.log 2>&1; then
    echo "  PASS          test/prob/beta_binomial"
  else
    echo "  FAIL          test/prob/beta_binomial"
    echo "     gtest failures: $(grep -c '^\[  FAILED  \]' bb_prob_beta_binomial.log)"
    echo "     compiler kills: $(grep -c 'Killed signal terminated program cc1plus\|Cannot allocate memory' bb_prob_beta_binomial.log)"
    grep -E '^\[  FAILED  \]|error:' bb_prob_beta_binomial.log | head -20
  fi
  echo "  ran $(grep -c '^\[       OK \]' bb_prob_beta_binomial.log) OK cases"
fi

echo ""
echo "== final header md5 (must equal stage patched_1):"
md5sum "$HDR"
echo "== done: $(date -Is)"
