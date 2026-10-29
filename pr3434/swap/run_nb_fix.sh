#!/bin/bash
#SBATCH --job-name=nb_fix
#SBATCH --output=nb_fix_%j.log
#SBATCH --time=06:00:00
#SBATCH --mem=96G
#SBATCH --cpus-per-task=16

# The negative binomial fixes after jachymb's review comment on branch
# stable-beta-binomial-lpmf (2026-10-12): the GLM partials without
# exp(theta), and lchoose(n + phi - 1, n) as -lbeta(phi, n + 1) - log(n + phi)
# in the GLM, neg_binomial_2_lpmf, neg_binomial_2_log_lpmf and
# neg_binomial_lpmf. Tree <tree>.
#
#   group 1  regression swap (STAN_MATH_DISTRIBUTIONS.md section 8) for
#            large_shapes_test: branch headers, then the 4 headers of the
#            branch head 468ff26dbb (the new rows must fail), then the
#            branch headers again
#   group 2  every unit test of the negative binomial family (prim, rev,
#            fwd, mix), large_shapes_test and the Laplace tests that use
#            neg_binomial_2_log
#   group 3  the generated test/prob tests of neg_binomial and
#            neg_binomial_2 (with neg_binomial_2_log), at -j$J3
#
# Deployed with this job into ~/nb_jobs: base_nb_headers.tar (the 468ff26dbb
# versions of the 4 changed CPU headers, paths relative to the tree).
# RUN_GROUPS selects groups, ':'-separated. Do not use `set -u` or `set -e`:
# /etc/profile reads unset variables, and nothing must abort.
set -o pipefail

source /etc/profile
module load triton/2025.1-gcc scicomp-r-env fftw/3.3.10

TREE="${TREE:?set TREE to the Stan Math tree}"
JOBDIR="$HOME/nb_jobs"
OUT="$JOBDIR/out_${SLURM_JOB_ID:-local}"
mkdir -p "$OUT"
cd "$TREE" || exit 1

J="${SLURM_CPUS_PER_TASK:-4}"
J3="${J3:-6}"
RUN_GROUPS="${RUN_GROUPS:-1:2:3}"
RUN_GROUPS="${RUN_GROUPS//:/ }"
echo "== host: $(hostname)   date: $(date -Is)   g++: $(g++ --version | head -1)"
echo "== groups: $RUN_GROUPS   J=$J  J3=$J3   out: $OUT"
[ -f make/local ] && echo "== WARNING: make/local exists"

BASE_TAR="$JOBDIR/base_nb_headers.tar"
[ -f "$BASE_TAR" ] || { echo "SETUP FAILED: no $BASE_TAR"; exit 1; }
mapfile -t HDRS < <(tar -tf "$BASE_TAR" | grep '\.hpp$')
echo "== ${#HDRS[@]} swapped headers"
for h in "${HDRS[@]}"; do
  [ -f "$h" ] || { echo "SETUP FAILED: no $h in the tree"; exit 1; }
done
echo "== md5 of the swapped headers (branch):"
md5sum "${HDRS[@]}"

SWAP_TESTS=(test/unit/math/rev/prob/large_shapes_test)

pass=0
fail=0
declare -a FAILED

run_one() {
  tag="$1"
  t="$2"
  base=$(basename "$t")
  rm -f "$t" "$t.o" "$t.d"
  if ! make -j"$J" "$t" > "$OUT/build_${tag}_${base}.log" 2>&1; then
    echo "  BUILD FAILED  $t"
    tail -25 "$OUT/build_${tag}_${base}.log"
    return 2
  fi
  if ./"$t" > "$OUT/run_${tag}_${base}.log" 2>&1; then
    n=$(grep -c '^\[       OK \]' "$OUT/run_${tag}_${base}.log")
    echo "  PASS          $t  ($n cases)"
    return 0
  fi
  nfail=$(grep -c '^\[  FAILED  \].*(' "$OUT/run_${tag}_${base}.log")
  npass=$(grep -c '^\[       OK \]' "$OUT/run_${tag}_${base}.log")
  nexp=$(grep -c 'Failure$' "$OUT/run_${tag}_${base}.log")
  echo "  FAIL          $t  ($nfail of $((nfail + npass)) cases, $nexp expectations)"
  grep '^\[  FAILED  \].*(' "$OUT/run_${tag}_${base}.log" | head -12
  return 1
}

count() {
  case "$1" in
    0) pass=$((pass + 1)) ;;
    *) fail=$((fail + 1)); FAILED+=("$2") ;;
  esac
}

if [[ " $RUN_GROUPS " == *" 1 "* ]]; then
  echo ""
  echo "######## group 1: regression swap ########"
  for t in "${SWAP_TESTS[@]}"; do
    [ -f "$t.cpp" ] || { echo "SETUP FAILED: missing $t.cpp"; exit 1; }
  done
  KEEP=$(mktemp -d)
  tar -cf "$KEEP/branch.tar" "${HDRS[@]}" || exit 1
  echo "-- stage branch_1"
  md5sum "${HDRS[@]}" | md5sum
  for t in "${SWAP_TESTS[@]}"; do run_one branch1 "$t"; count $? "swap branch_1 $t"; done
  echo "-- stage base (468ff26dbb)"
  tar -xf "$BASE_TAR" || exit 1
  md5sum "${HDRS[@]}" | md5sum
  for t in "${SWAP_TESTS[@]}"; do
    run_one base "$t"
    grep -E '^(NB|NB2|NB2L|GLM1|GLM5) ' "$OUT/run_base_$(basename "$t").log" \
      | sed -E 's/ (value|gradient [0-9]+)//' | sort | uniq -c
  done
  echo "-- stage branch_2"
  tar -xf "$KEEP/branch.tar" || exit 1
  md5sum "${HDRS[@]}" | md5sum
  for t in "${SWAP_TESTS[@]}"; do run_one branch2 "$t"; count $? "swap branch_2 $t"; done
  rm -rf "$KEEP"
fi

if [[ " $RUN_GROUPS " == *" 2 "* ]]; then
  echo ""
  echo "######## group 2: unit tests ########"
  LIST=$(mktemp)
  find test/unit/math -path '*/prob/*' -not -path '*/opencl/*' -name 'neg_binomial*_test.cpp' >> "$LIST"
  find test/unit/math/laplace -name 'laplace*neg_binomial*_test.cpp' >> "$LIST"
  echo test/unit/math/rev/prob/large_shapes_test.cpp >> "$LIST"
  # forward-mode GLM test of develop (#3418), after the merge of develop
  if [ -f test/unit/math/mix/prob/glm_forward_mode_tape_test.cpp ]; then
    echo test/unit/math/mix/prob/glm_forward_mode_tape_test.cpp >> "$LIST"
  fi
  sort -u "$LIST" -o "$LIST"
  echo "  $(wc -l < "$LIST") test sources"
  while IFS= read -r src; do
    run_one unit "${src%.cpp}"
    count $? "unit ${src%.cpp}"
  done < "$LIST"
  rm -f "$LIST"
fi

if [[ " $RUN_GROUPS " == *" 3 "* ]]; then
  echo ""
  echo "######## group 3: test/prob generated tests, -j$J3 ########"
  for fam in neg_binomial neg_binomial_2; do
    echo "  -- test/prob/$fam  $(date -Is)"
    if python3 runTests.py -j"$J3" "test/prob/$fam" > "$OUT/prob_${fam}.log" 2>&1; then
      echo "  PASS          test/prob/$fam ($(grep -c '^\[       OK \]' "$OUT/prob_${fam}.log") OK cases)"
      pass=$((pass + 1))
    else
      echo "  FAIL          test/prob/$fam"
      echo "     gtest failures: $(grep -c '^\[  FAILED  \]' "$OUT/prob_${fam}.log")"
      echo "     compiler kills: $(grep -c 'Killed signal terminated program cc1plus\|Cannot allocate memory' "$OUT/prob_${fam}.log")"
      grep -E '^\[  FAILED  \]|error:' "$OUT/prob_${fam}.log" | head -12
      fail=$((fail + 1))
      FAILED+=("prob $fam")
    fi
  done
fi

echo ""
echo "######## summary ########"
echo "passed: $pass   failed: $fail"
for f in "${FAILED[@]}"; do echo "  failed: $f"; done
echo "== final header set md5 (must equal stage branch_1):"
md5sum "${HDRS[@]}" | md5sum
echo "== done: $(date -Is)"
