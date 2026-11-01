#!/bin/bash
# Swap test: run the new and the existing initialize tests with the PR
# initialize.hpp, with the base (develop) version, and with the PR version
# again.  Run in a Stan clone at the PR branch with its submodules.
#   STAN=<stan clone> BASE=<base commit> ./swap_test.sh
set -e -o pipefail
: "${STAN:?set STAN to a Stan clone at the PR branch}"
BASE=${BASE:-a8ec3d18f}
H=src/stan/services/util/initialize.hpp
T="src/test/unit/services/util/initialize_partial_test.cpp src/test/unit/services/util/initialize_test.cpp"
cd "$STAN"
git diff --quiet -- "$H" || { echo "uncommitted changes in $H"; exit 1; }
run() {
  echo "== $1: $(md5sum "$H" | cut -c1-8)"
  rm -f test/unit/services/util/initialize_partial_test test/unit/services/util/initialize_test
  ./runTests.py $T 2>&1 | grep -E "^\[  (PASSED|FAILED) |^\[==========\]" || true
}
trap 'git checkout -- "$H"' EXIT
run "PR"
git show "$BASE:$H" > "$H"
run "develop $BASE"
git checkout -- "$H"
run "PR again"
