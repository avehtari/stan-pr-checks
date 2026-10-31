#!/bin/bash
# Swap test: run the new test with the develop versions of all changed
# prim headers at once, then with the PR headers again. Exit 0 only if the
# test passes with the PR headers, fails with the develop headers, and
# passes again after the restore.
#
# Usage: TREE=<clone at the PR head> ./swap/swap.sh
#
# Environment:
#   TREE      required; a stan-dev/math clone with the PR head checked out,
#             a built TBB, and no uncommitted changes in stan/ or test/
#   DEVELOP   optional; the merge base (default dc66941969)
#   PR_TREE   optional; the tree hash of the PR head (default 1741248cd7)
#   SUITE     optional; the test (default the new mix test)
#   J         optional; make -j (default 4)
# Logs go to out/ of this directory. On exit, also after an error or an
# interrupt, the script restores the PR headers with git checkout.
set -o pipefail

HERE="$(cd "$(dirname "$0")/.." && pwd)"
LOGS="$HERE/out"
mkdir -p "$LOGS"

: "${TREE:?set TREE to a clone at the PR head}"
DEVELOP="${DEVELOP:-dc66941969}"
PR_TREE="${PR_TREE:-1741248cd7}"
SUITE="${SUITE:-test/unit/math/mix/prob/discrete_cdfs_edge_elements_test}"
J="${J:-4}"
cd "$TREE" || exit 1

if [ -n "$(git status --short -- stan test)" ]; then
  echo "swap.sh: uncommitted changes in stan/ or test/; stop" >&2
  exit 1
fi
head_tree=$(git rev-parse --short=10 "HEAD^{tree}")
if [ "$head_tree" != "$PR_TREE" ]; then
  echo "swap.sh: HEAD has tree $head_tree, not the PR tree $PR_TREE; stop" >&2
  exit 1
fi
changed=$(git diff --name-only --diff-filter=M "$DEVELOP" HEAD -- stan/math/prim) \
  || exit 1
mapfile -t headers <<< "$changed"
if [ -z "$changed" ]; then
  echo "swap.sh: no changed prim headers found; stop" >&2
  exit 1
fi

trap 'git -C "$TREE" checkout HEAD -- stan' EXIT
trap 'exit 130' INT TERM

# Prints PASS or FAIL as the last word; BUILD_FAILED if the build fails
run_stage() {
  stage="$1"
  rm -f "$SUITE" "$SUITE.o" "$SUITE.d"
  if ! make -j"$J" "$SUITE" > "$LOGS/swap_build_$stage.log" 2>&1; then
    echo "$stage: BUILD_FAILED ($LOGS/swap_build_$stage.log)"
    return
  fi
  ./"$SUITE" > "$LOGS/swap_run_$stage.log" 2>&1
  status=$?
  ok=$(grep -c '^\[       OK \]' "$LOGS/swap_run_$stage.log")
  failed=$(grep -c '^\[  FAILED  \].*(' "$LOGS/swap_run_$stage.log")
  expectations=$(grep -c ': Failure' "$LOGS/swap_run_$stage.log")
  grep '^\[  FAILED  \].*(' "$LOGS/swap_run_$stage.log" | sed 's/^/  /' >&2
  echo "$stage: $ok passed, $failed failed, $expectations failed expectations" \
       "$([ $status -eq 0 ] && echo PASS || echo FAIL)"
}

r1=$(run_stage pr); echo "$r1"
for f in "${headers[@]}"; do
  git show "$DEVELOP:$f" > "$f" || exit 1
done
echo "develop versions of ${#headers[@]} headers in place"
r2=$(run_stage develop); echo "$r2"
git checkout HEAD -- stan || { echo "swap.sh: restore failed" >&2; exit 1; }
[ -z "$(git status --short -- stan)" ] && echo "PR headers restored"
r3=$(run_stage pr_again); echo "$r3"

if [[ "$r1" == *PASS && "$r2" == *FAIL && "$r3" == *PASS ]]; then
  echo "swap test: OK"
  exit 0
fi
echo "swap test: NOT OK"
exit 1
