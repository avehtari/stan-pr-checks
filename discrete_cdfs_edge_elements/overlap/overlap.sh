#!/bin/bash
# The merge-order note of the PR text: which headers of this PR the open
# PRs #3408, #3433, #3434, #3437, #3441 and #3445 also change.
#
# Usage: CLONE=<stan-dev/math clone> ./overlap/overlap.sh
#
# Environment:
#   CLONE     required; a clone with the remote origin = stan-dev/math
#   DEVELOP   optional; the merge base of this PR (default dc66941969)
#   PR        optional; the head of this PR (default 1315143b74)
# The script fetches the heads of the open PRs into refs/remotes/pr/<N>.
set -eo pipefail

: "${CLONE:?set CLONE to a stan-dev/math clone}"
DEVELOP="${DEVELOP:-dc66941969}"
PR="${PR:-1315143b74}"
OPEN_PRS="3408 3433 3434 3437 3441 3445"

for n in $OPEN_PRS; do
  git -C "$CLONE" fetch -q origin "+refs/pull/$n/head:refs/remotes/pr/$n"
done

mapfile -t headers < <(git -C "$CLONE" diff --name-only "$DEVELOP" "$PR" -- stan)
echo "headers changed by this PR: ${#headers[@]}"

declare -A any
for n in $OPEN_PRS; do
  ref="pr/$n"
  base=$(git -C "$CLONE" merge-base "$DEVELOP" "$ref")
  common=()
  for f in "${headers[@]}"; do
    if ! git -C "$CLONE" diff --quiet "$base" "$ref" -- "$f"; then
      common+=("$(basename "$f" .hpp)")
      any["$f"]=1
    fi
  done
  echo "#$n (head $(git -C "$CLONE" rev-parse --short "$ref")):" \
       "${#common[@]} common: ${common[*]}"
done
echo "headers that at least one open PR also changes: ${#any[@]}"

base_3434=$(git -C "$CLONE" merge-base "$DEVELOP" pr/3434)
total_3434=$(git -C "$CLONE" rev-list --count "$base_3434..pr/3434")
for n in 3437 3441 3445; do
  missing=$(git -C "$CLONE" rev-list --count "pr/3434" "^$base_3434" "^pr/$n")
  echo "#$n contains $((total_3434 - missing)) of the $total_3434 commits of #3434"
done
