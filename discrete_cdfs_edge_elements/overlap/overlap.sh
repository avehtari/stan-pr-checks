#!/bin/bash
# The merge-order note of the PR text: which headers of this PR the open
# PRs also change. Several of these PRs are stacked, so each PR is compared
# with its own base branch, not with develop.
#
# Usage: CLONE=<stan-dev/math clone> ./overlap/overlap.sh
#
# Environment:
#   CLONE     required; a clone with the remote origin = stan-dev/math
#   DEVELOP   optional; the merge base of this PR (default dc66941969)
#   PR        optional; the head of this PR (default 1315143b74)
# The script fetches the heads of the open PRs into refs/remotes/pr/<N>
# and their base branches into refs/remotes/origin/<branch>.
set -eo pipefail

: "${CLONE:?set CLONE to a stan-dev/math clone}"
DEVELOP="${DEVELOP:-dc66941969}"
PR="${PR:-1315143b74}"

# PR number and base branch, as on GitHub on 2026-10-31
OPEN_PRS="3408:develop 3409:develop 3489:stable-inc-beta
3490:stable-inc-beta-rev 3433:stable-inc-beta-cdfs 3434:develop
3437:stable-beta-binomial-lpmf 3441:stable-beta-binomial-lpmf
3445:fix-hypergeometric-3F2-termination"

mapfile -t headers < <(git -C "$CLONE" diff --name-only "$DEVELOP" "$PR" -- stan)
echo "headers changed by this PR: ${#headers[@]}"

declare -A any
for pair in $OPEN_PRS; do
  n="${pair%%:*}"
  base_branch="${pair#*:}"
  git -C "$CLONE" fetch -q origin "+refs/pull/$n/head:refs/remotes/pr/$n" \
    "+refs/heads/$base_branch:refs/remotes/origin/$base_branch"
  base=$(git -C "$CLONE" merge-base "origin/$base_branch" "pr/$n")
  common=()
  for f in "${headers[@]}"; do
    if ! git -C "$CLONE" diff --quiet "$base" "pr/$n" -- "$f"; then
      common+=("$(basename "$f" .hpp)")
      any["$f"]=1
    fi
  done
  echo "#$n (base $base_branch, head $(git -C "$CLONE" rev-parse --short "pr/$n")):" \
       "${#common[@]} common: ${common[*]}"
done
echo "headers that at least one open PR also changes: ${#any[@]}"
