#!/bin/bash
# Export stan/ and test/unit/ of the develop commit, the PR head and the
# head of stan-dev/math#3433 into out/trees/<name>/, for the probes and the
# swap test of this directory.
#
# Usage: CLONE=<stan-dev/math clone> ./make_trees.sh
#
# Environment:
#   CLONE     required; a stan-dev/math clone (remote origin = stan-dev/math)
#   DEVELOP   optional; the merge base (default dc66941969)
#   PR        optional; the PR head (default 1315143b74)
#   PR3433    optional; the head of #3433 (default 709a1a4d76)
# The script fetches the head of #3433 if the clone does not have it. The
# clone must have the PR branch.
set -eo pipefail
cd "$(dirname "$0")"

: "${CLONE:?set CLONE to a stan-dev/math clone}"
DEVELOP="${DEVELOP:-dc66941969}"
PR="${PR:-1315143b74}"
PR3433="${PR3433:-709a1a4d76}"

if ! git -C "$CLONE" cat-file -e "$PR3433^{commit}" 2>/dev/null; then
  git -C "$CLONE" fetch -q origin "+refs/pull/3433/head:refs/remotes/pr/3433"
fi

for pair in "develop:$DEVELOP" "pr:$PR" "pr3433:$PR3433"; do
  name="${pair%%:*}"
  commit="${pair#*:}"
  dir="out/trees/$name"
  rm -rf "$dir"
  mkdir -p "$dir"
  git -C "$CLONE" archive "$commit" stan test/unit | tar -x -C "$dir"
  echo "$name: $(git -C "$CLONE" rev-parse --short "$commit") -> $dir"
done
