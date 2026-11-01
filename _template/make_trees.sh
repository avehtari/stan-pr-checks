#!/bin/bash
# Export stan/ and test/unit/ of each named commit into out/trees/<name>/.
# The probes build against these trees (STAN_MATH=out/trees/<name>).
#
# Usage: CLONE=<stan-dev/math clone> ./make_trees.sh name:commit [name:commit ...]
#   for example: ./make_trees.sh develop:dc66941969 pr:1315143b74
#
# A commit of an open PR that the clone does not have: give it as
# name:pull/<N> and the script fetches refs/pull/<N>/head first.
# Stacked PRs: export the head of the stack that the claim needs, and say
# in the README which PRs it contains.
set -eo pipefail
cd "$(dirname "$0")"

: "${CLONE:?set CLONE to a stan-dev/math clone}"
[ $# -gt 0 ] || { echo "usage: CLONE=<clone> $0 name:commit ..." >&2; exit 2; }

for pair in "$@"; do
  name="${pair%%:*}"
  commit="${pair#*:}"
  case "$commit" in
    pull/*)
      n="${commit#pull/}"
      git -C "$CLONE" fetch -q origin "+refs/pull/$n/head:refs/remotes/pr/$n"
      commit="pr/$n" ;;
  esac
  dir="out/trees/$name"
  rm -rf "$dir"
  mkdir -p "$dir"
  git -C "$CLONE" archive "$commit" stan test/unit | tar -x -C "$dir"
  echo "$name: $(git -C "$CLONE" rev-parse --short "$commit") -> $dir"
done
