#!/bin/bash
# Pack the stan/ directories of the develop tree and the PR tree into
# speed/trees.tar.gz for run_cmp_beta_speed.sh.
#
#   WT     required; a clone or worktree of the PR branch (read only)
#   BASE   optional; the develop commit (default: merge base with
#          origin/develop)
# The PR tree is the working tree of WT, so uncommitted changes count.
set -euo pipefail

: "${WT:?set WT to a clone of the PR branch}"
HERE="$(cd "$(dirname "$0")" && pwd)"
G="git --no-optional-locks -C $WT"
BASE="${BASE:-$($G merge-base HEAD origin/develop)}"
OUT="$(mktemp -d)"
mkdir -p "$OUT/dev" "$OUT/pr"
$G archive "$BASE" stan | tar -x -C "$OUT/dev"
cp -r "$WT/stan" "$OUT/pr/"
{
  echo "== dev: $($G rev-parse --short "$BASE")"
  echo "== pr:  $($G rev-parse --short HEAD) plus the working tree;" \
       "changed files: $($G status --short -- stan | wc -l)"
} > "$OUT/commits.txt"
tar -czf "$HERE/trees.tar.gz" -C "$OUT" dev pr commits.txt
rm -rf "$OUT"
ls -la "$HERE/trees.tar.gz"
