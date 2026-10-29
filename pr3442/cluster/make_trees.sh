#!/bin/bash
# Pack stan/ of the merge base and of the PR head of stan-dev/math#3442 into
# base3442.tar.gz and head3442.tar.gz (the trees of every rerun).
#
# Usage: STAN_MATH=<clone with both commits> cluster/make_trees.sh [OUTDIR]
#   BASE  optional; default 43bc5684a3 (merge base, "develop" in the PR text)
#   HEAD  optional; default ad0b7460d4 (PR head)
# The clone is read only (git archive of the commits; the working tree is
# not used).
set -euo pipefail

: "${STAN_MATH:?set STAN_MATH to a clone of stan-dev/math with the PR branch}"
BASE="${BASE:-43bc5684a3}"
HEAD="${HEAD:-ad0b7460d4}"
OUTDIR="${1:-.}"
mkdir -p "$OUTDIR"
for x in "base3442:$BASE" "head3442:$HEAD"; do
  name="${x%%:*}"
  commit="${x#*:}"
  git -C "$STAN_MATH" archive --format=tar "$commit" stan \
    | gzip -n > "$OUTDIR/$name.tar.gz"
  echo "$name.tar.gz: $(git -C "$STAN_MATH" rev-parse --short "$commit")" \
       "$(md5sum < "$OUTDIR/$name.tar.gz" | cut -d' ' -f1)"
done
