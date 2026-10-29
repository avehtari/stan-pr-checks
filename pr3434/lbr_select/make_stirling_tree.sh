#!/bin/bash
# Copy a Stan Math stan/ tree and force internal::log_beta_ratio to the
# Stirling form for large shapes (no selection by term size).  Claim P3:
# "without the selection, points with N = 1e6 and one small shape lost up to
# 4e-11 where develop was correct".
#
# Usage: lbr_select/make_stirling_tree.sh <tree with stan/> <new tree>
set -eo pipefail

if [ $# -ne 2 ]; then
  echo "usage: $0 <tree with stan/> <new tree>" >&2
  exit 2
fi
SRC="$1"
DST="$2"
if [ -e "$DST" ]; then
  echo "$DST exists; choose a new directory" >&2
  exit 1
fi
mkdir -p "$DST"
cp -r "$SRC/stan" "$DST/stan"
F="$DST/stan/math/prim/fun/log_rising_factorial_ratio.hpp"
OLD='if (size_stirling <= size_lbeta) {'
NEW='if (true || size_stirling <= size_lbeta) {'
if [ "$(grep -cF "$OLD" "$F")" != 1 ]; then
  echo "expected one selection line in $F" >&2
  exit 1
fi
sed -i "s/if (size_stirling <= size_lbeta) {/$NEW/" "$F"
diff -r "$SRC/stan" "$DST/stan" || true
