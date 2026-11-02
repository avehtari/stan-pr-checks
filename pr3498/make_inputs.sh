#!/bin/bash
# Make the two inputs of run_branch_bench.sh from a stan-dev/math clone:
#   dev_normal_lpdf.hpp    the develop normal_lpdf, renamed normal_lpdf_dev
#   br_<commit>.tar.gz     stan/ of the branch commit, under tree/
#
# Usage: CLONE=<stan-dev/math clone> ./make_inputs.sh <develop commit> <branch commit>
# Example: CLONE=~/math ./make_inputs.sh dc66941969 760f6a1f65
set -eo pipefail

if [ $# -ne 2 ]; then
  echo "usage: CLONE=<clone> $0 <develop commit> <branch commit>" >&2
  exit 2
fi
: "${CLONE:?set CLONE to a stan-dev/math clone}"
DEV="$1"
BR="$2"

git -C "$CLONE" show "$DEV:stan/math/prim/prob/normal_lpdf.hpp" \
  | sed -e 's/STAN_MATH_PRIM_PROB_NORMAL_LPDF_HPP/NORMAL_LPDF_DEV_COPY_HPP/g' \
        -e 's/normal_lpdf(/normal_lpdf_dev(/g' \
        -e 's/normal_lpdf<false>(/normal_lpdf_dev<false>(/g' \
  > dev_normal_lpdf.hpp
git -C "$CLONE" archive --prefix=tree/ "$BR" stan | gzip > "br_$BR.tar.gz"
echo "wrote dev_normal_lpdf.hpp ($DEV) and br_$BR.tar.gz"
