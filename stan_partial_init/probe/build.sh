#!/bin/bash
# Build probe_init against a Stan tree, optionally with another
# initialize.hpp (the develop version) in front of the tree's src/.
#   STAN=<stan tree> STANC=<stanc> ./build.sh <out name> [<include dir>]
# <include dir> must contain stan/services/util/initialize.hpp.
# The Stan tree needs lib/stan_math with a built TBB in lib/stan_math/lib/tbb.
set -e -o pipefail
cd "$(dirname "$0")"
: "${STAN:?set STAN to a Stan tree}"
: "${STANC:?set STANC to a stanc binary}"
NAME=${1:?out name}
EXTRA=${2:-}
mkdir -p ../out
"$STANC" --o=../out/partial_init_dependent.hpp partial_init_dependent.stan
M=$STAN/lib/stan_math
INC=()
[ -n "$EXTRA" ] && INC+=(-I "$EXTRA")
g++ -std=c++17 -O2 -pthread -D_REENTRANT -DBOOST_DISABLE_ASSERTS \
  -Wno-ignored-attributes -Wno-deprecated-declarations \
  "${INC[@]}" -I ../out -I "$STAN/src" -I "$STAN/lib/rapidjson_1.1.0" -I "$M" \
  -isystem "$M/lib/eigen_5.0.1" -isystem "$M/lib/boost_1.87.0" \
  -isystem "$M/lib/sundials_6.1.1/include" -isystem "$M/lib/tbb_2020.3/include" \
  probe_init.cpp -o "../out/$NAME" \
  -L "$M/lib/tbb" -ltbb -Wl,-rpath,"$M/lib/tbb"
echo "built ../out/$NAME"
