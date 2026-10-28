#!/bin/bash
# Compile one C++ probe against a Stan Math tree.
# Used by every C++ probe of stan-dev/math#3409 in this directory.
#
# Usage: STAN_MATH=<tree> ./build.sh <source.cpp> <output binary>
#
# Environment:
#   STAN_MATH       required; a Stan Math tree that holds stan/
#   STAN_MATH_LIB   optional; the lib/ directory (default $STAN_MATH/lib).
#                   Set it when STAN_MATH holds stan/ only, for example a
#                   `git archive` of an older commit.
#   TBB_LIB         optional; a directory with a built libtbb.so
#                   (default $STAN_MATH_LIB/tbb)
#   CXX             optional; the compiler (default g++)
#   OPT             optional; the optimization flags (default -O3)
#   EXTRA_CXXFLAGS  optional; appended to the compile flags
set -eo pipefail

if [ $# -ne 2 ]; then
  echo "usage: STAN_MATH=<tree> $0 <source.cpp> <output binary>" >&2
  exit 2
fi
SRC="$1"
OUT="$2"

: "${STAN_MATH:?set STAN_MATH to a Stan Math tree}"
STAN_MATH_LIB="${STAN_MATH_LIB:-$STAN_MATH/lib}"
TBB_LIB="${TBB_LIB:-$STAN_MATH_LIB/tbb}"
CXX="${CXX:-g++}"
OPT="${OPT:--O3}"

for d in "$STAN_MATH/stan" "$STAN_MATH_LIB/eigen_5.0.1" \
         "$STAN_MATH_LIB/boost_1.87.0" "$STAN_MATH_LIB/sundials_6.1.1/include" \
         "$STAN_MATH_LIB/tbb_2020.3/include"; do
  if [ ! -d "$d" ]; then
    echo "build.sh: missing directory $d" >&2
    exit 1
  fi
done
if [ ! -e "$TBB_LIB/libtbb.so" ]; then
  echo "build.sh: no libtbb.so in $TBB_LIB" >&2
  echo "build.sh: build it with 'make -f make/standalone math-libs' in a" \
       "Stan Math clone and set TBB_LIB to its lib/tbb" >&2
  exit 1
fi

# -pthread is necessary: prim/fun/lgamma.hpp requires _REENTRANT.
# EXTRA_CXXFLAGS is not quoted on purpose: it can hold several flags.
# shellcheck disable=SC2086
"$CXX" -std=c++17 $OPT -pthread -DBOOST_DISABLE_ASSERTS $EXTRA_CXXFLAGS \
  -I "$STAN_MATH" \
  -I "$STAN_MATH_LIB/eigen_5.0.1" \
  -I "$STAN_MATH_LIB/boost_1.87.0" \
  -I "$STAN_MATH_LIB/sundials_6.1.1/include" \
  -I "$STAN_MATH_LIB/tbb_2020.3/include" \
  "$SRC" -o "$OUT" \
  -L "$TBB_LIB" -ltbb -Wl,-rpath,"$TBB_LIB"
echo "built $OUT from $SRC against $STAN_MATH"
