#!/bin/bash
# Which public functions pass an rvalue arena_matrix to ref_type_if, or an
# array var_value to as_array_or_scalar?  Copies stan/ of each tree, puts the
# develop helpers back with a static_assert, and compiles reach_probe.cpp.
# Usage (from pr3406/): TBB_LIB=<dir> varmat/reach.sh <name>=<tree> ...
#   example: varmat/reach.sh dev=$DEV head=$HEAD opened=$OPENED
set -u
OUT=out/varmat
for arg in "$@"; do
  tree=${arg%%=*}
  SM=${arg#*=}
  W=$OUT/instr_$tree
  rm -rf "$W"
  mkdir -p "$W"
  cp -r "$SM/stan" "$W/stan"
  ln -s "$(cd "$SM/lib" && pwd -P)" "$W/lib"
  python3 - "$W" <<'EOF'
import re, sys
w = sys.argv[1]
p = w + "/stan/math/prim/meta/ref_type.hpp"
s = open(p).read()
i = s.index("struct ref_type_if<Condition, T, require_arena_matrix_t<T>> {")
j = s.index("};", i)
s = s[:i] + """struct ref_type_if<Condition, T, require_arena_matrix_t<T>> {
  static_assert(!std::is_rvalue_reference<T>::value,
                "PROBE: ref_type_if of an rvalue arena_matrix");
  using type =
      typename ref_type_if<Condition, typename std::decay_t<T>::Base>::type;
""" + s[j:]
open(p, "w").write(s)
p = w + "/stan/math/rev/fun/as_array_or_scalar.hpp"
s = open(p).read()
i = s.index("template <typename T, require_var_matrix_t<T>* = nullptr>")
j = s.index("\n}\n", i) + 3
s = s[:i] + """template <typename T, require_var_matrix_t<T>* = nullptr>
inline auto as_array_or_scalar(T&& v) {
  static_assert(!is_eigen_array<value_type_t<T>>::value,
                "PROBE: as_array_or_scalar of an array var_value");
  return v.array();
}
""" + s[j:]
open(p, "w").write(s)
EOF
  log=$OUT/reach_$tree.log
  echo "== $tree"
  if STAN_MATH=$W EXTRA_CXXFLAGS="-fsyntax-only" OPT=-O0 \
      timeout 900 ./build.sh varmat/reach_probe.cpp $OUT/reach_$tree > "$log" 2>&1; then
    echo "compile: ok (no static_assert fired)"
  else
    echo "compile: FAILED"
    echo "static_asserts fired: $(grep -c 'PROBE:' "$log")"
    grep -E "required from (here|.*'auto stan::math::[a-z_A-Z]+(_lcdf|_lccdf|_cdf|_lpdf|Phi))" "$log" \
      | grep -oE "stan::math::[A-Za-z_]+(lcdf|lccdf|cdf|lpdf|Phi)[A-Za-z_]*" | sort | uniq -c
    grep -m 5 "error:" "$log" | cut -c1-250
  fi
done
