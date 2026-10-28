#!/usr/bin/env python3
"""Check the reference at the worst relative errors of the PR tree.

For each region, take the point and partial with the largest relative
error of the PR (where |r| / I >= 1e-3), compute both reference routes
again (refs/ref_mpmath.py and route A of refs/gen_test_refs.py), and print
their agreement and the error of the PR against the new reference.  If the
routes agree and the error does not change, the error is in the code.

Usage: python3 sweep/check_worst.py out/cmp_double_pr.txt \
         sweep/refs_beta_18000.txt
"""

import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "refs"))
sys.path.insert(0, HERE)
import ref_mpmath as R  # noqa: E402
import table_pr as T  # noqa: E402
from gen_test_refs import route_a  # noqa: E402

DBL_MIN = 2.2250738585072014e-308


def main():
  pr = T.parse(sys.argv[1])
  refs = T.load_refs(sys.argv[2])
  best = {}
  for reg, key, v, _ in pr:
    i, da, db = refs[key]
    if math.isnan(da) or i < DBL_MIN:
      continue
    ref = (da, db, da, db)
    for c in range(4):
      if abs(ref[c]) / i >= 1e-3 and abs(ref[c]) >= DBL_MIN:
        r = abs(v[c] - ref[c]) / abs(ref[c])
        if r > best.get(reg, (0.0,))[0]:
          best[reg] = (r, key, c, v[c])
  print(f"{'region':<10} {'a':>10} {'b':>10} {'z':>10} {'wrt':>3} "
        f"{'rel err':>9} {'routes':>7} {'rel err (new ref)':>17}")
  for reg, (r, key, c, v) in sorted(best.items(), key=lambda x: -x[1][0]):
    a, b, z = key
    wrt = "a" if c % 2 == 0 else "b"
    vb = (R.dIda if wrt == "a" else R.dIdb)(a, b, z)
    va = route_a(a, b, z, wrt)
    dig = float(R.rel_digits(va, vb))
    new = abs(v - float(vb)) / abs(float(vb))
    print(f"{reg[:10]:<10} {a:>10.6g} {b:>10.6g} {z:>10.6g} {wrt:>3} "
          f"{r:>9.2e} {dig:>7.1f} {new:>17.2e}")
  print("routes: agreeing digits of routes A and B (inf: at least 60)")
  return 0


if __name__ == "__main__":
  sys.exit(main())
