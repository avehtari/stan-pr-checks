#!/usr/bin/env python3
"""References for the digamma_diff unit test.

Usage: make_dd_test_refs.py points   write dd_test_points.txt for probe_dd
       make_dd_test_refs.py table    read out_dd_test.txt, print the errors
                                     of digamma_diff and of the plain
                                     difference, and C++ initializers

psi(x + d) - psi(x) at enough digits that the cancellation (log10(x) digits
for large x, log10(1/d) digits for small d) leaves 40.
"""

import math
import sys

from mpmath import mp, mpf

POINTS = [
    (1e-3, 0.5), (0.5, 1.0), (1.0, 3.0), (2.5, 1e-8), (9.75, 0.25),
    (10.0, 57.0), (10.5, 1e6), (37.25, 0.5), (1e3, 117.0), (1e6, 3.0),
    (1e8, 1e-8), (1e12, 57.0), (3e13, 400.0), (1e15, 1e6), (1e18, 1.0),
    (8.6963747076025050e17, 117.0), (5e-324 * 2 ** 1000, 2.0),
]


def ref(x, d):
  extra = (max(0, int(math.log10(max(x, 1.0))))
           + max(0, int(-math.log10(min(d, 1.0)))))
  with mp.workdps(60 + extra):
    return mp.digamma(mpf(x) + mpf(d)) - mp.digamma(mpf(x))


def write_points():
  with open("dd_test_points.txt", "w") as fh:
    for x, d in POINTS:
      fh.write(f"{float(x).hex()} {float(d).hex()}\n")
  print("wrote dd_test_points.txt")


def table():
  rows = {}
  for line in open("out_dd_test.txt"):
    if line.startswith("#"):
      continue
    left, right = line.split("=>")
    xs, ds = left.split()
    vals = right.split()
    rows[(float.fromhex(xs), float.fromhex(ds))] = (
        float.fromhex(vals[0]), float.fromhex(vals[1]))
  print(f"{'x':>12s} {'d':>8s} {'reference':>24s} {'dd rel':>8s} "
        f"{'plain rel':>9s}")
  for x, d in POINTS:
    r = ref(x, d)
    v, p = rows[(float(x), float(d))]
    ev = float(abs(mpf(v) - r) / abs(r))
    ep = float(abs(mpf(p) - r) / abs(r))
    print(f"{x:12.4g} {d:8.3g} {mp.nstr(r, 17):>24s} {ev:8.1e} {ep:9.1e}")
  print("\nC++ initializers {x, d, value}:")
  for x, d in POINTS:
    r = ref(x, d)
    print(f"    {{{float(x).hex()}, {float(d).hex()}, "
          f"{mp.nstr(r, 17, min_fixed=-1, max_fixed=-1, strip_zeros=False)}}},")


if __name__ == "__main__":
  if len(sys.argv) != 2 or sys.argv[1] not in ("points", "table"):
    print(__doc__)
    sys.exit(2)
  if sys.argv[1] == "points":
    write_points()
  else:
    table()
