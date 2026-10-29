#!/usr/bin/env python3
"""Errors of neg_binomial_lpmf at the points of stan-dev/math#1578.

Usage: score_1578.py <refs> <probe_sib output> [<probe_sib output> ...]
Prints per point the relative error of the value, d/dalpha and d/dbeta,
and the absolute error of alpha * d/dalpha (the log-scale gradient).
The references (ref_sib.py) have routes A and B; their worst disagreement
is printed first.
"""

import sys

from mpmath import mp, mpf

mp.dps = 60


def main():
  refs = {}
  worst = mpf(0)
  for line in open(sys.argv[1]):
    key, rest = line.rstrip("\n").split(" => ")
    a, b, _ = rest.split(" | ")
    a = [mpf(x) for x in a.split()]
    b = [mpf(x) for x in b.split()]
    for x, y in zip(a, b):
      worst = max(worst, abs(x - y) / max(abs(x), abs(y), mpf("1e-300")))
    refs[tuple(float.fromhex(h) for h in key.split()[2:])] = a
  print(f"route A against route B, worst relative: {mp.nstr(worst, 3)}")
  print(f"{'file':14s} {'n':>2s} {'alpha':>8s} {'value rel':>10s} "
        f"{'d/da rel':>10s} {'d/db rel':>10s} {'a*d/da abs':>10s} "
        f"{'d/dalpha':>20s}")
  for path in sys.argv[2:]:
    for line in open(path):
      left, right = line.split(" => ")
      args = tuple(float.fromhex(h) for h in left.split()[2:])
      v = [float.fromhex(x) for x in right.split()[:3]]
      r = refs[args]
      rels = [abs(mpf(x) - y) / abs(y) for x, y in zip(v, r)]
      ga = abs(mpf(args[1]) * (mpf(v[1]) - r[1]))
      name = path.rsplit("/", 1)[-1]
      print(f"{name:14s} {int(args[0]):2d} {args[1]:8.2g} "
            + " ".join(f"{float(e):10.1e}" for e in rels)
            + f" {float(ga):10.1e} {v[1]:20.16g}")
  return 0


if __name__ == "__main__":
  sys.exit(main())
