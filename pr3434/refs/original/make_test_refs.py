#!/usr/bin/env python3
"""References for the beta_binomial_lpmf unit tests (stan-dev/math#3153).

Usage: make_test_refs.py points   write test_points.txt for the probes
       make_test_refs.py table    read out_test_{develop,branch}.txt and
                                  print the references with both errors and
                                  C++ initializers

Values: 17 significant digits (route A of ref_mpmath, checked against
route B there). Gradients are given as log-shape gradients alpha * d/dalpha,
beta * d/dbeta, and for tied shapes d/dlc = s * (d/dalpha + d/dbeta).
"""

import math
import sys

from mpmath import mp, mpf

import ref_mpmath as R
import score_bb as S


def tied(lc):
  return math.exp(lc) / 2


POINTS = [
    # tag, n, N, alpha, beta
    ("REP32", 57, 117, tied(32.0), tied(32.0)),
    ("REP36", 57, 117, tied(36.0), tied(36.0)),
    ("REP42", 57, 117, tied(42.0), tied(42.0)),
    ("REP16", 57, 117, tied(16.0), tied(16.0)),
    ("REP18", 57, 117, tied(18.0), tied(18.0)),
    ("REP20", 57, 117, tied(20.0), tied(20.0)),
    ("REP26", 57, 117, tied(26.0), tied(26.0)),
    ("REP30", 57, 117, tied(30.0), tied(30.0)),
    ("UNEQ", 400, 1000, 3e13, 7e12),
    ("EDGE0", 0, 117, 1e17, 3e17),
    ("EDGEN", 117, 117, 1e17, 3e17),
    ("SEAM1", 3, 1000000, 1e6, 1e18),
    ("SEAM2", 0, 1000000, 5.85e17, 4.84),
    ("MIXED", 1000000, 1000000, 950000.0, 15.0),
    ("SMALL", 5, 20, 10.0, 25.0),
]


def write_points():
  with open("test_points.txt", "w") as fh:
    for tag, n, N, a, b in POINTS:
      fh.write(f"{tag} {n} {N} {float(a).hex()} {float(b).hex()}\n")
  print("wrote test_points.txt")


def cxx(x):
  return mp.nstr(x, 17, min_fixed=-1, max_fixed=-1, strip_zeros=False)


def table():
  dev = {r["tag"]: r for r in S.parse("out_test_develop.txt")}
  br = {r["tag"]: r for r in S.parse("out_test_branch.txt")}
  print(f"{'tag':6s} {'value':>24s} {'dev err':>8s} {'br err':>8s} "
        f"{'gc / d/dlc':>24s} {'dev err':>8s} {'br err':>8s}")
  for tag, n, N, a, b in POINTS:
    ref = R.ref_point(n, N, a, b)
    a_m, b_m = mpf(a), mpf(b)
    gc = ref["gc"]
    out = []
    for r in (dev[tag], br[tag]):
      ev = abs(mpf(r["v_var"]) - ref["value"])
      g = a_m * mpf(r["d_alpha"]) + b_m * mpf(r["d_beta"])
      eg = abs(g - gc)
      out.append((ev, eg))
    print(f"{tag:6s} {cxx(ref['value']):>24s} {float(out[0][0]):8.1e} "
          f"{float(out[1][0]):8.1e} {cxx(gc):>24s} {float(out[0][1]):8.1e} "
          f"{float(out[1][1]):8.1e}")
  print("\nC++ initializers {n, N, alpha, beta, value, ga, gb}:")
  for tag, n, N, a, b in POINTS:
    ref = R.ref_point(n, N, a, b)
    print(f"    {{{n}, {N}, {float(a).hex()}, {float(b).hex()}, "
          f"{cxx(ref['value'])}, {cxx(ref['ga'])}, {cxx(ref['gb'])}}},"
          f"  // {tag}")


if __name__ == "__main__":
  if len(sys.argv) != 2 or sys.argv[1] not in ("points", "table"):
    print(__doc__)
    sys.exit(2)
  if sys.argv[1] == "points":
    write_points()
  else:
    table()
