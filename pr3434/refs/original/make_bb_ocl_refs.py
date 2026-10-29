#!/usr/bin/env python3
"""References for the OpenCL beta_binomial_lpmf absolute test.

Uses ref_mpmath.ref_point() of beta_binomial_work (route A, checked against
route B there). Prints C++ rows {n, N, alpha, beta, value, ga, gb} with
ga = alpha * d/dalpha and gb = beta * d/dbeta.
"""

import os
import sys

from mpmath import mp

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                ".."))
import ref_mpmath as R  # noqa: E402

H = float.fromhex

POINTS = [
    # n, N, alpha, beta
    (57, 117, H("0x1.1f43fcc4b662cp+45"), H("0x1.1f43fcc4b662cp+45")),
    (400, 1000, H("0x1.b48eb57e00000p+44"), H("0x1.977420dc00000p+42")),
    (0, 117, H("0x1.6345785d8a000p+56"), H("0x1.0a741a4627800p+58")),
    (117, 117, H("0x1.6345785d8a000p+56"), H("0x1.0a741a4627800p+58")),
    # one shape below 10: the lbeta form inside log_beta_ratio
    (5, 117, 1e15, 4.84),
    (30, 117, 2.5, 1e14),
    (5, 20, 10.0, 25.0),
]


def cxx(x):
  return mp.nstr(x, 17, min_fixed=-1, max_fixed=-1, strip_zeros=False)


def main():
  for n, N, a, b in POINTS:
    ref = R.ref_point(n, N, a, b)
    print(f"    {{{n}, {N}, {float(a).hex()}, {float(b).hex()}, "
          f"{cxx(ref['value'])}, {cxx(ref['ga'])}, {cxx(ref['gb'])}}},")


if __name__ == "__main__":
  main()
