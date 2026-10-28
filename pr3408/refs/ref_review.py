#!/usr/bin/env python3
# Two-route references for the review probe (cdfs/probe_review.cpp).
# Supports the series-branch, #2006 fvar and gamma_lcdf(20000 | 800, 0.1)
# numbers.
"""Two-route mpmath references for the points added after the review.

  E2, E5  d/da Q(a, z) in the series branch of the old root
  E1      d/dt log Q(3 + t, (3 + t)(12 + t)) at t = 0, the fvar example
          with alpha = 3, beta = 3, y = 12 and all tangents 1
  E4      gamma_lcdf(20000 | 800, 0.1) = log P(800, 2000)

Route A: numerical differentiation of mpmath.gammainc (ref_mpmath.dQda_A).
Route B: an integral of the side that does not cancel (ref_mpmath.dQda_B,
probe_lccdf_tail.ref_quad).  Each line prints both routes and their
relative difference.

Usage: python3 refs/ref_review.py
"""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

from mpmath import mp, mpf  # noqa: E402

import ref_mpmath as rm  # noqa: E402
from probe_lccdf_tail import ref_quad  # noqa: E402

mp.dps = 80


def show(tag, a_val, b_val):
  d = abs(a_val - b_val) / abs(b_val)
  print(f"{tag:40s} A {mp.nstr(a_val, 20):>28s}  B {mp.nstr(b_val, 20):>28s}"
        f"  rel {mp.nstr(d, 3)}")


def main():
  # dQ/da, the quantity grad_reg_inc_gamma returns
  for a, z in [(0.5, 0.25), (8, 4), (170, 100), (165, 80), (150, 150)]:
    show(f"dQ/da ({a}, {z})", rm.dQda_A(a, z), rm.dQda_B(a, z))

  # E1: the path alpha = beta = 3 + t, y = 12 + t, so z = (3 + t)(12 + t)
  with mp.workdps(110):
    route_a = mp.diff(
        lambda t: mp.log(mp.gammainc(3 + t, (3 + t) * (12 + t), mp.inf,
                                     regularized=True)), mpf(0))
  q = mp.gammainc(3, 36, mp.inf, regularized=True)
  route_b = (rm.dQda_B(3, 36) + 15 * rm.dQdz(3, 36)) / q
  show("E1 d/dt log Q (3, 36), dz/dt = 15", route_a, route_b)
  lq_quad, _ = ref_quad(3, 36)
  show("E1 log Q(3, 36)", mp.log(q), lq_quad)

  # E4: log P(800, 2000) = log1p(-Q); Q from gammainc and from quadrature
  q_a = mp.gammainc(800, 2000, mp.inf, regularized=True)
  lq_b, _ = ref_quad(800, 2000)
  show("E4 log P(800, 2000)", mp.log1p(-q_a), mp.log1p(-mp.exp(lq_b)))


if __name__ == "__main__":
  main()
