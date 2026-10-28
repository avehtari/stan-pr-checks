#!/usr/bin/env python3
# References for gamma_lccdf_tail_reference.hpp (13 rows, second
# derivatives) in the #3408 tests.  Run: python3 refs/gen_refs.py
"""References for the gamma_lccdf deep-tail unit tests.

For gamma_lccdf(y | a, beta) with z = beta * y: log Q(a, z), d/da, d2/da2
and d/dy = -beta * z^(a-1) exp(-z) / (Gamma(a) Q).  Two routes for log Q,
d/da and d2/da2 (gammainc with mp.diff; scaled quadrature), compared where
gammainc converges.
"""
import mpmath as mp

from probe_lccdf_tail import ref_gammainc, ref_quad

mp.mp.dps = 50

POINTS = [  # (alpha, y, beta)
    (5, 6.5, 1), (5, 7.5, 1),
    (150, 1200, 1), (300, 1500, 1), (150, 600, 2),
    (10, 1000, 1), (20, 1000, 1), (200, 1600, 1), (1000, 3000, 1),
    (1e6, 1.1e6, 1),
    (700, 2100, 1), (1e4, 10100, 1), (200, 200, 1),
]


def d2_routes(a, z):
  q2 = mp.diff(lambda s: ref_quad(s, z)[1], a)
  try:
    zz = mp.mpf(z)
    lq = lambda s: mp.log(mp.gammainc(s, zz, mp.inf, regularized=True))
    g2 = mp.diff(lq, mp.mpf(a), 2)
  except mp.libmp.libhyper.NoConvergence:
    g2 = None
  return q2, g2


def main():
  for a, y, beta in POINTS:
    z = mp.mpf(y) * beta
    qv, qd = ref_quad(a, z)
    q2, g2 = d2_routes(a, z)
    try:
      gv, gd = ref_gammainc(a, z)
      agree = max(abs((qv - gv) / gv), abs((qd - gd) / gd),
                  abs((q2 - g2) / g2))
      agree_s = mp.nstr(agree, 2)
    except mp.libmp.libhyper.NoConvergence:
      agree_s = "quad only"
    dy = -beta * mp.exp((a - 1) * mp.log(z) - z - mp.loggamma(a) - qv)
    print(f"{{{a:g}, {y:g}, {beta:g}, {mp.nstr(qv, 17, min_fixed=-1, max_fixed=1)}, "
          f"{mp.nstr(qd, 17, min_fixed=-1, max_fixed=1)}, "
          f"{mp.nstr(q2, 17, min_fixed=-1, max_fixed=1)}, "
          f"{mp.nstr(dy, 17, min_fixed=-1, max_fixed=1)}}},  // {agree_s}")


if __name__ == "__main__":
  main()
