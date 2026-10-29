#!/usr/bin/env python3
"""Score probe_series output: the student_t_lpdf nu terms against mpmath.

Usage: score_series.py <probe output>

References at 60 digits, two routes each:
  C(nu)  A: loggamma((nu + 1) / 2) - loggamma(nu / 2) - log(nu) / 2
         B: log(rf(nu / 2, 1/2)) - log(nu) / 2 (mp.rf, the gamma ratio)
  D(nu)  A: digamma((nu + 1) / 2) - digamma(nu / 2)
         B: 2 * mp.diff of route A of C + 1 / nu (d/dnu of C is D / 2 - 1/(2 nu))
Errors in eps = 2^-52 relative to the reference.
"""

import sys

from mpmath import mp, mpf

mp.dps = 60
EPS = mpf(2) ** -52


def c_a(nu):
  return mp.loggamma((nu + 1) / 2) - mp.loggamma(nu / 2) - mp.log(nu) / 2


def main():
  worst = {"C": (mpf(0), None), "D": (mpf(0), None)}
  routes = {"C": mpf(0), "D": mpf(0)}
  n = 0
  for line in open(sys.argv[1]):
    left, right = line.split(" => ")
    nu = mpf(float.fromhex(left))
    c, d = (mpf(float.fromhex(t)) for t in right.split())
    ca = c_a(nu)
    cb = mp.log(mp.rf(nu / 2, mpf(1) / 2)) - mp.log(nu) / 2
    da = mp.digamma((nu + 1) / 2) - mp.digamma(nu / 2)
    with mp.workdps(100):
      db = 2 * (mp.diff(c_a, nu) + 1 / (2 * nu))
    routes["C"] = max(routes["C"], abs(ca - cb) / abs(ca))
    routes["D"] = max(routes["D"], abs(da - db) / abs(da))
    for k, got, ref in (("C", c, ca), ("D", d, da)):
      e = abs(got - ref) / abs(ref) / EPS
      if e > worst[k][0]:
        worst[k] = (e, float(nu) / 2)
    n += 1
  print(f"{n} points, h = nu / 2 from 20 to 1e15")
  print(f"route A against route B, worst relative: C {mp.nstr(routes['C'], 3)}"
        f", D {mp.nstr(routes['D'], 3)}")
  for k in ("C", "D"):
    e, h = worst[k]
    print(f"{k}: worst error {float(e):.2f} eps at h = {h:.6g}")
  return 0


if __name__ == "__main__":
  sys.exit(main())
