#!/usr/bin/env python3
# Prototype and two-route check of the log-space lower series (log P, d/da).
# Supports the lower-tail part of the log-space cdfs (Part 3 of #3408).
"""Prototype: log P(a, z) and d/da log P by the power series in log space.

P(a, z) = z^a e^-z / Gamma(a + 1) * S,  S = sum_k t_k,  t_0 = 1,
t_k = t_{k-1} z / (a + k).  d/da log t_k = r_k = r_{k-1} - 1 / (a + k).
d/da log P = log z - digamma(a + 1) + S' / S,  S' = sum_k t_k r_k.
Stop when a bound on the remaining terms is below eps of S and of S':
for q = z / (a + k + 1) < 1 the tail of t is at most t_k q / (1 - q).

Two reference routes: mpmath gammainc with mp.diff, and a quadrature of
the density on [0, z] scaled to 1 at t = z (valid for z < a - 1).
"""
import math
import mpmath as mp

mp.mp.dps = 50
EPS = 2.220446049250313e-16


def digamma(a):
  return float(mp.digamma(a))


def series(a, z, tol=EPS, max_steps=100000):
  t, r = 1.0, 0.0
  s, ds = 1.0, 0.0
  k = 0
  while k < max_steps:
    k += 1
    t *= z / (a + k)
    r -= 1.0 / (a + k)
    s += t
    ds += t * r
    q = z / (a + k + 1)
    if q < 1:
      tail = t * q / (1 - q)
      if tail <= tol * s and tail * (abs(r) + 1) <= tol * abs(ds):
        break
  lp = a * math.log(z) - z - math.lgamma(a + 1) + math.log(s)
  return lp, math.log(z) - digamma(a + 1) + ds / s, k


def ref_gammainc(a, z):
  a, z = mp.mpf(a), mp.mpf(z)
  f = lambda s: mp.log(mp.gammainc(s, 0, z, regularized=True))
  return f(a), mp.diff(f, a)


def ref_quad_lower(a, z):
  a, z = mp.mpf(a), mp.mpf(z)
  lw = lambda s: (a - 1) * mp.log1p(-s / z) + s
  dec = 1 / ((a - 1) / z - 1)
  nodes = [mp.mpf(0)] + [dec * 2**k for k in range(-4, 12) if dec * 2**k < z] + [z]
  i0 = mp.quad(lambda s: mp.exp(lw(s)), nodes)
  i1 = mp.quad(lambda s: mp.exp(lw(s)) * mp.log1p(-s / z), nodes)
  lp = (a - 1) * mp.log(z) - z - mp.loggamma(a) + mp.log(i0)
  return lp, mp.log(z) - mp.digamma(a) + i1 / i0


POINTS = [(1000, 100), (500, 30), (100, 1e-3), (200, 50), (50, 5), (0.5, 1e-300),
          (1e4, 6800), (1e6, 963000), (1e6, 950000), (171.6, 40), (3000, 2400)]


def main():
  print(f"{'a':>8} {'z':>8} {'logP':>10} {'val_err':>8} {'der_err':>8} {'steps':>6} "
        f"{'routes':>8}")
  for a, z in POINTS:
    lq, dq = ref_quad_lower(a, z)
    try:
      lg, dg = ref_gammainc(a, z)
      agree = float(max(abs((lq - lg) / lg), abs((dq - dg) / dg)))
    except mp.libmp.libhyper.NoConvergence:
      agree = float("nan")
    lp, dlp, k = series(a, z)
    print(f"{a:8g} {z:8g} {float(lq):10.1f} {float(abs((lp - lq) / lq)):8.1e} "
          f"{float(abs((dlp - dq) / dq)):8.1e} {k:6d} {agree:8.1e}")


if __name__ == "__main__":
  main()
