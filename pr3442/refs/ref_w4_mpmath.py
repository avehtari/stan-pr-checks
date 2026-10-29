#!/usr/bin/env python3
"""
ref_w4_mpmath.py -- mpmath references (50 digits) for the p4 and p5 points of
probe_w4.cpp.  Claims of stan-dev/math#3442: the reference of the 4- and
5-parameter accuracy numbers.

Usage: ref_w4_mpmath.py PROBE_W4_OUTPUT > REF_FILE

Reads the points from the probe output (any tree). Stan's functions are for
the upper boundary; the formulas below are for the lower boundary, used with
w -> 1 - w and v -> -v.

  standard density f(t | 0, 1, w): small-time series for t < 1 (k = -12..12),
      large-time series for t >= 1 (k = 1..59)
  density with drift variability (Blurton, Kesselmeier and Gondan, 2017):
      p(t) = exp((sv^2 a^2 w^2 - 2 a v w - v^2 t) / (2 (1 + sv^2 t)))
             / sqrt(1 + sv^2 t) * f(t / a^2 | 0, 1, w) / a^2
  sv = 0: P (absorption probability, closed form) and
      Q = P - F = (pi / a^2) e^{-v a w} sum_k k sin(k pi w) e^{-lam_k t} / lam_k,
      lam_k = v^2 / 2 + k^2 pi^2 / (2 a^2), k = 1..199; F = P - Q
  sv > 0: F = int_0^t p, P - F = int_t^inf p (mpmath.quad)

Output rows: <fn> <id> <value> <d/dy d/da d/dt0 d/dw d/dv d/dsv>, gradients
by mpmath.diff; "nan" where not computed (lcdf and lccdf with sv > 0).
"""
import sys

from mpmath import diff, exp, inf, log, mp, mpf, pi, quad, sin, sqrt

mp.dps = 50


def f_std(ts, w):
  if ts < 1:
    s = mpf(0)
    for k in range(-12, 13):
      wk = w + 2 * k
      s += wk * exp(-wk**2 / (2 * ts))
    return s / sqrt(2 * pi * ts**3)
  s = mpf(0)
  for k in range(1, 60):
    s += k * exp(-k**2 * pi**2 * ts / 2) * sin(k * pi * w)
  return pi * s


def dens_upper(t, a, w, v, sv):
  w = 1 - w
  v = -v
  q = 1 + sv**2 * t
  return (exp((sv**2 * a**2 * w**2 - 2 * a * v * w - v**2 * t) / (2 * q))
          / sqrt(q) * f_std(t / a**2, w) / a**2)


def p_q_upper(t, a, w, v):
  """(P, Q) for the upper boundary, sv = 0."""
  w = 1 - w
  v = -v
  if v != 0:
    P = (exp(-2 * v * a * w) - exp(-2 * v * a)) / (1 - exp(-2 * v * a))
  else:
    P = 1 - w
  c = pi / a**2 * exp(-v * a * w)
  Q = c * sum(k * sin(k * pi * w) * exp(-(v**2 / 2 + k**2 * pi**2
                                          / (2 * a**2)) * t)
              / (v**2 / 2 + k**2 * pi**2 / (2 * a**2)) for k in range(1, 200))
  return P, Q


def grads(fun, y, a, t0, w, v, sv, with_sv):
  x = [y, a, t0, w, v, sv]
  out = []
  for k in range(6):
    if k == 5 and not with_sv:
      out.append(mpf(0))
      continue
    def g(z, k=k):
      xx = list(x)
      xx[k] = z
      return fun(*xx)
    out.append(diff(g, x[k]))
  return out


def lpdf(y, a, t0, w, v, sv):
  return log(dens_upper(y - t0, a, w, v, sv))


def lcdf4(y, a, t0, w, v, sv):
  P, Q = p_q_upper(y - t0, a, w, v)
  return log(P - Q)


def lccdf4(y, a, t0, w, v, sv):
  return log(p_q_upper(y - t0, a, w, v)[1])


def main():
  seen = set()
  for line in open(sys.argv[1]):
    f = line.rstrip("\n").split("\t")
    fn, pid = f[0], f[2]
    if fn == "OLD" or (fn, pid) in seen:
      continue
    seen.add((fn, pid))
    y, a, t0, w, v, sv = (mpf(z) for z in f[3:9])
    has_sv = sv != 0
    if fn == "LPDF":
      val = lpdf(y, a, t0, w, v, sv)
      g = grads(lpdf, y, a, t0, w, v, sv, True)
    elif not has_sv:
      fun = lcdf4 if fn == "LCDF" else lccdf4
      val = fun(y, a, t0, w, v, sv)
      g = grads(fun, y, a, t0, w, v, sv, False)
    else:
      t = y - t0
      if fn == "LCDF":
        val = log(quad(lambda s: dens_upper(s, a, w, v, sv), [0, t]))
      else:
        val = log(quad(lambda s: dens_upper(s, a, w, v, sv), [t, t + 1, inf]))
      g = None
    print("%s\t%s\t%s\t%s" % (
        fn, pid, mp.nstr(val, 17),
        "\t".join("nan" if g is None else mp.nstr(z, 17)
                  for z in (g if g else range(6)))))
    sys.stdout.flush()


if __name__ == "__main__":
  main()
