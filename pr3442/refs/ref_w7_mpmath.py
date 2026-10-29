#!/usr/bin/env python3
"""
ref_w7_mpmath.py -- mpmath references for single points of the 7-parameter
lpdf, lcdf and lccdf (upper boundary).  Claim of stan-dev/math#3442: the second
route of the 7-parameter test references (the first is w7/probe_ref7.cpp).

Usage: ref_w7_mpmath.py FN ID y a t0 w v sv sw st0
  FN in {LPDF, LCDF, LCCDF}.  Prints one row:
  FN ID value d/dy d/da d/dt0 d/dw d/dv d/dsv d/dsw d/dst0 h_rel
Gradients by central differences of the log value with step
h = H_REL * max(|x|, X_MIN), at mp.dps = 30.  A parameter that is 0 among
sv, sw, st0 gets the gradient 0 (Stan's convention).

Model, with t = y - t0 and p5(s | w) the 5-parameter density (drift
variability in closed form, Blurton, Kesselmeier and Gondan, 2017):
  lpdf : (1/sw) int_{w - sw/2}^{w + sw/2} dw' (1/st0) int_{t - st0}^{t}
         p5(s | w') ds          (s > 0; sw = 0 or st0 = 0: no average)
  lcdf : (1/sw) int dw' (1/st0) int_0^t p5(s | w') min(t - s, st0) ds
  lccdf: (1/sw) int dw' (1/st0) int_0^inf p5(s | w')
         clamp(s - t + st0, 0, st0) ds
  (st0 = 0: the weights become the indicators s <= t and s > t.)
Each integrand is divided by p5 at a node near its peak, so that mp.quad,
whose stop rule is absolute, sees values near 1.
"""
import sys

from mpmath import exp, inf, log, mp, mpf, pi, quad, sin, sqrt

mp.dps = 30
H_REL = mpf("1e-9")
X_MIN = mpf("1e-3")


def f_std(ts, w):
  if ts < 1:
    s = mpf(0)
    for k in range(-15, 16):
      wk = w + 2 * k
      s += wk * exp(-wk**2 / (2 * ts))
    return s / sqrt(2 * pi * ts**3)
  s = mpf(0)
  for k in range(1, 80):
    s += k * exp(-k**2 * pi**2 * ts / 2) * sin(k * pi * w)
  return pi * s


def p5(t, a, w, v, sv):
  """Density of the upper boundary at time t > 0."""
  if t <= 0:
    return mpf(0)
  w = 1 - w
  v = -v
  q = 1 + sv**2 * t
  return (exp((sv**2 * a**2 * w**2 - 2 * a * v * w - v**2 * t) / (2 * q))
          / sqrt(q) * f_std(t / a**2, w) / a**2)


def time_integral(fn, t, a, w, v, sv, st0, scale):
  """The integral over s for one w, divided by scale."""
  f = lambda s: p5(s, a, w, v, sv) / scale
  if fn == "LPDF":
    if st0 == 0:
      return f(t)
    lo = max(t - st0, mpf(0))
    return quad(f, [lo, (lo + t) / 2, t]) / st0
  if fn == "LCDF":
    if st0 == 0:
      pts = [mpf(0), t / 4, t / 2, 3 * t / 4, t]
      return quad(f, pts)
    kink = t - st0
    g = lambda s: f(s) * min(t - s, st0) / st0
    pts = sorted(set([mpf(0), t / 4, t / 2, 3 * t / 4, t]
                     + ([kink] if kink > 0 else [])))
    return quad(g, pts)
  # LCCDF
  tail = [t, t + mpf("0.01"), t + mpf("0.1"), t + 1, t + 10, inf]
  if st0 == 0:
    return quad(f, tail)
  lo = max(t - st0, mpf(0))
  g = lambda s: f(s) * min(max(s - t + st0, 0), st0) / st0
  head = quad(g, [lo, (lo + t) / 2, t]) if t > lo else mpf(0)
  return head + quad(f, tail)


def value(fn, y, a, t0, w, v, sv, sw, st0):
  t = y - t0
  # scale: p5 at the end of the time range, or at its peak for an lcdf
  if fn == "LCDF":
    ss = [t * k / 64 for k in range(1, 65)]
    scale = max(p5(s, a, w, v, sv) for s in ss)
  else:
    scale = p5(t, a, w, v, sv)
    if scale == 0:
      scale = p5(t + st0, a, w, v, sv)
  if sw == 0:
    return log(scale) + log(time_integral(fn, t, a, w, v, sv, st0, scale))
  inner = lambda ww: time_integral(fn, t, a, ww, v, sv, st0, scale)
  return log(scale) + log(quad(inner, [w - sw / 2, w, w + sw / 2]) / sw)


def main():
  fn, pid = sys.argv[1], sys.argv[2]
  x = [mpf(z) for z in sys.argv[3:11]]
  val = value(fn, *x)
  grads = []
  for k in range(8):
    if k >= 5 and x[k] == 0:
      grads.append(mpf(0))
      continue
    h = H_REL * max(abs(x[k]), X_MIN)
    xp, xm = list(x), list(x)
    xp[k] += h
    xm[k] -= h
    grads.append((value(fn, *xp) - value(fn, *xm)) / (2 * h))
  print("%s\t%s\t%s\t%s\t%s" % (fn, pid, mp.nstr(val, 20),
                                "\t".join(mp.nstr(g, 15) for g in grads),
                                mp.nstr(H_REL, 3)))
  sys.stdout.flush()


if __name__ == "__main__":
  main()
