#!/usr/bin/env python3
"""
ref_w5grad_mpmath.py -- mpmath value and gradients of the 5-parameter Wiener
log density (upper boundary), compared with the output of probe_w5grad.cpp.
Claim of stan-dev/math#3442: the mpmath references of the 5-parameter test points.

Usage: ref_w5grad_mpmath.py PROBE_OUTPUT

Density (lower boundary; Blurton, Kesselmeier and Gondan, 2017):
  p(t | a, v, w, sv) = 1 / sqrt(1 + sv^2 t)
      * exp((sv^2 a^2 w^2 - 2 a v w - v^2 t) / (2 (1 + sv^2 t)))
      * (pi / a^2) sum_k k sin(k pi w) exp(-k^2 pi^2 t / (2 a^2))
The upper boundary is the lower one with w -> 1 - w and v -> -v.
Large-time series, 400 terms, 50 digits; gradients by mpmath.diff.
Gradient order: y, a, t0, w, v, sv (as wiener_lpdf).
"""
import sys

from mpmath import diff, exp, log, mp, mpf, pi, sin, sqrt

mp.dps = 50


def log_upper(t, a, w, v, sv, n_terms=400):
  w = 1 - w
  v = -v
  s = sum(k * sin(k * pi * w) * exp(-k**2 * pi**2 * t / (2 * a**2))
          for k in range(1, n_terms))
  q = 1 + sv**2 * t
  return (-log(q) / 2 + (sv**2 * a**2 * w**2 - 2 * a * v * w - v**2 * t)
          / (2 * q) + log(pi / a**2 * s))


def ref(t, a, w, v, sv):
  x = [mpf(t), mpf(a), mpf(w), mpf(v), mpf(sv)]
  val = log_upper(*x)
  g = []
  for k in range(5):
    def f(z, k=k):
      y = list(x)
      y[k] = z
      return log_upper(*y)
    g.append(diff(f, x[k]))
  # y, a, t0, w, v, sv
  return [val, g[0], g[1], -g[0], g[2], g[3], g[4]]


def main():
  rows = {}
  for line in open(sys.argv[1]):
    f = line.rstrip("\n").split("\t")
    if f[0] in ("PUB", "PUB8", "REL"):
      rows.setdefault(f[1], {})[f[0]] = [float(z) for z in f[2:]]
  names = ["value", "y", "a", "t0", "w", "v", "sv"]
  print("%-6s %-5s %-9s %-9s %-4s" % ("id", "tag", "|dvalue|", "max grel", "at"))
  for pid, d in rows.items():
    p = d["PUB"][:5]
    r = ref(*p)
    for tag in ("PUB", "PUB8", "REL"):
      if tag not in d:
        continue
      z = d[tag][5:]
      dv = abs(z[0] - float(r[0]))
      ge, gn = 0.0, "-"
      for k in range(1, 7):
        rk = float(r[k])
        e = abs(z[k] - rk) / max(abs(rk), 1.0)
        if e > ge:
          ge, gn = e, names[k]
      print("%-6s %-5s %-9.1e %-9.1e %-4s" % (pid, tag, dv, ge, gn))
    print("%-6s ref   value %.15g grads %s" % (
        pid, float(r[0]), " ".join("%.10g" % float(q) for q in r[1:])))


if __name__ == "__main__":
  main()
