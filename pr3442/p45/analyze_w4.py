#!/usr/bin/env python3
"""
analyze_w4.py -- accuracy and time of the 4- and 5-parameter Wiener functions
(probe_w4.cpp) against the mpmath references (refs/ref_w4_mpmath.py).  Claims
of stan-dev/math#3442: the 4- and 5-parameter rows of the accuracy table.

Usage: analyze_w4.py REF_FILE DIR TREE...

Per function, point set and tree: value error (absolute, log scale) and
largest gradient error per point (relative to max(|ref|, 1); gradients
y a t0 w v sv), max and 99th percentile, number of points with value error
> 1e-6 (inf: of these, the points with a value of -inf), and the sum of the per-point median times (ms per pass over the
300 points). OLD (the older 4-parameter wiener_lpdf) is compared with the
LPDF reference of the same point.
"""
import math
import os
import sys

NAMES = ["y", "a", "t0", "w", "v", "sv"]


def pct(xs, q):
  xs = sorted(xs)
  return xs[min(len(xs) - 1, int(q * len(xs)))] if xs else float("nan")


def main():
  ref = {}
  for line in open(sys.argv[1]):
    f = line.rstrip("\n").split("\t")
    ref[(f[0], f[1])] = [float(x) for x in f[2:9]]
  d = sys.argv[2]
  print("%-6s %-3s %-9s %-8s %-8s %-6s %-4s %-8s %-8s %-3s %-5s %-8s" % (
      "fn", "set", "tree", "v max", "v p99", "n>1e-6", "inf", "g max", "g p99",
      "at", "throw", "ms/pass"))
  for fn, st in (("OLD", "p4"), ("LPDF", "p4"), ("LPDF", "p5"),
                 ("LCDF", "p4"), ("LCCDF", "p4"), ("LCDF", "p5"),
                 ("LCCDF", "p5")):
    for t in sys.argv[3:]:
      ve, ge, tm, thr, worst = [], [], 0.0, 0, (0.0, "-")
      for line in open(os.path.join(d, "probe_w4_%s.out" % t)):
        f = line.rstrip("\n").split("\t")
        if f[0] != fn or f[1] != st:
          continue
        if f[9] == "THROW":
          thr += 1
          continue
        z = [float(x) for x in f[9:16]]
        tm += float(f[16])
        r = ref.get(("LPDF" if fn == "OLD" else fn, f[2]))
        if r is None:
          continue
        ve.append(abs(z[0] - r[0]))
        if not math.isnan(r[1]):
          g = 0.0
          for k in range(6):
            if fn == "OLD" and k == 5:
              continue
            e = abs(z[k + 1] - r[k + 1]) / max(abs(r[k + 1]), 1.0)
            if e > g:
              g = e
              if e > worst[0]:
                worst = (e, NAMES[k])
          ge.append(g)
      if not ve:
        continue
      print("%-6s %-3s %-9s %-8.1e %-8.1e %-6d %-4d %-8s %-8s %-3s %-5d %-8.2f" % (
          fn, st, t, max(ve), pct(ve, 0.99), sum(1 for e in ve if e > 1e-6),
          sum(1 for e in ve if math.isinf(e)),
          "%.1e" % max(ge) if ge else "-", "%.1e" % pct(ge, 0.99) if ge else "-",
          worst[1], thr, tm / 1000))


if __name__ == "__main__":
  main()
