#!/usr/bin/env python3
"""
analyze_tol.py -- accuracy and time of the 7-parameter wiener_lpdf at the 300
D points against the references of probe_ref7.cpp (REF_D=1).  Claim of
stan-dev/math#3442: "7-parameter lpdf, ordinary" in the accuracy table.

Usage: analyze_tol.py REF_FILE DACC_DIR TREE...

DACC_DIR holds probe_dacc_<tree>.out. Errors: value absolute; gradients
relative to max(|ref|, 1). Gradient order: y a t0 w v sv sw st0.
Columns: max and 99th percentile of the value error and of the largest
gradient error per point, and the summed median time per call over the
300 points (ms; the time of one pass over all points).
"""
import os
import sys

NAMES = ["y", "a", "t0", "w", "v", "sv", "sw", "st0"]


def pct(xs, q):
  xs = sorted(xs)
  return xs[min(len(xs) - 1, int(q * len(xs)))]


def main():
  refs = {}
  for line in open(sys.argv[1]):
    f = line.rstrip("\n").split("\t")
    refs[f[1]] = [float(x) for x in f[10:19]]
  print("max ref error estimate: %s" % max(
      float(line.rstrip("\n").split("\t")[21]) for line in open(sys.argv[1])))
  d = sys.argv[2]
  print("%-8s %-6s %-9s %-9s %-9s %-9s %-4s %-8s" % (
      "tree", "prec", "v max", "v p99", "g max", "g p99", "at", "ms/pass"))
  for t in sys.argv[3:]:
    rows = {}
    for line in open(os.path.join(d, "probe_dacc_%s.out" % t)):
      f = line.rstrip("\n").split("\t")
      rows.setdefault(f[10], []).append(
          (f[1], [float(x) for x in f[11:20]], float(f[20])))
    for prec, lst in rows.items():
      ve, ge, tm = [], [], 0.0
      worst = (0.0, "-")
      for pid, z, us in lst:
        r = refs[pid]
        ve.append(abs(z[0] - r[0]))
        g = 0.0
        for k in range(1, 9):
          e = abs(z[k] - r[k]) / max(abs(r[k]), 1.0)
          if e > g:
            g = e
            if e > worst[0]:
              worst = (e, NAMES[k - 1])
        ge.append(g)
        tm += us
      print("%-8s %-6s %-9.1e %-9.1e %-9.1e %-9.1e %-4s %-8.2f" % (
          t, prec, max(ve), pct(ve, 0.99), max(ge), pct(ge, 0.99), worst[1],
          tm / 1000))


if __name__ == "__main__":
  main()
