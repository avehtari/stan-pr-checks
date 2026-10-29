#!/usr/bin/env python3
"""
analyze_cdf_fd.py -- gradients of the 7-parameter lcdf and lccdf against
finite differences of the references of probe_ref7.cpp (REF_CDF_FD, REF_P5).
Claims of stan-dev/math#3442: the 7-parameter lcdf and lccdf gradient errors.

Usage:
  analyze_cdf_fd.py d REF_FILES... -- DIR TREE[@PREC]...   (D points;
      probe_cdfacc_<tree>.out, rows at precision PREC, default 1e-08)
  analyze_cdf_fd.py p5 REF_FILE -- DIR TREE...            (p5 points;
      probe_w4_<tree>.out, rows LCDF/LCCDF of set p5)

Error per point: largest |g - ref| / max(|ref|, 1) over the gradients.
The FD noise is estimated by |FD d/dy + FD d/dt0| / max(|FD d/dy|, 1),
which is 0 for exact derivatives. Points whose value differs from the
reference by more than 1e-2 are counted separately ("far").
"""
import math
import os
import sys

NAMES = ["y", "a", "t0", "w", "v", "sv", "sw", "st0"]


def pct(xs, q):
  xs = sorted(xs)
  return xs[min(len(xs) - 1, int(q * len(xs)))] if xs else float("nan")


def load_refs(files):
  refs = {}
  for path in files:
    for line in open(path):
      f = line.rstrip("\n").split("\t")
      if f[0] != "RFD":
        continue
      g = [float(x) for x in f[4:20]]
      refs[("LCDF", f[1])] = (float(f[2]), g[:8])
      refs[("LCCDF", f[1])] = (float(f[3]), g[8:])
  return refs


def main():
  mode = sys.argv[1]
  sep = sys.argv.index("--")
  refs = load_refs(sys.argv[2:sep])
  d = sys.argv[sep + 1]
  trees = sys.argv[sep + 2:]
  noise = [abs(g[0] + g[2]) / max(abs(g[0]), 1.0) for _, g in refs.values()]
  print("references: %d rows; FD noise max %.1e, p99 %.1e" %
        (len(refs), max(noise), pct(noise, 0.99)))
  print("%-6s %-14s %-4s %-8s %-8s %-3s %-4s %-8s %-8s" % (
      "fn", "tree", "n", "g max", "g p99", "at", "far", "v max", "v p99"))
  for fn in ("LCDF", "LCCDF"):
    for spec in trees:
      tree, prec = (spec.split("@") + ["1e-08"])[:2]
      ge, ve, far, worst = [], [], 0, (0.0, "-")
      if mode == "d":
        path = os.path.join(d, "probe_cdfacc_%s.out" % tree)
      else:
        path = os.path.join(d, "probe_w4_%s.out" % tree)
      for line in open(path):
        f = line.rstrip("\n").split("\t")
        if mode == "d":
          if f[0] != fn or f[10] != prec or f[11] == "THROW":
            continue
          pid, z = f[1], [float(x) for x in f[11:20]]
        else:
          if f[0] != fn or f[1] != "p5" or f[9] == "THROW":
            continue
          pid, z = f[2], [float(x) for x in f[9:16]]
        r = refs.get((fn, pid))
        if r is None:
          continue
        if not math.isfinite(z[0]) or abs(z[0] - r[0]) > 1e-2:
          far += 1
          continue
        ve.append(abs(z[0] - r[0]))
        e = 0.0
        for k in range(len(z) - 1):
          if math.isnan(r[1][k]):
            continue
          ek = abs(z[k + 1] - r[1][k]) / max(abs(r[1][k]), 1.0)
          if not math.isfinite(ek):
            ek = math.inf
          if ek > e:
            e = ek
            if ek > worst[0]:
              worst = (ek, NAMES[k])
        ge.append(e)
      print("%-6s %-14s %-4d %-8.1e %-8.1e %-3s %-4d %-8.1e %-8.1e" % (
          fn, spec, len(ge), max(ge), pct(ge, 0.99), worst[1], far,
          max(ve), pct(ve, 0.99)))


if __name__ == "__main__":
  main()
