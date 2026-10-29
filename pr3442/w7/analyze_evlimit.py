#!/usr/bin/env python3
"""
analyze_evlimit.py -- error and time of the 7-parameter lcdf and lccdf at the
D points, by number of integration dimensions (sv, sw, st0 non-zero).  Claims
of stan-dev/math#3442: the rows "1-2 integration dims"; 3 dims, 6000, 4.3e-5.

Usage: analyze_evlimit.py REF_FILE[,REF_FILE...] DIR TREE[@PREC]...

REF_FILE rows: probe_ref7.cpp REF_CDF ("RCDF id lcdf lccdf", values only)
or REF_CDF_FD ("RFD id lcdf lccdf <8 lcdf grads> <8 lccdf grads>").
Value error is absolute (log scale).  Gradient error per point is the
largest |g - ref| / max(|ref|, 1) (only with RFD rows; nan reference
gradients are skipped).  Time = sum of per-point median times (ms).
"""
import math
import os
import sys


def pct(xs, q):
  xs = sorted(xs)
  return xs[min(len(xs) - 1, int(q * len(xs)))] if xs else float("nan")


def load_refs(spec):
  ref = {}
  for path in spec.split(","):
    for line in open(path):
      f = line.rstrip("\n").split("\t")
      if f[0] not in ("RCDF", "RFD"):
        continue
      g = [float(x) for x in f[4:20]] if f[0] == "RFD" else None
      ref[("LCDF", f[1])] = (float(f[2]), g[:8] if g else None)
      ref[("LCCDF", f[1])] = (float(f[3]), g[8:] if g else None)
  return ref


def main():
  ref = load_refs(sys.argv[1])
  print("%-6s %-12s %-3s %-4s %-8s %-8s %-8s %-8s %-9s" % (
      "fn", "tree", "dim", "n", "v max", "v p99", "g max", "g p99",
      "ms/pass"))
  for fn in ("LCDF", "LCCDF"):
    for spec in sys.argv[3:]:
      tree, prec = (spec.split("@") + ["1e-08"])[:2]
      by_dim = {}
      for line in open(os.path.join(sys.argv[2], "probe_cdfacc_%s.out" % tree)):
        f = line.rstrip("\n").split("\t")
        if f[0] != fn or f[10] != prec or f[11] == "THROW":
          continue
        r = ref.get((fn, f[1]))
        if r is None:
          continue
        p = [float(x) for x in f[2:10]]
        z = [float(x) for x in f[11:20]]
        dim = (p[5] != 0) + (p[6] != 0) + (p[7] != 0)
        e = abs(z[0] - r[0]) if math.isfinite(z[0]) else math.inf
        ge = float("nan")
        if r[1] is not None:
          ge = 0.0
          for k in range(8):
            if math.isnan(r[1][k]):
              continue
            ek = abs(z[k + 1] - r[1][k]) / max(abs(r[1][k]), 1.0)
            ge = max(ge, ek if math.isfinite(ek) else math.inf)
        by_dim.setdefault(dim, []).append((e, ge, float(f[20])))
      for dim in sorted(by_dim):
        es = [x[0] for x in by_dim[dim]]
        gs = [x[1] for x in by_dim[dim] if not math.isnan(x[1])]
        print("%-6s %-12s %-3d %-4d %-8.1e %-8.1e %-8.1e %-8.1e %-9.1f" % (
            fn, spec, dim, len(es), max(es), pct(es, 0.99),
            max(gs) if gs else float("nan"), pct(gs, 0.99),
            sum(x[2] for x in by_dim[dim]) / 1000))


if __name__ == "__main__":
  main()
