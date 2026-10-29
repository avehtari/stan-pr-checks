#!/usr/bin/env python3
# analyze_oldw4.py -- score probe_oldw4.cpp (the older 4-parameter wiener_lpdf
# and the 5-parameter lpdf with sv = 0) against ref_w4_mpmath.py.
# Claim: "returns -inf at short times where the log density is below about
# -708 (35 of 600 grid points); where it is finite, it is accurate".
# Usage, from pr3442/:
#   python3 p45/analyze_oldw4.py results/probe_oldw4_head3442.out \
#       results/ref_oldw4.out
import collections
import math
import sys


def main():
  probe, refs = sys.argv[1], sys.argv[2]
  ref = {}
  for line in open(refs):
    f = line.rstrip("\n").split("\t")
    ref[f[1]] = [float(x) for x in f[2:8]]
  for fn in ("OLD", "LPDF"):
    bad, ev, eg, min_fin = [], 0.0, 0.0, 0.0
    for line in open(probe):
      f = line.rstrip("\n").split("\t")
      if f[0] != fn:
        continue
      p = [float(x) for x in f[3:9]]
      r = ref[f[2]]
      if f[9] == "THROW" or not math.isfinite(float(f[9])):
        bad.append((p[0] - p[2], p[1], r[0]))
        continue
      o = [float(x) for x in f[9:15]]
      ev = max(ev, abs(o[0] - r[0]))
      eg = max(eg, max(abs(a - b) / max(1.0, abs(b))
                       for a, b in zip(o[1:6], r[1:6])))
      min_fin = min(min_fin, r[0])
    print("%-4s non-finite %d; finite: value %.2e, gradients %.2e, "
          "lowest log density %.4g" % (fn, len(bad), ev, eg, min_fin))
    if bad:
      cnt = collections.Counter((round(t, 6), a) for t, a, _ in bad)
      print("     non-finite at (t, a): %s" % sorted(cnt.items()))
      print("     true log density there: %.4g to %.4g"
            % (max(b[2] for b in bad), min(b[2] for b in bad)))


if __name__ == "__main__":
  main()
