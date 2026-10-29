#!/usr/bin/env python3
"""
analyze_swst0.py -- d/dsw and d/dst0 of the 7-parameter wiener_lpdf at small
sw and st0 against the score-form references (probe_ref7.cpp REF_SCORE).
Claims of stan-dev/math#3442: d/dsw 1.1 -> 3.3e-9; end-point error 3e-10/sw.

Usage: analyze_swst0.py REF_FILE[,REF_FILE...] DIR TREE...

Error = |g - ref| / max(|ref|, 1).  Rows: series (w: sw varies; t: st0
varies; a: the other variability 0, b: non-zero), size index j (0..4 =
1e-8, 1e-6, 1e-4, 1e-2, 0.1 or 0.5 d), tree and precision; columns: max
error of the varied gradient, max error of all other gradients, max value
error, and the summed time.  First: the reference check |score - endpoint|
for the large sizes (j = 3, 4), where the endpoint form is well conditioned.
"""
import math
import os
import re
import sys

SIZES = ["1e-8", "1e-6", "1e-4", "1e-2", "big"]


def load_refs(spec):
  refs = {}
  for path in spec.split(","):
    for line in open(path):
      f = line.rstrip("\n").split("\t")
      if f[0] != "REF" or len(f) < 24 or f[10] == "FAIL":
        continue
      refs[f[1]] = ([float(x) for x in f[10:19]], float(f[21]),
                    float(f[22]), float(f[23]))
  return refs


def rel(x, r):
  if not math.isfinite(x):
    return math.inf
  return abs(x - r) / max(abs(r), 1.0)


def main():
  refs = load_refs(sys.argv[1])
  d = sys.argv[2]
  trees = sys.argv[3:]
  print("references: %d" % len(refs))
  chk = {}
  for pid, (r, err, g6e, g7e) in refs.items():
    m = re.match(r"k\d+([wt])(\d)([ab])", pid)
    s, j = m.group(1), int(m.group(2))
    k = 6 if s == "w" else 7
    end = g6e if s == "w" else g7e
    key = (s, j)
    chk[key] = max(chk.get(key, 0.0), rel(end, r[k + 1]))
  print("reference check |score - endpoint| / max(|ref|, 1) by series, size:")
  for key in sorted(chk):
    print("  %s %-5s %.1e" % (key[0], SIZES[key[1]], chk[key]))
  print("%-3s %-5s %-6s %-6s %-8s %-8s %-8s %-8s %-4s" % (
      "ser", "size", "tree", "prec", "g var", "g other", "value", "ms", "n"))
  for series in ("wa", "wb", "ta", "tb"):
    for j in range(5):
      for t in trees:
        stats = {}
        for line in open(os.path.join(d, "probe_swst0_%s.out" % t)):
          f = line.rstrip("\n").split("\t")
          m = re.match(r"k\d+([wt])(\d)([ab])", f[1])
          if m.group(1) + m.group(3) != series or int(m.group(2)) != j:
            continue
          if f[1] not in refs:
            continue
          prec = f[10]
          st = stats.setdefault(prec, [0.0, 0.0, 0.0, 0.0, 0])
          r = refs[f[1]][0]
          if f[11] == "THROW":
            st[0] = st[1] = st[2] = math.inf
            st[4] += 1
            continue
          z = [float(x) for x in f[11:20]]
          k = 6 if series[0] == "w" else 7
          st[0] = max(st[0], rel(z[k + 1], r[k + 1]))
          for m2 in range(8):
            if m2 != k:
              st[1] = max(st[1], rel(z[m2 + 1], r[m2 + 1]))
          st[2] = max(st[2], abs(z[0] - r[0]) if math.isfinite(z[0])
                      else math.inf)
          st[3] += float(f[20]) / 1000
          st[4] += 1
        for prec in sorted(stats, key=float, reverse=True):
          st = stats[prec]
          print("%-3s %-5s %-6s %-6s %-8.1e %-8.1e %-8.1e %-8.2f %-4d" % (
              series, SIZES[j], t, prec, st[0], st[1], st[2], st[3], st[4]))


if __name__ == "__main__":
  main()
