#!/usr/bin/env python3
"""Summarize RESULT lines of run_branch_bench.sh logs.

Usage: summarize.py LOG

Per (case, propto, N): ns per element of develop and branch (median over
rounds of the per-round medians), the ratio branch / develop of these
medians, and the smallest and largest per-round ratio.
"""
import statistics
import sys
from collections import defaultdict

t = defaultdict(dict)  # (case, propto, N) -> {(tree, round): ns}
for line in open(sys.argv[1]):
    f = line.split()
    if len(f) < 9 or f[0] != "RESULT" or f[2] == "MISSING":
        continue
    rnd, tree, case, propto, n, med = f[1], f[2], f[3], int(f[4]), int(f[5]), float(f[6])
    t[(case, propto, n)][(tree, rnd)] = med

print(f"{'case':4} {'propto':6} {'N':>6} {'dev ns/el':>9} {'br ns/el':>9} "
      f"{'ratio':>6} {'min':>6} {'max':>6}")
for key in sorted(t, key=lambda k: (k[0], -k[1], k[2])):
    case, propto, n = key
    d = t[key]
    rounds = sorted({r for (_, r) in d})
    dev = [d[("dev", r)] for r in rounds if ("dev", r) in d]
    br = [d[("br", r)] for r in rounds if ("br", r) in d]
    ratios = [d[("br", r)] / d[("dev", r)] for r in rounds
              if ("br", r) in d and ("dev", r) in d]
    md, mb = statistics.median(dev), statistics.median(br)
    print(f"{case:4} {propto:6d} {n:6d} {md / n:9.2f} {mb / n:9.2f} "
          f"{mb / md:6.2f} {min(ratios):6.2f} {max(ratios):6.2f}")
