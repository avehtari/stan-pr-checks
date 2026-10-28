#!/usr/bin/env python3
# Bins the error of one tree near z = a by decade of a (cached references).
# Supports "3e-11 for a < 1 000, 2.9e-7 for a < 10 000, 5.2e-6 at a = z = 1e4".
"""Worst and median relative error of grad_reg_inc_gamma near z = a.

Usage: python3 bin_near_diag.py PROBE_OUTPUT [CACHE]

PROBE_OUTPUT is probe_cmp_roots output in mode double.  CACHE defaults to
refs_cache_49990.txt next to this script.  Rows: region C4 ("z near a")
binned by decade of a, and the points of region C10 ("corners") with a
and z both above 9 000.  Points with a subnormal or zero reference are
left out, as in score_cmp_roots.py.
"""
import math
import os
import sys

DBL_MIN = 2.2250738585072014e-308
HERE = os.path.dirname(os.path.abspath(__file__))


def main():
  path = sys.argv[1]
  cache = (sys.argv[2] if len(sys.argv) > 2
           else os.path.join(HERE, "refs_cache_49990.txt"))
  refs = {}
  for line in open(cache):
    a, z, r = line.split()
    refs[(float.fromhex(a), float.fromhex(z))] = float.fromhex(r)
  bins = {}
  region = None
  for line in open(path):
    if line.startswith("# probe") or line.startswith("# evaluated"):
      continue
    if line.startswith("#"):
      region = line[2:].strip()
      continue
    t = line.split()
    a, z, up = (float.fromhex(x) for x in t[:3])
    r = refs[(a, z)]
    if math.isnan(r) or abs(r) < DBL_MIN:
      continue
    e = abs(up - r) / abs(r) if math.isfinite(up) else math.inf
    if region.startswith("C4 "):
      lo = 10 ** math.floor(math.log10(a))
      key = f"C4  [{lo:g}, {10 * lo:g})"
    elif region.startswith("C10 ") and a > 9000 and z > 9000:
      key = "C10 a, z > 9 000"
    else:
      continue
    bins.setdefault(key, []).append(e)
  print(f"{'bin':<24} {'n':>5} {'worst':>9} {'median':>9}")
  for key in sorted(bins, key=lambda k: (k[:3], float(k.split('[')[-1].split(',')[0])
                                         if '[' in k else 0.0)):
    v = sorted(bins[key])
    print(f"{key:<24} {len(v):>5} {v[-1]:>9.2e} {v[len(v) // 2]:>9.2e}")


if __name__ == "__main__":
  main()
