#!/usr/bin/env python3
"""Check every reference of refs_bb.txt against two more routes.

refs_bb.txt holds route A of ref_bb.py (80 digits: loggamma value,
closed-form digamma differences). At each point this script computes
  route B  the value at 130 digits and mp.diff of the value at 130 digits
  route C  the gradients as exact finite sums (N <= 2000 only)
and prints the worst disagreement per region of points_bb.txt.

Usage: check_bb_routes.py <points_bb.txt> <refs_bb.txt> [processes] [limit]
`limit` checks only the first points of each region (for a quick run).
"""

import multiprocessing
import sys

from mpmath import mp, mpf

import ref_bb as R


def load(path):
  refs = {}
  for line in open(path):
    p = line.split()
    key = (int(p[0]), int(p[1]), float.fromhex(p[2]), float.fromhex(p[3]))
    refs[key] = p[4:]
  return refs


def rel(a, b):
  s = max(abs(a), abs(b))
  return abs(a - b) / s if s > mpf("1e-60") else abs(a - b)


def one(item):
  region, key, ref = item
  n, N, a, b = key
  mp.dps = R.DPS
  va, da, db = (mpf(x) for x in ref[1:4])
  with mp.workdps(R.DPS_CHECK):
    vb = R.lp_bb(n, N, a, b)
  gb = R.grad_bb_diff(n, N, a, b)
  mp.dps = 60
  dv = abs(va - vb) / (1 + abs(va))
  dg = max(rel(da, gb[0]), rel(db, gb[1]))
  dc = None
  if N <= R.NMAX_SUM:
    mp.dps = R.DPS
    gc = R.grad_bb_sum(n, N, a, b)
    mp.dps = 60
    dc = max(rel(da, gc[0]), rel(db, gc[1]))
  return region, float(dv), float(dg), None if dc is None else float(dc)


def main():
  if len(sys.argv) < 3:
    print(__doc__)
    return 2
  procs = int(sys.argv[3]) if len(sys.argv) > 3 else 8
  limit = int(sys.argv[4]) if len(sys.argv) > 4 else None
  refs = load(sys.argv[2])
  items = []
  per = {}
  region = None
  for line in open(sys.argv[1]):
    if line.startswith("#") or not line.strip():
      continue
    tag, n, N, a, b = line.split()
    key = (int(n), int(N), float.fromhex(a), float.fromhex(b))
    per[tag] = per.get(tag, 0) + 1
    if limit is not None and per[tag] > limit:
      continue
    items.append((tag, key, refs[key]))
  with multiprocessing.Pool(procs) as pool:
    res = pool.map(one, items, chunksize=16)
  worst = {}
  for region, dv, dg, dc in res:
    w = worst.setdefault(region, [0, 0.0, 0.0, None])
    w[0] += 1
    w[1] = max(w[1], dv)
    w[2] = max(w[2], dg)
    if dc is not None:
      w[3] = dc if w[3] is None else max(w[3], dc)
  print("route A (80 digits) against B (130 digits, mp.diff) and C (sums)")
  print(f"{'region':8s} {'points':>6s} {'value A-B':>10s} {'grad A-B':>10s} "
        f"{'grad A-C':>10s}")
  for region, w in worst.items():
    c = "-" if w[3] is None else f"{w[3]:.1e}"
    print(f"{region:8s} {w[0]:6d} {w[1]:10.1e} {w[2]:10.1e} {c:>10s}")
  return 0


if __name__ == "__main__":
  sys.exit(main())
