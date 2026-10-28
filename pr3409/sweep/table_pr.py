#!/usr/bin/env python3
"""The accuracy table of the PR text, from two probe_cmp_beta outputs.

Usage: python3 table_pr.py <develop output> <PR output> [refs file]

Both outputs come from `probe_cmp_beta <points> double`.  The reference is
refs_beta_18000.txt (I, d/da, d/db per point).

Per region and path (A = inc_beta_dda / inc_beta_ddb, B = the 9-argument
grad_reg_inc_beta), a point counts as
  non-finite  if either partial is NaN, infinite, or threw;
  wrong       if it is finite and either partial has err_log > 1e-2,
where err_log = |v - r| / I is the absolute error of d/dtheta log I.
Points where the double I is 0 are left out.

For the PR tree: the worst err_log and the worst relative error
(where |r| / I >= 1e-3) over both partials and both paths, with the points
where I is subnormal left out, and the number of those points.
"""

import math
import sys

DBL_MIN = 2.2250738585072014e-308


def parse(path):
  rows = []
  region = "?"
  for line in open(path):
    line = line.strip()
    if not line:
      continue
    if line.startswith("#"):
      if not line.startswith(("# probe", "# evaluated")):
        region = line[2:]
      continue
    f = line.split()
    if not f[0].startswith(("0x", "-0x")):
      continue
    key = tuple(float.fromhex(x) for x in f[:3])
    vals = [float.fromhex(x) for x in f[3:7]]
    rows.append((region, key, vals, f[7]))
  return rows


def load_refs(path):
  refs = {}
  for line in open(path):
    f = line.split()
    refs[tuple(float.fromhex(x) for x in f[:3])] = tuple(
        float.fromhex(x) for x in f[3:6])
  return refs


def bad(vals, status, cols):
  return any(status[c] == "T" or not math.isfinite(vals[c]) for c in cols)


def main():
  dev = parse(sys.argv[1])
  pr = parse(sys.argv[2])
  refs = load_refs(sys.argv[3] if len(sys.argv) > 3 else "refs_beta_18000.txt")
  assert [r[1] for r in dev] == [r[1] for r in pr], "point files differ"
  regions = []
  for r in dev:
    if r[0] not in regions:
      regions.append(r[0])
  print(f"{'region':<26} {'n':>5} | {'A nf':>5} {'A >1e-2':>7} "
        f"{'B nf':>5} {'B >1e-2':>7} | {'PR nf':>5} {'worst':>8} "
        f"{'rel*':>8} {'sub':>4}")
  tot = [0] * 5
  for reg in regions:
    n = 0
    cnt = [0, 0, 0, 0]  # A nf, A wrong, B nf, B wrong
    pr_nf = 0
    worst = 0.0
    relw = 0.0
    sub = 0
    for (rg, key, dv, ds), (_, _, pv, ps) in zip(dev, pr):
      if rg != reg:
        continue
      i, rda, rdb = refs[key]
      if math.isnan(rda) or i <= 0.0:
        continue
      n += 1
      ref = (rda, rdb, rda, rdb)
      for p, cols in ((0, (0, 1)), (1, (2, 3))):
        if bad(dv, ds, cols):
          cnt[2 * p] += 1
        elif any(abs(dv[c] - ref[c]) / i > 1e-2 for c in cols):
          cnt[2 * p + 1] += 1
      if bad(pv, ps, range(4)):
        pr_nf += 1
        continue
      if i < DBL_MIN:
        sub += 1
        continue
      for c in range(4):
        worst = max(worst, abs(pv[c] - ref[c]) / i)
        if abs(ref[c]) / i >= 1e-3 and abs(ref[c]) >= DBL_MIN:
          relw = max(relw, abs(pv[c] - ref[c]) / abs(ref[c]))
    for k in range(4):
      tot[k] += cnt[k]
    tot[4] += pr_nf
    print(f"{reg[:26]:<26} {n:>5d} | {cnt[0]:>5d} {cnt[1]:>7d} {cnt[2]:>5d} "
          f"{cnt[3]:>7d} | {pr_nf:>5d} {worst:>8.1e} {relw:>8.1e} {sub:>4d}")
  print(f"{'total':<26} {'':>5} | {tot[0]:>5d} {tot[1]:>7d} {tot[2]:>5d} "
        f"{tot[3]:>7d} | {tot[4]:>5d}")
  print("columns: nf = non-finite or threw; >1e-2 = finite with err_log > "
        "1e-2; worst = worst err_log; rel* = worst relative error where "
        "|r| / I >= 1e-3; sub = points with subnormal I, left out of worst "
        "and rel*")
  return 0


if __name__ == "__main__":
  sys.exit(main())
