#!/usr/bin/env python3
"""Value errors of beta_binomial_lpmf with and without the log_beta_ratio
selection (claim P3).

Usage (from pr3434/):
  python3 lbr_select/score_select.py out/bb_dev.txt out/bb_pr.txt \
      out/bb_stirling.txt

The three files are probe_bb outputs on the same points: develop, the PR
head, and the PR head with the Stirling form forced
(make_stirling_tree.sh). References: refs/refs_bb.txt (route A of ref_bb.py;
check_bb_routes.py checks it against routes B and C).
Prints per region and per N the worst absolute value error of each tree,
and the points where the Stirling-only tree is worse than both develop and
the PR by more than a factor of 4 and by more than 1e-13.
"""

import os
import sys
from collections import OrderedDict

from mpmath import mp, mpf

mp.dps = 60
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def load_refs():
  refs = {}
  for line in open(os.path.join(ROOT, "refs", "refs_bb.txt")):
    p = line.split()
    key = (int(p[0]), int(p[1]), float.fromhex(p[2]), float.fromhex(p[3]))
    refs[key] = mpf(p[5])
  return refs


def load(path, refs):
  out = OrderedDict()
  for line in open(path):
    if line.startswith("#") or not line.strip():
      continue
    left, right = line.split("=>")
    tag, n, N, a, b = left.split()
    key = (int(n), int(N), float.fromhex(a), float.fromhex(b))
    v = right.split()[1]
    try:
      err = abs(mpf(float.fromhex(v)) - refs[key])
    except ValueError:
      err = mpf("inf")
    out[(tag, key)] = err
  return out


def fmt(x):
  return f"{float(x):.1e}"


def main():
  if len(sys.argv) != 4:
    print(__doc__)
    return 2
  refs = load_refs()
  dev, pr, st = (load(p, refs) for p in sys.argv[1:4])
  groups = OrderedDict()
  for k in dev:
    tag, (n, N, a, b) = k
    groups.setdefault((tag, N), []).append(k)
  print(f"{'region':6s} {'N':>8s} {'pts':>5s} {'develop':>8s} {'PR':>8s} "
        f"{'Stirling':>8s} {'#worse':>6s}")
  bad = []
  for (tag, N), keys in groups.items():
    worse = 0
    for k in keys:
      lim = 4 * max(dev[k], pr[k])
      if st[k] > lim and st[k] > mpf("1e-13"):
        worse += 1
        bad.append((st[k], k))
    print(f"{tag:6s} {N:8d} {len(keys):5d} "
          f"{fmt(max(dev[k] for k in keys)):>8s} "
          f"{fmt(max(pr[k] for k in keys)):>8s} "
          f"{fmt(max(st[k] for k in keys)):>8s} {worse:6d}")
  bad.sort(reverse=True)
  print(f"\n{len(bad)} points where the Stirling-only tree is worse than "
        "4 x max(develop, PR) and above 1e-13; the worst 10:")
  for e, (tag, (n, N, a, b)) in bad[:10]:
    k = (tag, (n, N, a, b))
    print(f"  {tag:5s} n {n:7d} N {N:7d} alpha {a:10.4g} beta {b:10.4g}  "
          f"Stirling {fmt(e)}  develop {fmt(dev[k])}  PR {fmt(pr[k])}")
  return 0


if __name__ == "__main__":
  sys.exit(main())
