#!/usr/bin/env python3
"""Count the outcomes of beta_neg_binomial_cdf on the grid (claim S10).

Usage: count_bnb.py <probe_bnb output, develop> <probe_bnb output, PR>
The outputs come from `probe_bnb <points> <limit> 1` (cdf only).
Outcome per point: finite (finished, value and gradients finite), nonfinite
(finished with a NaN or inf), threw, limit (time limit), died.
Prints the counts per tree, the transitions between the trees, and the
points where develop is finite and the PR is not.
"""

import math
import sys
from collections import Counter


def load(path):
  out = {}
  for line in open(path):
    left, right = line.split(" => ")
    t = right.split()
    vals = [float.fromhex(x) if "nan" not in x else math.nan for x in t[1:5]]
    status = t[5]
    if status == ".":
      kind = "finite" if all(math.isfinite(v) for v in vals) else "nonfinite"
    else:
      kind = {"T": "threw", "X": "limit", "C": "died"}[status]
    key = tuple(float.fromhex(h) for h in left.split()[2:])
    out[key] = (kind, vals, float(t[6]))
  return out


def main():
  dev, pr = load(sys.argv[1]), load(sys.argv[2])
  for name, d in (("develop", dev), ("PR", pr)):
    c = Counter(v[0] for v in d.values())
    print(f"{name:8s} {len(d)} points: " + ", ".join(
        f"{k} {c[k]}" for k in ("finite", "nonfinite", "threw", "limit",
                                "died")))
  trans = Counter((dev[k][0], pr[k][0]) for k in dev)
  print("\ntransitions develop -> PR (changed only):")
  for (a, b), n in sorted(trans.items()):
    if a != b:
      print(f"  {a:9s} -> {b:9s} {n}")
  print("\ndevelop finite, PR not finite (n, r, alpha, beta: develop value,"
        " PR value and gradients):")
  for k in dev:
    if dev[k][0] == "finite" and pr[k][0] != "finite":
      args = ", ".join(f"{x:g}" for x in k)
      print(f"  ({args}): {dev[k][1][0]:.3g}; {pr[k][0]} "
            + " ".join(f"{v:.3g}" for v in pr[k][1]))
  return 0


if __name__ == "__main__":
  sys.exit(main())
