#!/usr/bin/env python3
# Writes points.txt, the 1 800-point sweep in nine regions (seed fixed).
# Supports the nine-region accuracy table of #3408 Part 1.
"""Generate the (a, z) sweep point file for the incomplete gamma roots.

Format, the same convention as erfcx_libcerf_work/bench/score_mpmath.py:

    # <region name>
    <a in hex> <z in hex>
    ...

Hex is exact in both directions, so the C++ reads exactly the points that
the Python scorer will use.

The regions follow the branch structure of develop `grad_reg_inc_gamma` and
`grad_reg_lower_inc_gamma`, not a uniform grid. A uniform grid hides a
branch that is wrong on a small set.
"""

import math
import random
import sys

SEED = 20260918
N_DEFAULT = 200


def log_uniform(rng, lo, hi):
  return math.exp(rng.uniform(math.log(lo), math.log(hi)))


def safe_z_max(a, hard_cap):
  """Keep exp(a * log z) representable, so the series branch cannot
  overflow. a * log z < 709 gives z < exp(709 / a)."""
  if a <= 0:
    return hard_cap
  limit = math.exp(min(709.0 / a, 700.0))
  return max(1.0001, min(hard_cap, limit))


def regions(n, rng):
  out = []

  # R1: series, a < 1
  pts = []
  for _ in range(n):
    a = log_uniform(rng, 0.01, 1.0)
    pts.append((a, log_uniform(rng, 1e-3, min(8.0, safe_z_max(a, 8.0)))))
  out.append(("R1 series, a<1, z<8", pts))

  # R2: series, moderate a, z < min(a, 8)
  pts = []
  for _ in range(n):
    a = log_uniform(rng, 1.0, 10.0)
    hi = min(a, 8.0)
    pts.append((a, log_uniform(rng, 0.1, max(0.11, hi))))
  out.append(("R2 series, 1<=a<=10, z<min(a,8)", pts))

  # R3: series, large a, z < a. Capped so exp(a log z) stays finite.
  pts = []
  for _ in range(n):
    a = log_uniform(rng, 10.0, 170.0)
    hi = min(a, safe_z_max(a, 60.0))
    pts.append((a, log_uniform(rng, 1.0, max(1.01, hi))))
  out.append(("R3 series, 10<a<=170, z<a", pts))

  # R4: asymptotic branch where it is not valid
  pts = []
  for _ in range(n):
    a = log_uniform(rng, 8.0, 170.0)
    pts.append((a, a * rng.uniform(1.0, 2.0)))
  out.append(("R4 asymptotic, a<=z<=2a", pts))

  # R5: asymptotic branch where it is valid
  pts = []
  for _ in range(n):
    a = log_uniform(rng, 0.1, 60.0)
    lo = max(8.0, 3.0 * a)
    pts.append((a, min(600.0, lo * rng.uniform(1.0, 4.0))))
  out.append(("R5 asymptotic, z>=3a", pts))

  # R6: the seam at z = 8
  pts = [(log_uniform(rng, 1.0, 8.0), rng.uniform(7.5, 8.5))
         for _ in range(n)]
  out.append(("R6 seam z=8", pts))

  # R7: integer a, the discrete-count case
  pts = []
  for _ in range(n):
    a = float(rng.randint(1, 40))
    hi = min(3.0 * a, safe_z_max(a, 600.0))
    pts.append((a, rng.uniform(0.5 * a, max(0.51 * a, hi))))
  out.append(("R7 integer a, 0.5a<z<3a", pts))

  # R8: above the tgamma overflow limit
  pts = []
  for _ in range(n):
    a = log_uniform(rng, 150.0, 400.0)
    pts.append((a, a * rng.uniform(1.0, 1.5)))
  out.append(("R8 a>171.6 (tgamma overflow)", pts))

  # R9: far upper tail, where grad_reg_lower_inc_gamma cancels
  pts = []
  for _ in range(n):
    a = log_uniform(rng, 12.0, 120.0)
    pts.append((a, min(600.0, a * rng.uniform(3.0, 6.0))))
  out.append(("R9 z>=3a, a>=12 (lower-tail cancel)", pts))

  return out


def main():
  n = int(sys.argv[1]) if len(sys.argv) > 1 else N_DEFAULT
  path = sys.argv[2] if len(sys.argv) > 2 else "points.txt"
  rng = random.Random(SEED)
  total = 0
  with open(path, "w") as fh:
    for name, pts in regions(n, rng):
      fh.write(f"# {name}\n")
      for a, z in pts:
        fh.write(f"{a.hex()} {z.hex()}\n")
        total += 1
  print(f"wrote {total} points to {path}")
  return 0


if __name__ == "__main__":
  sys.exit(main())
