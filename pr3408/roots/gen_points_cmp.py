#!/usr/bin/env python3
# Writes points_cmp.txt, the 49 990-point comparison (seed fixed; n = 5000).
# Supports the 6 379 / 0 non-finite counts and the z = a limit of #3408.
"""Point file for the develop-against-branch comparison.

Usage: python3 gen_points_cmp.py [n_per_region] [outfile]

Ranges agreed with the user on 2026-09-19: a and z both in [1e-3, 1e4],
log-spaced, with denser sampling on the structures the code branches on.

Format: `# <region>` header lines, then one `a z` pair of hex doubles per
line. Hex is exact in both directions, so the C++ reads exactly the points
the Python scorer uses.
"""

import math
import random
import sys

SEED = 20260919

LO = 1e-3
HI = 1e4

# a * log z = 709 is where exp(a * log z) overflows in the develop series.
OVERFLOW_PRODUCT = 709.782712893384
# tgamma(a) overflows just above this.
TGAMMA_LIMIT = 171.6243769563027


def log_uniform(rng, lo, hi):
  return math.exp(rng.uniform(math.log(lo), math.log(hi)))


def jitter(rng, x, rel):
  """A point near x, within a relative band."""
  return x * math.exp(rng.uniform(-rel, rel))


def regions(n, rng):
  out = []

  def add(name, pts):
    out.append((name, [(a, z) for a, z in pts if a > 0 and z > 0]))

  # C1: the whole box, log-uniform in both
  add("C1 whole box", [(log_uniform(rng, LO, HI), log_uniform(rng, LO, HI))
                       for _ in range(n)])

  # C2: small shape
  add("C2 a<1", [(log_uniform(rng, LO, 1.0), log_uniform(rng, LO, HI))
                 for _ in range(n)])

  # C3: z below the mode, the deep left tail
  add("C3 z<<a", [(lambda a: (a, log_uniform(rng, LO, max(LO * 1.01, a / 3))))(
      log_uniform(rng, 1.0, HI)) for _ in range(n)])

  # C4: z near a, where develop takes the invalid asymptotic branch
  add("C4 z near a", [(lambda a: (a, jitter(rng, a, 0.7)))(
      log_uniform(rng, 1.0, HI)) for _ in range(n)])

  # C5: z well above a, where the develop Gautschi branch cancels
  add("C5 z>=3a", [(lambda a: (a, min(HI, a * rng.uniform(3.0, 30.0))))(
      log_uniform(rng, LO, HI / 3)) for _ in range(n)])

  # C6: the develop branch seam at z = 8
  add("C6 seam z=8", [(log_uniform(rng, LO, HI), jitter(rng, 8.0, 0.08))
                      for _ in range(n)])

  # C7: integer shape, the discrete-count case
  add("C7 integer a", [(lambda a: (a, jitter(rng, a, 1.5)))(
      float(rng.randint(1, 400))) for _ in range(n)])

  # C8: across the tgamma overflow limit
  add("C8 a near 171.62", [(jitter(rng, TGAMMA_LIMIT, 0.05),
                            log_uniform(rng, LO, HI)) for _ in range(n)])

  # C9: across the series overflow limit, a * log z = 709
  pts = []
  for _ in range(n):
    z = log_uniform(rng, 1.05, HI)
    a = jitter(rng, OVERFLOW_PRODUCT / math.log(z), 0.05)
    if LO <= a <= HI:
      pts.append((a, z))
  add("C9 a*log z near 709", pts)

  # C10: the far corners
  corners = [(LO, LO), (LO, HI), (HI, LO), (HI, HI), (LO, 1.0), (HI, 1.0),
             (1.0, LO), (1.0, HI), (TGAMMA_LIMIT, TGAMMA_LIMIT)]
  pts = []
  for a, z in corners:
    pts.append((a, z))
    for _ in range(max(1, n // len(corners))):
      pts.append((jitter(rng, a, 0.02), jitter(rng, z, 0.02)))
  add("C10 corners", pts)

  return out


def main():
  n = int(sys.argv[1]) if len(sys.argv) > 1 else 200
  path = sys.argv[2] if len(sys.argv) > 2 else "points_cmp.txt"
  rng = random.Random(SEED)
  total = 0
  with open(path, "w") as fh:
    for name, pts in regions(n, rng):
      fh.write(f"# {name}\n")
      for a, z in pts:
        fh.write(f"{a.hex()} {z.hex()}\n")
        total += 1
  print(f"wrote {total} points in 10 regions to {path}")
  return 0


if __name__ == "__main__":
  sys.exit(main())
