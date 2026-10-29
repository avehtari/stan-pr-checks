#!/usr/bin/env python3
"""Grid of beta_neg_binomial cdf points for probe_bnb.cpp and ref_bnb.py.

Output line: BNB grid <n> <r> <alpha> <beta>, every number a hex double, so
C++ and Python read the same values.

Usage: gen_points_bnb.py <points file>
"""

import sys

ALPHA = [0.25, 0.5, 1.0, 2.0, 4.0, 10.0, 100.0, 1e4]
SHAPE = [0.5, 3.0, 30.0, 300.0, 3e3, 3e4, 1e6]
COUNT = [0, 2, 20, 200, 2000, 20000]


def main():
  if len(sys.argv) != 2:
    print(__doc__)
    return 2
  lines = []
  for a in ALPHA:
    for r in SHAPE:
      for b in SHAPE:
        for n in COUNT:
          hexes = " ".join(float(x).hex() for x in (n, r, a, b))
          lines.append(f"BNB grid {hexes}")
  with open(sys.argv[1], "w") as fh:
    fh.write("\n".join(lines) + "\n")
  print(f"{len(lines)} points")
  return 0


if __name__ == "__main__":
  sys.exit(main())
