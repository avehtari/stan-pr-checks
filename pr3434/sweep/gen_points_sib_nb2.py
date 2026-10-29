#!/usr/bin/env python3
"""Dense NB2 / NB2L phi sample: points_nb2.txt.

The develop NB2 phi partial is round(round(small - psi(phi)) + psi(n + phi)).
Its error is 0 or about one ulp of psi(phi), depending on rounding. A sparse
sample can miss the one-ulp case, so sample 40 phi values per decade.
"""

import math


def main():
  out = []
  for s, name in [(1e10, "1e10"), (1e12, "1e12"), (1e15, "1e15")]:
    for k in range(40):
      phi = s * (1.0 + 0.0731 * k)
      for n, mu in [(1, 3.0), (57, 50.0)]:
        for tag, m in [("NB2", mu), ("NB2L", math.log(mu))]:
          out.append(f"{tag} dense:{name} " + " ".join(
              float(x).hex() for x in [n, m, phi]))
  with open("points_nb2.txt", "w") as fh:
    fh.write("\n".join(out) + "\n")
  print(f"points_nb2.txt: {len(out)} lines")


if __name__ == "__main__":
  main()
