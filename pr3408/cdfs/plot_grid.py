#!/usr/bin/env python3
# Scores grid_probe on the 266-cell grid and writes the CSV for plot_grid.R.
# Supports the 266-cell grid claim (all finite, 3.4e-12) of #3408.
"""Relative error of d/dalpha gamma_lccdf(y | alpha, 1) on a grid.

Same grid as the plot in the stan-dev/math#3408 comment by mjcarter95:
alpha from 5 to 1000, y / alpha from 1 to 100. Reference: mpmath, two
routes (probe_lccdf_tail.ref). Writes a CSV for plot_grid.R.

Usage: plot_grid.py PROBE_BINARY OUT_CSV
"""
import math
import subprocess
import os
import sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                "..", "refs"))

from probe_lccdf_tail import ref  # noqa: E402

ALPHAS = [5, 10, 20, 30, 50, 70, 100, 150, 200, 300, 500, 700, 850, 1000]
RATIOS = [1, 1.05, 1.1, 1.2, 1.3, 1.5, 2, 3, 4, 5, 6, 8, 10, 15, 20, 30, 50,
          75, 100]


def main():
  probe, out = sys.argv[1], sys.argv[2]
  pts = [(a, r) for r in RATIOS for a in ALPHAS]
  stdin = "".join(f"{a!r} {r * a!r}\n" for a, r in pts)
  lines = subprocess.run([probe], input=stdin, capture_output=True, text=True,
                         check=True).stdout.strip().split("\n")
  worst = 0.0
  n_bad = 0
  with open(out, "w") as f:
    f.write("alpha,ratio,log10_rel_err\n")
    for (a, r), line in zip(pts, lines):
      g = float.fromhex(line.split()[2])
      if not math.isfinite(g):
        n_bad += 1
        f.write(f"{a},{r},NA\n")
        continue
      _, d_ref, _ = ref(a, r * a)
      e = float(abs((g - d_ref) / d_ref))
      worst = max(worst, e)
      f.write(f"{a},{r},{math.log10(max(e, 1e-300))}\n")
  print(f"worst {worst:.2e} non-finite {n_bad} cells {len(pts)}")


if __name__ == "__main__":
  main()
