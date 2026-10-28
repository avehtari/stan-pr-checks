#!/usr/bin/env python3
# Scores the fvar<fvar<double>> second derivative that probe_cpp writes to
# stderr.  Supports "second derivatives within 1.9e-10 (7e-9 at alpha = 5)".
"""Relative error of d2/dalpha2 gamma_lccdf(y | alpha, 1), fvar<fvar<double>>.

Usage: python3 score_cpp_ffd2.py LABEL=FILE [LABEL=FILE ...]

FILE is the stderr of probe_cpp (alpha, y, d2 in hex, one line per point).
Reference: mp.diff of the quadrature route of refs/probe_lccdf_tail.py, the
same reference as the d2 column of score_cpp.py.
"""
import math
import os
import sys

import mpmath as mp

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                "..", "refs"))
from probe_lccdf_tail import ref_quad  # noqa: E402


def read(path):
  return [[float.fromhex(t) for t in line.split()] for line in open(path)]


def main():
  files = [arg.split("=", 1) for arg in sys.argv[1:]]
  data = {label: read(path) for label, path in files}
  first = data[files[0][0]]
  print(f"{'alpha':>7} {'y':>8} {'reference':>12}"
        + "".join(f" {label[:10]:>12}" for label, _ in files))
  worst = {label: 0.0 for label, _ in files}
  for k, row in enumerate(first):
    a, y = row[0], row[1]
    d2 = mp.diff(lambda s: ref_quad(s, y)[1], a)
    line = f"{a:7g} {y:8g} {float(d2):12.5g}"
    for label, _ in files:
      v = data[label][k][2]
      if not math.isfinite(v):
        line += f" {'nonfin':>12}"
        worst[label] = math.inf
        continue
      e = float(abs((v - d2) / d2))
      worst[label] = max(worst[label], e)
      line += f" {e:12.1e}"
    print(line)
  print("worst: " + ", ".join(f"{lb} {w:.2e}" for lb, w in worst.items()))


if __name__ == "__main__":
  main()
