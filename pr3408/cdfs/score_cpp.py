#!/usr/bin/env python3
# Scores probe_cpp output per autodiff mode (refs/probe_lccdf_tail.py).
# Supports 5.6e-12 / 1.9e-10 / 7e-9 of #3408 Part 3.
"""Score probe_cpp output files: relative error of d/dalpha per mode.

Usage: score_cpp.py LABEL=FILE [LABEL=FILE ...]
"""
import mpmath as mp
import os
import sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                "..", "refs"))

from probe_lccdf_tail import ref, ref_quad  # noqa: E402

MODES = ["var", "fvar", "fvar_var", "fvar_fvar"]


def read(path):
  rows = []
  for line in open(path):
    if line.startswith("#"):
      continue
    rows.append([float.fromhex(t) for t in line.split()])
  return rows


def rel(x, r):
  if x != x or x in (float("inf"), float("-inf")):
    return "nonfin"
  return f"{float(abs((x - r) / r)):.1e}"


def main():
  files = [arg.split("=", 1) for arg in sys.argv[1:]]
  data = {label: read(path) for label, path in files}
  first = data[files[0][0]]
  head = f"{'alpha':>7} {'y':>8} {'logQ':>9}"
  for label, _ in files:
    head += "".join(f" {label[:4] + ':' + m[:6]:>12}" for m in MODES)
  print(head)
  for k, row in enumerate(first):
    a, y = row[0], row[1]
    lq, dlq, _ = ref(a, y)
    line = f"{a:7g} {y:8g} {float(lq):9.1f}"
    for label, _ in files:
      vals = data[label][k][3:7]
      line += "".join(f" {rel(v, dlq):>12}" for v in vals)
    if len(row) > 7:
      d2 = mp.diff(lambda s: ref_quad(s, y)[1], a)
      line += "  d2:" + " ".join(rel(data[lb][k][7], d2) for lb, _ in files)
    print(line)


if __name__ == "__main__":
  main()
