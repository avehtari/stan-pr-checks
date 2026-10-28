#!/usr/bin/env python3
# Scores ref324_probe output of one or two trees against the 324 rows.
# Supports "199 of 324 non-finite", "1.6e+93" and "values change only
# where they were -inf or lost digits" (#3408 summary table, Part 3).
"""Score ref324_probe output against the reference rows it prints.

Usage: python3 ref324_score.py LABEL=FILE [LABEL=FILE]

Per file: the count of rows with a throw or a non-finite value or gradient
component, the worst relative value and gradient error over the finite
results, and every row whose value is off by more than 1e-12 (the
tolerance of the reference test).  Relative error uses max(|ref|, 1e-300)
in the denominator, as the reference test does.

With two files (first = develop, second = PR): every row where the two
values differ by more than 1e-12 relative, with the develop value
classified as non-finite, finite and off by more than 1e-12, or within
1e-12 of the reference.  The claim "values change only where they were
-inf or lost digits" holds if the last class is empty.
"""
import math
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
HEADER = os.path.join(HERE, "incomplete_gamma_cdfs_reference_324.hpp")
TOL_VALUE = 1e-12
ROW_RE = re.compile(r'^\s*\{"(\w+)", ([^,]+), ([^,]+), ([^,]+), .*//\s*(\w+)')


def read_args():
  """(function, y, p1, p2, class) of each header row, in order."""
  out = []
  for line in open(HEADER):
    m = ROW_RE.match(line)
    if m:
      out.append((m.group(1), m.group(2), m.group(3), m.group(4), m.group(5)))
  return out


def read(path):
  rows = []
  for line in open(path):
    t = line.split()
    if not t:
      continue
    idx, name, n_grad = int(t[0]), t[1], int(t[2])
    ref = [float.fromhex(x) for x in t[3:7]]
    if t[7] == "THROW":
      got = None
    else:
      got = [float.fromhex(x) for x in t[7:11]]
    rows.append((idx, name, n_grad, ref[0], ref[1:1 + n_grad],
                 None if got is None else got[0],
                 None if got is None else got[1:1 + n_grad]))
  return rows


def rel(x, r):
  if not math.isfinite(x):
    return math.inf
  return abs(x - r) / max(abs(r), 1e-300)


def where(args, idx):
  name, y, p1, p2, cls = args[idx]
  return f"{name}({y} | {p1}, {p2}) [{cls}]"


def score(label, rows, args):
  print(f"== {label}: {len(rows)} rows")
  n_throw = sum(1 for r in rows if r[5] is None)
  n_nonfin = 0
  worst_v, worst_g = (0.0, None), (0.0, None)
  off = []
  for idx, name, ng, rv, rg, gv, gg in rows:
    if gv is None:
      continue
    if not (math.isfinite(gv) and all(math.isfinite(g) for g in gg)):
      n_nonfin += 1
    ev = rel(gv, rv)
    if math.isfinite(ev) and ev > worst_v[0]:
      worst_v = (ev, idx)
    if ev > TOL_VALUE:
      off.append((idx, gv, rv, ev))
    for k, (g, r) in enumerate(zip(gg, rg)):
      eg = rel(g, r)
      if math.isfinite(eg) and eg > worst_g[0]:
        worst_g = (eg, (idx, k, g, r))
  print(f"   throws {n_throw}; rows with a non-finite value or gradient "
        f"{n_nonfin}; total {n_throw + n_nonfin} of {len(rows)}")
  if worst_v[1] is not None:
    print(f"   worst finite value error    {worst_v[0]:.2e}  "
          f"{where(args, worst_v[1])}")
  if worst_g[1] is not None:
    idx, k, g, r = worst_g[1]
    print(f"   worst finite gradient error {worst_g[0]:.2e}  "
          f"{where(args, idx)} component {k}: got {g:.6g} ref {r:.6g}")
  print(f"   values off by more than {TOL_VALUE:g}: {len(off)}")
  for idx, gv, rv, ev in off:
    print(f"     {where(args, idx):60s} got {gv:.6g} ref {rv:.6g} "
          f"rel {ev:.2e}")
  return {r[0]: r for r in rows}


def compare(la, a, lb, b, args):
  print(f"== value changes from {la} to {lb} (relative difference > "
        f"{TOL_VALUE:g})")
  counts = {"non-finite": 0, "off": 0, "within": 0}
  for idx in sorted(a):
    va, vb, rv = a[idx][5], b[idx][5], a[idx][3]
    if va is None or vb is None:
      continue
    if math.isfinite(va) and math.isfinite(vb) and rel(vb, va) <= TOL_VALUE:
      continue
    if va == vb:
      continue
    ea = rel(va, rv)
    cls = ("non-finite" if not math.isfinite(va)
           else "off" if ea > TOL_VALUE else "within")
    counts[cls] += 1
    print(f"   {where(args, idx):60s} {la} {va:.6g} (rel {ea:.1e}) "
          f"{lb} {vb:.6g} ref {rv:.6g}  {cls}")
  print(f"   {la} value non-finite: {counts['non-finite']}; finite and off by "
        f"more than {TOL_VALUE:g}: {counts['off']}; within {TOL_VALUE:g} of "
        f"the reference: {counts['within']}")


def main():
  args = read_args()
  files = [a.split("=", 1) for a in sys.argv[1:]]
  data = []
  for label, path in files:
    rows = read(path)
    assert len(rows) == len(args), f"{path}: {len(rows)} rows, header {len(args)}"
    data.append((label, score(label, rows, args)))
  if len(data) == 2:
    compare(data[0][0], data[0][1], data[1][0], data[1][1], args)


if __name__ == "__main__":
  main()
