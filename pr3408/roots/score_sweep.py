#!/usr/bin/env python3
# Scores sweep_roots output per region against mpmath (refs/ref_mpmath.py).
# Supports the nine-region table of #3408 Part 1.
"""Score sweep_values.txt against the mpmath reference, per region.

Usage: python3 score_sweep.py [sweep_values.txt] [n_processes]

Columns in the input, all hex:
  a  z  grad_reg_inc_gamma  grad_reg_lower_inc_gamma  Eigen igamma_der_a  status

Sign convention: column 3 is dQ/da, columns 4 and 5 are dP/da = -dQ/da.
Every column is therefore scored against the reference with its own sign.

Reported per region: median, 99th percentile and worst **relative** error,
and the count of non-finite results. Relative error is used rather than ulp
because the errors here span 1e-16 to 1e+133, and an ulp figure of 1e+150
carries no more information than the relative one.
"""

import math
import multiprocessing as mputil
import os
import sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                "..", "refs"))

import ref_mpmath as R  # noqa: E402


def rel(value, reference):
  ref = float(reference)
  if not math.isfinite(value):
    return math.inf
  if ref == 0.0:
    return 0.0 if value == 0.0 else math.inf
  return abs(value - ref) / abs(ref)


def score_point(rec):
  a, z, up, lo, ei = rec
  try:
    ref_q = float(R.dQda(a, z))
  except Exception:
    return None
  ref_p = -ref_q
  return (a, z, rel(up, ref_q), rel(lo, ref_p), rel(ei, ref_p),
          math.isfinite(up), math.isfinite(lo), math.isfinite(ei))


def summarize(name, rows):
  out = [f"{name:<36}"]
  for idx in (2, 3, 4):
    vals = sorted(r[idx] for r in rows)
    finite = [v for v in vals if math.isfinite(v)]
    med = finite[len(finite) // 2] if finite else math.inf
    worst = vals[-1]
    out.append(f"{med:>9.2e} {worst:>9.2e}")
  nf = [sum(0 if r[i] else 1 for r in rows) for i in (5, 6, 7)]
  out.append(f"  {nf[0]:>4d} {nf[1]:>4d} {nf[2]:>4d}")
  return " ".join(out)


def main():
  path = sys.argv[1] if len(sys.argv) > 1 else "sweep_values.txt"
  nproc = int(sys.argv[2]) if len(sys.argv) > 2 else 6

  regions = []
  current_name, current = None, []
  for line in open(path):
    line = line.strip()
    if not line:
      continue
    if line.startswith("#"):
      if line.startswith("# evaluated"):
        continue
      if current_name is not None:
        regions.append((current_name, current))
      current_name, current = line[2:], []
      continue
    parts = line.split()
    if len(parts) < 5:
      continue
    vals = [float.fromhex(p) for p in parts[:5]]
    current.append(tuple(vals))
  if current_name is not None:
    regions.append((current_name, current))

  print(f"# {path} against mpmath, relative error per region")
  print("# 'up' = grad_reg_inc_gamma, 'lo' = grad_reg_lower_inc_gamma, "
        "'ei' = Eigen numext::igamma_der_a")
  print("# the last three columns count non-finite results in each column\n")
  header = (f"{'region':<36} {'up med':>9} {'up wrst':>9} {'lo med':>9} "
            f"{'lo wrst':>9} {'ei med':>9} {'ei wrst':>9}   "
            f"{'nfU':>4} {'nfL':>4} {'nfE':>4}")
  print(header)

  with mputil.Pool(nproc) as pool:
    for name, recs in regions:
      rows = [r for r in pool.map(score_point, recs) if r is not None]
      if rows:
        print(summarize(name, rows))
  return 0


if __name__ == "__main__":
  sys.exit(main())
