#!/usr/bin/env python3
"""Two-route mpmath references for the points of sweep/points_st.txt.

Route A: the value and the analytic gradients of ref_sweep.py at 50 digits
(loggamma for the values, closed-form digamma differences).
Route B: the value at 110 digits and mp.diff of the value at 110 digits.
MSTF (multi_student_t_lpdf with Sigma = I) has the reference of MST
(multi_student_t_cholesky_lpdf with L = I).

Usage: ref_st.py <points file> <refs file> [processes]
Writes one line per distinct point:
  <TAG> <hex args> => <A value> <A grads> | <B value> <B grads>
and prints the worst disagreement of the two routes per tag: relative where
the larger magnitude is at least 1e-40, absolute below.
"""

import multiprocessing
import sys

from mpmath import mp, mpf

import ref_sweep as R

NSTR = 45


def ref_tag(tag):
  return "MST" if tag == "MSTF" else tag


def one(item):
  tag, hexes = item
  args = [float.fromhex(h) for h in hexes]
  t = ref_tag(tag)
  mp.dps = R.DPS
  va = R.LP[t](*args)
  ga = R.GRAD[t](*args)
  with mp.workdps(R.DPS_DIFF):
    vb = R.LP[t](*args)
  gb = R.grad_route_b(t, args)
  mp.dps = R.DPS
  a = " ".join(mp.nstr(x, NSTR) for x in [va] + ga)
  b = " ".join(mp.nstr(x, NSTR) for x in [vb] + gb)
  return f"{tag} {' '.join(hexes)} => {a} | {b}"


def disagreement(x, y):
  # absolute below 1e-40: route A at 50 digits leaves rounding of about
  # 1e-50 where the exact gradient is 0 (ST d/dsigma at |y - mu| = sigma)
  scale = max(abs(x), abs(y))
  if scale < mpf("1e-40"):
    return abs(x - y)
  return abs(x - y) / scale


def main():
  if len(sys.argv) not in (3, 4):
    print(__doc__)
    return 2
  procs = int(sys.argv[3]) if len(sys.argv) == 4 else 8
  items = []
  seen = set()
  for line in open(sys.argv[1]):
    if line.startswith("#") or not line.strip():
      continue
    parts = line.split()
    key = (parts[0], tuple(parts[1:]))
    if key not in seen:
      seen.add(key)
      items.append(key)
  with multiprocessing.Pool(procs) as pool:
    lines = pool.map(one, items, chunksize=8)
  with open(sys.argv[2], "w") as fh:
    fh.write("\n".join(lines) + "\n")
  mp.dps = 60
  worst = {}
  for line in lines:
    left, right = line.split(" => ")
    tag = left.split()[0]
    a, b = right.split(" | ")
    w = worst.setdefault(tag, [mpf(0), mpf(0), 0])
    av = [mpf(x) for x in a.split()]
    bv = [mpf(x) for x in b.split()]
    w[0] = max(w[0], disagreement(av[0], bv[0]))
    for x, y in zip(av[1:], bv[1:]):
      w[1] = max(w[1], disagreement(x, y))
    w[2] += 1
  print("route A (50 digits, analytic) against route B (110 digits, mp.diff)")
  print(f"{'tag':6s} {'points':>6s} {'value A-B':>10s} {'grad A-B':>10s}")
  for tag, w in worst.items():
    print(f"{tag:6s} {w[2]:6d} {float(w[0]):10.1e} {float(w[1]):10.1e}")
  return 0


if __name__ == "__main__":
  sys.exit(main())
