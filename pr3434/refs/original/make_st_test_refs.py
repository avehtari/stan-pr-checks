#!/usr/bin/env python3
"""Choose test rows for ST, MST, YS and LBETA and print C++ rows.

Usage: make_st_test_refs.py   (run in beta_binomial_work/st)

Inputs: out_dev.txt and out_br.txt, outputs of density_sweep_work/
probe_sweep.cpp on points_st.txt, built against develop and the branch.
References: density_sweep_work/ref_mpmath.py (lp_* and grad_*, 50 digits,
route A checked against mp.diff at 110 digits there).

Per tag up to N_FAIL rows where develop fails by at least FAIL_RATIO times
the tolerance and the branch is within the tolerance / 10, plus one guard
row where both pass. Gradients of positive parameters are compared as
x * d/dx. Output rows: {"TAG", {args as hex}, value, {raw gradients}}.
"""

import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "density_sweep_work"))

from mpmath import mp, mpf  # noqa: E402

import ref_mpmath as R  # noqa: E402

TAGS = ["ST", "MST", "YS", "LBETA"]
LOG_ARGS = {"ST": [-1, 1, -1, 3], "MST": [-1, -1, 2], "YS": [1],
            "LBETA": [0, 1]}
TV, TG = 1e-12, 1e-11
N_FAIL = 3
FAIL_RATIO = 100.0


def parse(path):
  out = {}
  for line in open(path):
    if line.startswith("#") or "=>" not in line:
      continue
    left, right = line.split("=>")
    f = left.split()
    if f[0] not in TAGS:
      continue
    args = tuple(float.fromhex(t) for t in f[1:])
    vals = []
    for t in right.split()[:-1]:
      t = t.replace("-nan", "nan")
      vals.append(float(t) if t in ("nan", "inf", "-inf") else float.fromhex(t))
    out[(f[0], args)] = vals
  return out


def ref(tag, args):
  v = R.LP[tag](*args)
  g = R.GRAD[tag](*args)
  return [v] + list(g)


def cols(tag, args, vals):
  out = [mpf(vals[0])]
  for j, la in enumerate(LOG_ARGS[tag]):
    s = mpf(args[la]) if la >= 0 else mpf(1)
    out.append(s * mpf(vals[1 + j]))
  return out


def ratio(tag, args, got, r):
  cg, cr = cols(tag, args, got), cols(tag, args, r)
  worst = 0.0
  for j, (a, b) in enumerate(zip(cg, cr)):
    if not mp.isfinite(a):
      return math.inf
    tol = (TV if j == 0 else TG) * max(1, abs(b))
    worst = max(worst, float(abs(a - b) / tol))
  return worst


def cxx(x):
  return mp.nstr(mpf(x), 17, min_fixed=-1, max_fixed=-1, strip_zeros=False)


def main():
  dev, br = parse("out_dev.txt"), parse("out_br.txt")
  rows = []
  for tag in TAGS:
    cands = []
    for k in br:
      if k[0] != tag or k not in dev:
        continue
      r = ref(tag, k[1])
      cands.append((ratio(tag, k[1], dev[k], r), ratio(tag, k[1], br[k], r),
                    k, r))
    fails = sorted((c for c in cands if c[1] <= 0.1 and c[0] >= FAIL_RATIO),
                   key=lambda c: -c[0])[:N_FAIL]
    guards = sorted((c for c in cands if c[1] <= 0.1 and c[0] <= 0.1),
                    key=lambda c: max(abs(a) for a in c[2][1]))[:1]
    for rd, rb, k, r in fails + guards:
      kind = "fail" if rd >= FAIL_RATIO else "guard"
      print(f"# {tag:6s} dev/tol {rd:9.2e} br/tol {rb:9.2e} {kind}")
      rows.append((k, r, kind))
  for k, r, kind in rows:
    args = ", ".join(float(a).hex() for a in k[1])
    grads = ", ".join(cxx(g) for g in r[1:])
    print(f'    {{"{k[0]}", {{{args}}}, {cxx(r[0])}, {{{grads}}}}},')


if __name__ == "__main__":
  main()
