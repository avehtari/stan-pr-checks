#!/usr/bin/env python3
"""Choose test points for the root-cause families and print C++ rows.

Usage: make_family_test_refs.py   (run in beta_binomial_work/sib)

Inputs: refs.txt (certified references), out_develop*.txt and
out_branch*.txt (probe outputs on develop and on the branch).

Per tag, up to N_FAIL points where develop fails by at least FAIL_RATIO
times the tolerance and the branch is within the tolerance / 10 in every
column, plus one guard point where both pass. Columns are compared as
log-shape gradients x * d/dx for the shape-type parameters (SHAPE of
score_sib.py) and raw for the others.

Output: a C++ initializer list, one row per point:
  {"TAG", {args as hex}, value, {raw gradients}},
"""

import math

from mpmath import mp, mpf

from score_sib import SHAPE, OTHER

mp.dps = 60
TOL = {  # (value, gradient) tolerance factors
    "NB": (1e-12, 1e-11), "NB2": (1e-12, 1e-11), "NB2L": (1e-12, 1e-11),
    "GLM1": (1e-12, 1e-11), "GLM5": (1e-12, 1e-11), "BNB": (1e-12, 1e-11),
    "DM": (1e-12, 1e-11), "LKJ": (1e-12, 1e-11), "LKJC": (1e-12, 1e-11),
    "YSC": (1e-12, 1e-11), "YSLC": (1e-12, 1e-11), "YSLCC": (1e-12, 1e-11),
    "BBC": (1e-6, 1e-4), "BBLC": (1e-6, 1e-4), "BBLCC": (1e-6, 1e-4),
    "BNBC": (1e-6, 1e-4), "BNBLCC": (1e-6, 1e-4),
}
N_FAIL = 3
FAIL_RATIO = 100.0


def parse(t):
  t = t.replace("-nan", "nan")
  return float(t) if t in ("nan", "inf", "-inf") else float.fromhex(t)


def key(f):
  return (f[0], f[1]) + tuple(float.fromhex(t) for t in f[2:])


def load(paths, is_ref=False):
  out = {}
  for p in paths:
    for line in open(p):
      f = line.split()
      if not f or f[0] not in TOL or "=>" not in line:
        continue
      left, right = line.split("=>")
      r = right.split("|")[0].split()
      out[key(left.split())] = ([mpf(t) for t in r] if is_ref
                                else [parse(t) for t in r[:-2]])
  return out


def columns(tag, args, vals):
  """value, then one entry per gradient: log-shape or raw."""
  out = [vals[0]]
  shape_idx = dict(SHAPE.get(tag, []))
  for g in range(len(vals) - 1):
    v = vals[1 + g]
    if g in shape_idx:
      out.append(mpf(args[shape_idx[g]]) * mpf(v))
    else:
      out.append(mpf(v))
  return out


def ratios(tag, args, got, ref):
  """worst error / tolerance over all columns (inf if non-finite)."""
  tv, tg = TOL[tag]
  cg = columns(tag, args, got)
  cr = columns(tag, args, ref)
  worst = 0.0
  for j, (a, b) in enumerate(zip(cg, cr)):
    if not mp.isfinite(a):
      return math.inf
    tol = (tv if j == 0 else tg) * max(1, abs(b))
    worst = max(worst, float(abs(a - b) / tol))
  return worst


def cxx(x):
  return mp.nstr(mpf(x), 17, min_fixed=-1, max_fixed=-1, strip_zeros=False)


def main():
  ref = load(["refs.txt"], True)
  dev = load(["out_develop.txt", "out_develop_extra.txt"])
  br = load(["out_branch.txt", "out_branch_extra.txt"])
  rows = []
  print(f"{'tag':7s} {'label':10s} {'dev/tol':>9s} {'br/tol':>9s}  kind")
  for tag in TOL:
    cands = []
    for k in br:
      if k[0] != tag or k not in dev or k not in ref:
        continue
      args = k[2:]
      rb = ratios(tag, args, br[k], ref[k])
      rd = ratios(tag, args, dev[k], ref[k])
      cands.append((rd, rb, k))
    fails = sorted((c for c in cands if c[1] <= 0.1 and c[0] >= FAIL_RATIO),
                   key=lambda c: -c[0])[:N_FAIL]
    guards = sorted((c for c in cands if c[1] <= 0.1 and c[0] <= 0.1),
                    key=lambda c: max(c[2][2:]))[:1]
    for rd, rb, k in fails + guards:
      kind = "fail" if rd >= FAIL_RATIO else "guard"
      print(f"{tag:7s} {k[1]:10s} {rd:9.2e} {rb:9.2e}  {kind}")
      rows.append((k, ref[k], kind))
  print()
  for k, r, kind in rows:
    args = ", ".join(float(a).hex() for a in k[2:])
    grads = ", ".join(cxx(g) for g in r[1:])
    print(f'    {{"{k[0]}", {{{args}}}, {cxx(r[0])}, {{{grads}}}}},'
          f'  // {k[1]} {kind}')


if __name__ == "__main__":
  main()
