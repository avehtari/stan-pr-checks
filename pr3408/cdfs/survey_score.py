#!/usr/bin/env python3
# Scores survey_probe output against survey_refs.txt; run from cdfs/.
# Supports the 2e-13 / 8e-13 / 2.7e-10 claim of #3408 Part 3.
"""Score survey_probe output against survey_refs.txt.

Usage: survey_score.py OUT [--detail]
Prints, per function and point class, the worst relative error of the
value and of the gradient (worst over its components), with markers:
  inf   a non-finite result where the reference is finite
  THR   the function threw
Relative error uses max(|ref|, 1e-300) in the denominator.
"""
import math
import sys
from collections import defaultdict

cases = [l.split() for l in open("survey_cases.txt")]
refs = [[float(t) for t in l.split()] for l in open("survey_refs.txt")]
outs = [l.split() for l in open(sys.argv[1])]
detail = "--detail" in sys.argv


def rel(x, r):
  if not math.isfinite(x):
    return math.inf
  return abs(x - r) / max(abs(r), 1e-300)


table = defaultdict(lambda: [0.0, 0.0])
for c, r, o in zip(cases, refs, outs):
  name, cls, a, z = c[0], c[4], c[5], c[6]
  key = (name, cls)
  if o == ["THROW"]:
    table[key][0] = table[key][1] = math.nan
    if detail:
      print(f"{name:28} a={a:>6} z={z:>6} THROW")
    continue
  vals = [float.fromhex(t) for t in o]
  ev = rel(vals[0], r[0])
  eg = max(rel(v, rr) for v, rr in zip(vals[1:], r[1:]))
  t = table[key]
  t[0] = max(t[0], ev) if not math.isnan(t[0]) else t[0]
  t[1] = max(t[1], eg) if not math.isnan(t[1]) else t[1]
  if detail:
    gs = " ".join(f"{rel(v, rr):8.1e}" for v, rr in zip(vals[1:], r[1:]))
    print(f"{name:28} a={a:>6} z={z:>6} val {ev:8.1e} ref {r[0]:10.4g} grad {gs}")


def fmt(x):
  if math.isnan(x):
    return "THR"
  if math.isinf(x):
    return "inf"
  return f"{x:.0e}"


classes = ["upper_deep", "upper_tail", "bulk", "bulk_big", "lower_tail",
           "lower_deep"]
names = list(dict.fromkeys(c[0] for c in cases))
print(f"{'function':28}" + "".join(f" {c[:10]:>11}" for c in classes))
for n in names:
  row = f"{n:28}"
  for c in classes:
    v, g = table.get((n, c), [math.nan, math.nan])
    row += f" {fmt(v):>5}/{fmt(g):>5}"
  print(row)
print("cell = worst value error / worst gradient error")
