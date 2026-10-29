#!/usr/bin/env python3
"""Worst errors at large shapes over ALL regions, two trees (claim P1).

Usage (from pr3434/):
  python3 sweep/score_p1.py <dev label> <pr label> [--regions]

Reads, for each label L in out/:
  sib_L.txt sib_extra_L.txt sib_slow_L.txt sib_nb2_L.txt   (probe_sib)
  st_L.txt                                                 (probe_st)
  bb_L.txt                                                 (probe_bb)
and the references refs/refs_sib.txt, refs/refs_st.txt, refs/refs_bb.txt.

Metrics, the worst over every point of every region of a function:
  value  |stan - ref|, absolute
  grad   log-scale gradient error |x (d/dx stan - d/dx ref)| over the
         positive shape-type arguments x (sib: the SHAPE table of
         score_sib.py; st: ST nu, MST/MSTF nu, YS alpha, LBETA x;
         bb: alpha and beta)
At the worst point the relative error is printed too, because some values
and gradients are large (lkj_corr at eta = 3e15: x d/dx is about 1e15).
Counts per function and tree: points, throws, finished calls with a
non-finite value or gradient (NaN and inf separately).
"""

import math
import os
import sys
from collections import OrderedDict

from mpmath import mp, mpf

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(ROOT, "refs"))

import score_sib as S  # noqa: E402

mp.dps = 60
INF = mpf("inf")

# rows of the worst-error table of the PR text: (row, [function keys])
ROWS = [
    ("beta_binomial_lpmf", ["BB"]),
    ("beta_neg_binomial_lpmf", ["BNB"]),
    ("neg_binomial_lpmf", ["NB"]),
    ("neg_binomial_2_lpmf", ["NB2"]),
    ("neg_binomial_2_log_lpmf", ["NB2L"]),
    ("neg_binomial_2_log_glm_lpmf", ["GLM1", "GLM5"]),
    ("dirichlet_multinomial_lpmf", ["DM"]),
    ("lkj_corr_lpdf", ["LKJ"]),
    ("lkj_corr_cholesky_lpdf", ["LKJC"]),
    ("yule_simon_lccdf", ["YSLCC"]),
    ("yule_simon_cdf", ["YSC"]),
    ("yule_simon_lcdf", ["YSLC"]),
    ("beta_binomial_cdf", ["BBC"]),
    ("beta_binomial_lcdf", ["BBLC"]),
    ("beta_binomial_lccdf", ["BBLCC"]),
    ("beta_neg_binomial_cdf", ["BNBC"]),
    ("beta_neg_binomial_lcdf", ["BNBLC"]),
    ("beta_neg_binomial_lccdf", ["BNBLCC"]),
    ("student_t_lpdf", ["ST"]),
    ("multi_student_t_cholesky_lpdf", ["MST"]),
    ("multi_student_t_lpdf", ["MSTF"]),
    ("yule_simon_lpmf", ["YS"]),
    ("lbeta (var), d/dx", ["LBETA"]),
]

# st tags: (gradient index, argument index) of the log-scale gradients
ST_SHAPE = {"ST": [(1, 1)], "MST": [(2, 2)], "MSTF": [(2, 2)],
            "YS": [(0, 1)], "LBETA": [(1, 1)]}


def hexf(tok):
  try:
    return float.fromhex(tok)
  except ValueError:
    return math.nan


def finite(x):
  return not (math.isnan(x) or math.isinf(x))


def new_point(fn, region, args, status, value, grads, verr, vrel, gerr, grel):
  return {"fn": fn, "region": region, "args": args, "status": status,
          "value": value, "grads": grads, "verr": verr, "vrel": vrel,
          "gerr": gerr, "grel": grel}


# ── sib: 13 families, references with routes A, B, C ─────────────────────

def load_sib(label, refs):
  pts = []
  for f in ["sib", "sib_extra", "sib_slow", "sib_nb2"]:
    path = os.path.join(ROOT, "out", f"{f}_{label}.txt")
    for r in S.load_probe(path):
      ref = refs[r["key"]]
      st = "T" if r["status"] == "T" else "."
      e = S.errors(r, ref)
      if e is None:
        pts.append(new_point(r["tag"], r["label"], r["args"], "T", r["value"],
                             r["grads"], None, None, None, None))
        continue
      vr = ref["A"][0]
      pts.append(new_point(r["tag"], r["label"], r["args"], st, r["value"],
                           r["grads"], e["val_abs"], e["val_rel"],
                           e["lsg_abs"], e["lsg_rel"]))
  return pts


# ── st: student_t family, yule_simon_lpmf, lbeta ─────────────────────────

def norm(tag, hexes):
  return (tag,) + tuple(float.fromhex(h) for h in hexes)


def load_st_refs(path):
  refs = {}
  for line in open(path):
    left, right = line.rstrip("\n").split(" => ")
    parts = left.split()
    a, b = right.split(" | ")
    refs[norm(parts[0], parts[1:])] = [mpf(x) for x in a.split()]
  return refs


def load_st(label, refs):
  pts = []
  region = "(none)"
  for line in open(os.path.join(ROOT, "out", f"st_{label}.txt")):
    line = line.rstrip("\n")
    if line.startswith("# B"):
      region = line[2:].split()[0]
      continue
    if line.startswith("#") or "=>" not in line:
      continue
    left, right = line.split(" => ")
    parts = left.split()
    tag = parts[0]
    args = [float.fromhex(h) for h in parts[1:]]
    out = right.split()
    status = out[-1]
    nums = [hexf(t) for t in out[:-1]]
    if status == "T":
      pts.append(new_point(tag, region, args, "T", math.nan, nums[1:], None,
                           None, None, None))
      continue
    ref = refs[norm(tag, parts[1:])]
    value, grads = nums[0], nums[1:]
    if finite(value):
      verr = abs(mpf(value) - ref[0])
      vrel = verr / abs(ref[0]) if ref[0] != 0 else verr
    else:
      verr = vrel = INF
    gerr = mpf(0)
    grel = mpf(0)
    for gi, ai in ST_SHAPE[tag]:
      x = mpf(args[ai])
      if not finite(grads[gi]):
        gerr = grel = INF
        continue
      d = abs(x * (mpf(grads[gi]) - ref[1 + gi]))
      gerr = max(gerr, d)
      t = abs(x * ref[1 + gi])
      if t > mpf("1e-6"):
        grel = max(grel, d / t)
    if tag == "LBETA":
      verr = vrel = None  # the row is about the derivative only
    pts.append(new_point(tag, region, args, status, value, grads, verr, vrel,
                         gerr, grel))
  return pts


# ── bb: beta_binomial_lpmf, 13 177 points ────────────────────────────────

def load_bb_refs(path):
  refs = {}
  for line in open(path):
    parts = line.split()
    key = (int(parts[0]), int(parts[1]), float.fromhex(parts[2]),
           float.fromhex(parts[3]))
    refs[key] = [mpf(x) for x in parts[4:]]  # lchoose value d_a d_b
  return refs


def load_bb(label, refs):
  import ref_bb as RB
  pts = []
  for line in open(os.path.join(ROOT, "out", f"bb_{label}.txt")):
    if line.startswith("#") or not line.strip():
      continue
    left, right = line.split("=>")
    tag, n, N, a, b = left.split()
    vals = right.split()
    status = vals[-1]
    v_dbl, v_var, d_a, d_b = (hexf(t) for t in vals[:4])
    key = (int(n), int(N), float.fromhex(a), float.fromhex(b))
    if key not in refs:
      r = RB.ref_point(*key)
      refs[key] = [r["lchoose"], r["value"], r["d_alpha"], r["d_beta"]]
    ref = refs[key]
    args = [key[0], key[1], key[2], key[3]]
    if "T" in status:
      pts.append(new_point("BB", tag, args, "T", v_var, [d_a, d_b], None,
                           None, None, None))
      continue
    if finite(v_var):
      verr = abs(mpf(v_var) - ref[1])
      vrel = verr / abs(ref[1]) if ref[1] != 0 else verr
    else:
      verr = vrel = INF
    gerr = mpf(0)
    grel = mpf(0)
    for g, x, gr in [(d_a, key[2], ref[2]), (d_b, key[3], ref[3])]:
      if not finite(g):
        gerr = grel = INF
        continue
      d = abs(mpf(x) * (mpf(g) - gr))
      gerr = max(gerr, d)
      t = abs(mpf(x) * gr)
      if t > mpf("1e-6"):
        grel = max(grel, d / t)
    pts.append(new_point("BB", tag, args, status, v_var, [d_a, d_b], verr,
                         vrel, gerr, grel))
  return pts


# ── summary ──────────────────────────────────────────────────────────────

def fmt(x):
  if x is None:
    return "-"
  x = float(x)
  if math.isinf(x):
    return "inf"
  if math.isnan(x):
    return "nan"
  if x == 0:
    return "0"
  return f"{x:.1e}"


def counts(pts):
  thr = sum(1 for p in pts if p["status"] == "T")
  nan = inf = 0
  for p in pts:
    if p["status"] == "T":
      continue
    nums = [p["value"]] + list(p["grads"])
    if any(math.isnan(x) for x in nums):
      nan += 1
    elif any(math.isinf(x) for x in nums):
      inf += 1
  return thr, nan, inf


def worst(pts, key):
  best = None
  for p in pts:
    if p[key] is None:
      continue
    if best is None or p[key] > best[key]:
      best = p
  return best


def point_str(p):
  if p is None:
    return "-"
  return f"{p['region']} (" + ", ".join(f"{a:.4g}" for a in p["args"]) + ")"


def summary(dev, pr, labels):
  print("Columns: worst absolute error (relative error at that point);")
  print("max rel: the worst relative error over all points.\n")
  print("| function | tree | points | throw | NaN | inf | worst value "
        "(rel) | worst log-scale grad (rel) | max rel value | max rel grad |")
  print("|---|---|---|---|---|---|---|---|---|---|")
  for row, keys in ROWS:
    for side, pts in zip(labels, (dev, pr)):
      sub = [p for p in pts if p["fn"] in keys]
      thr, nan, inf = counts(sub)
      wv = worst(sub, "verr")
      wg = worst(sub, "gerr")
      v = "-" if wv is None else f"{fmt(wv['verr'])} ({fmt(wv['vrel'])})"
      g = "-" if wg is None else f"{fmt(wg['gerr'])} ({fmt(wg['grel'])})"
      rv = worst(sub, "vrel")
      rg = worst(sub, "grel")
      rv = "-" if rv is None else fmt(rv["vrel"])
      rg = "-" if rg is None else fmt(rg["grel"])
      print(f"| {row} | {side} | {len(sub)} | {thr} | {nan} | {inf} | {v} "
            f"| {g} | {rv} | {rg} |")
  print("\nWorst points (region, arguments):\n")
  for row, keys in ROWS:
    for side, pts in zip(labels, (dev, pr)):
      sub = [p for p in pts if p["fn"] in keys]
      print(f"- {row}, {side}: value at {point_str(worst(sub, 'verr'))}; "
            f"grad at {point_str(worst(sub, 'gerr'))}")


def per_region(dev, pr, labels):
  print("\nPer region: worst value and log-scale gradient error\n")
  print("| function | region | n | value dev | value pr | grad dev "
        "| grad pr | thr/nonfin dev | thr/nonfin pr |")
  print("|---|---|---|---|---|---|---|---|---|")
  for row, keys in ROWS:
    regions = OrderedDict()
    for p in dev:
      if p["fn"] in keys:
        regions.setdefault(p["region"], None)
    for reg in regions:
      cells = []
      subs = [[p for p in pts if p["fn"] in keys and p["region"] == reg]
              for pts in (dev, pr)]
      for k in ("verr", "gerr"):
        for sub in subs:
          w = worst(sub, k)
          cells.append(fmt(None if w is None else w[k]))
      nf = []
      for sub in subs:
        thr, nan, inf = counts(sub)
        nf.append(f"{thr}/{nan + inf}")
      print(f"| {row} | {reg} | {len(subs[0])} | " + " | ".join(cells)
            + " | " + " | ".join(nf) + " |")


def main():
  if len(sys.argv) < 3:
    print(__doc__)
    return 2
  labels = sys.argv[1:3]
  sib_refs = S.load_refs(os.path.join(ROOT, "refs", "refs_sib.txt"))
  st_refs = load_st_refs(os.path.join(ROOT, "refs", "refs_st.txt"))
  bb_refs = load_bb_refs(os.path.join(ROOT, "refs", "refs_bb.txt"))
  sides = []
  for label in labels:
    sides.append(load_bb(label, bb_refs) + load_sib(label, sib_refs)
                 + load_st(label, st_refs))
  summary(sides[0], sides[1], labels)
  if "--regions" in sys.argv:
    per_region(sides[0], sides[1], labels)
  return 0


if __name__ == "__main__":
  sys.exit(main())
