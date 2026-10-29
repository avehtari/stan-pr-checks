#!/usr/bin/env python3
"""Score probe_bb output against the mpmath reference (stan-dev/math#3153).

Plan: ../2026-10-06_PLAN_fix_issue3153_beta_binomial_lpmf.md, steps 5.3, 5.4.

Usage: score_bb.py <probe output> [<second probe output>]

With two files, also print the largest change between them per region
(plan section 4, third item: no change where develop is already accurate).

References are cached in refs_bb.txt (keyed by n, N, alpha hex, beta hex).

Metrics per region (all absolute unless named rel):
  val      |value - ref|
  val/sc   |value - ref| / (1 + |lchoose| + |ref|): the unavoidable rounding
           scale; lchoose(1e6, 5e5) is 7e5, so 1e-12 absolute is not
           reachable there
  dA rel   max over d/dalpha and d/dbeta of the relative error
  gab      max over alpha * d/dalpha and beta * d/dbeta (log-shape gradients)
  gc       alpha * d/dalpha + beta * d/dbeta, the log-concentration gradient
           (d/dlc of the issue for tied shapes)
  gc rel   relative error of gc where |gc ref| > 1e-6 (plan target 1e-8)
"""

import math
import os
import sys

from mpmath import mp, mpf

sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(
    os.path.abspath(__file__))), "refs"))
import ref_bb as R  # noqa: E402

CACHE = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(
    __file__))), "refs", "refs_bb.txt")
NSTR = 45


def load_cache():
  cache = {}
  if not os.path.exists(CACHE):
    return cache
  with open(CACHE) as fh:
    for line in fh:
      parts = line.split()
      key = tuple(parts[:4])
      vals = [mpf(x) for x in parts[4:]]
      cache[key] = dict(zip(["lchoose", "value", "d_alpha", "d_beta"], vals))
  return cache


def save_cache(cache, new_keys):
  if not new_keys:
    return
  with open(CACHE, "a") as fh:
    for key in new_keys:
      r = cache[key]
      vals = [mp.nstr(r[k], NSTR) for k in ["lchoose", "value", "d_alpha",
                                            "d_beta"]]
      fh.write(" ".join(list(key) + vals) + "\n")


def parse(path):
  rows = []
  with open(path) as fh:
    for line in fh:
      if line.startswith("#") or not line.strip():
        continue
      left, right = line.split("=>")
      tag, n, N, a, b = left.split()
      vals = right.split()
      status = vals[-1]
      nums = [float.fromhex(x) if x not in ("nan", "-nan") else math.nan
              for x in vals[:-1]]
      rows.append({"tag": tag, "key": (n, N, a, b), "n": int(n),
                   "N": int(N), "alpha": float.fromhex(a),
                   "beta": float.fromhex(b), "v_dbl": nums[0],
                   "v_var": nums[1], "d_alpha": nums[2], "d_beta": nums[3],
                   "status": status})
  return rows


def refs_for(rows, cache):
  new_keys = []
  for r in rows:
    if r["key"] in cache:
      continue
    ref = R.ref_point(r["n"], r["N"], r["alpha"], r["beta"])
    cache[r["key"]] = {k: ref[k] for k in ["lchoose", "value", "d_alpha",
                                           "d_beta"]}
    new_keys.append(r["key"])
  save_cache(cache, new_keys)


def _dpsi(x, k):
  return mp.digamma(x + k) - mp.digamma(x)


def errors(r, ref):
  a, b = mpf(r["alpha"]), mpf(r["beta"])
  n, N = r["n"], r["N"]
  out = {}
  out["val"] = abs(mpf(r["v_var"]) - ref["value"])
  out["val_dbl"] = abs(mpf(r["v_dbl"]) - ref["value"])
  out["val/sc"] = out["val"] / (1 + abs(ref["lchoose"]) + abs(ref["value"]))
  rels = []
  for k in ["d_alpha", "d_beta"]:
    e = abs(mpf(r[k]) - ref[k])
    rels.append(e / abs(ref[k]) if ref[k] != 0 else e)
  out["dA rel"] = max(rels)
  # scale of the terms that cancel inside each partial: an error of a few
  # eps times this is the input-rounding limit (backward stable)
  tail = abs(_dpsi(a + b, N))
  t_a = abs(_dpsi(a, n)) + tail
  t_b = abs(_dpsi(b, N - n)) + tail
  e_a = abs(mpf(r["d_alpha"]) - ref["d_alpha"])
  e_b = abs(mpf(r["d_beta"]) - ref["d_beta"])
  out["dA/sc"] = max(e_a / t_a, e_b / t_b)
  ga_ref, gb_ref = a * ref["d_alpha"], b * ref["d_beta"]
  ga, gb = a * mpf(r["d_alpha"]), b * mpf(r["d_beta"])
  out["gab"] = max(abs(ga - ga_ref), abs(gb - gb_ref))
  gc_ref = ga_ref + gb_ref
  out["gc"] = abs((ga + gb) - gc_ref)
  out["gc/sc"] = out["gc"] / (a * t_a + b * t_b)
  out["gc rel"] = (out["gc"] / abs(gc_ref)) if abs(gc_ref) > mpf("1e-6") \
      else mpf(0)
  out["gc_ref"] = gc_ref
  return out


def fmt(x):
  x = float(x)
  if math.isnan(x):
    return "nan"
  if math.isinf(x):
    return "inf"
  return f"{x:8.1e}"


COLS = ["val", "val/sc", "dA/sc", "gab", "gc", "gc/sc", "gc rel"]


def by_n(rows):
  """Worst absolute value and gc error per N, all regions except REP."""
  print("\nper N (all regions except REP):")
  print(f"{'N':>8s} {'val':>8s} {'val/sc':>8s} {'gc':>8s} {'gc/sc':>8s}")
  for N in sorted({r["N"] for r in rows}):
    sub = [r["err"] for r in rows if r["N"] == N and r["tag"] != "REP"]
    if not sub:
      continue
    print(f"{N:8d} " + " ".join(
        f"{fmt(max(e[c] for e in sub)):>8s}"
        for c in ["val", "val/sc", "gc", "gc/sc"]))


def score(rows, cache, label):
  print(f"\n== {label}: {len(rows)} points, "
        f"{sum(1 for r in rows if 'T' in r['status'])} threw")
  regions = []
  for r in rows:
    if r["tag"] not in regions:
      regions.append(r["tag"])
  head = f"{'region':6s} " + " ".join(f"{c:>8s}" for c in COLS) + \
      f" {'#v>1e-12':>8s} {'#g>1e-10':>8s}"
  print(head)
  worst = {}
  for tag in regions:
    sub = [r for r in rows if r["tag"] == tag]
    mx = {c: mpf(0) for c in COLS}
    arg = {c: None for c in COLS}
    nbad_val = 0
    nbad_gc = 0
    for r in sub:
      e = errors(r, cache[r["key"]])
      r["err"] = e
      if mp.isnan(e["val"]) or e["val"] > mpf("1e-12"):
        nbad_val += 1
      if mp.isnan(e["gc"]) or e["gc"] > mpf("1e-10"):
        nbad_gc += 1
      for c in COLS:
        if mp.isnan(e[c]) or e[c] > mx[c]:
          if not mp.isnan(mx[c]):
            mx[c] = e[c]
            arg[c] = r
    worst[tag] = arg
    print(f"{tag:6s} " + " ".join(f"{fmt(mx[c]):>8s}" for c in COLS)
          + f" {nbad_val:8d} {nbad_gc:8d}")
  return worst


def print_worst(worst, cols=("val", "gc")):
  print("\nworst points (n, N, alpha, beta):")
  for tag, arg in worst.items():
    parts = []
    for c in cols:
      r = arg[c]
      if r is None:
        continue
      parts.append(f"{c}: ({r['n']}, {r['N']}, {r['alpha']:.3g}, "
                   f"{r['beta']:.3g})")
    print(f"  {tag:6s} " + "  ".join(parts))


def rep_table(rows, cache):
  """The reporter's table: difference to the binomial limit and d/dlc."""
  sub = [r for r in rows if r["tag"] == "REP"]
  if not sub:
    return
  print("\nREP: n = 57, N = 117, alpha = beta = exp(lc) / 2")
  print(f"{'lc':>5s} {'stan diff':>12s} {'true diff':>12s} "
        f"{'stan d/dlc':>12s} {'true d/dlc':>12s} {'gc rel':>8s}")
  log2 = mp.log(2)
  for i, r in enumerate(sub):
    lc = 0.25 * i
    if lc < 16 or (lc % 2) != 0:
      continue
    ref = cache[r["key"]]
    limit = ref["lchoose"] - 117 * log2
    s = mpf(r["alpha"])
    stan_gc = s * (mpf(r["d_alpha"]) + mpf(r["d_beta"]))
    true_gc = s * (ref["d_alpha"] + ref["d_beta"])
    rel = abs(stan_gc - true_gc) / abs(true_gc)
    print(f"{lc:5.1f} {float(mpf(r['v_var']) - limit):12.4e} "
          f"{float(ref['value'] - limit):12.4e} {float(stan_gc):12.4e} "
          f"{float(true_gc):12.4e} {float(rel):8.1e}")
  # bounds over every REP point, lc = 0, 0.25, ..., 44
  worst_v = (mpf(0), None)
  worst_g = (mpf(0), None)
  for i, r in enumerate(sub):
    ref = cache[r["key"]]
    s = mpf(r["alpha"])
    ev = abs(mpf(r["v_var"]) - ref["value"])
    eg = abs(s * (mpf(r["d_alpha"]) + mpf(r["d_beta"]))
             - s * (ref["d_alpha"] + ref["d_beta"]))
    if mp.isnan(ev) or ev > worst_v[0]:
      worst_v = (ev, 0.25 * i)
    if mp.isnan(eg) or eg > worst_g[0]:
      worst_g = (eg, 0.25 * i)
  print(f"REP bounds over {len(sub)} points, lc 0 to 44: value error "
        f"{float(worst_v[0]):.2e} (lc {worst_v[1]}), d/dlc error "
        f"{float(worst_g[0]):.2e} (lc {worst_g[1]})")


def compare(rows1, rows2):
  """Largest change between two probe outputs, per region."""
  print("\nchange between the two outputs (value abs, partial rel):")
  by_key = {(r["tag"], r["key"]): r for r in rows2}
  regions = []
  for r in rows1:
    if r["tag"] not in regions:
      regions.append(r["tag"])
  for tag in regions:
    dv = 0.0
    dg = 0.0
    for r in rows1:
      if r["tag"] != tag:
        continue
      s = by_key[(tag, r["key"])]
      dv = max(dv, abs(r["v_var"] - s["v_var"]))
      for k in ["d_alpha", "d_beta"]:
        den = max(abs(r[k]), 1e-300)
        dg = max(dg, abs(r[k] - s[k]) / den)
    print(f"  {tag:6s} value {dv:8.1e}  partial rel {dg:8.1e}")


def main():
  if len(sys.argv) < 2:
    print(__doc__)
    return 2
  cache = load_cache()
  outputs = []
  for path in sys.argv[1:]:
    rows = parse(path)
    refs_for(rows, cache)
    outputs.append(rows)
    worst = score(rows, cache, os.path.basename(path))
    print_worst(worst)
    by_n(rows)
    rep_table(rows, cache)
  if len(outputs) == 2:
    compare(outputs[0], outputs[1])
  return 0


if __name__ == "__main__":
  sys.exit(main())
