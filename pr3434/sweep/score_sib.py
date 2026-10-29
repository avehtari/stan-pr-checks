#!/usr/bin/env python3
"""Score probe_sib output against refs.txt (ref_sib.py).

Usage: score_sib.py <refs file> <probe output> [<probe output> ...]

Prints markdown: the route certification table, one table per tag (rows =
config:decade), and the per-tag summary.

Metrics per row (worst over the points of the row):
  val abs   |stan - ref|
  val rel   |stan - ref| / |ref|
  lsg abs   max over shape-type x of |x * (d/dx stan - d/dx ref)|
            (log-shape gradient; DM also the kappa direction sum a_i d/da_i)
  lsg rel   the same divided by |x * d/dx ref|, only where that is > 1e-6
  oth       other gradients (mu, eta, GLM alpha and beta):
            |stan - ref| / max(1, |ref|)
  thr       calls that threw; a non-finite result without a throw counts as
            an infinite error
"""

import math
import sys
from collections import OrderedDict

from mpmath import mp, mpf

mp.dps = 60

# (gradient index, argument index) of each shape-type parameter
SHAPE = {
    "NB": [(0, 1), (1, 2)], "NB2": [(1, 2)], "NB2L": [(1, 2)],
    "GLM1": [(2, 4)], "GLM5": [(2, 12)],
    "BNB": [(0, 1), (1, 2), (2, 3)], "DM": [(0, 3), (1, 4), (2, 5)],
    "LKJ": [(0, 3)], "LKJC": [(0, 5)],
    "BBC": [(0, 2), (1, 3)], "BBLC": [(0, 2), (1, 3)],
    "BBLCC": [(0, 2), (1, 3)],
    "YSC": [(0, 1)], "YSLC": [(0, 1)], "YSLCC": [(0, 1)],
    "BNBC": [(0, 1), (1, 2), (2, 3)], "BNBLC": [(0, 1), (1, 2), (2, 3)],
    "BNBLCC": [(0, 1), (1, 2), (2, 3)],
}
OTHER = {"NB2": [0], "NB2L": [0], "GLM1": [0, 1], "GLM5": [0, 1]}

TITLE = OrderedDict([
    ("NB", "1. neg_binomial_lpmf(n | alpha, beta), beta = alpha / mu"),
    ("NB2", "2. neg_binomial_2_lpmf(n | mu, phi)"),
    ("NB2L", "3. neg_binomial_2_log_lpmf(n | log(mu), phi)"),
    ("GLM1", "4a. neg_binomial_2_log_glm_lpmf, one observation"),
    ("GLM5", "4b. neg_binomial_2_log_glm_lpmf, five observations"),
    ("BNB", "5. beta_neg_binomial_lpmf(n | r, alpha, beta)"),
    ("DM", "6. dirichlet_multinomial_lpmf(ns | kappa * (0.2, 0.3, 0.5))"),
    ("LKJ", "7a. lkj_corr_lpdf(Omega | eta), K = 3"),
    ("LKJC", "7b. lkj_corr_cholesky_lpdf(L | eta), K = 3"),
    ("BBC", "8a. beta_binomial_cdf(n | 117, alpha, beta)"),
    ("BBLC", "8b. beta_binomial_lcdf(n | 117, alpha, beta)"),
    ("BBLCC", "8c. beta_binomial_lccdf(n | 117, alpha, beta)"),
    ("YSC", "9a. yule_simon_cdf(n | alpha)"),
    ("YSLC", "9b. yule_simon_lcdf(n | alpha)"),
    ("YSLCC", "9c. yule_simon_lccdf(n | alpha)"),
    ("BNBC", "10a. beta_neg_binomial_cdf(n | r, alpha, beta)"),
    ("BNBLC", "10b. beta_neg_binomial_lcdf(n | r, alpha, beta)"),
    ("BNBLCC", "10c. beta_neg_binomial_lccdf(n | r, alpha, beta)"),
])

DEFECT = mpf("1e-8")


def num(tok):
  if tok in ("none",):
    return None
  return mpf(tok)


def norm_key(text):
  """Key by values: C printf %a and Python float.hex spell hex differently."""
  parts = text.split()
  return tuple(parts[:2]) + tuple(float.fromhex(h) for h in parts[2:])


def load_refs(path):
  refs = {}
  with open(path) as fh:
    for line in fh:
      key, rest = line.rstrip("\n").split(" => ")
      key = norm_key(key)
      a, b, c = rest.split(" | ")
      refs[key] = {"A": [num(x) for x in a.split()],
                   "B": [num(x) for x in b.split()], "C": num(c)}
  return refs


def hexnum(tok):
  try:
    return float.fromhex(tok)
  except ValueError:
    return math.nan


def load_probe(path):
  rows = []
  with open(path) as fh:
    for line in fh:
      left, right = line.split(" => ")
      parts = left.split()
      vals = right.split()
      rows.append({"key": norm_key(left), "tag": parts[0], "label": parts[1],
                   "args": [float.fromhex(h) for h in parts[2:]],
                   "value": hexnum(vals[0]),
                   "grads": [hexnum(v) for v in vals[1:-2]],
                   "status": vals[-2], "ms": float(vals[-1])})
  return rows


def rel(a, b):
  """Relative difference; absolute where both are below 1e-60.

  mp.diff at 140 digits cannot resolve a derivative far below the value
  (YSC at ccdf 1e-800, an exact zero of d/dmu at n = mu), so a gradient
  below 1e-60 is compared in absolute terms.
  """
  scale = max(abs(a), abs(b))
  if scale < mpf("1e-60"):
    return abs(a - b)
  return abs(a - b) / scale


def certify(refs):
  """Worst route A vs B (value, gradients) and A vs C (value) per tag."""
  out = OrderedDict()
  for key, r in refs.items():
    tag = key[0]
    w = out.setdefault(tag, [mpf(0), mpf(0), None, 0])
    w[0] = max(w[0], rel(r["A"][0], r["B"][0]))
    for x, y in zip(r["A"][1:], r["B"][1:]):
      w[1] = max(w[1], rel(x, y))
    if r["C"] is not None:
      d = rel(r["A"][0], r["C"])
      w[2] = d if w[2] is None else max(w[2], d)
      w[3] += 1
  return out


def finite(x):
  return not (math.isnan(x) or math.isinf(x))


def errors(row, ref):
  """Error dictionary of one point, or None if the call threw."""
  if row["status"] == "T":
    return None
  inf = mpf("inf")
  e = {}
  vr = ref["A"][0]
  if finite(row["value"]):
    e["val_abs"] = abs(mpf(row["value"]) - vr)
    e["val_rel"] = e["val_abs"] / abs(vr) if vr != 0 else e["val_abs"]
  else:
    e["val_abs"] = e["val_rel"] = inf
  gref = ref["A"][1:]
  lsg_abs = mpf(0)
  lsg_rel = mpf(0)
  pairs = []
  for gi, ai in SHAPE[row["tag"]]:
    x = mpf(row["args"][ai])
    pairs.append((x * mpf(row["grads"][gi]) if finite(row["grads"][gi])
                  else None, x * gref[gi]))
  if row["tag"] == "DM":
    if all(p[0] is not None for p in pairs):
      pairs.append((sum(p[0] for p in pairs), sum(p[1] for p in pairs)))
    else:
      pairs.append((None, sum(p[1] for p in pairs)))
  for s, t in pairs:
    if s is None:
      lsg_abs = lsg_rel = inf
      continue
    d = abs(s - t)
    lsg_abs = max(lsg_abs, d)
    if abs(t) > mpf("1e-6"):
      lsg_rel = max(lsg_rel, d / abs(t))
  e["lsg_abs"] = lsg_abs
  e["lsg_rel"] = lsg_rel
  oth = None
  for gi in OTHER.get(row["tag"], []):
    g = row["grads"][gi]
    d = (abs(mpf(g) - gref[gi]) / max(1, abs(gref[gi])) if finite(g)
         else inf)
    oth = d if oth is None else max(oth, d)
  e["oth"] = oth
  return e


def fmt(x):
  if x is None:
    return "-"
  x = float(x)
  if math.isinf(x):
    return "inf"
  if x == 0:
    return "0"
  return f"{x:.1e}"


COLS = ["val_abs", "val_rel", "lsg_abs", "lsg_rel", "oth"]
HEAD = ["val abs", "val rel", "lsg abs", "lsg rel", "oth"]


def decade_of(label):
  return label.split(":")[1]


def main():
  if len(sys.argv) < 3:
    print(__doc__)
    return 2
  refs = load_refs(sys.argv[1])
  rows = []
  for path in sys.argv[2:]:
    rows += load_probe(path)

  print("## Reference certification\n")
  print("Worst relative disagreement per tag: route A (analytic, 90 digits)")
  print("against route B (140 digits, mp.diff), and A against C (second")
  print("formula, value only; count of points in parentheses).\n")
  print("| tag | value A-B | grad A-B | value A-C |")
  print("|---|---|---|---|")
  for tag, w in certify(refs).items():
    c = "-" if w[2] is None else f"{fmt(w[2])} ({w[3]})"
    print(f"| {tag} | {fmt(w[0])} | {fmt(w[1])} | {c} |")

  summary = []
  for tag, title in TITLE.items():
    sub = [r for r in rows if r["tag"] == tag]
    if not sub:
      continue
    print(f"\n## {title}\n")
    print("| config:decade | pts | thr | " + " | ".join(HEAD)
          + " | ms max |")
    print("|---|---|---|" + "---|" * len(HEAD) + "---|")
    groups = OrderedDict()
    cfg_order = []
    for r in sub:
      cfg = r["label"].split(":")[0]
      if cfg not in cfg_order:
        cfg_order.append(cfg)
    for r in sorted(sub, key=lambda r: (
        cfg_order.index(r["label"].split(":")[0]),
        float(decade_of(r["label"])))):
      groups.setdefault(r["label"], []).append(r)
    val_from = None
    grad_from = None
    worst_val = (mpf(0), None)
    worst_lsg = (mpf(0), None)
    n_thr = 0
    for label, grp in groups.items():
      mx = {c: None for c in COLS}
      thr = 0
      for r in grp:
        e = errors(r, refs[r["key"]])
        if e is None:
          thr += 1
          continue
        for c in COLS:
          if e[c] is not None and (mx[c] is None or e[c] > mx[c]):
            mx[c] = e[c]
        if e["val_abs"] > worst_val[0]:
          worst_val = (e["val_abs"], r)
        if e["lsg_abs"] > worst_lsg[0]:
          worst_lsg = (e["lsg_abs"], r)
        dec = decade_of(label)
        if (e["val_abs"] > DEFECT and e["val_rel"] > DEFECT):
          if val_from is None or float(dec) < float(val_from):
            val_from = dec
        if e["lsg_abs"] > DEFECT and (e["lsg_rel"] > DEFECT
                                      or e["lsg_rel"] == 0):
          if grad_from is None or float(dec) < float(grad_from):
            grad_from = dec
      n_thr += thr
      ms = max(r["ms"] for r in grp)
      print(f"| {label} | {len(grp)} | {thr} | "
            + " | ".join(fmt(mx[c]) for c in COLS) + f" | {ms:.0f} |")
    summary.append((tag, val_from, grad_from, worst_val, worst_lsg, n_thr,
                    len(sub)))

  print("\n## Summary (defect: abs and rel error both > 1e-8)\n")
  print("| tag | value defect from | grad defect from | worst val abs "
        "| worst lsg abs | thr/pts |")
  print("|---|---|---|---|---|---|")
  for tag, vf, gf, wv, wg, thr, n in summary:
    print(f"| {tag} | {vf or '-'} | {gf or '-'} | {fmt(wv[0])} "
          f"| {fmt(wg[0])} | {thr}/{n} |")
  print("\n## Worst points\n")
  for tag, vf, gf, wv, wg, thr, n in summary:
    for name, w in [("val", wv), ("lsg", wg)]:
      if w[1] is None:
        continue
      args = ", ".join(f"{a:.4g}" for a in w[1]["args"])
      print(f"- {tag} {name} {fmt(w[0])} at {w[1]['label']} ({args})")
  return 0


if __name__ == "__main__":
  sys.exit(main())
