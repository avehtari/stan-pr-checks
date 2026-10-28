#!/usr/bin/env python3
"""Score the Level 0 develop-against-branch comparison for the beta roots.

Usage:
  python3 score_cmp_beta.py double cmp_beta_double_develop.txt \\
                                   cmp_beta_double_branch.txt refs_beta_18000.txt
  python3 score_cmp_beta.py all    cmp_beta_all_develop.txt \\
                                   cmp_beta_all_branch.txt refs_beta_18000.txt

Columns of the probe output: mode=double has A_da A_db B_g1 B_g2; mode=all
has the same four, each in the four configurations double, var,
fvar<double>, fvar<fvar<double>>. The status field marks a throw with 'T'.

Metrics per region, see summarize_beta.py: nan (non-finite or throw),
e>1e-6 and e>1e-2 (err_log = |v - r| / I), e worst, rel* (relative error
where |r| / I >= 1e-3). In mode=all a final table gives the largest spread
between the four configurations per region, normal-range points only.
"""

import math
import sys

DBL_MIN = 2.2250738585072014e-308
CONFIGS = ["double", "var", "fvar", "ffvar"]
PATHS = ["A_da", "A_db", "B_g1", "B_g2"]


def parse(path):
  rows = []
  region = "?"
  for line in open(path):
    line = line.strip()
    if not line:
      continue
    if line.startswith("#"):
      if line.startswith("# probe") or line.startswith("# evaluated"):
        continue
      region = line[2:]
      continue
    f = line.split()
    if not f[0].startswith(("0x", "-0x")):
      continue  # the /usr/bin/time line, or anything else that is not data
    key = tuple(float.fromhex(x) for x in f[:3])
    status = f[-1]
    vals = [float.fromhex(x) for x in f[3:-1]]
    rows.append((region, key, vals, status))
  return rows


def load_refs(path):
  refs = {}
  for line in open(path):
    f = line.split()
    refs[tuple(float.fromhex(x) for x in f[:3])] = tuple(
        float.fromhex(x) for x in f[3:6])
  return refs


def table(title, rows, refs, col, ri):
  print(f"\n## {title}\n")
  print(f"{'region':<26} {'n':>5} {'nan':>5} {'e>1e-6':>6} {'e>1e-2':>6} "
        f"{'e worst':>8} {'rel*':>8}  worst err_log at")
  by = {}
  for reg, key, vals, status in rows:
    by.setdefault(reg, []).append((key, vals, status))
  for reg, items in by.items():
    n = nan = big6 = big2 = 0
    eworst = (0.0, None)
    relw = 0.0
    for key, vals, status in items:
      i, r = refs[key][0], refs[key][ri]
      if math.isnan(r) or i <= 0.0:
        continue
      n += 1
      v = vals[col]
      if not math.isfinite(v) or status[col] == "T":
        nan += 1
        continue
      e = abs(v - r) / i
      if e > 1e-6:
        big6 += 1
      if e > 1e-2:
        big2 += 1
      if e > eworst[0]:
        eworst = (e, (key, v, r, i))
      if abs(r) / i >= 1e-3 and abs(r) >= DBL_MIN:
        relw = max(relw, abs(v - r) / abs(r))
    w = ""
    if eworst[1] is not None:
      key, v, r, i = eworst[1]
      w = (f"({key[0]:.4g}, {key[1]:.4g}, {key[2]:.4g}) got {v:.3g} "
           f"ref {r:.3g} I {i:.2g}")
    print(f"{reg:<26} {n:>5d} {nan:>5d} {big6:>6d} {big2:>6d} "
          f"{eworst[0]:>8.1e} {relw:>8.1e}  {w}")


def spread_table(title, rows, refs, path_index, nconf):
  print(f"\n## {title}: largest spread between configurations\n")
  print(f"{'region':<26} {'n':>5} " + " ".join(f"{p:>9}" for p in PATHS))
  by = {}
  for reg, key, vals, status in rows:
    by.setdefault(reg, []).append((key, vals, status))
  for reg, items in by.items():
    out = []
    n = 0
    for p in range(4):
      worst = 0.0
      for key, vals, status in items:
        i = refs[key][0]
        if math.isnan(i) or i <= 0.0:
          continue
        if p == 0:
          n += 1
        cols = vals[p * nconf:(p + 1) * nconf]
        fin = [v for v in cols if math.isfinite(v)]
        if len(fin) < nconf:
          if fin:
            worst = math.inf
          continue
        m = max(abs(v) for v in fin)
        if m == 0.0:
          continue
        worst = max(worst, (max(fin) - min(fin)) / m)
      out.append(f"{worst:>9.1e}")
    print(f"{reg:<26} {n:>5d} " + " ".join(out))


def main():
  mode = sys.argv[1]
  dev = parse(sys.argv[2])
  bra = parse(sys.argv[3])
  refs = load_refs(sys.argv[4])
  assert len(dev) == len(bra), "row counts differ"
  for d, b in zip(dev, bra):
    assert d[1] == b[1], "point files are not aligned"
  nconf = 4 if mode == "all" else 1
  print(f"# Level 0, develop against branch, mode={mode}, {len(dev)} points")
  print("# reference: refs_beta_18000.txt (density-form quadrature, 80 "
        "digits, two routes)")
  for p, pname in enumerate(PATHS):
    ri = 1 if pname.endswith("da") or pname.endswith("g1") else 2
    for c in range(nconf):
      col = p * nconf + c
      for tag, rows in (("develop", dev), ("branch", bra)):
        table(f"{pname} {CONFIGS[c]}, {tag}", rows, refs, col, ri)
  if nconf == 4:
    spread_table("develop", dev, refs, 0, nconf)
    spread_table("branch", bra, refs, 0, nconf)
  return 0


if __name__ == "__main__":
  sys.exit(main())
