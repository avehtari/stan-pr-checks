#!/usr/bin/env python3
# Scores probe_cmp_roots output of two trees against mpmath (cached refs).
# Supports 6 379 / 0 non-finite, 3e-11, 2.9e-7, 5.2e-6 and the AD spread.
"""Score the Level 0 develop-against-branch comparison.

Usage:
  python3 score_cmp_roots.py double cmp_roots_double_develop.txt \\
                                    cmp_roots_double_branch.txt [nproc]
  python3 score_cmp_roots.py all    cmp_roots_all_develop.txt \\
                                    cmp_roots_all_branch.txt [nproc]

The mpmath references are computed once per point file and cached in
`refs_cache_<n>.txt` (a, z, dQ/da in hex), so a rescoring is instant.
Any refs_cache_*.txt file in the working directory or next to this
script that holds every point is used, so the 2 000-point subset reads
its references from refs_cache_49990.txt.

Three classes of point are reported separately, because they answer
different questions:

  normal     |ref| >= 2.2e-308. The ulp question lives here.
  subnormal  4.9e-324 <= |ref| < 2.2e-308. A result of 0 here is flush to
             zero, not an accuracy defect; reported as a count and as the
             largest |ref| that was flushed, which is the underflow
             threshold of axis (a).
  zero       |ref| below the double range. Only "did it return 0" matters.
"""

import glob
import math
import multiprocessing as mputil
import os
import sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                "..", "refs"))

import ref_mpmath as R  # noqa: E402

UPPER = "grad_reg_inc_gamma (dQ/da)"
LOWER = "grad_reg_lower_inc_gamma (dP/da)"
DBL_MIN = 2.2250738585072014e-308
DBL_TRUE_MIN = 4.9406564584124654e-324


def rel(value, reference):
  if not math.isfinite(value):
    return math.inf
  if reference == 0.0:
    return 0.0 if value == 0.0 else math.inf
  return abs(value - reference) / abs(reference)


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
    parts = line.split()
    if len(parts) < 4:
      continue
    nums = [float.fromhex(p) for p in parts[:-1]]
    rows.append((region, nums[0], nums[1], nums[2:], parts[-1]))
  return rows


def ref_one(az):
  a, z = az
  try:
    return float(R.dQda(a, z))
  except Exception:
    return math.nan


def load_or_compute_refs(points, nproc):
  cache = f"refs_cache_{len(points)}.txt"
  here = os.path.dirname(os.path.abspath(__file__))
  candidates = [cache] + sorted(glob.glob("refs_cache_*.txt")) + sorted(
      glob.glob(os.path.join(here, "refs_cache_*.txt")))
  for cand in candidates:
    if not os.path.exists(cand):
      continue
    refs = {}
    for line in open(cand):
      a, z, r = line.split()
      refs[(float.fromhex(a), float.fromhex(z))] = float.fromhex(r)
    if all(p in refs for p in points):
      print(f"# references loaded from {cand}")
      return [refs[p] for p in points]
  with mputil.Pool(nproc) as pool:
    vals = pool.map(ref_one, points, chunksize=64)
  with open(cache, "w") as fh:
    for (a, z), r in zip(points, vals):
      fh.write(f"{a.hex()} {z.hex()} {r.hex()}\n")
  print(f"# references computed and cached in {cache}")
  return vals


def summarize(vals):
  finite = sorted(v for v in vals if math.isfinite(v))
  if not finite:
    return math.nan, math.nan
  return finite[len(finite) // 2], max(vals) if vals else math.nan


def main():
  mode = sys.argv[1]
  dev = parse(sys.argv[2])
  bra = parse(sys.argv[3])
  nproc = int(sys.argv[4]) if len(sys.argv) > 4 else 6
  assert len(dev) == len(bra), "row counts differ"
  points = []
  for d, b in zip(dev, bra):
    assert d[1] == b[1] and d[2] == b[2], "point files are not aligned"
    points.append((d[1], d[2]))
  refs_q = load_or_compute_refs(points, nproc)

  ncol = len(dev[0][3])
  half = ncol // 2
  config_names = (["double", "var", "fvar", "ffvar"] if half == 4
                  else ["double"])

  print(f"# Level 0, develop against branch, mode={mode}, {len(points)} points")
  print("# reference: mpmath dQ/da, two routes agreeing to 57 digits; "
        "dP/da = -dQ/da\n")

  for fn_index, fn_name in ((0, UPPER), (half, LOWER)):
    sign = 1.0 if fn_index == 0 else -1.0
    for ci, cname in enumerate(config_names):
      col = fn_index + ci
      print(f"\n## {fn_name}, configuration {cname}\n")
      print(f"{'region':<24} {'n':>5} {'dev med':>9} {'dev worst':>9} "
            f"{'bra med':>9} {'bra worst':>9} {'dev nf':>6} {'bra nf':>6} "
            f"{'subn':>5} {'dev 0':>5} {'bra 0':>5}")
      regions = {}
      for d, b, rq in zip(dev, bra, refs_q):
        regions.setdefault(d[0], []).append((d, b, sign * rq))
      flush_dev, flush_bra = 0.0, 0.0
      worst_dev, worst_bra = (0.0, None), (0.0, None)
      for region, rows in regions.items():
        dev_rel, bra_rel = [], []
        dnf = bnf = subn = d0 = b0 = 0
        for d, b, r in rows:
          vd, vb = d[3][col], b[3][col]
          if math.isnan(r):
            continue
          if abs(r) < DBL_MIN:
            subn += 1
            if vd == 0.0:
              d0 += 1
              flush_dev = max(flush_dev, abs(r))
            if vb == 0.0:
              b0 += 1
              flush_bra = max(flush_bra, abs(r))
            continue
          if not math.isfinite(vd):
            dnf += 1
          if not math.isfinite(vb):
            bnf += 1
          e_d, e_b = rel(vd, r), rel(vb, r)
          dev_rel.append(e_d)
          bra_rel.append(e_b)
          if math.isfinite(e_d) and e_d > worst_dev[0]:
            worst_dev = (e_d, (d[1], d[2], vd, r))
          if math.isfinite(e_b) and e_b > worst_bra[0]:
            worst_bra = (e_b, (d[1], d[2], vb, r))
        dm, dw = summarize(dev_rel)
        bm, bw = summarize(bra_rel)
        print(f"{region:<24} {len(rows):>5d} {dm:>9.2e} {dw:>9.2e} "
              f"{bm:>9.2e} {bw:>9.2e} {dnf:>6d} {bnf:>6d} {subn:>5d} "
              f"{d0:>5d} {b0:>5d}")
      print(f"  flush-to-zero threshold, largest |ref| returned as 0: "
            f"develop {flush_dev:.3e}, branch {flush_bra:.3e}")
      for tag, w in (("develop", worst_dev), ("branch", worst_bra)):
        if w[1] is not None:
          a, z, v, r = w[1]
          print(f"  worst finite {tag:<8} rel={w[0]:.3e} at a={a.hex()} "
                f"z={z.hex()}  (a={a:.8g} z={z:.8g} got={v:.6g} "
                f"ref={r:.6g})")

  if half == 4:
    print("\n\n## spread between the four autodiff configurations "
          "(normal region only)\n")
    print(f"{'region':<24} {'dev up':>10} {'bra up':>10} {'dev lo':>10} "
          f"{'bra lo':>10}")
    regions = {}
    for d, b, rq in zip(dev, bra, refs_q):
      regions.setdefault(d[0], []).append((d, b, rq))

    def spread(vals):
      fin = [v for v in vals if math.isfinite(v)]
      if len(fin) < 2:
        return math.inf if len(vals) > len(fin) else 0.0
      m = max(abs(v) for v in fin)
      return (max(fin) - min(fin)) / m if m > 0 else 0.0

    for region, rows in regions.items():
      rows = [x for x in rows if not math.isnan(x[2]) and abs(x[2]) >= DBL_MIN]
      if not rows:
        continue
      du = max(spread(x[0][3][0:4]) for x in rows)
      bu = max(spread(x[1][3][0:4]) for x in rows)
      dl = max(spread(x[0][3][4:8]) for x in rows)
      bl = max(spread(x[1][3][4:8]) for x in rows)
      print(f"{region:<24} {du:>10.2e} {bu:>10.2e} {dl:>10.2e} {bl:>10.2e}")
  return 0


if __name__ == "__main__":
  sys.exit(main())
