#!/usr/bin/env python3
"""Incomplete beta gradient roots: path A against path B, and both against
the mpmath reference, per region.

  path A  inc_beta_dda, inc_beta_ddb          (beta_cdf, neg_binomial_*_cdf)
  path B  grad_reg_inc_beta through grad_2F1  (beta_lcdf, beta_lccdf, ...)

Usage:
  python3 sweep_beta.py gen   [n_per_region] [points.txt]
  python3 sweep_beta.py eval  points.txt values.txt [nproc]
  python3 sweep_beta.py score points.txt values.txt [nproc] [name=extra.txt ...]

`gen` writes the point file (hex, exact in both directions). `eval` runs the
binary64 transcriptions of both paths at every point, writes the values file
and prints the A-against-B spread. `score` computes the mpmath reference at
every point (cached in refs_beta_<n>.txt) and prints the per-region error
tables for each path. `eval` runs anywhere; `score` runs on Triton.

Regions follow the branch structure of inc_beta_dda (the reflection seams)
and the consumers (neg_binomial with integer b), not a uniform grid.

Error columns:
  rel      |v - r| / |r|, the usual relative error, for |r| >= DBL_MIN
  err_log  |v - r| / I_z(a, b), the absolute error of d/dtheta log I, which
           is what beta_lcdf and beta_cdf (after the division by Pn) use
"""

import math
import multiprocessing as mputil
import os
import random
import sys

# ref_mpmath.py is in ../refs
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                "..", "refs"))

SEED = 20260920
DBL_MIN = 2.2250738585072014e-308


def log_uniform(rng, lo, hi):
  return math.exp(rng.uniform(math.log(lo), math.log(hi)))


# ── Regions ──────────────────────────────────────────────────────────────
def regions(n, rng):
  out = []

  def add(name, pts):
    out.append((name, pts))

  add("R1 a,b<1",
      [(log_uniform(rng, 0.01, 1.0), log_uniform(rng, 0.01, 1.0),
        rng.uniform(0.01, 0.99)) for _ in range(n)])
  add("R2 1<=a,b<=20 z<0.75",
      [(log_uniform(rng, 1.0, 20.0), log_uniform(rng, 1.0, 20.0),
        rng.uniform(0.01, 0.75)) for _ in range(n)])
  add("R3 seam z>0.75 a<500",
      [(log_uniform(rng, 1.0, 500.0), log_uniform(rng, 1.0, 500.0),
        rng.uniform(0.75, 0.999)) for _ in range(n)])
  pts = []
  for _ in range(n):
    b = log_uniform(rng, 500.0, 12500.0)
    a = log_uniform(rng, 0.5, b)
    pts.append((a, b, rng.uniform(0.1, 0.75)))
  add("R4 seam b>a b>500 z<=0.75", pts)
  add("R5 b>>a z small",
      [(log_uniform(rng, 0.1, 5.0), log_uniform(rng, 100.0, 12500.0),
        log_uniform(rng, 1e-4, 0.1)) for _ in range(n)])
  add("R6 large a,b z uniform",
      [(log_uniform(rng, 100.0, 12500.0), log_uniform(rng, 100.0, 12500.0),
        rng.uniform(0.2, 0.8)) for _ in range(n)])
  pts = []
  for _ in range(n):
    a = log_uniform(rng, 100.0, 12500.0)
    b = log_uniform(rng, 100.0, 12500.0)
    mode = a / (a + b)
    sd = math.sqrt(a * b / ((a + b) ** 2 * (a + b + 1)))
    z = min(0.999999, max(1e-6, mode + rng.uniform(-4.0, 4.0) * sd))
    pts.append((a, b, z))
  add("R7 large a,b z near mode", pts)
  add("R8 negbin n<=200",
      [(log_uniform(rng, 0.1, 1000.0), float(rng.randint(0, 200) + 1),
        rng.uniform(0.05, 0.95)) for _ in range(n)])
  add("R9 negbin 200<n<=5000",
      [(log_uniform(rng, 0.1, 1000.0), float(rng.randint(201, 5000) + 1),
        rng.uniform(0.05, 0.95)) for _ in range(n)])
  return out


def gen(n, path):
  rng = random.Random(SEED)
  total = 0
  with open(path, "w") as fh:
    for name, pts in regions(n, rng):
      fh.write(f"# {name}\n")
      for a, b, z in pts:
        fh.write(f"{a.hex()} {b.hex()} {z.hex()}\n")
        total += 1
  print(f"wrote {total} points to {path}")


def read_points(path):
  rows = []
  region = "?"
  for line in open(path):
    line = line.strip()
    if not line:
      continue
    if line.startswith("#"):
      region = line[2:]
      continue
    a, b, z = (float.fromhex(p) for p in line.split()[:3])
    rows.append((region, a, b, z))
  return rows


# ── eval: both transcriptions ────────────────────────────────────────────
def eval_one(abz):
  import transcript as T
  a, b, z = abz
  out = []
  for fn in (T.inc_beta_dda, T.inc_beta_ddb):
    try:
      v, k, br = fn(a, b, z)
    except Exception:
      v, k, br = math.nan, -1, "throw"
    out.append((v, k, br))
  try:
    g1, g2, k, br = T.grad_reg_inc_beta(a, b, z)
  except Exception:
    g1, g2, k, br = math.nan, math.nan, -1, "throw"
  out.append((g1, k, br))
  out.append((g2, k, br))
  return out


def spread(x, y):
  if not (math.isfinite(x) and math.isfinite(y)):
    return math.inf
  m = max(abs(x), abs(y))
  return 0.0 if m == 0.0 else abs(x - y) / m


def run_eval(points_path, values_path, nproc):
  pts = read_points(points_path)
  with mputil.Pool(nproc) as pool:
    vals = pool.map(eval_one, [(a, b, z) for _, a, b, z in pts], chunksize=8)
  with open(values_path, "w") as fh:
    region = None
    for (reg, a, b, z), v in zip(pts, vals):
      if reg != region:
        fh.write(f"# {reg}\n")
        region = reg
      fh.write(f"{a.hex()} {b.hex()} {z.hex()} "
               f"{v[0][0].hex()} {v[1][0].hex()} {v[2][0].hex()} {v[3][0].hex()} "
               f"{v[0][1]} {v[1][1]} {v[2][1]} {v[0][2]} {v[1][2]} {v[2][2]}\n")
  print(f"wrote {len(pts)} rows to {values_path}\n")

  print("# path A (inc_beta_dda/ddb) against path B (grad_reg_inc_beta)")
  print("# spread = |A - B| / max(|A|, |B|)")
  print("# nfA, nfB = non-finite results per path; C0 = path B skipped "
        "grad_2F1 because z^a (1-z)^b underflowed;")
  print("# bAB0 = beta(a, b) underflowed, path B divides by 0\n")
  print(f"{'region':<28} {'n':>5} {'da med':>9} {'da worst':>9} "
        f"{'db med':>9} {'db worst':>9} {'>1e-6':>6} {'nfA':>4} {'nfB':>4} "
        f"{'C0':>4} {'bAB0':>4} "
        f"{'kA med':>7} {'kA max':>7} {'kB med':>7} {'kB max':>7}")
  by_region = {}
  for (reg, a, b, z), v in zip(pts, vals):
    by_region.setdefault(reg, []).append(((a, b, z), v))
  for reg, rows in by_region.items():
    sa = sorted(spread(v[0][0], v[2][0]) for _, v in rows)
    sb = sorted(spread(v[1][0], v[3][0]) for _, v in rows)
    big = sum(1 for x in sa + sb if x > 1e-6)
    nfa = sum(1 for _, v in rows
              if not all(math.isfinite(v[i][0]) for i in range(2)))
    nfb = sum(1 for _, v in rows
              if not all(math.isfinite(v[i][0]) for i in range(2, 4)))
    cu = sum(1 for _, v in rows if v[2][2] == "C_underflow")
    bu = sum(1 for _, v in rows if v[2][2] == "betaAB_underflow")
    ka = sorted(max(v[0][1], v[1][1]) for _, v in rows)
    kb = sorted(v[2][1] for _, v in rows)
    n = len(rows)
    print(f"{reg:<28} {n:>5d} {sa[n // 2]:>9.2e} {sa[-1]:>9.2e} "
          f"{sb[n // 2]:>9.2e} {sb[-1]:>9.2e} {big:>6d} {nfa:>4d} {nfb:>4d} "
          f"{cu:>4d} {bu:>4d} "
          f"{ka[n // 2]:>7d} {ka[-1]:>7d} {kb[n // 2]:>7d} {kb[-1]:>7d}")
    worst = max(rows, key=lambda r: spread(r[1][0][0], r[1][2][0]))
    (a, b, z), v = worst
    print(f"    worst da: a={a:.6g} b={b:.6g} z={z:.6g}  A={v[0][0]:.9g} "
          f"({v[0][2]}, k={v[0][1]})  B={v[2][0]:.9g} ({v[2][2]}, k={v[2][1]})")


# ── score: mpmath reference ──────────────────────────────────────────────
def ref_one(abz):
  import ref_mpmath as R
  from mpmath import mp, mpf
  a, b, z = abz
  try:
    # betainc_reg, not mp.betainc: the latter raises for a + b above a few
    # thousand near the mode, which left 576 reference NaNs in job 20359555
    with mp.workdps(40):
      i = float(R.betainc_reg(a, b, z))
    da = float(R.dIda(a, b, z))
    db = float(R.dIdb(a, b, z))
  except Exception:
    return (math.nan, math.nan, math.nan)
  return (i, da, db)


def load_or_compute_refs(pts, nproc):
  cache = f"refs_beta_{len(pts)}.txt"
  keys = [(a, b, z) for _, a, b, z in pts]
  refs = {}
  if os.path.exists(cache):
    for line in open(cache):
      f = line.split()
      refs[tuple(float.fromhex(x) for x in f[:3])] = tuple(
          float.fromhex(x) for x in f[3:6])
  # a NaN entry is a failed reference; recompute it, keep the rest
  todo = [k for k in keys if k not in refs or math.isnan(refs[k][1])]
  if not todo:
    print(f"# references loaded from {cache}")
    return [refs[k] for k in keys]
  print(f"# computing {len(todo)} references ({len(keys) - len(todo)} "
        f"loaded from {cache})")
  # mpmath's internal caches grow with the number of distinct high-precision
  # gamma and digamma arguments; job 20358916 reached 1.3 GB per worker
  # and was throttled, then OOM-killed, at 20 workers. Recycle each worker
  # after a bounded number of points.
  with mputil.Pool(nproc, maxtasksperchild=32) as pool:
    vals = pool.map(ref_one, todo, chunksize=8)
  for k, v in zip(todo, vals):
    refs[k] = v
  with open(cache, "w") as fh:
    for (a, b, z) in keys:
      i, da, db = refs[(a, b, z)]
      fh.write(f"{a.hex()} {b.hex()} {z.hex()} {i.hex()} {da.hex()} "
               f"{db.hex()}\n")
  print(f"# references cached in {cache}")
  return [refs[k] for k in keys]


def read_values(path):
  rows = []
  for line in open(path):
    line = line.strip()
    if not line or line.startswith("#"):
      continue
    f = line.split()
    rows.append((tuple(float.fromhex(x) for x in f[:3]),
                 tuple(float.fromhex(x) for x in f[3:7]),
                 tuple(int(x) for x in f[7:10]), tuple(f[10:13])))
  return rows


def rel(v, r):
  if not math.isfinite(v):
    return math.inf
  if r == 0.0:
    return 0.0 if v == 0.0 else math.inf
  return abs(v - r) / abs(r)


def read_extra_values(path):
  """A third-party values file: a b z v_da v_db per line, hex, 'nan' for a
  non-finite result. Returns rows in the read_values() layout."""
  rows = []
  for line in open(path):
    line = line.strip()
    if not line or line.startswith("#"):
      continue
    f = line.split()
    v = tuple(math.nan if x == "nan" else float.fromhex(x) for x in f[3:5])
    rows.append((tuple(float.fromhex(x) for x in f[:3]), v, (0, 0, 0),
                 ("-", "-", "-")))
  return rows


def run_score(points_path, values_path, nproc, extra=None):
  pts = read_points(points_path)
  vals = read_values(values_path)
  assert len(pts) == len(vals), "row counts differ"
  for p, v in zip(pts, vals):
    assert (p[1], p[2], p[3]) == v[0], "point files are not aligned"
  refs = load_or_compute_refs(pts, nproc)

  # series: (name, rows, column offset of d/da in rows[i][1])
  series = [("A inc_beta_dda/ddb", vals, 0), ("B grad_reg_inc_beta", vals, 2)]
  for name, path in (extra or []):
    ev = read_extra_values(path)
    assert len(ev) == len(pts), f"row count differs in {path}"
    for p, v in zip(pts, ev):
      assert (p[1], p[2], p[3]) == v[0], f"{path} is not aligned"
    series.append((name, ev, 0))

  print("# reference: ref_mpmath.dIda / dIdb, density-form quadrature at "
        "80 digits")
  cols = (("d/da", 0, 1), ("d/db", 1, 2))
  for cname, ci, ri in cols:
    for pname, vals, off in series:
      col = off + ci
      print(f"\n## {cname}, path {pname}\n")
      print(f"{'region':<28} {'n':>5} {'rel med':>9} {'rel p99':>9} "
            f"{'rel worst':>9} {'err_log':>9} {'>1e-6':>6} {'nf':>4} "
            f"{'sub':>4} {'ref nan':>7}")
      by_region = {}
      for p, v, r in zip(pts, vals, refs):
        by_region.setdefault(p[0], []).append((p, v, r))
      for reg, rows in by_region.items():
        rels, elog = [], []
        nf = sub = rnan = 0
        worst = (0.0, None)
        worst_log = (0.0, None)
        for p, v, r in rows:
          i, rr = r[0], r[ri]
          val = v[1][col]
          if math.isnan(rr):
            rnan += 1
            continue
          if abs(rr) < DBL_MIN:
            sub += 1
            continue
          if not math.isfinite(val):
            nf += 1
          e = rel(val, rr)
          rels.append(e)
          if math.isfinite(e) and e > worst[0]:
            worst = (e, (p, val, rr, v[2], v[3]))
          if i > 0 and math.isfinite(val):
            el = abs(val - rr) / i
            elog.append(el)
            if el > worst_log[0]:
              worst_log = (el, (p, val, rr, v[2], v[3]))
        fin = sorted(x for x in rels if math.isfinite(x))
        n = len(rows)
        med = fin[len(fin) // 2] if fin else math.nan
        p99 = fin[int(0.99 * (len(fin) - 1))] if fin else math.nan
        wl = max(elog) if elog else math.nan
        big = sum(1 for x in rels if x > 1e-6)
        print(f"{reg:<28} {n:>5d} {med:>9.2e} {p99:>9.2e} "
              f"{worst[0]:>9.2e} {wl:>9.2e} {big:>6d} {nf:>4d} {sub:>4d} "
              f"{rnan:>7d}")
        if worst[1] is not None:
          p, val, rr, ks, brs = worst[1]
          print(f"    worst rel at a={p[1]:.8g} b={p[2]:.8g} z={p[3]:.8g} "
                f"got={val:.6g} ref={rr:.6g} k={ks} br={brs}")
        if worst_log[1] is not None and worst_log[1] is not worst[1]:
          p, val, rr, ks, brs = worst_log[1]
          print(f"    worst err_log at a={p[1]:.8g} b={p[2]:.8g} "
                f"z={p[3]:.8g} got={val:.6g} ref={rr:.6g} k={ks} br={brs}")


def main():
  mode = sys.argv[1] if len(sys.argv) > 1 else ""
  if mode == "gen":
    n = int(sys.argv[2]) if len(sys.argv) > 2 else 2000
    path = sys.argv[3] if len(sys.argv) > 3 else "points_beta.txt"
    gen(n, path)
  elif mode == "eval":
    nproc = int(sys.argv[4]) if len(sys.argv) > 4 else 6
    run_eval(sys.argv[2], sys.argv[3], nproc)
  elif mode == "score":
    nproc = int(sys.argv[4]) if len(sys.argv) > 4 else 6
    extra = [tuple(x.split("=", 1)) for x in sys.argv[5:]]
    run_score(sys.argv[2], sys.argv[3], nproc, extra)
  else:
    print(__doc__)
    return 2
  return 0


if __name__ == "__main__":
  sys.exit(main())
