#!/usr/bin/env python3
"""Points for, and scoring of, digamma_diff(x, d) = psi(x + d) - psi(x).

Usage:
  score_dd.py gen <points file>     write the points
  score_dd.py score <probe output>  score probe_dd output against mpmath

Regions (x log-uniform in each band; d per region):
  INT_<band>   d an integer from D_INT
  REAL_<band>  d log-uniform in [1e-6, 1e7]
  HALF_<band>  d = 0.5 (student_t), 1 and 2.5
  EDGE         x and d at the seam x = 10, x tiny, d tiny, d huge

Errors in ulp of the reference (|err| / (ulp of |ref|)). The derivative
d/dx is psi'(x + d) - psi'(x); its error is reported relative.
"""

import math
import os
import random
import sys

from mpmath import mp, mpf

mp.dps = 60
ROUTE_B = os.environ.get("DD_ROUTE_B", "1") != "0"

BANDS = {"a": (1e-3, 10.0), "b": (10.0, 1e4), "c": (1e4, 1e10),
         "d": (1e10, 1e18)}
D_INT = [1, 2, 3, 5, 10, 57, 117, 1000, 10 ** 6]
PER_REGION = 600


def log_uniform(rng, lo, hi):
  return math.exp(rng.uniform(math.log(lo), math.log(hi)))


def gen(path):
  rng = random.Random(20261007)
  out = []
  for band, (lo, hi) in BANDS.items():
    rows = [(log_uniform(rng, lo, hi), float(rng.choice(D_INT)))
            for _ in range(PER_REGION)]
    out.append(f"# region INT_{band} {len(rows)}")
    out += [f"{x.hex()} {d.hex()}" for x, d in rows]
    rows = [(log_uniform(rng, lo, hi), log_uniform(rng, 1e-6, 1e7))
            for _ in range(PER_REGION)]
    out.append(f"# region REAL_{band} {len(rows)}")
    out += [f"{x.hex()} {d.hex()}" for x, d in rows]
    rows = [(log_uniform(rng, lo, hi), rng.choice([0.5, 1.0, 2.5]))
            for _ in range(PER_REGION)]
    out.append(f"# region HALF_{band} {len(rows)}")
    out += [f"{x.hex()} {d.hex()}" for x, d in rows]
  edge = []
  for x in [9.999999999999998, 10.0, 10.000000000000002, 1e-300, 5e-324,
            1.0, 0.5, 1e18, 1e300]:
    for d in [5e-324, 1e-300, 1e-10, 0.5, 1.0, 1e6, 1e18, 1e300, 1.7e308]:
      edge.append((x, d))
  out.append(f"# region EDGE {len(edge)}")
  out += [f"{float(x).hex()} {float(d).hex()}" for x, d in edge]
  with open(path, "w") as fh:
    fh.write("\n".join(out) + "\n")
  print(f"wrote {path}")


def ulp(x):
  x = abs(float(x))
  if x == 0 or math.isinf(x):
    return 5e-324
  return math.ulp(x)


def score(path):
  region = None
  stats = {}
  order = []
  with open(path) as fh:
    for line in fh:
      if line.startswith("# region"):
        region = line.split()[2]
        order.append(region)
        stats[region] = {"n": 0, "dd": 0.0, "plain": 0.0, "dx": 0.0,
                         "thrown": 0, "nonfinite": 0, "worst": None,
                         "routes": 0.0, "n_b": 0}
        continue
      if line.startswith("#") or not line.strip():
        continue
      left, right = line.split("=>")
      xs, ds = left.split()
      vals = right.split()
      s = stats[region]
      s["n"] += 1
      if vals[-1] == "T":
        s["thrown"] += 1
        continue
      x, d = float.fromhex(xs), float.fromhex(ds)
      v, plain, dx = (float.fromhex(t) if "nan" not in t else math.nan
                      for t in vals[:3])
      mx, md = mpf(x), mpf(d)
      # psi(x + d) - psi(x) is about d / x: the cancellation costs
      # log10(x) digits for large x and log10(1 / d) digits for small d
      extra = (max(0, int(math.log10(max(x, 1.0))))
               + max(0, int(-math.log10(min(d, 1.0)))))
      with mp.workdps(80 + extra):
        ref = mp.digamma(mx + md) - mp.digamma(mx)
        dref = mp.polygamma(1, mx + md) - mp.polygamma(1, mx)
      # route B: mp.diff in x of loggamma(x + d) - loggamma(x), with 40
      # more digits; skipped where x + d overflows a double or x > 1e100
      if ROUTE_B and math.isfinite(x + d) and x < 1e100:
        # the step is relative to x (x - h > 0 for tiny x); loggamma(x + d)
        # costs log10(d) more digits
        with mp.workdps(120 + extra + max(0, int(math.log10(max(d, 1.0))))):
          ref_b = mp.diff(lambda t: mp.loggamma(t + md) - mp.loggamma(t), mx,
                          h=mx * mpf(10) ** -50)
        s["routes"] = max(s["routes"],
                          float(abs(ref_b - ref) / abs(ref)) if ref != 0
                          else float(abs(ref_b)))
        s["n_b"] += 1
      if not math.isfinite(v):
        s["nonfinite"] += 1
        continue
      e = float(abs(mpf(v) - ref)) / ulp(ref)
      ep = float(abs(mpf(plain) - ref)) / ulp(ref) if math.isfinite(plain) \
          else math.inf
      edx = float(abs(mpf(dx) - dref) / abs(dref)) if dref != 0 else 0.0
      if e > s["dd"]:
        s["dd"] = e
        s["worst"] = (x, d)
      s["plain"] = max(s["plain"], ep)
      s["dx"] = max(s["dx"], edx)
  print("routes: worst relative disagreement of route A (digamma difference)")
  print("and route B (mp.diff of the loggamma difference), over n_B points")
  print(f"{'region':8s} {'n':>5s} {'dd ulp':>9s} {'plain ulp':>10s} "
        f"{'d/dx rel':>9s} {'thrown':>6s} {'nonfin':>6s} {'routes':>8s} "
        f"{'n_B':>5s}  worst (x, d)")
  total = {"n": 0, "dd": 0.0, "plain": 0.0}
  for r in order:
    s = stats[r]
    w = s["worst"]
    ws = f"({w[0]:.4g}, {w[1]:.4g})" if w else ""
    print(f"{r:8s} {s['n']:5d} {s['dd']:9.2f} {s['plain']:10.3g} "
          f"{s['dx']:9.1e} {s['thrown']:6d} {s['nonfinite']:6d} "
          f"{s['routes']:8.1e} {s['n_b']:5d}  {ws}")
    total["n"] += s["n"]
    total["dd"] = max(total["dd"], s["dd"])
    total["plain"] = max(total["plain"], s["plain"])
  print(f"all regions: {total['n']} points, worst {total['dd']:.2f} ulp, "
        f"plain difference {total['plain']:.3g} ulp")


if __name__ == "__main__":
  if len(sys.argv) != 3 or sys.argv[1] not in ("gen", "score"):
    print(__doc__)
    sys.exit(2)
  if sys.argv[1] == "gen":
    gen(sys.argv[2])
  else:
    score(sys.argv[2])
