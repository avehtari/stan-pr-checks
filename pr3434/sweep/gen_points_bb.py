#!/usr/bin/env python3
"""Sweep points for beta_binomial_lpmf (stan-dev/math#3153).

Plan: ../2026-10-06_PLAN_fix_issue3153_beta_binomial_lpmf.md, step 5.2.

Output line:  <REGION> <n> <N> <alpha hex> <beta hex>
Block header: # region <REGION> <description> <count>

Hex (float.hex) is exact in both directions, so the probe and the scorer
see the same doubles.

Regions, POINTS_PER_REGION points each (REP has its own grid):
  EQ_a..EQ_d  alpha = beta = s, s log-uniform in a band
  UN_a..UN_d  alpha, beta independent, log-uniform in the band
  MC_a..MC_d  mean / concentration: mu ~ U(0.01, 0.99), kappa log-uniform
              in twice the band, alpha = mu * kappa, beta = (1 - mu) * kappa
  SEAM        one shape in [1, 100], the other in [1e6, 1e18]
  REP         the reporter: n = 57, N = 117, alpha = beta = exp(lc) / 2,
              lc = 0, 0.25, ..., 44

Bands: a [1e-2, 10], b [10, 1e6], c [1e6, 1e12], d [1e12, 1e18].
N is drawn from N_SET; n is 0, N, round(N * alpha / (alpha + beta)) or
uniform on 0..N, with equal probability.
"""

import math
import random
import sys

SEED = 20261007
POINTS_PER_REGION = 1000
N_SET = [1, 10, 117, 1000, 1000000]
BANDS = {
    "a": (1e-2, 10.0),
    "b": (10.0, 1e6),
    "c": (1e6, 1e12),
    "d": (1e12, 1e18),
}


def log_uniform(rng, lo, hi):
  return math.exp(rng.uniform(math.log(lo), math.log(hi)))


def draw_counts(rng, alpha, beta):
  N = rng.choice(N_SET)
  kind = rng.randrange(4)
  if kind == 0:
    n = 0
  elif kind == 1:
    n = N
  elif kind == 2:
    n = int(round(N * alpha / (alpha + beta)))
  else:
    n = rng.randint(0, N)
  return n, N


def line(tag, n, N, alpha, beta):
  return f"{tag} {n} {N} {float(alpha).hex()} {float(beta).hex()}"


def region_eq(rng, band):
  lo, hi = BANDS[band]
  s = log_uniform(rng, lo, hi)
  return s, s


def region_un(rng, band):
  lo, hi = BANDS[band]
  return log_uniform(rng, lo, hi), log_uniform(rng, lo, hi)


def region_mc(rng, band):
  lo, hi = BANDS[band]
  mu = rng.uniform(0.01, 0.99)
  kappa = log_uniform(rng, 2 * lo, 2 * hi)
  return mu * kappa, (1 - mu) * kappa


def region_seam(rng):
  small = log_uniform(rng, 1.0, 100.0)
  large = log_uniform(rng, 1e6, 1e18)
  if rng.random() < 0.5:
    return small, large
  return large, small


def main():
  rng = random.Random(SEED)
  out = []
  counts = {}

  def emit_region(tag, desc, draw):
    rows = []
    for _ in range(POINTS_PER_REGION):
      alpha, beta = draw()
      n, N = draw_counts(rng, alpha, beta)
      rows.append(line(tag, n, N, alpha, beta))
    out.append(f"# region {tag} {desc} {len(rows)}")
    out.extend(rows)
    counts[tag] = len(rows)

  for band, (lo, hi) in BANDS.items():
    emit_region(f"EQ_{band}", f"alpha=beta_in_[{lo:g},{hi:g}]",
                lambda band=band: region_eq(rng, band))
  for band, (lo, hi) in BANDS.items():
    emit_region(f"UN_{band}", f"alpha,beta_independent_in_[{lo:g},{hi:g}]",
                lambda band=band: region_un(rng, band))
  for band, (lo, hi) in BANDS.items():
    emit_region(f"MC_{band}", f"mu*kappa,kappa_in_[{2 * lo:g},{2 * hi:g}]",
                lambda band=band: region_mc(rng, band))
  emit_region("SEAM", "one_in_[1,100]_other_in_[1e6,1e18]",
              lambda: region_seam(rng))

  rows = []
  for i in range(0, 177):
    lc = 0.25 * i
    s = math.exp(lc) / 2
    rows.append(line("REP", 57, 117, s, s))
  out.append(f"# region REP n=57,N=117,alpha=beta=exp(lc)/2,lc=0..44 "
             f"{len(rows)}")
  out.extend(rows)
  counts["REP"] = len(rows)

  path = sys.argv[1] if len(sys.argv) > 1 else "points_bb.txt"
  with open(path, "w") as fh:
    fh.write("\n".join(out) + "\n")
  total = sum(counts.values())
  for tag, c in counts.items():
    print(f"{tag:6s} {c:6d}")
  print(f"total  {total:6d}  -> {path}")


if __name__ == "__main__":
  main()
