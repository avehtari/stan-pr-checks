#!/usr/bin/env python3
"""Prototype (double precision) of the exp_mod_normal lcdf value for z < 0
through the erfcx ratio, against the mpmath reference of score_emg_phi.py.

  F = Phi(z) (1 - R),  R = erfcx(u + d) / erfcx(u),  u = -z/sqrt2, d = v/sqrt2
  (the z^2/2 terms of b - a cancel exactly).
  A: log R = log erfcx(u + d) - log erfcx(u)
  B: log R = -log1p(d / u) + g(u) - g(u + d),  g(t) = log(sqrt(pi) t erfcx(t))
  log F = log Phi(z) + log(-expm1(log R)); log Phi(z) from mpmath rounded to
  double (the PR's kernel is accurate there).

Prints worst ulp of the value per z, for A, B and the PR (out_pr.txt).
Usage: proto_emg_left_tail.py out_pr.txt
"""
import math
import sys

import mpmath
from scipy.special import erfcx

sys.path.insert(0, ".")
from score_emg_phi import EPS, emg_log, load  # noqa: E402

SQRT_PI = math.sqrt(math.pi)


def g(t):
    return math.log(SQRT_PI * t * erfcx(t))


def proto(z, v, form):
    u = -z / math.sqrt(2)
    d = v / math.sqrt(2)
    if form == "A":
        log_r = math.log(erfcx(u + d)) - math.log(erfcx(u))
    else:
        log_r = -math.log1p(d / u) + (g(u) - g(u + d))
    log_phi = float(mpmath.log(mpmath.ncdf(z)))
    return log_phi + math.log(-math.expm1(log_r))


def main(path):
    tree, _ = load(path)
    zs = sorted({k[1] for k in tree if k[1] < 0})
    vs = sorted({k[2] for k in tree})
    print("lcdf value, worst ulp over v in %s" % vs)
    print("%8s %10s %10s %10s" % ("z", "PR", "A", "B"))
    for z in zs:
        w = {"PR": 0.0, "A": 0.0, "B": 0.0}
        for v in vs:
            ref = emg_log(False, mpmath.mpf(z), 0, 1, mpmath.mpf(v))
            got = {"PR": tree[("lcdf", z, v)][0], "A": proto(z, v, "A"),
                   "B": proto(z, v, "B")}
            for k, x in got.items():
                if not math.isfinite(x):
                    w[k] = float("inf")
                else:
                    w[k] = max(w[k], float(abs((x - ref) / ref) / EPS))
        print("%8g %10.1e %10.1e %10.1e" % (z, w["PR"], w["A"], w["B"]))


if __name__ == "__main__":
    main(sys.argv[1])
