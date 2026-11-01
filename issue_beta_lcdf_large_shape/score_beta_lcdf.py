#!/usr/bin/env python3
# Score the output of probe_beta_lcdf against mpmath.
# Route A: log of the lower regularized incomplete beta I_x(a, b).
# Route B: log1p(-I_{1-x}(b, a)), the complement of the upper tail.
# Columns: relative error of the Stan value, d/dx, d/da, d/db, of the Boost
# ibeta value, and the agreement of the routes (largest over the four).
# Usage: python3 score_beta_lcdf.py < probe output
import sys
import mpmath as mp

mp.mp.dps = 40


def route_a(x, a, b):
    return mp.log(mp.betainc(a, b, 0, x, regularized=True))


def route_b(x, a, b):
    return mp.log1p(-mp.betainc(b, a, 0, 1 - x, regularized=True))


def rel(s, r):
    return abs(mp.mpf(s) - r) / abs(r) if r != 0 else abs(mp.mpf(s))


def refs(route, x, a, b):
    return [route(x, a, b),
            mp.diff(lambda v: route(v, a, b), x),
            mp.diff(lambda v: route(x, v, b), a),
            mp.diff(lambda v: route(x, a, v), b)]


print(f"{'a':>4s} {'b':>11s} {'b*x':>5s} {'value':>9s} {'d/dx':>9s} {'d/da':>9s}"
      f" {'d/db':>9s} {'ibeta':>9s} {'A vs B':>8s}")
for line in sys.stdin:
    f = line.split()
    if len(f) != 8:
        continue
    a, b, x = (mp.mpf(v) for v in f[:3])
    stan = [float(v) for v in f[3:7]]
    ib = float(f[7])
    ra = refs(route_a, x, a, b)
    rb = refs(route_b, x, a, b)
    agree = max(rel(p, q) for p, q in zip(ra, rb))
    errs = [rel(s, r) for s, r in zip(stan, ra)]
    e_ib = rel(ib, mp.exp(ra[0]))
    print(f"{mp.nstr(a, 3):>4s} {mp.nstr(b, 4):>11s} {mp.nstr(b * x, 2):>5s} "
          + " ".join(f"{mp.nstr(e, 2):>9s}" for e in errs + [e_ib])
          + f" {mp.nstr(agree, 2):>8s}")
    sys.stdout.flush()
