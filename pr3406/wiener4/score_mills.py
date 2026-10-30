#!/usr/bin/env python3
"""Score the log Mills ratio forms of probe_mills.cpp against mpmath.

Claim (PR comment): std_normal_lcdf(-x) - std_normal_lpdf(x) loses accuracy
without bound, also with correctly rounded inputs; the erfcx form does not.
Usage: score_mills.py <probe output> [<probe output> ...]
ulp = relative error / 2^-53.  Route A: log(ncdf(-x)) - log(npdf(x)) at
60 digits.  Route B: log(sqrt(pi/2) * Gamma(1/2, x^2/2) exp(x^2/2) / sqrt(pi)).
"""
import sys

import mpmath as mp

mp.mp.dps = 60
EPS = mp.mpf(2) ** -53


def route_a(x):
    return mp.log(mp.ncdf(-x)) - mp.log(mp.npdf(x))


def route_b(x):
    u = x / mp.sqrt(2)
    erfcx = mp.gammainc(mp.mpf(1) / 2, u * u, mp.inf) * mp.exp(u * u)
    return mp.log(mp.sqrt(mp.pi / 2) * erfcx / mp.sqrt(mp.pi))


def ulp(v, ref):
    return float(abs(mp.mpf(v) - ref) / abs(ref) / EPS)


def score(path):
    rows = []
    with open(path) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            rows.append([float.fromhex(s) for s in line.split()])
    print(f"# {path}")
    print(f"{'x':>10} {'difference':>11} {'CR inputs':>10} {'erfcx form':>10}"
          f" {'routes agree':>12}")
    worst = [0.0, 0.0, 0.0]
    for k, (x, lc, lp, diff, ex) in enumerate(rows):
        xm = mp.mpf(x)
        ra, rb = route_a(xm), route_b(xm)
        agree = abs(ra - rb) / abs(ra)
        # correctly rounded lcdf(-x) and lpdf(x), subtracted in double
        cr = float(mp.log(mp.ncdf(-xm))) - float(mp.log(mp.npdf(xm)))
        u = [ulp(diff, ra), ulp(cr, ra), ulp(ex, ra)]
        worst = [max(w, v) for w, v in zip(worst, u)]
        if k < 4:
            print(f"{x:10.4g} {u[0]:11.1f} {u[1]:10.1f} {u[2]:10.2f} "
                  f"{mp.nstr(agree, 2):>12}")
        elif k == 4:
            print("# log grid 1 to 1e4, 20 points per decade")
        if k >= 4 and (k - 4) % 10 == 0:
            print(f"{x:10.4g} {u[0]:11.1f} {u[1]:10.1f} {u[2]:10.2f} "
                  f"{mp.nstr(agree, 2):>12}")
    print(f"{'worst':>10} {worst[0]:11.1f} {worst[1]:10.1f} {worst[2]:10.2f}")


if __name__ == "__main__":
    for p in sys.argv[1:]:
        score(p)
