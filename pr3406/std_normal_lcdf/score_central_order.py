#!/usr/bin/env python3
"""Claim N2 (comment on 5f33a10d31): 9.6 ulp on a denser grid of 20001 points.
Scores probe_central_order output (three summation orders and the slope) or
probe_dedup output (columns value, slope) against mpmath at 50 digits.

Two routes: A erfc / log1p(-erfc) as in score_dedup.py, B log(ncdf(z)) at
70 digits.  Prints the worst ulp (relative error / 2^-53) for z < 0 and
z >= 0 per column, and the largest relative difference between the routes.
Usage: score_central_order.py <file> [--lo LO --hi HI]
"""
import argparse

import mpmath as mp

mp.mp.dps = 50
EPS = mp.mpf(2) ** -53


def ref_a(z):
    if z < 0:
        q = mp.erfc(-z / mp.sqrt(2)) / 2
        return mp.log(q), mp.npdf(z) / q
    q = mp.erfc(z / mp.sqrt(2)) / 2
    return mp.log1p(-q), mp.npdf(z) / (1 - q)


def ref_b(z):
    with mp.workdps(70):
        p = mp.ncdf(z)
        return +mp.log(p), +(mp.npdf(z) / p)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("file")
    ap.add_argument("--lo", type=float, default=-0.66291)
    ap.add_argument("--hi", type=float, default=0.66291)
    a = ap.parse_args()
    rows = []
    for line in open(a.file):
        if line.startswith("#"):
            continue
        rows.append([float.fromhex(t) for t in line.split()])
    rows = [r for r in rows if a.lo <= r[0] <= a.hi]
    ncol = len(rows[0]) - 1
    # probe_central_order: a b c slope; probe_dedup: value slope
    kinds = ["value"] * (ncol - 1) + ["slope"]
    # a = LOG_HALF, h = z^2 / 2, L = log(erfcx(-z / sqrt 2))
    names = (["(a-h)+L", "(a+L)-h", "a+(L-h)", "slope"] if ncol == 4
             else ["value", "slope"])
    worst = {(s, j): (0.0, None) for s in ("z<0", "z>=0") for j in range(ncol)}
    route_diff = 0.0
    seen = set()
    for r in rows:
        z = mp.mpf(r[0])
        if r[0] in seen:
            continue
        seen.add(r[0])
        va, sa = ref_a(z)
        vb, sb = ref_b(z)
        route_diff = max(route_diff, float(abs(vb / va - 1)), float(abs(sb / sa - 1)))
        side = "z<0" if r[0] < 0 else "z>=0"
        for j in range(ncol):
            ref = va if kinds[j] == "value" else sa
            u = float(abs((mp.mpf(r[j + 1]) - ref) / ref) / EPS)
            if u > worst[(side, j)][0]:
                worst[(side, j)] = (u, r[0])
    print("file %s, %d points in [%g, %g]" % (a.file, len(seen), a.lo, a.hi))
    print("%-6s" % "side" + "".join("%26s" % n for n in names))
    for side in ("z<0", "z>=0"):
        cells = ["%8.3g @ %-14.9g" % worst[(side, j)] if worst[(side, j)][1] is not None
                 else "-" for j in range(ncol)]
        print("%-6s" % side + "".join("%26s" % c for c in cells))
    print("largest relative difference of routes A and B: %.2g" % route_diff)


if __name__ == "__main__":
    main()
