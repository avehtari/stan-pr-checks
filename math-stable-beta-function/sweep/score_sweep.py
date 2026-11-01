#!/usr/bin/env python3
"""Score beta(a, b) of two trees (develop, branch) on the sweep grid.

Reference: mpmath at 50 digits by two routes, mp.beta(a, b) and
exp(loggamma(a) + loggamma(b) - loggamma(a + b)); the script reports their
largest relative difference.  Regions: small (both arguments below 10),
mixed (one below 10), large (both at least 10).  Points whose reference
is below the smallest normal double (2.2e-308) form their own row
("underflow"); there the error is the absolute difference in units of
the smallest subnormal (4.9e-324).

For each region: the worst relative error of each side and its point; the
number of points where the branch error is above 4 times the develop
error and above 1e-15 ("worse"); the non-finite results of each side; and
the number of points where the two sides differ bitwise.

Usage: python3 score_sweep.py DEVELOP_OUT BRANCH_OUT [NPROC]
"""
import math
import sys
from multiprocessing import Pool

import mpmath as mp

DBL_MIN = 2.2250738585072014e-308
EPS = 2.220446049250313e-16
REGIONS = ["small", "mixed", "large", "underflow"]


def parse(path):
    grid, edges = [], []
    for line in open(path):
        tag, a, b, v = line.split()
        a, b = float.fromhex(a), float.fromhex(b)
        v = None if v == "THROW" else (float(v) if "n" in v.lower()
                                       else float.fromhex(v))
        (grid if tag == "G" else edges).append((a, b, v, line.split()[3]))
    return grid, edges


def ref(point):
    a, b = point
    with mp.workdps(50):
        ma, mb = mp.mpf(a), mp.mpf(b)
        r1 = mp.beta(ma, mb)
        r2 = mp.exp(mp.loggamma(ma) + mp.loggamma(mb) - mp.loggamma(ma + mb))
        agree = abs(r1 - r2) / r1
        return r1, float(agree)


def err(v, r):
    """Relative error, or for an underflow reference the absolute error in
    units of 4.9e-324; inf for a non-finite v."""
    if v is None or not math.isfinite(v):
        return math.inf
    with mp.workdps(50):
        if r < DBL_MIN:
            return float(abs(mp.mpf(v) - r) / mp.mpf(5e-324))
        return float(abs(mp.mpf(v) - r) / r)


def region(a, b, r):
    if r < DBL_MIN:
        return "underflow"
    return REGIONS[(a >= 10) + (b >= 10)]


def main():
    dev_grid, dev_edges = parse(sys.argv[1])
    br_grid, br_edges = parse(sys.argv[2])
    nproc = int(sys.argv[3]) if len(sys.argv) > 3 else 1
    assert len(dev_grid) == len(br_grid)
    points = [(a, b) for a, b, _, _ in dev_grid]
    for (a1, b1, _, _), (a2, b2, _, _) in zip(dev_grid, br_grid):
        assert a1 == a2 and b1 == b2
    with Pool(nproc) as pool:
        refs = pool.map(ref, points, chunksize=2000)

    stats = {
        r: dict(n=0, wd=0.0, pd=None, wb=0.0, pb=None, worse=0, nfd=0,
                nfb=0, diff=0, nd=0.0, nb=0.0)
        for r in REGIONS
    }
    worst_agree = 0.0
    worse_points = []
    for (a, b, vd, hd), (_, _, vb, hb), (r, agree) in zip(dev_grid, br_grid,
                                                          refs):
        worst_agree = max(worst_agree, agree)
        reg = region(a, b, r)
        s = stats[reg]
        s["n"] += 1
        ed, eb = err(vd, r), err(vb, r)
        if vd is None or not math.isfinite(vd):
            s["nfd"] += 1
        if vb is None or not math.isfinite(vb):
            s["nfb"] += 1
        if hd != hb:
            s["diff"] += 1
        if math.isfinite(ed) and ed > s["wd"]:
            s["wd"], s["pd"] = ed, (a, b)
        if math.isfinite(eb) and eb > s["wb"]:
            s["wb"], s["pb"] = eb, (a, b)
        scale = 1.0
        if reg != "underflow":
            # exp() multiplies an absolute error of log B by B; scale the
            # error by (1 + |log B|) eps, the rounding level of exp(log B)
            with mp.workdps(50):
                scale = float((1 + abs(mp.log(r))) * EPS)
            s["nd"] = max(s["nd"], ed / scale)
            s["nb"] = max(s["nb"], eb / scale)
        limit = 1e-15 if reg != "underflow" else 1.0
        if eb > 4 * ed and eb > limit:
            s["worse"] += 1
            worse_points.append((reg, a, b, ed, eb, eb / scale))

    print("# two mpmath routes: largest relative difference %.2e"
          % worst_agree)
    print("# region     points  develop_worst  at (a, b)               "
          "branch_worst  at (a, b)               worse  nonfinite_d/b  "
          "bitwise_diff  scaled_d  scaled_b")
    for reg in REGIONS:
        s = stats[reg]

        def at(p):
            return "(%.4g, %.4g)" % p if p else "-"

        print("%-10s %7d  %12.2e  %-22s  %12.2e  %-22s  %5d  %6d/%-6d  %6d"
              "  %8.3g  %8.3g"
              % (reg, s["n"], s["wd"], at(s["pd"]), s["wb"], at(s["pb"]),
                 s["worse"], s["nfd"], s["nfb"], s["diff"], s["nd"],
                 s["nb"]))
    print("# underflow errors are in units of 4.9e-324, not relative")
    print("# scaled: the worst relative error / ((1 + |log B|) eps)")
    for reg, a, b, ed, eb, sb in worse_points:
        print("WORSE %s a=%.17g b=%.17g develop=%.3e branch=%.3e "
              "branch_scaled=%.3g" % (reg, a, b, ed, eb, sb))

    print("# edge inputs: develop -> branch, where they differ")
    changed = 0
    for (a, b, vd, hd), (_, _, vb, hb) in zip(dev_edges, br_edges):
        if hd != hb:
            changed += 1
            print("EDGE a=%r b=%r develop=%s branch=%s" % (a, b, hd, hb))
    print("# edge inputs: %d of %d differ" % (changed, len(dev_edges)))


if __name__ == "__main__":
    main()
