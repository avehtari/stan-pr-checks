#!/usr/bin/env python3
"""Claim N4 (comment on 0d4ac36e34): score probe_cdfs outputs against mpmath.

Reference at the double inputs the probe used (y = 0.3 + 1.7 z rounded, and
exp of that for lognormal), at 100 digits.  Partials by mpmath.diff (route A)
and, for d/dy, by the analytic density (route B); the script prints the
largest relative difference of the two routes.  A component whose reference
is below DBL_MIN counts as correct when the result is too.  Prints the worst
ulp (relative error / 2^-53) per z bin and component for every file, and the
count of non-finite results.

Usage: score_cdfs.py out_A.txt [out_B.txt ...]
"""
import math
import sys

import mpmath

mpmath.mp.dps = 100
EPS = mpmath.mpf(2) ** -53
TINY = mpmath.mpf(2) ** -1022
MU, SIGMA, LAMBDA = 0.3, 1.7, 1.0
BINS = [(-38, -30), (-30, -10), (-10, -3), (-3, 3), (3, 8.01)]
COMP = ["value", "dy", "dmu", "dsigma", "dlambda"]


def cdf(name, y, mu, s, lam):
    if name == "std_normal":
        return mpmath.ncdf(y)
    if name == "normal":
        return mpmath.ncdf((y - mu) / s)
    if name == "lognormal":
        return mpmath.ncdf((mpmath.log(y) - mu) / s)
    z = (y - mu) / s
    v = lam * s
    return mpmath.ncdf(z) - mpmath.exp(v * v / 2 - lam * (y - mu)) * mpmath.ncdf(z - v)


def dy_density(name, y, mu, s):
    """Route B for d/dy: the density, no numerical differentiation."""
    if name == "std_normal":
        return mpmath.npdf(y)
    if name == "normal":
        return mpmath.npdf((y - mu) / s) / s
    if name == "lognormal":
        return mpmath.npdf((mpmath.log(y) - mu) / s) / (s * y)
    return None


def reference(name, z):
    y = MU + SIGMA * z
    if name == "std_normal":
        y = z
    elif name == "lognormal":
        y = math.exp(y)
    args = [mpmath.mpf(y), mpmath.mpf(MU), mpmath.mpf(SIGMA), mpmath.mpf(LAMBDA)]
    out = [cdf(name, *args)]
    for k in range(4):
        def f(t, k=k):
            a = list(args)
            a[k] = t
            return cdf(name, *a)
        out.append(mpmath.diff(f, args[k]))
    b = dy_density(name, *args[:3])
    diff = 0.0
    if b is not None and b != 0:
        diff = float(abs(out[1] / b - 1))
    return out, diff


def used(name):
    return {"std_normal": 2, "normal": 4, "lognormal": 4}.get(name, 5)


def ulp(got, exact):
    if not math.isfinite(got):
        return float("inf")
    if abs(exact) < TINY:
        return 0.0 if abs(got) < 2 ** -1022 else float("inf")
    return float(abs((mpmath.mpf(got) - exact) / exact) / EPS)


def load(path):
    rows = {}
    for line in open(path):
        p = line.split()
        vals = [float("nan") if "nan" in t else float.fromhex(t) for t in p[1:]]
        rows[(p[0], vals[0])] = vals[1:]
    return rows


def main(paths):
    trees = [load(p) for p in paths]
    tags = [p.rsplit("/", 1)[-1] for p in paths]
    print("files: " + " / ".join(tags))
    names = ["std_normal", "normal", "lognormal", "exp_mod_normal"]
    for name in names:
        nc = used(name)
        w = [{b: [0.0] * nc for b in BINS} for _ in trees]
        bad = [0] * len(trees)
        route = 0.0
        for (n, z) in sorted(k for k in trees[0] if k[0] == name):
            ref, diff = reference(name, z)
            route = max(route, diff)
            for t, tree in enumerate(trees):
                got = tree[(n, z)]
                bad[t] += sum(1 for g in got[:nc] if not math.isfinite(g))
                for b in BINS:
                    if b[0] <= z < b[1]:
                        for i in range(nc):
                            w[t][b][i] = max(w[t][b][i], ulp(got[i], ref[i]))
        print("\n%s: worst ulp per file; non-finite %s; d/dy routes differ by %.2g"
              % (name, " / ".join(str(c) for c in bad), route))
        print("%-12s" % "z" + "".join("%30s" % c for c in COMP[:nc]))
        for b in BINS:
            cells = [" /".join("%8.2g" % w[t][b][i] for t in range(len(trees)))
                     for i in range(nc)]
            print("[%4g,%4g)  " % (b[0], round(b[1])) + "".join("%30s" % c for c in cells))


if __name__ == "__main__":
    main(sys.argv[1:])
