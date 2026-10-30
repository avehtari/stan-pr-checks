#!/usr/bin/env python3
"""Claims M1a and M2 of README_emg.md (reference() used by score_emg_grid.py):
score probe_emg_phi output (develop and the PR) against mpmath.

exp_mod_normal with y = z, mu = 0, sigma = 1, lambda = v:
  F = Phi(z) - exp(v^2/2 - v z) Phi(z - v),  lcdf = log F, lccdf = log(1 - F).
Gradients by mpmath.diff on the 4-parameter form (100 digits); dy is
cross-checked against the analytic derivative.

Prints, per kind (lcdf, lccdf) and tree, a z x v table of the worst ulp over
value and the 4 gradients ("x" = non-finite or throw where the reference is
finite), and the worst-ulp field.  Then Phi relative error per x bin.

Usage: score_emg_phi.py out_dev.txt out_pr.txt
"""
import math
import sys

import mpmath

mpmath.mp.dps = 100
EPS = mpmath.mpf(2) ** -53


def emg_log(upper, y, mu, s, lam):
    # F and 1 - F are both formed directly; the log of the larger one is
    # taken as log1p of minus the smaller one, so neither side rounds to 1.
    z = (y - mu) / s
    v = lam * s
    b = mpmath.exp(v * v / 2 - lam * (y - mu)) * mpmath.ncdf(z - v)
    f_low = mpmath.ncdf(z) - b
    f_up = mpmath.ncdf(-z) + b
    if upper:
        return mpmath.log(f_up) if f_up < 0.5 else mpmath.log1p(-f_low)
    return mpmath.log(f_low) if f_low < 0.5 else mpmath.log1p(-f_up)


def reference(upper, z, v):
    args = [mpmath.mpf(z), mpmath.mpf(0), mpmath.mpf(1), mpmath.mpf(v)]
    val = emg_log(upper, *args)
    grads = []
    for k in range(4):
        def f(t, k=k):
            a = list(args)
            a[k] = t
            return emg_log(upper, *a)
        grads.append(mpmath.diff(f, args[k]))
    # second route for d/dy: (f(z) - (-v) b ... ) analytic
    zz, vv = args[0], args[3]
    b = mpmath.exp(vv * vv / 2 - vv * zz) * mpmath.ncdf(zz - vv)
    if upper:
        dy = (-mpmath.npdf(zz) - vv * b + mpmath.exp(vv * vv / 2 - vv * zz)
              * mpmath.npdf(zz - vv)) / mpmath.exp(val)
    else:
        dy = (mpmath.npdf(zz) + vv * b - mpmath.exp(vv * vv / 2 - vv * zz)
              * mpmath.npdf(zz - vv)) / mpmath.exp(val)
    check = abs(dy - grads[0]) / abs(dy) if dy != 0 else 0
    return [val] + grads, check


def ulp(got, exact):
    if math.isnan(got) or math.isinf(got):
        return float("inf")
    if exact == 0:
        return 0.0 if got == 0 else float("inf")
    return float(abs((mpmath.mpf(got) - exact) / exact) / EPS)


def fmt(u):
    if u == float("inf"):
        return "x"
    if u < 10:
        return "%.0f" % u
    return "%.0e" % u


def load(path):
    emg, phi = {}, []
    for line in open(path):
        p = line.split()
        if p[0] == "E":
            vals = [float.fromhex(t) if t != "nan" and "nan" not in t
                    else float("nan") for t in p[2:]]
            emg[(p[1], vals[0], vals[1])] = vals[2:]
        elif p[0] == "P":
            phi.append((float.fromhex(p[1]), float.fromhex(p[2])))
    return emg, phi


def main(p_dev, p_pr):
    dev, phi_dev = load(p_dev)
    pr, phi_pr = load(p_pr)
    zs = sorted({k[1] for k in dev})
    vs = sorted({k[2] for k in dev})
    worst_check = (0.0, None)
    refs = {}
    for key in dev:
        refs[key], c = reference(key[0] == "lccdf", key[1], key[2])
        if float(c) > worst_check[0]:
            worst_check = (float(c), key)
    print("max rel. difference analytic vs mpmath.diff for d/dy: %.1e at %s"
          % worst_check)
    names = ["value", "dy", "dmu", "dsigma", "dlambda"]
    for kind in ("lcdf", "lccdf"):
        for label, tree in (("develop", dev), ("PR", pr)):
            print("\n%s %s: worst ulp over value and gradients" % (kind, label))
            print("%8s " % "z \\ v" + " ".join("%6g" % v for v in vs))
            for z in zs:
                cells = []
                for v in vs:
                    got = tree[(kind, z, v)]
                    ref = refs[(kind, z, v)]
                    u = max(ulp(g, r) for g, r in zip(got, ref)
                            if mpmath.isfinite(r))
                    cells.append(fmt(u))
                print("%8g " % z + " ".join("%6s" % c for c in cells))
    # Phi
    print("\nPhi: worst relative error in ulp per x bin")
    bins = [(-38.5, -30), (-30, -20), (-20, -10), (-10, -5), (-5, 0),
            (0, 5), (5, 8.5)]
    print("%-14s %10s %10s" % ("x bin", "develop", "PR"))
    for lo, hi in bins:
        wd = wp = 0.0
        for (x, a), (_, b) in zip(phi_dev, phi_pr):
            if lo <= x < hi:
                r = mpmath.ncdf(mpmath.mpf(x))
                if r < mpmath.mpf(2) ** -1022:
                    continue
                wd = max(wd, ulp(a, r))
                wp = max(wp, ulp(b, r))
        print("[%5.1f, %5.1f) %10.1f %10.1f" % (lo, hi, wd, wp))


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
