#!/usr/bin/env python3
"""Claim M5 of README_emg.md: the 6 deep-tail reference rows of the
exp_mod_normal lcdf / lccdf unit test (exp_mod_normal_lcdf_reference.hpp).

Same row format as ref_emg_lcdf_tests.py:
  {y, mu, sigma, lambda, lcdf, dy, dmu, dsigma, dlambda,
                         lccdf, dy, dmu, dsigma, dlambda},
at the exact double inputs, mu = 0.5, sigma = 2.

Two routes, both at 800 digits:
  direct:  F = Phi(z) - exp(v^2 / 2 - v z) Phi(z - v)
  erfcx:   F = exp(-z^2 / 2) / 2 (E(u) - E(t)), u = -z / sqrt 2,
           t = u + v / sqrt 2, E = erfcx from its asymptotic series
Partials by mp.diff of log F and log(1 - F).  The script prints the largest
relative difference between the routes where both apply.

Usage: ref_emg_lcdf_deep_tail.py > rows.txt
"""
import sys

import mpmath as mp

mp.mp.dps = 800

MU = 0.5
SIGMA = 2.0

# (z, v, comment); log F stays finite in double for |z| <= 1.8e154
POINTS = [
    (-1e5, 1e-300, "deep tail, Delta subnormal"),
    (-1e20, 1e-300, "deep tail, Delta subnormal"),
    (-1e100, 1e-300, "deep tail, 1 - r underflows"),
    (-1e150, 1e-10, "deep tail, erfcx' near underflow"),
    (-1e154, 0.1, "deep tail, erfcx' underflows"),
    (-1e154, 1.0, "deep tail, erfcx' underflows"),
]


def erfcx_asym(x):
    """erfcx(x) for x >= 1e4 from the asymptotic series."""
    total, term, k = mp.mpf(1), mp.mpf(1), 1
    eps = mp.mpf(10) ** (-mp.mp.dps)
    while True:
        term = -term * (2 * k - 1) / (2 * x * x)
        if abs(term) < eps:
            break
        total += term
        k += 1
    return total / (x * mp.sqrt(mp.pi))


def log_f_erfcx(y, mu, sigma, lam):
    z = (y - mu) / sigma
    v = lam * sigma
    u = -z / mp.sqrt(2)
    t = u + v / mp.sqrt(2)
    delta = erfcx_asym(u) - erfcx_asym(t)
    return -z * z / 2 + mp.log(delta / 2), delta * mp.exp(-z * z / 2) / 2


def log_f_direct(y, mu, sigma, lam):
    z = (y - mu) / sigma
    v = lam * sigma
    f = mp.ncdf(z) - mp.exp(v * v / 2 - v * z) * mp.ncdf(z - v)
    return mp.log(f), f


def row(route, y, mu, sigma, lam):
    args = [mp.mpf(a) for a in (y, mu, sigma, lam)]
    out = []
    for upper in (False, True):
        def logf(*a, upper=upper):
            lf, f = route(*a)
            return mp.log1p(-f) if upper else lf
        out.append(logf(*args))
        for k in range(4):
            def g(x, k=k):
                a = list(args)
                a[k] = x
                return logf(*a)
            out.append(mp.diff(g, args[k]))
    return out


def fmt(x):
    return repr(float(x))


def main():
    worst = mp.mpf(0)
    for z, v, comment in POINTS:
        y = MU + SIGMA * z
        lam = v / SIGMA
        vals = row(log_f_erfcx, y, MU, SIGMA, lam)
        # The direct route also works beyond |z| = 1e20 at 800 digits: the
        # cancellation in F loses about -log10(v / |z|) <= 400 digits here.
        vals_d = row(log_f_direct, y, MU, SIGMA, lam)
        row_worst = mp.mpf(0)
        for a, b in zip(vals, vals_d):
            if a != 0:
                row_worst = max(row_worst, abs((a - b) / a))
        worst = max(worst, row_worst)
        print(f"z = {z:g}, v = {v:g}: largest relative difference of the "
              f"routes {mp.nstr(row_worst, 3)}", file=sys.stderr)
        nums = [fmt(a) for a in (y, MU, SIGMA, lam)] + [fmt(a) for a in vals]
        print(f"    // {comment}: z = {z:g}, v = {v:g}")
        print("    {" + ", ".join(nums) + "},")
    print(f"largest relative difference of the two routes: {mp.nstr(worst, 3)}",
          file=sys.stderr)


if __name__ == "__main__":
    main()
