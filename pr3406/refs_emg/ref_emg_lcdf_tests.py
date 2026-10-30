#!/usr/bin/env python3
"""Claim M5 of README_emg.md: the 13 reference rows of the exp_mod_normal
lcdf / lccdf unit test (exp_mod_normal_lcdf_reference.hpp), two routes.

F = Phi(z) - exp(v^2 / 2 - v z) Phi(z - v), z = (y - mu) / sigma,
v = lambda sigma, in mpmath at 80 digits.  Value and the partials with
respect to y, mu, sigma and lambda of log F and log(1 - F), at the exact
double inputs.  One or more points per branch of exp_mod_normal_lcdf_kernel.

Route A (the rows): the direct form above; partials by mp.diff.
Route B (check): F = g (E(u) - E(t)), 1 - F = g (E(-u) + E(t)) with
g = exp(-z^2 / 2) / 2, u = -z / sqrt 2, t = u + v / sqrt 2,
E(x) = exp(x^2) erfc(x); partials from the closed form
F_z = v b, F_v = phi(z) - (v - z) b, b = exp(v^2 / 2 - v z) Phi(z - v).
Route B runs at 120 digits.  The script prints to stderr, per row, the
number of digits to which the routes agree (minimum over the 10 numbers).

Prints C++ initializer rows:
  {y, mu, sigma, lambda, lcdf, dy, dmu, dsigma, dlambda,
                         lccdf, dy, dmu, dsigma, dlambda},
Usage: ref_emg_lcdf_tests.py > rows.txt
"""
import sys

import mpmath as mp

mp.mp.dps = 80

MU = 0.5
SIGMA = 2.0

# (z, v, branch)
POINTS = [
    (-1.0, 1.7, "z <= 0, direct difference"),
    (-1.0, 1e-6, "z <= 0, midpoint series"),
    (-30.0, 1e-6, "z <= 0, Cody tail form"),
    (-100.0, 1.0, "far left tail"),
    (-5.0, 0.1, "z <= 0, moderate"),
    (0.5, 0.01, "z > 0, F < 1/4, t >= 0"),
    (3.0, 1e-3, "z > 0, F < 1/4, t < 0"),
    (0.5, 2.0, "z > 0, F >= 1/4, t >= 0"),
    (3.0, 1.7, "z > 0, F >= 1/4, t < 0"),
    (20.0, 1.7, "z > 0, P = b only"),
    (50.0, 45.0, "lccdf, P near underflow, t < 0"),
    (40.0, 60.0, "lccdf, P near underflow, t >= 0"),
    (0.5, 100.0, "large v"),
]


def cdf(y, mu, sigma, lam):
    z = (y - mu) / sigma
    v = lam * sigma
    return mp.ncdf(z) - mp.exp(v * v / 2 - v * z) * mp.ncdf(z - v)


def ccdf(y, mu, sigma, lam):
    z = (y - mu) / sigma
    v = lam * sigma
    return mp.ncdf(-z) + mp.exp(v * v / 2 - v * z) * mp.ncdf(z - v)


def log_tail(fun, other, *a):
    """log(fun) from the smaller of fun and other = 1 - fun."""
    x = fun(*a)
    return mp.log(x) if x < 0.5 else mp.log1p(-other(*a))


def row(y, mu, sigma, lam):
    args = [mp.mpf(a) for a in (y, mu, sigma, lam)]
    out = []
    for fun, other in ((cdf, ccdf), (ccdf, cdf)):
        def logf(*a, fun=fun, other=other):
            return log_tail(fun, other, *a)
        out.append(logf(*args))
        for k in range(4):
            def g(x, k=k):
                a = list(args)
                a[k] = x
                return logf(*a)
            out.append(mp.diff(g, args[k]))
    return out


def erfcx(x):
    return mp.exp(x * x) * mp.erfc(x)


def row_b(y, mu, sigma, lam):
    """Route B: erfcx form for the values, closed form for the partials."""
    with mp.workdps(120):
        y, mu, sigma, lam = (mp.mpf(a) for a in (y, mu, sigma, lam))
        z = (y - mu) / sigma
        v = lam * sigma
        u = -z / mp.sqrt(2)
        t = u + v / mp.sqrt(2)
        g = mp.exp(-z * z / 2) / 2
        f = g * (erfcx(u) - erfcx(t))
        p = g * (erfcx(-u) + erfcx(t))
        b = mp.exp(v * v / 2 - v * z) * mp.ncdf(z - v)
        f_z = v * b
        f_v = mp.npdf(z) - (v - z) * b
        # dF/dy, dF/dmu, dF/dsigma, dF/dlambda
        df = [f_z / sigma, -f_z / sigma, -z / sigma * f_z + lam * f_v,
              sigma * f_v]
        # log of the side near 1 from log1p of the other side
        log_f = mp.log(f) if f < 0.5 else mp.log1p(-p)
        log_p = mp.log(p) if p < 0.5 else mp.log1p(-f)
        out = [log_f] + [d / f for d in df]
        out += [log_p] + [-d / p for d in df]
        return [+x for x in out]


def digits(a, b):
    if a == b:
        return float("inf")
    scale = max(abs(a), abs(b))
    return float(-mp.log10(abs(a - b) / scale))


def fmt(x):
    return repr(float(x))


def main():
    worst = float("inf")
    for z, v, branch in POINTS:
        y = MU + SIGMA * z
        lam = v / SIGMA
        vals = row(y, MU, SIGMA, lam)
        agree = min(digits(a, b) for a, b in zip(vals, row_b(y, MU, SIGMA, lam)))
        worst = min(worst, agree)
        print("z = %g, v = %g: routes agree to %.1f digits" % (z, v, agree),
              file=sys.stderr)
        nums = [fmt(a) for a in (y, MU, SIGMA, lam)] + [fmt(a) for a in vals]
        print(f"    // {branch}: z = {z}, v = {v}")
        print("    {" + ", ".join(nums) + "},")
    print("smallest agreement of the two routes: %.1f digits" % worst,
          file=sys.stderr)


if __name__ == "__main__":
    main()
