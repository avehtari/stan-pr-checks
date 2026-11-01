#!/usr/bin/env python3
"""Score the output of probe_caller.cpp against mpmath at 50 digits.

beta_lcdf(x | a, b) = log I_x(a, b); its x derivative is the beta density
divided by I_x(a, b).  Route 1 uses mp.betainc.  Route 2 uses the binomial
sum for integer a and b: I_x(a, b) = P(Binomial(a + b - 1, x) >= a), here
1 - (1 - x)^(b + 1) - (b + 1) x (1 - x)^b for a = 2.

Usage: python3 score_caller.py < probe_output.txt
"""
import sys

import mpmath as mp


def refs(a, b, x):
    with mp.workdps(50):
        a, b, x = mp.mpf(a), mp.mpf(b), mp.mpf(x)
        dens = x ** (a - 1) * (1 - x) ** (b - 1) / mp.beta(a, b)
        i1 = mp.betainc(a, b, 0, x, regularized=True)
        assert a == 2
        i2 = 1 - (1 - x) ** (b + 1) - (b + 1) * x * (1 - x) ** b
        agree = abs(i1 - i2) / i1
        return mp.log(i1), dens / i1, agree


def main():
    print("# a b b*x  rel_err_lcdf  rel_err_d_dx  routes_agree")
    worst = [0.0, 0.0, 0.0]
    for line in sys.stdin:
        if line.startswith("#"):
            continue
        a, b, x, lcdf, ddx = (float.fromhex(t) for t in line.split())
        r_lcdf, r_ddx, agree = refs(a, b, x)
        with mp.workdps(50):
            e1 = float(abs((mp.mpf(lcdf) - r_lcdf) / r_lcdf))
            e2 = float(abs((mp.mpf(ddx) - r_ddx) / r_ddx))
        worst = [max(worst[0], e1), max(worst[1], e2),
                 max(worst[2], float(agree))]
        print("%g %.10g %g  %.2e  %.2e  %.1e" % (a, b, b * x, e1, e2,
                                                 float(agree)))
    print("# worst: lcdf %.2e  d_dx %.2e  routes %.1e" % tuple(worst))


if __name__ == "__main__":
    main()
