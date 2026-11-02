#!/usr/bin/env python3
"""Reference values for the mpmath test in test/unit/math/mix/prob/normal_test.cpp
(normal_lpdf_test::check_reference).

Value (propto = false), d/dmu and d/dsigma of sum_i log N(y_i | mu_i, sigma_i)
at the double inputs of the test, with 40 digits; printed with 17
significant digits (enough for a correctly rounded double).
d/dy = -d/dmu.  The propto = true value adds N log(sqrt(2 pi)).
"""
import mpmath as mp

mp.mp.dps = 40


def f(s):
    """Exact value of the double nearest to the decimal string s."""
    return mp.mpf(float(s))


Y = [f("-0.3"), f("1.2"), f("2.5"), f("0.0")]
MU = [f("0.1"), f("-0.4"), f("2.0"), f("0.7")]
SIGMA_VEC = [f("0.5"), f("1.3"), f("2.2"), f("0.9")]
SIGMA_SCALAR = f("1.7")
C = mp.log(2 * mp.pi) / 2


def show(tag, sigma):
    z = [(a - b) / s for a, b, s in zip(Y, MU, sigma)]
    value = sum(-q * q / 2 - mp.log(s) - C for q, s in zip(z, sigma))
    d_mu = [q / s for q, s in zip(z, sigma)]
    d_sigma = [(q * q - 1) / s for q, s in zip(z, sigma)]
    print(f"{tag} value {mp.nstr(value, 17)}")
    print(f"{tag} d_mu {', '.join(mp.nstr(v, 17) for v in d_mu)}")
    print(f"{tag} d_sigma {', '.join(mp.nstr(v, 17) for v in d_sigma)}")
    print(f"{tag} d_sigma_sum {mp.nstr(sum(d_sigma), 17)}")


print(f"log_sqrt_two_pi {mp.nstr(C, 20)}")
show("scalar_sigma", [SIGMA_SCALAR] * 4)
show("vector_sigma", SIGMA_VEC)
