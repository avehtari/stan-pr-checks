#!/usr/bin/env python3
# Score the output of probe_nb2_lcdf_large_n against mpmath.
# neg_binomial_2_lcdf(n | mu, phi) = log I_x(phi, n + 1), x = phi / (mu + phi).
# Route A: the lower incomplete beta; where its series does not converge
# (n much larger than mu), log1p(-S) with S the direct sum of the pmf over
# k > n. Route B: log1p(-I_{1-x}(n + 1, phi)), the upper incomplete beta.
# For the values at phi = 2 the routes are independent: route B is then a
# polynomial of two terms (hyp2f1 with 1 - phi = -1), route A a series or a
# direct sum. The phi derivative (mp.diff) also evaluates route B at
# phi +- h, where it is a general hyp2f1.
# Usage: python3 score_large_n.py < probe output
import sys
import mpmath as mp

mp.mp.dps = 60


def tail_sum(n, mu, phi):
    # The term ratio r_k = (k + phi) / (k + 1) (1 - q) decreases to 1 - q
    # for phi > 1, so the remainder after a term is at most
    # term r_k / (1 - r_k). Stop when that bound is below 1e-55 of the sum.
    q = phi / (mu + phi)
    k = n + 1
    term = mp.exp(phi * mp.log(q) - mp.loggamma(phi) + mp.loggamma(k + phi)
                  - mp.loggamma(k + 1) + k * mp.log(1 - q))
    total = mp.mpf(0)
    while True:
        total += term
        r = (k + phi) / (k + 1) * (1 - q)
        term *= r
        k += 1
        if term * r / (1 - r) < total * mp.mpf(10) ** -55:
            return total + term


def route_a(n, mu, phi, use_sum):
    if use_sum:
        return mp.log1p(-tail_sum(n, mu, phi))
    return mp.log(mp.betainc(phi, n + 1, 0, phi / (mu + phi), regularized=True))


def route_b(n, mu, phi):
    return mp.log1p(-mp.betainc(n + 1, phi, 0, mu / (mu + phi), regularized=True))


def rel(x, ref):
    if ref == 0:
        return abs(x)
    return abs(x - ref) / abs(ref)


print(f"{'n':>11s} {'mu':>11s}  {'quantity':8s} {'Stan':>23s} {'mpmath (B)':>23s}"
      f" {'abs err':>9s} {'rel err':>9s} {'A vs B':>9s}")
for line in sys.stdin:
    f = line.split()
    if len(f) != 6:
        continue
    n = int(f[0])
    mu, phi = mp.mpf(f[1]), mp.mpf(f[2])
    stan = [float(x) for x in f[3:6]]
    use_sum = n > 100 * mu
    ref_b = [route_b(n, mu, phi),
             mp.diff(lambda m: route_b(n, m, phi), mu),
             mp.diff(lambda p: route_b(n, mu, p), phi)]
    ref_a = [route_a(n, mu, phi, use_sum),
             mp.diff(lambda m: route_a(n, m, phi, use_sum), mu),
             mp.diff(lambda p: route_a(n, mu, p, use_sum), phi)]
    for name, s, a, b in zip(("value", "d/dmu", "d/dphi"), stan, ref_a, ref_b):
        print(f"{n:11d} {mp.nstr(mu, 6):>11s}  {name:8s} {s:23.15g}"
              f" {mp.nstr(b, 15):>23s} {mp.nstr(abs(mp.mpf(s) - b), 3):>9s}"
              f" {mp.nstr(rel(mp.mpf(s), b), 3):>9s} {mp.nstr(rel(a, b), 2):>9s}")
