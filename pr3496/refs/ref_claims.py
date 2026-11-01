#!/usr/bin/env python3
# mpmath references for the numbers that the PR text quotes, two routes each:
# neg_binomial_lccdf(3 | 2.5, 1.5) and neg_binomial_2_lcdf(3 | 4, 2) with its
# mu and phi derivatives; poisson_binomial_lccdf(3 | 0.5, 0.2, 0.7).
import mpmath as mp

mp.mp.dps = 40


def nb_pmf(k, alpha, p):
    # Negative binomial with shape alpha and success probability p
    return mp.binomial(k + alpha - 1, k) * p**alpha * (1 - p) ** k


def report(name, a, b):
    rel = abs(a - b) / abs(b) if b != 0 else abs(a - b)
    print(f"{name:52s} route A {mp.nstr(a, 15):>22s}  route B "
          f"{mp.nstr(b, 15):>22s}  rel diff {mp.nstr(rel, 3)}")


# neg_binomial_lccdf(3 | alpha = 2.5, beta = 1.5): p = beta / (1 + beta)
alpha, beta, n = mp.mpf("2.5"), mp.mpf("1.5"), 3
p = beta / (1 + beta)
a = mp.log(1 - mp.fsum(nb_pmf(k, alpha, p) for k in range(n + 1)))
b = mp.log(mp.betainc(n + 1, alpha, 0, 1 - p, regularized=True))
report("neg_binomial_lccdf(3 | 2.5, 1.5)", a, b)


# neg_binomial_2_lcdf(3 | mu = 4, phi = 2): p = phi / (mu + phi)
def nb2_lcdf_sum(mu, phi):
    q = phi / (mu + phi)
    return mp.log(mp.fsum(nb_pmf(k, phi, q) for k in range(n + 1)))


def nb2_lcdf_beta(mu, phi):
    return mp.log(mp.betainc(phi, n + 1, 0, phi / (mu + phi), regularized=True))


mu, phi = mp.mpf(4), mp.mpf(2)
report("neg_binomial_2_lcdf(3 | 4, 2)", nb2_lcdf_sum(mu, phi),
       nb2_lcdf_beta(mu, phi))
report("d/dmu neg_binomial_2_lcdf(3 | 4, 2)",
       mp.diff(lambda m: nb2_lcdf_sum(m, phi), mu),
       mp.diff(lambda m: nb2_lcdf_beta(m, phi), mu))
report("d/dphi neg_binomial_2_lcdf(3 | 4, 2)",
       mp.diff(lambda f: nb2_lcdf_sum(mu, f), phi),
       mp.diff(lambda f: nb2_lcdf_beta(mu, f), phi))

# poisson_binomial_lccdf(3 | 0.5, 0.2, 0.7): with N = 3 trials, Y <= 3, so
# P(Y <= 3) = 1 and the lccdf at y = 3 is log(0) = -inf. Route A: the
# probabilities of Y = 0..3 by the recursion over the trials. Route B: by
# enumeration of the 8 outcomes.
theta = [mp.mpf("0.5"), mp.mpf("0.2"), mp.mpf("0.7")]
probs = [mp.mpf(1)]
for t in theta:
    nxt = [mp.mpf(0)] * (len(probs) + 1)
    for k, pk in enumerate(probs):
        nxt[k] += pk * (1 - t)
        nxt[k + 1] += pk * t
    probs = nxt
enum = [mp.mpf(0)] * 4
for b0 in (0, 1):
    for b1 in (0, 1):
        for b2 in (0, 1):
            enum[b0 + b1 + b2] += ((theta[0] if b0 else 1 - theta[0])
                                   * (theta[1] if b1 else 1 - theta[1])
                                   * (theta[2] if b2 else 1 - theta[2]))
for k in range(4):
    report(f"P(Y = {k}), poisson_binomial(0.5, 0.2, 0.7)", probs[k], enum[k])
report("P(Y <= 3)", mp.fsum(probs), mp.fsum(enum))
print("poisson_binomial_lccdf(3 | 0.5, 0.2, 0.7) = log(1 - P(Y <= 3)) = log(0) = -inf")
print(f"log(2^-53) = {mp.nstr(mp.log(mp.mpf(2) ** -53), 17)}")
