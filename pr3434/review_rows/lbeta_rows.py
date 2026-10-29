# References for the large-argument lbeta derivative test (two routes).
from mpmath import mp, mpf, loggamma, digamma, diff
mp.dps = 120
def fmt(x):
    return "%.17g" % float(x)
rows = [(0.5, 1e14), (1.5, 1e14), (0.5, 1e8), (12.5, 1e10)]
for a, b in rows:
    A, B = mpf(a), mpf(b)
    lb = lambda u, v: loggamma(u) + loggamma(v) - loggamma(u + v)
    v = lb(A, B)
    da = digamma(A) - digamma(A + B); db = digamma(B) - digamma(A + B)
    da2 = diff(lambda u: lb(u, B), A); db2 = diff(lambda w: lb(A, w), B)
    ag = max(abs(da / da2 - 1), abs(db / db2 - 1))
    print(f"    {{{float(a).hex()}, {float(b).hex()}, {fmt(v)}, {fmt(da)}, {fmt(db)}}},  // agree {float(ag):.0e}")
