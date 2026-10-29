from mpmath import mp, mpf, loggamma, digamma, log1p, exp, log, factorial, diff
mp.dps = 90
def fmt(x):
    s = "%.17g" % float(x)
    if "e" not in s and "." not in s: s += ".0"
    return s
def ys(tag, n, a):
    a = mpf(a)
    cA = exp(log(n) + loggamma(n) + loggamma(a + 1) - loggamma(n + a + 1))
    dA = cA * (digamma(a + 1) - digamma(n + a + 1))
    cB = factorial(n); s = mpf(0)
    for j in range(1, n + 1):
        cB /= (a + j); s += 1 / (a + j)
    dB = -cB * s
    ag = max(abs(cA / cB - 1), abs(dA / dB - 1))
    if tag == "YSC": return 1 - cB, -dB, ag
    return log1p(-cB), -dB / (1 - cB), ag
for tag in ["YSC", "YSLC"]:
    for n, a in [(1, 1e12), (2, 1e8)]:
        v, d, ag = ys(tag, n, a)
        print(f'    {{"{tag}", {{{float(n).hex()}, {float(a).hex()}}}, {fmt(v)}, {{{fmt(d)}}}}},  // routes agree to {float(ag):.0e}')
# NB: n = 0, alpha = 1e10 + 1, beta = 0.8
a0, b0 = mpf(1e10 + 1), mpf(float.fromhex("0x1.999999999999ap-1"))
def lp(a, b, n=0):
    return loggamma(n + a) - loggamma(a) - loggamma(n + 1) + a * log(b / (1 + b)) - n * log1p(b)
vA = a0 * log(b0 / (1 + b0)); daA = log(b0 / (1 + b0)); dbA = a0 / (b0 * (1 + b0))
mp.dps = 140
vB = lp(a0, b0); daB = diff(lambda a: lp(a, b0), a0); dbB = diff(lambda b: lp(a0, b), b0)
ag = max(abs(vA / vB - 1), abs(daA / daB - 1), abs(dbA / dbB - 1))
print(f'    {{"NB", {{0x0.0p+0, {float(a0).hex()}, {float(b0).hex()}}}, {fmt(vA)}, {{{fmt(daA)}, {fmt(dbA)}}}}},  // routes agree to {float(ag):.0e}')
