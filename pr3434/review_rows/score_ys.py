from mpmath import mp, mpf, loggamma, digamma, log1p, exp, log, factorial
mp.dps = 80
def ref(tag, n, a):
    a = mpf(a)
    lA = log(n) + loggamma(n) + loggamma(a + 1) - loggamma(n + a + 1)
    cA = exp(lA); dA = cA * (digamma(a + 1) - digamma(n + a + 1))
    cB = factorial(n); s = mpf(0)
    for j in range(1, n + 1):
        cB /= (a + j); s += 1 / (a + j)
    dB = -cB * s
    agree = max(abs(cA / cB - 1), abs(dA / dB - 1))
    if tag == "YSC":
        return 1 - cB, -dB * a, agree
    return log1p(-cB), (-dB / (1 - cB)) * a, agree
rows = {}
for c in ["510d263078", "64407093ba"]:
    for l in open(f"out/review_rows_{c}.txt"):
        p = l.split()
        if p[0] == "NB": continue
        rows.setdefault((p[0], int(p[1]), float.fromhex(p[2])), {})[c] = (float.fromhex(p[3]), float.fromhex(p[4]))
print("tag   n  alpha   value      g=a*d/da   route   dev v.err/tol dev g.err/tol  PR v/tol  PR g/tol  dev g.rel")
for (tag, n, a), d in rows.items():
    v, g, agree = ref(tag, n, a)
    tv, tg = 1e-12 * max(1, abs(v)), 1e-11 * max(1, abs(g))
    e = {c: (abs(d[c][0] - v) / tv, abs(d[c][1] - g) / tg) for c in d}
    rel = abs(d["510d263078"][1] - g) / abs(g)
    print(f"{tag:5s} {n:2d} {a:7.0e} {float(v):10.3e} {float(g):10.3e} {float(agree):7.0e} {float(e['510d263078'][0]):10.2e} {float(e['510d263078'][1]):10.2e} {float(e['64407093ba'][0]):9.2e} {float(e['64407093ba'][1]):9.2e} {float(rel):9.2e}")
for c in ["510d263078", "64407093ba"]:
    print(c, [l for l in open(f"out/review_rows_{c}.txt") if l.startswith("NB")][0].strip())
