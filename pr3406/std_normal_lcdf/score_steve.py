"""Claim N5 (about 2.7 x 2^-52 after 1f8b6ba94d): score the ERFCX and LCDF
lines of probe_steve outputs against mpmath.  One route per regime; the lcdf
reference is checked by three routes in score_lcdf_grid.py.

Usage: python3 score_steve.py out_probe_A.txt [out_probe_B.txt ...]
Relative error in units of 2^-52 (ulp-like), worst and median per section
and per z region; also the count of non-finite results where the reference
is finite.
"""
import math
import sys

import mpmath as mp

mp.mp.dps = 60


def ref_erfcx(x):
    x = mp.mpf(x)
    if x > 1e6:
        # asymptotic series; the terms fall by 1/(2 x^2) <= 5e-13 each
        total, term, k = mp.mpf(1), mp.mpf(1), 1
        while True:
            term = -term * (2 * k - 1) / (2 * x * x)
            if abs(term) < mp.mpf(10) ** -70:
                break
            total += term
            k += 1
        return total / (x * mp.sqrt(mp.pi))
    return mp.exp(x * x) * mp.erfc(x)


def ref_lcdf(z):
    z = mp.mpf(z)
    if z < 0:
        # Phi(z) = erfcx(-z / sqrt 2) exp(-z^2 / 2) / 2
        e = ref_erfcx(-z / mp.sqrt(2))
        value = -z * z / 2 + mp.log(e / 2)
        slope = mp.sqrt(2 / mp.pi) / e
        return value, slope
    # upper side: 1 - Phi(z) = Phi(-z) from erfcx, so no digits are lost
    q = ref_erfcx(z / mp.sqrt(2)) * mp.exp(-z * z / 2) / 2
    value = mp.log1p(-q)
    slope = mp.npdf(z) / (1 - q)
    return value, slope


def rel_ulp(got, ref):
    if ref == 0:
        return 0.0 if got == 0 else math.inf
    if not math.isfinite(got):
        return math.inf
    return float(abs((mp.mpf(got) - ref) / ref)) / 2.0**-52


def region(z):
    if z <= -4 * math.sqrt(2):
        return "z<=-5.66"
    if z < 0:
        return "-5.66<z<0"
    return "z>=0"


def score(path, cache):
    rows = {}
    for line in open(path):
        parts = line.split()
        if parts[0] == "ERFCX":
            x, v = float(parts[1]), float(parts[2])
            key = ("erfcx", x)
            if key not in cache:
                cache[key] = ref_erfcx(x)
            rows.setdefault("erfcx", []).append(rel_ulp(v, cache[key]))
        elif parts[0] == "LCDF":
            z, v, s = (float(p) for p in parts[1:4])
            key = ("lcdf", z)
            if key not in cache:
                cache[key] = ref_lcdf(z)
            rv, rs = cache[key]
            reg = region(z)
            rows.setdefault("value " + reg, []).append(rel_ulp(v, rv))
            rows.setdefault("slope " + reg, []).append(rel_ulp(s, rs))
    return rows


def main():
    cache = {}
    results = [(p, score(p, cache)) for p in sys.argv[1:]]
    keys = sorted(results[0][1])
    print("%-22s" % "section" + "".join("%28s" % p.split("_")[-1] for p, _ in results))
    for k in keys:
        cells = []
        for _, rows in results:
            vals = rows[k]
            finite = [v for v in vals if math.isfinite(v)]
            nonfin = len(vals) - len(finite)
            worst = max(finite) if finite else float("nan")
            med = sorted(finite)[len(finite) // 2] if finite else float("nan")
            cells.append("%9.3g / %7.3g / nf %d" % (worst, med, nonfin))
        print("%-22s" % k + "".join("%28s" % c for c in cells))


if __name__ == "__main__":
    main()
