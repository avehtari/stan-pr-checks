#!/usr/bin/env python3
"""References for the range-end unit tests (mu = 0, sigma = 1, lambda = 1,
alpha = 1).  Reads probe_test_points output files; prints the reference
value and d/dy by two routes of log Phi, their agreement, and the relative
error of each file.

Usage: python3 ref_test_points.py out_dev.txt out_head.txt
"""
import sys
import mpmath as mp

# 700 digits: the gradient subtracts logs of size up to 1.5e616
mp.mp.dps = 700


def log_erfcx_cf(x):
    # Laplace continued fraction: sqrt(pi) erfcx(x) = 1/(x+ (1/2)/(x+ 1/(x+ ...)))
    t = x
    for k in range(400, 0, -1):
        t = x + (mp.mpf(k) / 2) / t
    return -mp.log(t) - mp.log(mp.pi) / 2


def log_phi_a(z):
    if z < -1e50:
        x = -z / mp.sqrt(2)
        return log_erfcx_cf(x) - x * x - mp.log(2)
    if z > 1e50:
        x = z / mp.sqrt(2)
        return mp.log1p(-mp.exp(log_erfcx_cf(x) - x * x) / 2)
    return mp.log(mp.erfc(-z / mp.sqrt(2)) / 2)


def log_phi_b(z):
    # asymptotic series below -60 (its smallest term, about exp(-z^2/2), is
    # then below 10^-700), quadrature of the density above
    if z <= -60:
        s, term, k = mp.mpf(1), mp.mpf(1), 1
        while k < 5000:
            nxt = term * (-(2 * k - 1) / (z * z))
            if abs(nxt) < mp.mpf(10) ** (-mp.mp.dps) or abs(nxt) > abs(term):
                break
            term = nxt
            s += term
            k += 1
        return -z * z / 2 - mp.log(-z) - mp.log(2 * mp.pi) / 2 + mp.log(s)
    if z <= -8:
        # continued fraction; mp.quad stops on an absolute tolerance and is
        # not accurate for a density this small
        x = -z / mp.sqrt(2)
        return log_erfcx_cf(x) - x * x - mp.log(2)
    if z >= 8:
        return mp.log1p(-mp.exp(log_phi_b(-z)))
    with mp.workdps(60):
        if z <= 0:
            return +mp.log(mp.quad(mp.npdf, [-mp.inf, z]))
        return +mp.log1p(-mp.quad(mp.npdf, [z, mp.inf]))


def log_npdf(z):
    return -z * z / 2 - mp.log(2 * mp.pi) / 2


def refs(name, y, lphi):
    y = mp.mpf(y)
    if name == "normal_lcdf":
        v = lphi(y)
        return v, mp.exp(log_npdf(y) - v)
    if name == "normal_lccdf":
        v = lphi(-y)
        return v, -mp.exp(log_npdf(y) - v)
    if name == "lognormal_lcdf":
        z = mp.log(y)
        v = lphi(z)
        return v, mp.exp(log_npdf(z) - v) / y
    if name == "lognormal_lccdf":
        z = mp.log(y)
        v = lphi(-z)
        return v, -mp.exp(log_npdf(z) - v) / y
    if name == "skew_normal_lpdf":
        v = mp.log(2) + log_npdf(y) + lphi(y)
        return v, -y + mp.exp(log_npdf(y) - lphi(y))
    # exp_mod_normal, lambda = sigma = 1: log f = 1/2 - y + log Phi(y - 1)
    with mp.workdps(700):
        log_f = mp.mpf(1) / 2 - y + lphi(y - 1)
        log_b = log_f  # b = exp(-y + 1/2) Phi(y - 1) equals f when lambda = 1
        if name == "exp_mod_normal_lpdf":
            return log_f, -1 + mp.exp(log_npdf(y - 1) - lphi(y - 1))
        log_cdf = lphi(y) + mp.log1p(-mp.exp(log_b - lphi(y)))
        log_ccdf = mp.log(mp.exp(lphi(-y)) + mp.exp(log_b))
        if name == "exp_mod_normal_lcdf":
            return log_cdf, mp.exp(log_f - log_cdf)
        if name == "exp_mod_normal_lccdf":
            return log_ccdf, -mp.exp(log_f - log_ccdf)
        if name == "exp_mod_normal_cdf":
            return mp.exp(log_cdf), mp.exp(log_f)
    raise ValueError(name)


def rel(x, r):
    if mp.isnan(x) or mp.isinf(x):
        return "nonfinite"
    if r == 0:
        return f"{float(abs(x)):.1e}abs"
    return f"{float(abs(x - r) / abs(r)):.1e}"


def main():
    files = [open(f).read().split("\n") for f in sys.argv[1:]]
    for lines in zip(*files):
        if not lines[0].strip():
            continue
        name, yh = lines[0].split()[:2]
        y = float.fromhex(yh)
        va, da = refs(name, y, log_phi_a)
        vb, db = refs(name, y, log_phi_b)
        agree_v = float(abs(va - vb) / abs(va)) if va != 0 else 0.0
        agree_d = float(abs(da - db) / abs(da)) if da != 0 else 0.0
        print(f"{name} y={y!r}")
        print(f"  ref value {mp.nstr(va, 20)}  dy {mp.nstr(da, 20)}")
        print(f"  routes agree: value {agree_v:.1e}  dy {agree_d:.1e}")
        for f, line in zip(sys.argv[1:], lines):
            p = line.split()
            v, d = float.fromhex(p[2]), float.fromhex(p[3])
            print(f"  {f}: value {v!r} ({rel(mp.mpf(v), va)})"
                  f"  dy {d!r} ({rel(mp.mpf(d), da)})")


if __name__ == "__main__":
    main()
