#!/usr/bin/env python3
"""Claim N6: recompute the reference constants of the #3406 tests by two
independent mpmath routes and compare them with the constants in the files.
Files (read from TREE): std_normal_lcdf_impl_test.cpp, Phi_test.cpp,
normal_lcdf_tail_test_helpers.hpp (expect_branch_accuracy).

Route A: erfc, log(erfc(-z/sqrt2)/2) or log1p(-erfc(z/sqrt2)/2), at 60 digits
plus the digits that log1p cancels.  Route B: quadrature of the density on
the small side, Phi(-a) = phi(a) M, M = int_0^inf exp(-a u - u^2/2) du, at
45 digits (asymptotic series of M for a > 1e6).  Slope: route A
exp(lpdf - lcdf) with extra digits, route B 1 / M (z < 0).
Columns: the constant, route A rounded to double, |B/A - 1| (route agreement),
|constant/A - 1| (test constant against the reference).
Usage: gen_refs.py TREE
"""
import os
import re
import sys

import mpmath as mp

DPS_A = 60
DPS_B = 45


def small_tail_quad(a):
    a = abs(a)
    step = 1 / (a + 1)
    m = mp.quad(lambda u: mp.exp(-a * u - u * u / 2),
                [0, step, 4 * step, 16 * step, mp.inf])
    return mp.npdf(a) * m


def phi_routes(x):
    """Phi(x) by route A and route B; for x > 0 also 1 - Phi(x)."""
    x = mp.mpf(x)
    with mp.workdps(DPS_A):
        a = mp.erfc(-x / mp.sqrt(2)) / 2
    with mp.workdps(DPS_B):
        q = small_tail_quad(x)
        b = q if x < 0 else 1 - q
    return a, b


def lcdf_routes(z):
    """(value, slope) of log Phi(z) by route A and by route B.

    Both routes form the small tail t = Phi(-|z|) (no cancellation; mpmath
    numbers do not underflow), then value log(t) or log1p(-t) and slope
    phi(z) / Phi(z), without exp(lpdf - lcdf), which cancels for large |z|.
    """
    z = mp.mpf(z)
    a = abs(z)
    # Extra digits for exp(lpdf - lcdf) in route A when |z| is large.
    extra = int(2 * mp.log10(a + 1)) + 10
    with mp.workdps(DPS_A + extra):
        lpdf = -z * z / 2 - mp.log(mp.sqrt(2 * mp.pi))
        t = mp.erfc(a / mp.sqrt(2)) / 2
        va = mp.log(t) if z < 0 else mp.log1p(-t)
        sa = mp.exp(lpdf - va)
        route_a = (+va, +sa)
    with mp.workdps(DPS_B):
        # Phi(-a) = phi(a) M: log value and slope 1 / M need no subtraction.
        if a > 1e6:
            # Asymptotic series of the Mills ratio; terms fall by 1 / a^2.
            m, term, k = mp.mpf(1), mp.mpf(1), 1
            while abs(term) > mp.eps:
                term = -term * (2 * k - 1) / (a * a)
                m += term
                k += 1
            m /= a
        else:
            step = 1 / (a + 1)
            m = mp.quad(lambda u: mp.exp(-a * u - u * u / 2),
                        [0, step, 4 * step, 16 * step, mp.inf])
        if z < 0:
            route_b = (lpdf + mp.log(m), 1 / m)
        else:
            t = mp.exp(lpdf) * m
            route_b = (mp.log1p(-t), mp.exp(lpdf) / (1 - t))
    # Keep the full precision: unary + here would round to the default 15 digits.
    return [route_a, route_b]


def rel(a, b):
    with mp.workdps(DPS_A):
        if b == 0:
            return mp.mpf(0) if a == 0 else mp.inf
        return abs(mp.mpf(a) / b - 1)


def row(label, const, ref_a, ref_b):
    print("%-26s %24r %24r %9.1e %9.1e" % (
        label, const, float(ref_a), float(rel(ref_b, ref_a)), float(rel(const, ref_a))))


def header(title):
    print("\n" + title)
    print("%-26s %24s %24s %9s %9s" % ("input", "constant", "route A", "B vs A",
                                        "const vs A"))


NUM = r"[-+]?[0-9.]+(?:e[-+]?[0-9]+)?"


def lcdf_impl(tree):
    path = os.path.join(tree, "test/unit/math/prim/prob/std_normal_lcdf_impl_test.cpp")
    text = open(path).read()
    block = text[text.index("scalar_tail_kernel"):]
    block = block[:block.index("};")]
    rows = re.findall(r"\{(%s),\s*(%s),\s*(%s)\}" % (NUM, NUM, NUM), block)
    header("std_normal_lcdf_impl_test.cpp scalar_tail_kernel (%d rows): "
           "value, then slope" % len(rows))
    for z, v, s in rows:
        (va, sa), (vb, sb) = lcdf_routes(float(z))
        row("lcdf(%s)" % z, float(v), va, vb)
        row("slope(%s)" % z, float(s), sa, sb)


def phi_test(tree):
    path = os.path.join(tree, "test/unit/math/prim/fun/Phi_test.cpp")
    text = open(path).read()
    block = text[text.index("TEST(MathFunctions, Phi)"):]
    block = block[:block.index("\n}\n")]
    pairs = [(x, c) for c, x in re.findall(
        r"EXPECT_NEAR\((%s),\s*stan::math::Phi\((%s)\)\s*," % (NUM, NUM), block)]
    pairs += re.findall(r"Phi\((%s)\)\s*/\s*(%s)" % (NUM, NUM), block)
    header("Phi_test.cpp TEST(MathFunctions, Phi) (%d rows)" % len(pairs))
    for x, c in pairs:
        a, b = phi_routes(float(x))
        row("Phi(%s)" % x, float(c), a, b)


def branch_accuracy(tree):
    path = os.path.join(tree, "test/unit/math/mix/prob/normal_lcdf_tail_test_helpers.hpp")
    text = open(path).read()
    block = text[text.index("void expect_branch_accuracy"):]
    block = block[block.index("cases[]"):]
    block = block[:block.index("};")]
    rows = re.findall(r"\{(%s),\s*(%s)\}" % (NUM, NUM), block)
    header("normal_lcdf_tail_test_helpers.hpp expect_branch_accuracy (%d rows): "
           "d/dy log Phi(y)" % len(rows))
    for y, d in rows:
        (_, sa), (_, sb) = lcdf_routes(float(y))
        row("slope(%s)" % y, float(d), sa, sb)


def main():
    tree = sys.argv[1]
    print("tree: %s" % tree)
    lcdf_impl(tree)
    phi_test(tree)
    branch_accuracy(tree)


if __name__ == "__main__":
    main()
