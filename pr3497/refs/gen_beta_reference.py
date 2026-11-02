#!/usr/bin/env python3
"""Reference rows of beta(a, b) at large arguments for the unit tests.

Routes: mp.beta(a, b), and exp(loggamma(a) + loggamma(b) - loggamma(a + b)).
For a row with the argument 2, also the closed form 1 / (b (b + 1)).  The
precision grows with the arguments (60 + 2 log10(max(a, b)) digits):
loggamma(1.7e308) is about 1.2e311, so the loggamma route needs more than
311 digits there, and mp.beta at 50 digits returns 1.77 for
beta(0.5, 1.7e308).  The script prints the largest relative difference of
the routes to stderr.

Usage:
  python3 gen_beta_reference.py                 # print the C++ header
  python3 gen_beta_reference.py --check FILE    # compare the rows of FILE
"""
import math
import re
import sys

import mpmath as mp

A_VALUES = [0.5, 2.0, 10.0]
B_VALUES = [1e6, 1e9, 2.0**31]
HUGE = [(0.5, 1.7e308)]

HEADER = """\
#ifndef TEST_UNIT_MATH_PRIM_FUN_BETA_REFERENCE_HPP
#define TEST_UNIT_MATH_PRIM_FUN_BETA_REFERENCE_HPP

#include <array>
#include <vector>

namespace stan {
namespace test {
// Reference values of beta(a, b) from mpmath (mp.beta, checked against
// exp(loggamma(a) + loggamma(b) - loggamma(a + b)) and, for the argument 2,
// against 1 / (b (b + 1))). Columns: a, b, beta(a, b). One argument is
// large, in both argument orders.
inline const std::vector<std::array<double, 3>> beta_large_reference{{
"""

FOOTER = """\
}};
}  // namespace test
}  // namespace stan

#endif
"""


def points():
    for a in A_VALUES:
        for b in B_VALUES:
            yield a, b
            yield b, a
    for a, b in HUGE:
        yield a, b
        yield b, a


def reference(a, b):
    dps = 60 + 2 * int(math.log10(max(a, b, 1.0)))
    with mp.workdps(dps):
        ma, mb = mp.mpf(a), mp.mpf(b)
        r1 = mp.beta(ma, mb)
        r2 = mp.exp(mp.loggamma(ma) + mp.loggamma(mb) - mp.loggamma(ma + mb))
        rel = abs(r1 - r2) / abs(r1)
        if 2.0 in (a, b):
            other = mb if a == 2.0 else ma
            r3 = 1 / (other * (other + 1))
            rel = max(rel, abs(r1 - r3) / abs(r1))
        return float(r1), float(rel)


def main():
    worst = 0.0
    rows = []
    for a, b in points():
        v, rel = reference(a, b)
        worst = max(worst, rel)
        rows.append((a, b, v))
    sys.stderr.write("largest relative difference of the routes: %.3g\n"
                     % worst)
    if len(sys.argv) == 3 and sys.argv[1] == "--check":
        text = open(sys.argv[2]).read()
        found = re.findall(r"\{([-+0-9.e]+),\s*([-+0-9.e]+),\s*([-+0-9.e]+)\}",
                           text)
        got = [tuple(float(x) for x in t) for t in found]
        bad = [r for r in rows if r not in got]
        print("rows in %s: %d; generated rows: %d; generated rows not found: %d"
              % (sys.argv[2], len(got), len(rows), len(bad)))
        for r in bad:
            print("missing", r)
        sys.exit(1 if bad else 0)
    out = [HEADER]
    for a, b, v in rows:
        out.append("    {%r, %r, %r},\n" % (a, b, v))
    out.append(FOOTER)
    sys.stdout.write("".join(out))


if __name__ == "__main__":
    main()
