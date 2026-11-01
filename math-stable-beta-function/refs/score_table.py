#!/usr/bin/env python3
"""Relative error of beta(a, b) of two probe_table.cpp outputs.

The reference is reference() of gen_beta_reference.py (two or three mpmath
routes; the precision grows with the arguments).

Usage: python3 score_table.py DEVELOP_OUT BRANCH_OUT
"""
import math
import sys

import mpmath as mp

from gen_beta_reference import reference


def num(t):
    try:
        return float.fromhex(t)
    except ValueError:
        return float(t)


def read(path):
    return [tuple(num(t) for t in l.split()) for l in open(path)]


def rel(v, r):
    if not math.isfinite(v):
        return str(v)
    with mp.workdps(50):
        return "%.2g" % float(abs(mp.mpf(v) - mp.mpf(r)) / mp.mpf(r))


def main():
    d, b = read(sys.argv[1]), read(sys.argv[2])
    print("| a | b | develop | branch |")
    print("|---|---|---|---|")
    for (a, bb, vd), (_, _, vb) in zip(d, b):
        r, _ = reference(a, bb)
        print("| %g | %g | %s | %s |" % (a, bb, rel(vd, r), rel(vb, r)))


if __name__ == "__main__":
    main()
