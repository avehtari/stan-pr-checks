#!/usr/bin/env python3
"""Compare the beta(a, b) column of two probe_edges.cpp outputs.

Counts the identical pairs, the pairs whose class (finite, 0, inf, NaN)
changes, and the pairs with another finite value; lists the class changes
and the finite changes with both arguments below 10.  Also checks that
the second file follows its rule column.

Usage: python3 compare_edges.py EDGES_DEVELOP EDGES_BRANCH
"""
import math
import sys


def rows(path):
    return [l.split() for l in open(path) if not l.startswith("#")]


def cls(s):
    v = float(s)
    if math.isnan(v):
        return "NaN"
    if math.isinf(v):
        return "inf" if v > 0 else "-inf"
    return "0" if v == 0 else "finite"


def main():
    d, b = rows(sys.argv[1]), rows(sys.argv[2])
    assert len(d) == len(b)
    same, cls_chg, val_chg, small_chg, rule_mis = 0, [], 0, [], 0
    for x, y in zip(d, b):
        assert x[:2] == y[:2]
        if y[2] != y[4]:
            rule_mis += 1
        if x[2] == y[2]:
            same += 1
        elif cls(x[2]) != cls(y[2]):
            cls_chg.append((x[0], x[1], x[2], y[2]))
        else:
            val_chg += 1
            if float(x[0]) < 10 and float(x[1]) < 10:
                small_chg.append((x[0], x[1]))
    print("pairs %d; identical %d; class changes %d; other finite value "
          "changes %d; of these with both arguments below 10: %d; branch "
          "differs from its rule column: %d"
          % (len(d), same, len(cls_chg), val_chg, len(small_chg), rule_mis))
    for c in cls_chg:
        print("CLASS a=%s b=%s develop=%s branch=%s" % c)
    for c in small_chg:
        print("SMALL a=%s b=%s" % c)


if __name__ == "__main__":
    main()
