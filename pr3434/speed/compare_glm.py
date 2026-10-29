#!/usr/bin/env python3
"""Compare two glm_compare outputs, line by line.

Usage: compare_glm.py out_a.txt out_b.txt

For each label kind (the part after the case tag) prints the number of
values, the number of bitwise-equal values, the largest relative
difference (|a - b| / max(|a|, |b|)) and the case where it occurs.
Non-finite values must match exactly; a mismatch is listed.
"""
import collections
import math
import sys


def read(path):
    out = collections.OrderedDict()
    for line in open(path):
        p = line.split()
        out[p[0]] = [float.fromhex(v) for v in p[1:]]
    return out


def kind(label):
    for k in ("_pt0_", "_pt1_"):
        if k in label:
            i = label.index(k)
            return label[i + 1:]
    return label


def main():
    a = read(sys.argv[1])
    b = read(sys.argv[2])
    if list(a) != list(b):
        print("LABEL MISMATCH")
    stats = collections.OrderedDict()
    bad = []
    for lab in a:
        if lab not in b:
            continue
        s = stats.setdefault(kind(lab), {"n": 0, "eq": 0, "max": 0.0,
                                         "at": ""})
        for i, (x, y) in enumerate(zip(a[lab], b[lab])):
            s["n"] += 1
            if x == y or (math.isnan(x) and math.isnan(y)):
                s["eq"] += 1
                continue
            if not (math.isfinite(x) and math.isfinite(y)):
                bad.append((lab, i, x, y))
                continue
            r = abs(x - y) / max(abs(x), abs(y))
            if r > s["max"]:
                s["max"] = r
                s["at"] = "%s[%d] %r %r" % (lab, i, x, y)
    for k, s in stats.items():
        print("%-22s %5d values  %5d bitwise  max rel %.2e  %s"
              % (k, s["n"], s["eq"], s["max"], s["at"]))
    print("non-finite mismatches: %d" % len(bad))
    for x in bad[:20]:
        print("  ", x)


if __name__ == "__main__":
    main()
