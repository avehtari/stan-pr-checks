#!/usr/bin/env python3
"""Compare propto_probe outputs of develop and the branch.

Usage: compare_propto.py out_dev.txt out_branch.txt

A value matches if |branch - develop| <= 1e-12 or the relative difference is
<= 1e-13 (the develop forms are accurate at the probe points). Prints one
line per (function, propto) with the count of mismatches and the distinct
differences, so that a dropped constant shows as one repeated number.
"""
import collections
import sys


def read(path):
    vals = {}
    for line in open(path):
        fn, mask, propto, v = line.split()
        vals[(fn, int(mask), int(propto))] = float(v)
    return vals


def main():
    dev = read(sys.argv[1])
    br = read(sys.argv[2])
    if set(dev) != set(br):
        print("KEY MISMATCH: %d only in develop, %d only in branch"
              % (len(set(dev) - set(br)), len(set(br) - set(dev))))
    stats = collections.OrderedDict()
    for key in dev:
        if key not in br:
            continue
        fn, mask, propto = key
        d, b = dev[key], br[key]
        diff = b - d
        ok = abs(diff) <= 1e-12 or abs(diff) <= 1e-13 * abs(d)
        s = stats.setdefault((fn, propto), {"n": 0, "bad": 0, "diffs": [],
                                            "masks": []})
        s["n"] += 1
        if not ok:
            s["bad"] += 1
            s["diffs"].append(diff)
            s["masks"].append(mask)
    total_bad = 0
    for (fn, propto), s in stats.items():
        total_bad += s["bad"]
        if s["bad"] == 0:
            line = "ok"
        else:
            distinct = sorted(set("%.12g" % x for x in s["diffs"]))
            line = "MISMATCH masks %s diffs %s" % (
                ",".join(str(m) for m in s["masks"]), " ".join(distinct[:6]))
        print("%-28s propto=%d  %3d values  %3d bad  %s"
              % (fn, propto, s["n"], s["bad"], line))
    print("TOTAL mismatches: %d" % total_bad)


if __name__ == "__main__":
    main()
