#!/usr/bin/env python3
"""
compare_ref7.py -- probe_shift outputs (E and D rows) against the references
of probe_ref7.cpp.  Claims of stan-dev/math#3442: the extreme grid (finite
points; d/dsw error 2e3 with sw > 0, st0 = 0).

Usage: compare_ref7.py REF_FILE RESULTS_DIR [TREE ...]

For each point with a finite reference, print the value error and the
largest relative gradient error (relative to max(|ref|, 1)) per tree, and
the name of the worst gradient. Gradient order: y a t0 w v sv sw st0.
"""
import math
import os
import sys

NAMES = ["y", "a", "t0", "w", "v", "sv", "sw", "st0"]


def load_ref(path):
    refs = {}
    for line in open(path):
        f = line.rstrip("\n").split("\t")
        if len(f) < 19 or f[10] == "FAIL":
            continue
        refs[f[1]] = ([float(x) for x in f[2:10]], [float(x) for x in f[10:19]])
    return refs


def load_tree(path):
    rows = {}
    for line in open(path):
        f = line.rstrip("\n").split("\t")
        if f[0] not in ("E", "D"):
            continue
        prec = f[10]
        if f[0] == "D" and prec != "0.0001":
            continue
        rows[f[1]] = "THROW" if f[11] == "THROW" else [float(x) for x in f[11:20]]
    return rows


def err(r, ref):
    if r == "THROW":
        return "THROW"
    if not math.isfinite(r[0]):
        return "%-8s" % ("%g" % r[0])
    ve = abs(r[0] - ref[0])
    ge, gn = 0.0, "-"
    for k in range(1, 9):
        if not math.isfinite(r[k]):
            ge, gn = math.inf, NAMES[k - 1]
            break
        e = abs(r[k] - ref[k]) / max(abs(ref[k]), 1.0)
        if e > ge:
            ge, gn = e, NAMES[k - 1]
    return "%-8.1e %-8.1e %-3s" % (ve, ge, gn)


def main():
    refs = load_ref(sys.argv[1])
    d = sys.argv[2]
    trees = sys.argv[3:] or ["base3442", "head3442"]
    data = {t: load_tree(os.path.join(d, "probe_shift_%s.out" % t)) for t in trees}
    print("%-5s %-11s %s" % ("id", "ref value", "  ".join(
        "%-21s" % ("%s: |dv| gerr name" % t) for t in trees)))
    for pid, (p, ref) in refs.items():
        if not all(math.isfinite(x) for x in ref):
            continue
        cells = []
        for t in trees:
            r = data[t].get(pid)
            cells.append("%-21s" % ("-" if r is None else err(r, ref)))
        print("%-5s %-11.5g %s" % (pid, ref[0], "  ".join(cells)))


if __name__ == "__main__":
    main()
