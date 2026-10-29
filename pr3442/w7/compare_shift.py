#!/usr/bin/env python3
"""
compare_shift.py -- probe_shift outputs of two trees: the #3329 point and
scan, the D points, the extreme grid.  Claims of stan-dev/math#3442: the
#3329 value row (develop -inf) and the finite count of the extreme grid.

Usage: compare_shift.py RESULTS_DIR [BASE_TREE NEW_TREE]
  Reads probe_shift_<tree>.out; default trees base3442 head3442.

Sections (see probe_shift.cpp): P and S against the R rows of the new tree
(log(F5 / st0) and p5 / F5 at sw -> 0, by Boost Gauss-Kronrod on the
5-parameter density); E finite count, throws, and the new tree against the
base tree where both are finite.  The R rows leave out sw = 1e-8: at w = 1e-4
this changes the log density by about sw^2 / (24 w^2) = 4e-10, so the P rows
with sw > 0 show about 4e-10 even if exact.  For the D points and the E
references use analyze_tol.py and compare_ref7.py.
"""
import math
import os
import sys

NAMES = ["y", "a", "t0", "w", "v", "sv", "sw", "st0"]


def load(path):
    rows, refs = {}, {}
    for line in open(path):
        f = line.rstrip("\n").split("\t")
        if f[0] == "R":
            refs[f[1]] = (float(f[2]), float(f[3]), float(f[4]))
            continue
        key = (f[0], f[1], f[10])
        if f[11] == "THROW":
            rows[key] = "THROW"
        else:
            rows[key] = [float(x) for x in f[11:20]]
    return rows, refs


def rel(x, r, floor=1.0):
    if not (math.isfinite(x) and math.isfinite(r)):
        return math.inf
    return abs(x - r) / max(abs(r), floor)


def main():
    d = sys.argv[1]
    base, new = sys.argv[2:4] if len(sys.argv) > 3 else ("base3442", "head3442")
    trees = [base, new]
    data = {t: load(os.path.join(d, "probe_shift_%s.out" % t)) for t in trees}
    refs = data[new][1]

    print("== P and S: value vs log(F5/st0) or the 5-parameter value; "
          "d/dy vs p5/F5; d/dst0 vs -1/st0")
    print("%-14s %-9s %-12s %-12s %-12s %-12s" %
          ("id", "tree", "value", "|v err|", "d/dy rel", "d/dst0"))
    pts = [("P", "both", "point"), ("P", "st0_only", "point"),
           ("P", "sw_only", "five_param")] + \
          [("S", "s%d" % k, "s%d" % k) for k in range(16)]
    for sec, idx, rid in pts:
        rv, rdy = refs[rid][1], refs[rid][2]
        for t in trees:
            r = data[t][0].get((sec, idx, "1e-08"))
            if r is None or r == "THROW":
                print("%-14s %-9s %s" % (idx, t, r))
                continue
            dy = "-" if rid == "five_param" else "%.1e" % rel(r[1], rdy)
            dst0 = "-" if rid == "five_param" else "%.10g" % r[8]
            print("%-14s %-9s %-12.6f %-12.1e %-12s %-12s" %
                  (idx, t, r[0], abs(r[0] - rv), dy, dst0))

    print("\n== E: extreme grid")
    ids = sorted({k[1] for k in data[base][0] if k[0] == "E"},
                 key=lambda s: int(s[1:]))
    for t in trees:
        rows = [data[t][0].get(("E", i, "1e-06")) for i in ids]
        fin = sum(1 for r in rows if r not in (None, "THROW")
                  and all(math.isfinite(x) for x in r))
        thr = sum(1 for r in rows if r == "THROW")
        ninf = sum(1 for r in rows if r not in (None, "THROW")
                   and r[0] == -math.inf)
        print("%-9s finite %d of %d, value -inf %d, throws %d" %
              (t, fin, len(ids), ninf, thr))
    vmax, gmax, nb = 0.0, 0.0, 0
    for i in ids:
        a = data[new][0].get(("E", i, "1e-06"))
        b = data[base][0].get(("E", i, "1e-06"))
        if a in (None, "THROW") or b in (None, "THROW"):
            continue
        if not all(math.isfinite(x) for x in a + b):
            continue
        nb += 1
        vmax = max(vmax, abs(a[0] - b[0]))
        gmax = max(gmax, max(rel(a[k], b[k]) for k in range(1, 9)))
    print("%s vs %s where both finite (%d points): value %.1e, "
          "gradients %.1e" % (new, base, nb, vmax, gmax))


if __name__ == "__main__":
    main()
