#!/usr/bin/env python3
"""
analyze_scan.py -- the 7-parameter wiener_lpdf on the #3329 scan against the
score-form references of probe_ref7.cpp.  Claims of stan-dev/math#3442: the
#3329 row of the accuracy table (value 4.8e-11, d/dsw -1e8 -> 5.3e-10).

Usage: analyze_scan.py REF_FILE DIR TREE...
  REF_FILE: probe_ref7.cpp REF_SCORE=1 with REF_LIST = the 16 scan points
  (ids l0 .. l15 = scan points s0 .. s15 of probe_shift.cpp).
  DIR holds probe_shift_<tree>.out.

Per scan point and tree: value (absolute error), d/dy (relative error), and
d/dsw (absolute error, as the true d/dsw is 1e-4 to 1e-2).
"""
import math
import os
import sys


def load_refs(path):
    refs = {}
    for line in open(path):
        f = line.rstrip("\n").split("\t")
        if f[0] != "REF" or len(f) < 19 or f[10] == "FAIL":
            continue
        refs["s" + f[1][1:]] = [float(x) for x in f[10:19]]
    return refs


def main():
    refs = load_refs(sys.argv[1])
    d = sys.argv[2]
    trees = sys.argv[3:]
    rows = {}
    for t in trees:
        for line in open(os.path.join(d, "probe_shift_%s.out" % t)):
            f = line.rstrip("\n").split("\t")
            if f[0] != "S":
                continue
            rows[(t, f[1])] = None if f[11] == "THROW" else \
                [float(x) for x in f[11:20]]
    print("%-4s %-9s %-11s %-9s %-9s %-11s %-9s" % (
        "id", "tree", "value", "|dv|", "d/dy rel", "d/dsw", "|d sw|"))
    worst = {t: [0.0, 0.0, 0.0, 0] for t in trees}
    for k in range(16):
        pid = "s%d" % k
        r = refs.get(pid)
        if r is None:
            continue
        for t in trees:
            z = rows.get((t, pid))
            w = worst[t]
            if z is None or not math.isfinite(z[0]):
                w[3] += 1
                print("%-4s %-9s %s" % (pid, t, "THROW" if z is None else z[0]))
                continue
            dv = abs(z[0] - r[0])
            dy = abs(z[1] - r[1]) / max(abs(r[1]), 1.0)
            dsw = abs(z[7] - r[7])
            w[0], w[1], w[2] = max(w[0], dv), max(w[1], dy), max(w[2], dsw)
            print("%-4s %-9s %-11.5g %-9.1e %-9.1e %-11.4g %-9.1e" % (
                pid, t, z[0], dv, dy, z[7], dsw))
    print("worst over the finite points:")
    for t in trees:
        w = worst[t]
        print("  %-9s value %.1e  d/dy %.1e  d/dsw %.1e  not finite %d" % (
            t, w[0], w[1], w[2], w[3]))


if __name__ == "__main__":
    main()
