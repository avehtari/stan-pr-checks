#!/usr/bin/env python3
"""Claim N2 (comment on 5f33a10d31): the before/after table.
Score two probe_dedup outputs (first "PR" = before, second "dedup" = after)
against mpmath.  The reference route is erfc; refs_normal/gen_refs.py checks
it against a second route.

ulp = relative error / 2^-53.  Per range: worst ulp of value log Phi(z) and
slope phi(z)/Phi(z) for PR, prototype ("dedup") and, where valid, the
central erfcx_small variant ("cx").  Also counts bitwise-equal results and
points where the prototype error is above 4 x the PR error and above 2 ulp.

Usage: score_dedup.py out_pr.txt out_dedup.txt
"""
import sys

import mpmath

mpmath.mp.dps = 50
EPS = mpmath.mpf(2) ** -53
CENTRAL_MAX = 0.46875 * 2 ** 0.5


def reference(z):
    if z < 0:
        q = mpmath.erfc(-z / mpmath.sqrt(2)) / 2
        return mpmath.log(q), mpmath.npdf(z) / q
    q = mpmath.erfc(z / mpmath.sqrt(2)) / 2
    return mpmath.log1p(-q), mpmath.npdf(z) / (1 - q)


def ulp(got, exact):
    if exact == 0:
        return 0.0 if got == 0 else float("inf")
    return float(abs((mpmath.mpf(got) - exact) / exact) / EPS)


def load(path):
    blocks = []
    with open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                _, name, lo, hi, _n = line.split()
                blocks.append((name, []))
            else:
                blocks[-1][1].append(line.split())
    return blocks


def main(p_pr, p_dd):
    pr, dd = load(p_pr), load(p_dd)
    hdr = "%-11s %15s %15s %15s %6s %6s" % (
        "range", "PR val/slope", "dedup val/slope", "cx val/slope", "same", "worse")
    print(hdr)
    for (name, rows_a), (_, rows_b) in zip(pr, dd):
        w = {"pr": [0, 0], "dd": [0, 0], "cx": [0, 0]}
        same = worse = 0
        has_cx = False
        for a, b in zip(rows_a, rows_b):
            z = float.fromhex(a[0])
            v, s = reference(mpmath.mpf(z))
            ua = (ulp(float.fromhex(a[1]), v), ulp(float.fromhex(a[2]), s))
            ub = (ulp(float.fromhex(b[1]), v), ulp(float.fromhex(b[2]), s))
            for k in range(2):
                w["pr"][k] = max(w["pr"][k], ua[k])
                w["dd"][k] = max(w["dd"][k], ub[k])
                if ub[k] > 4 * ua[k] and ub[k] > 2:
                    worse += 1
            same += a[1] == b[1] and a[2] == b[2]
            if abs(z) < CENTRAL_MAX and len(b) == 5:
                has_cx = True
                uc = (ulp(float.fromhex(b[3]), v), ulp(float.fromhex(b[4]), s))
                for k in range(2):
                    w["cx"][k] = max(w["cx"][k], uc[k])
        cx = "%6.2f / %6.2f" % tuple(w["cx"]) if has_cx else "-"
        print("%-11s %6.2f / %6.2f %6.2f / %6.2f %15s %6d %6d" % (
            name, w["pr"][0], w["pr"][1], w["dd"][0], w["dd"][1], cx, same, worse))


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
