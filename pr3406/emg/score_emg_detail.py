#!/usr/bin/env python3
"""Claim M1a of README_emg.md (original scorer of the failure table; the
same table is in score_emg_grid.py): per-component errors of
exp_mod_normal lcdf/lccdf on one tree.

Uses the reference of score_emg_phi.py.  A component whose reference is
below DBL_MIN in magnitude counts as correct when the result is also below
DBL_MIN (underflow, correctly rounded to 0 or a subnormal).

Prints, per kind, the worst ulp of each component over regime groups, the
number of non-finite results, and the 8 worst cells with all components.

Usage: score_emg_detail.py out_pr.txt
"""
import math
import sys

import mpmath

sys.path.insert(0, ".")
from score_emg_phi import EPS, load, reference  # noqa: E402

TINY = mpmath.mpf(2) ** -1022
NAMES = ["value", "dy", "dmu", "dsigma", "dlambda"]


def ulp(got, exact):
    if math.isnan(got) or math.isinf(got):
        return float("inf")
    if abs(exact) < TINY:
        return 0.0 if abs(got) < 2 ** -1022 else float("inf")
    return float(abs((mpmath.mpf(got) - exact) / exact) / EPS)


def regime(z, v):
    if v <= 1e-4:
        return "small v (<= 1e-4)"
    if z <= -30:
        return "left tail (z <= -30)"
    if v >= 10:
        return "large v (>= 10)"
    return "rest"


def main(path):
    tree, _ = load(path)
    for kind in ("lcdf", "lccdf"):
        groups = {}
        cells = []
        nonfinite = 0
        for (k, z, v), got in tree.items():
            if k != kind:
                continue
            ref, _ = reference(kind == "lccdf", z, v)
            us = [ulp(g, r) for g, r in zip(got, ref)]
            nonfinite += sum(1 for g in got if not math.isfinite(g))
            g = groups.setdefault(regime(z, v), [0.0] * 5)
            for i, u in enumerate(us):
                g[i] = max(g[i], u)
            cells.append((max(us), z, v, us))
        print("\n%s (#3406): worst ulp per component; non-finite results: %d"
              % (kind, nonfinite))
        print("%-22s" % "regime" + "".join("%9s" % n for n in NAMES))
        for name in ("rest", "small v (<= 1e-4)", "left tail (z <= -30)",
                     "large v (>= 10)"):
            g = groups.get(name)
            if g:
                print("%-22s" % name + "".join("%9.1e" % u for u in g))
        print("8 worst cells (z, v): " + " ".join(NAMES))
        for w, z, v, us in sorted(cells, reverse=True)[:8]:
            print("  (%g, %g): " % (z, v) + " ".join("%.1e" % u for u in us))


if __name__ == "__main__":
    main(sys.argv[1])
