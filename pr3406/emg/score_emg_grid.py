#!/usr/bin/env python3
"""Claims M1a and M2 of README_emg.md: score probe_emg_phi / probe_emg_dense
outputs against mpmath (100 digits, reference() of score_emg_phi.py).

Per file and kind (lcdf, lccdf): worst ulp per component and regime (the
regimes of score_emg_detail.py), the non-finite results, the worst ulp over
the value and the four partials, the same without dsigma, the worst dsigma,
and the 5 worst points.  The ulp rule is the one of score_emg_detail.py: a
component whose reference is below DBL_MIN counts as correct when the result
is also below DBL_MIN.  Route check: d/dy from mpmath.diff against the
analytic derivative.

The reference depends only on the (z, v) points, so it is cached in a file.

Usage: score_emg_grid.py <ref cache> <probe output> [<probe output> ...]
"""
import json
import math
import multiprocessing
import os
import sys

import mpmath

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from score_emg_phi import EPS, load, reference  # noqa: E402

TINY = mpmath.mpf(2) ** -1022
NAMES = ["value", "dy", "dmu", "dsigma", "dlambda"]
REGIMES = ["rest", "small v (<= 1e-4)", "left tail (z <= -30)",
           "large v (>= 10)"]


def ulp(got, exact):
    if math.isnan(got) or math.isinf(got):
        return float("inf")
    if abs(exact) < TINY:
        return 0.0 if abs(got) < 2 ** -1022 else float("inf")
    return float(abs((mpmath.mpf(got) - exact) / exact) / EPS)


def regime(z, v):
    if v <= 1e-4:
        return REGIMES[1]
    if z <= -30:
        return REGIMES[2]
    if v >= 10:
        return REGIMES[3]
    return REGIMES[0]


def ref_job(key):
    kind, z, v = key
    vals, check = reference(kind == "lccdf", z, v)
    return key, [mpmath.nstr(x, 40) for x in vals], float(check)


def get_refs(cache_path, keys):
    cache = {}
    if os.path.exists(cache_path):
        with open(cache_path) as f:
            for k, vals, chk in json.load(f):
                cache[(k[0], float.fromhex(k[1]), float.fromhex(k[2]))] = (
                    vals, chk)
    todo = [k for k in keys if k not in cache]
    if todo:
        with multiprocessing.Pool() as pool:
            for key, vals, chk in pool.imap_unordered(ref_job, todo, 8):
                cache[key] = (vals, chk)
        with open(cache_path, "w") as f:
            json.dump([[[k[0], k[1].hex(), k[2].hex()], v[0], v[1]]
                       for k, v in cache.items()], f)
    return {k: ([mpmath.mpf(x) for x in cache[k][0]], cache[k][1])
            for k in keys}


def score(path, refs):
    tree, _ = load(path)
    print("\n== %s (%d lines)" % (path, len(tree)))
    worst_chk = max(refs[k][1] for k in tree)
    print("route check, max rel. difference of d/dy (diff vs analytic): %.1e"
          % worst_chk)
    for kind in ("lcdf", "lccdf"):
        groups, cells = {}, []
        nf_comp = nf_pts = 0
        for (k, z, v), got in tree.items():
            if k != kind:
                continue
            ref = refs[(k, z, v)][0]
            us = [ulp(g, r) for g, r in zip(got, ref)]
            n = sum(1 for g in got if not math.isfinite(g))
            nf_comp += n
            nf_pts += n > 0
            g = groups.setdefault(regime(z, v), [0.0] * 5)
            for i, u in enumerate(us):
                g[i] = max(g[i], u)
            cells.append((max(us), z, v, us))
        print("\n%s: %d points; non-finite components %d at %d points"
              % (kind, len(cells), nf_comp, nf_pts))
        print("%-22s" % "regime" + "".join("%9s" % n for n in NAMES))
        for name in REGIMES:
            g = groups.get(name)
            if g:
                print("%-22s" % name + "".join("%9.1e" % u for u in g))
        w_all = max(cells)
        w_nos = max((max(us[:3] + us[4:]), z, v) for _, z, v, us in cells)
        w_sig = max((us[3], z, v) for _, z, v, us in cells)
        print("worst over value and 4 partials: %.0f ulp at (%g, %g)"
              % w_all[:3])
        print("worst without dsigma:            %.0f ulp at (%g, %g)" % w_nos)
        print("worst dsigma:                    %.0f ulp at (%g, %g)" % w_sig)
        print("5 worst points (z, v): " + " ".join(NAMES))
        for w, z, v, us in sorted(cells, reverse=True)[:5]:
            print("  (%g, %g): " % (z, v) + " ".join("%.1e" % u for u in us))


def main(cache_path, paths):
    keys = set()
    for p in paths:
        keys.update(load(p)[0].keys())
    refs = get_refs(cache_path, sorted(keys))
    for p in paths:
        score(p, refs)


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2:])
