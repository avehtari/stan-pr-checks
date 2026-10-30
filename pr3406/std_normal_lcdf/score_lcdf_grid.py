#!/usr/bin/env python3
"""Claim N1 (comment of 2026-09-18): score probe_lcdf_grid outputs against mpmath.

Value log Phi(y) by three routes: A log(erfc(-y/sqrt2)/2), B log(ncdf(y)),
C quadrature of the small tail, phi(|y|) * int_0^inf exp(-|y| u - u^2/2) du
(peak 1 at u = 0), then log or log1p.  mpmath ncdf calls erfc for y < 0, so
B is not independent of A; C is.  Gradient exp(lpdf - lcdf) for each route.
Route A is the reference; the script prints the largest relative difference
of B and C from A per range (value and gradient).

Usage: score_lcdf_grid.py out_A.txt [out_B.txt ...]
Error units: rel = |got - ref| / |ref|; ulp = rel / 2^-53; ulp52 = rel / 2^-52.
"""
import math
import multiprocessing
import sys

import mpmath as mp

DPS = 50
DPS_QUAD = 40
EPS53 = mp.mpf(2) ** -53


def small_tail_quad(y):
    """Phi(-|y|) by quadrature of the density, scaled to a peak of 1."""
    a = abs(y)
    step = 1 / (a + 1)
    m = mp.quad(lambda u: mp.exp(-a * u - u * u / 2),
                [0, step, 4 * step, 16 * step, mp.inf])
    return mp.npdf(a) * m


def lcdf_routes(y):
    """Return log Phi(y) and its gradient by routes A, B, C."""
    y = mp.mpf(y)
    # For y > 0, Phi is near 1: add the digits that log(Phi) cancels.
    extra = 0
    if y > 0:
        with mp.workdps(20):
            q = mp.erfc(y / mp.sqrt(2)) / 2
            extra = int(-mp.log10(q)) + 5
    with mp.workdps(DPS + extra):
        a = mp.log(mp.erfc(-y / mp.sqrt(2)) / 2)
        b = mp.log(mp.ncdf(y))
    with mp.workdps(DPS_QUAD):
        g = small_tail_quad(y)
        c = mp.log(g) if y < 0 else mp.log1p(-g)
    with mp.workdps(DPS + extra):
        lpdf = -y * y / 2 - mp.log(mp.sqrt(2 * mp.pi))
        vals = (a, b, c)
        grads = tuple(mp.exp(lpdf - v) for v in vals)
        return [+v for v in vals], [+d for d in grads]


def _routes_of(y):
    return y, lcdf_routes(y)


def rel(got, ref):
    if not math.isfinite(got):
        return math.inf
    return float(abs((mp.mpf(got) - ref) / ref))


def load(path):
    blocks = []
    with open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                _, name, _lo, _hi, _n = line.split()
                blocks.append((name, []))
            else:
                blocks[-1][1].append([float.fromhex(t) for t in line.split()])
    return blocks


def main(paths):
    data = [(p, load(p)) for p in paths]
    blocks0 = data[0][1]
    ys = sorted({row[0] for _, blocks in data for _, rows in blocks for row in rows})
    with multiprocessing.Pool() as pool:
        cache = dict(pool.map(_routes_of, ys, chunksize=64))
    agree = {}
    for name, rows in blocks0:
        worst_b = worst_c = 0.0
        for row in rows:
            v, d = cache[row[0]]
            for k in (0, 1):
                ref = (v, d)[k][0]
                worst_b = max(worst_b, float(abs(((v, d)[k][1] - ref) / ref)))
                worst_c = max(worst_c, float(abs(((v, d)[k][2] - ref) / ref)))
        agree[name] = (worst_b, worst_c)
    for path, blocks in data:
        print("\n== %s" % path)
        print("std_normal_lcdf: value, d/dy; normal_lcdf(y | 0, 1): d/dy (rel only)")
        print("%-16s %5s %21s %33s %21s %3s %8s %8s" % (
            "range", "n", "value ulp @ y", "grad rel / ulp / ulp52 @ y",
            "normal_lcdf rel @ y", "nf", "B-A", "C-A"))
        for (name, rows), (name0, _) in zip(blocks, blocks0):
            assert name == name0
            if name == "test_points":
                continue
            wv, yv, wg, yg, wn, yn, nf = 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0
            for row in rows:
                y = row[0]
                v, d = cache[y]
                ev = rel(row[1], v[0])
                eg = rel(row[3], d[0])
                en = rel(row[4], d[0])
                nf += sum(not math.isfinite(t) for t in (row[1], row[3], row[4]))
                nf += row[1] != row[2]  # double and var values differ
                if ev > wv:
                    wv, yv = ev, y
                if eg > wg:
                    wg, yg = eg, y
                if en > wn:
                    wn, yn = en, y
            print("%-16s %5d %8.3g @ %-10.6g %9.3g / %8.3g / %8.3g @ %-10.6g "
                  "%8.3g @ %-10.6g %3d %8.2g %8.2g" % (
                      name, len(rows), wv / float(EPS53), yv,
                      wg, wg / float(EPS53), wg / float(EPS53) / 2, yg,
                      wn, yn, nf, agree[name][0], agree[name][1]))
        for name, rows in blocks:
            if name != "test_points":
                continue
            print("test points of expect_branch_accuracy: y, grad rel error of "
                  "std_normal_lcdf, of normal_lcdf(y | 0, 1)")
            for row in rows:
                d = cache[row[0]][1]
                with mp.workdps(60):
                    db = float(abs(d[1] / d[0] - 1))
                    dc = float(abs(d[2] / d[0] - 1))
                print("  %-20r %10.3g %10.3g  (routes B-A %.1g, C-A %.1g)" % (
                    row[0], rel(row[3], d[0]), rel(row[4], d[0]), db, dc))


if __name__ == "__main__":
    main(sys.argv[1:])
