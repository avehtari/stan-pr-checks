#!/usr/bin/env python3
"""Claim N3: score Phi(x) columns against mpmath (50 digits), per range of x.
--table uses the rows of the posted tables.  A second route (quadrature of
the density) checks the reference at each worst point.

Input: lines "x v1 v2 ..." in hex (%a), as printed by probe_phi_*.cpp.
Error unit: ulp = |v - ref| / |ref| / 2^-53.  Prints the worst error and
its x for every column and range; --worst-in LO HI prints the worst point
of each column inside [LO, HI].
"""
import argparse
import mpmath as mp

mp.mp.dps = 50
EPS = mp.mpf(2) ** -53
RANGES = [(-38.5, -26.5), (-26.5, -5.0), (-5.0, -0.67), (-0.67, 0.67),
          (0.67, 5.0), (5.0, 8.5)]
# --table: the rows of the posted tables of 2026-10-19 and 2026-10-20
TABLE_RANGES = [(-38.5, -30.0), (-30.0, -20.0), (-20.0, -10.0), (-10.0, -5.0),
                (-5.0, -1.0), (-1.0, -0.6), (-0.6, 0.0), (0.0, 1.0), (1.0, 5.0),
                (5.0, 8.5)]


def ref_phi(x):
    # erfc form: no cancellation for x < 0, 50 digits for x > 0
    return mp.erfc(-mp.mpf(x) / mp.sqrt(2)) / 2


def ref_phi_quad(x):
    # Second route: quadrature of the density on the small side, scaled to
    # a peak of 1: Phi(-a) = phi(a) int_0^inf exp(-a u - u^2 / 2) du.
    a = abs(mp.mpf(x))
    step = 1 / (a + 1)
    m = mp.quad(lambda u: mp.exp(-a * u - u * u / 2),
                [0, step, 4 * step, 16 * step, mp.inf])
    q = mp.npdf(a) * m
    return q if x < 0 else 1 - q


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("file")
    ap.add_argument("--names", default="")
    ap.add_argument("--worst-in", nargs=2, type=float)
    ap.add_argument("--table", action="store_true",
                    help="use the ranges of the posted tables")
    ap.add_argument("--normal-only", action="store_true",
                    help="skip points where Phi(x) < 2^-1022 (as the posted tables)")
    a = ap.parse_args()
    rows = []
    with open(a.file) as f:
        for line in f:
            p = line.split()
            if not p or p[0].startswith("#"):
                continue
            rows.append([float.fromhex(s) for s in p])
    ncol = len(rows[0]) - 1
    names = a.names.split(",") if a.names else [f"c{j}" for j in range(ncol)]
    ranges = ([tuple(a.worst_in)] if a.worst_in
              else TABLE_RANGES if a.table else RANGES)
    worst = {(r, j): (0.0, None) for r in ranges for j in range(ncol)}
    for row in rows:
        x = row[0]
        rr = [r for r in ranges if r[0] <= x < r[1]]
        if not rr:
            continue
        ref = ref_phi(x)
        if ref == 0 or (a.normal_only and ref < mp.mpf(2) ** -1022):
            continue
        for j in range(ncol):
            u = float(abs(mp.mpf(row[j + 1]) - ref) / abs(ref) / EPS)
            if u > worst[(rr[0], j)][0]:
                worst[(rr[0], j)] = (u, x)
    print("range," + ",".join(names))
    for r in ranges:
        cells = []
        for j in range(ncol):
            u, x = worst[(r, j)]
            cells.append(f"{u:.3g}@{x!r}" if x is not None else "0")
        print(f"[{r[0]},{r[1]})," + ",".join(cells))
    # Second route at every worst point: relative difference of the routes.
    xs = {x for (u, x) in worst.values() if x is not None}
    diff = max((float(abs(ref_phi_quad(x) / ref_phi(x) - 1)) for x in xs
                if ref_phi(x) != 0), default=0.0)
    print(f"routes erfc and quadrature at the {len(xs)} worst points: "
          f"largest relative difference {diff:.2g}")


if __name__ == "__main__":
    main()
