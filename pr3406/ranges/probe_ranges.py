#!/usr/bin/env python3
"""Find the valid input ranges of the stan-dev/math#3406 table on two trees.

Claim: "The valid ranges of inputs (before over/underflow) across all
functions are expanded" (PR description, 12-row table).
Usage: probe_ranges.py --dev <probe binary> --head <probe binary>
"""
import argparse
import math
import struct
import subprocess
import sys

import mpmath as mp

from refs_normal import REF_A, REF_B, reference, rel_diff

DBL_MAX = sys.float_info.max
TRUE_MIN = 5e-324
FLOOR = 2.0 ** -1074          # absolute floor: one subnormal ulp
TOL = {"T": 1e-2, "A": 1e-8, "C": 1e-8}

# ── Rows of the PR table: label, functions, axis, PR cells ──
# PR cells: (lo, hi) per (dev value, head value, dev grad, head grad).
# "L" means the representable limit; ("<", x) and (">", x) one-sided.
ROWS = [
    ("normal_lcdf, std_normal_lcdf", ["normal_lcdf", "std_normal_lcdf"], "z",
     [(-5e30, 38.5), (-1e154, 38.5), ("L", 38.6), ("L", 38.6)]),
    ("normal_lccdf, std_normal_lccdf", ["normal_lccdf", "std_normal_lccdf"],
     "z", [(-38.5, 5e30), (-38.5, 1e154), (-38.6, "L"), (-38.6, "L")]),
    ("normal_cdf, std_normal_cdf", ["normal_cdf", "std_normal_cdf"], "z",
     [(-37.5, 8.24), (-38.5, 8.28), (-37.5, 38.6), (-38.5, 38.6)]),
    ("Phi", ["Phi"], "z",
     [(-37.5, 8.24), (-38.5, 8.28), (-38.6, 38.6), (-38.6, 38.6)]),
    ("lognormal_lcdf", ["lognormal_lcdf"], "logy",
     [(-38.5, 8.28), (("<", -500), 38.5), (-38.5, 37.6),
      (("<", -500), 37.6)]),
    ("lognormal_lccdf", ["lognormal_lccdf"], "logy",
     [(-8.28, 38.5), (-38.5, (">", 500)), (-38.6, 38.5),
      (-38.6, (">", 500))]),
    ("lognormal_cdf", ["lognormal_cdf"], "logy",
     [(-38.5, 8.28), (-38.5, 8.28), (-27.2, 37.6), (-38.5, 37.6)]),
    ("exp_mod_normal_lpdf", ["exp_mod_normal_lpdf"], "z",
     [(-37.5, "L"), (-1e154, "L"), (-37.5, "L"), ("L", "L")]),
    ("exp_mod_normal_lcdf", ["exp_mod_normal_lcdf"], "z",
     [(-7.2, 37.9), (-2e5, 745), (-8.2, 745), (-2e5, 745)]),
    ("exp_mod_normal_lccdf", ["exp_mod_normal_lccdf"], "z",
     [(-7.2, 745), (-38.4, "L"), (-8.2, 745), (-38.6, "L")]),
    ("exp_mod_normal_cdf", ["exp_mod_normal_cdf"], "z",
     [(-7.2, 37.9), (-38.4, 37.9), (-8.2, 745), (-38.4, 745)]),
    ("skew_normal_lpdf", ["skew_normal_lpdf"], "z",
     [(-38.5, 1e154), (-1e154, 1e154), (-38.5, "L"), (-5e307, "L")]),
]


def saturation(name, part):
    """Values that mean over/underflow to a limit of the function."""
    if part == 1:
        return {0.0}
    if name.endswith("_cdf") or name == "Phi":
        return {0.0, 1.0}
    if name.endswith("lcdf") or name.endswith("lccdf"):
        return {0.0}
    return set()


# ── Probe binaries ──
class Probe:
    def __init__(self, path):
        self.p = subprocess.Popen([path], stdin=subprocess.PIPE,
                                  stdout=subprocess.PIPE, text=True,
                                  bufsize=1)
        self.cache = {}

    def __call__(self, name, x):
        key = (name, x)
        if key not in self.cache:
            self.p.stdin.write(f"{name} {x.hex()}\n")
            self.p.stdin.flush()
            out = self.p.stdout.readline().split()
            self.cache[key] = tuple(parse(s) for s in out)
        return self.cache[key]

    def batch(self, name, xs):
        todo = [x for x in xs if (name, x) not in self.cache]
        if todo:
            inp = "".join(f"{name} {x.hex()}\n" for x in todo)
            # a new process per batch avoids a pipe deadlock
            out = subprocess.run([self.p.args[0]], input=inp, text=True,
                                 capture_output=True, check=True).stdout
            for x, line in zip(todo, out.splitlines()):
                self.cache[(name, x)] = tuple(parse(s) for s in line.split())


def parse(s):
    try:
        return float.fromhex(s)
    except ValueError:
        return float(s)


REF_CACHE = {}


def ref(name, x):
    key = (name, x)
    if key not in REF_CACHE:
        REF_CACHE[key] = reference(name, x)
    return REF_CACHE[key]


def check(name, part, v, x):
    """Return (ok per criterion, relative error, reason)."""
    r = ref(name, x)
    if r is None:
        return {"T": False, "A": False, "C": False}, None, "domain"
    rv = r[part]
    if math.isnan(v):
        return {"T": False, "A": False, "C": False}, None, "nan"
    if abs(rv) > DBL_MAX:
        return {"T": False, "A": False, "C": False}, None, "ref-overflow"
    if math.isinf(v):
        return {"T": False, "A": False, "C": False}, None, "inf"
    err = abs(mp.mpf(v) - rv)
    rel = float(err / abs(rv)) if rv != 0 else (0.0 if v == 0 else math.inf)
    sat = v in saturation(name, part) and rv != v
    ok = {}
    for c, tol in TOL.items():
        acc = err <= tol * abs(rv) + FLOOR
        ok[c] = bool(acc and (c == "C" or not sat))
    reason = "sat" if sat else ("inaccurate" if not ok["A"] else "ok")
    return ok, rel, reason


# ── Axes ──
def scan_points(axis, side):
    """Scan points from the center outward; side is +1 or -1."""
    if axis == "z":
        mags = [10 ** (k / 40) for k in range(-80, 121)]
        mags += [10 ** (k / 4) for k in range(13, 1233)]
        mags = [m for m in mags if m < DBL_MAX] + [DBL_MAX]
        return [0.0] + [side * m for m in mags]
    zs = [k * 0.05 for k in range(0, 801)] + list(range(41, 746))
    xs = []
    for z in zs:
        x = float(mp.exp(side * mp.mpf(z)))
        if 0 < x < DBL_MAX and x not in xs:
            xs.append(x)
    xs.append(TRUE_MIN if side < 0 else DBL_MAX)
    return xs


def zcoord(axis, x):
    if axis == "z":
        return x
    return float(mp.log(x)) if x > 0 else -math.inf


def d2i(x):
    return struct.unpack("<q", struct.pack("<d", abs(x)))[0]


def i2d(i, sign):
    return math.copysign(struct.unpack("<d", struct.pack("<q", i))[0], sign)


def bisect(probe, name, part, crit, good, bad):
    sign = math.copysign(1.0, bad) if bad != 0 else 1.0
    ig, ib = d2i(good), d2i(bad)
    while abs(ib - ig) > 1:
        im = (ig + ib) // 2
        xm = i2d(im, sign)
        g, b = abs(i2d(ig, sign)), abs(i2d(ib, sign))
        if g > 0 and (b - g) <= 1e-7 * g:
            break
        ok, _, _ = check(name, part, probe(name, xm)[part], xm)
        if ok[crit]:
            ig = im
        else:
            ib = im
    return i2d(ig, sign), i2d(ib, sign)


def find_end(probe, name, part, crit, axis, side):
    pts = scan_points(axis, side)
    probe.batch(name, pts)
    oks = [check(name, part, probe(name, x)[part], x)[0][crit] for x in pts]
    if not oks[0]:
        return {"end": None, "note": "invalid at the center"}
    try:
        i = oks.index(False)
    except ValueError:
        return {"end": pts[-1], "limit": True, "note": ""}
    good, bad = bisect(probe, name, part, crit, pts[i - 1], pts[i])
    vb = probe(name, bad)[part]
    ok, rel, reason = check(name, part, vb, bad)
    r = ref(name, bad)
    note = f"first bad {zcoord(axis, bad):.6g}: {vb:.4g} ({reason}"
    if r is not None and abs(r[part]) <= DBL_MAX:
        note += f", ref {mp.nstr(r[part], 4)}"
    note += ")"
    last = max(j for j, o in enumerate(oks) if o)
    if last > i:
        note += (f"; valid again at {sum(oks[i:])} scan points, outermost "
                 f"{zcoord(axis, pts[last]):.3g}")
    lim = reason == "ref-overflow"
    return {"end": good, "limit": lim, "note": note}


def interior_worst(probe, name, part, axis, lo, hi):
    """Worst relative error and A-failures at scan points inside [lo, hi]."""
    worst, nbad, at = 0.0, 0, None
    for side in (-1, 1):
        for x in scan_points(axis, side):
            if lo is not None and hi is not None and lo <= x <= hi:
                ok, rel, _ = check(name, part, probe(name, x)[part], x)
                if rel is not None and rel > worst:
                    worst, at = rel, x
                nbad += not ok["A"]
    return worst, nbad, at


def fmt_end(e, axis, side):
    if e["end"] is None:
        return "none"
    z = zcoord(axis, e["end"])
    s = f"{z:.3g}"
    if e.get("limit"):
        s += "L"
    return s


def fmt_range(lo, hi, axis):
    return f"[{fmt_end(lo, axis, -1)}, {fmt_end(hi, axis, 1)}]"


def matches(pr, e, axis):
    if e["end"] is None:
        return False
    z = zcoord(axis, e["end"])
    if pr == "L":
        return bool(e.get("limit")) or abs(z) > 1e307
    if isinstance(pr, tuple):
        return z <= pr[1] if pr[0] == "<" else z >= pr[1]
    return abs(z - pr) <= 0.02 * abs(pr) + 0.05


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dev", required=True)
    ap.add_argument("--head", required=True)
    ap.add_argument("--crit", default="T,A,C")
    a = ap.parse_args()
    probes = {"dev": Probe(a.dev), "head": Probe(a.head)}
    results = {}
    print("# criteria: T finite, not saturated, rel err <= 1e-2;"
          " A as T with 1e-8; C rel err <= 1e-8, saturation allowed")
    print("# abs floor 2^-1074 in every criterion; 'L' = representable limit")
    for label, names, axis, pr in ROWS:
        for name in names:
            for tree, probe in probes.items():
                for part in (0, 1):
                    for crit in a.crit.split(","):
                        lo = find_end(probe, name, part, crit, axis, -1)
                        hi = find_end(probe, name, part, crit, axis, 1)
                        results[(name, tree, part, crit)] = (lo, hi)
                        pn = ["value", "grad"][part]
                        print(f"{name:22s} {tree:4s} {pn:5s} {crit} "
                              f"{fmt_range(lo, hi, axis)}")
                        for e, s in ((lo, "lo"), (hi, "hi")):
                            if e.get("note"):
                                print(f"    {s}: {e['note']}")
                    lo, hi = results[(name, tree, part, "T")]
                    if lo["end"] is not None and hi["end"] is not None:
                        w, nb, at = interior_worst(probe, name, part, axis,
                                                   lo["end"], hi["end"])
                        print(f"    inside T range: worst rel err {w:.2e} at "
                              f"z={zcoord(axis, at) if at else 0:.4g}, "
                              f"{nb} scan points fail A")
            sys.stdout.flush()
    # ── Tables in the PR layout ──
    for crit in a.crit.split(","):
        print(f"\n## Table, criterion {crit}  (* = differs from the PR cell)")
        print("| function | dev value | head value | dev grad | head grad |")
        print("|---|---|---|---|---|")
        for label, names, axis, pr in ROWS:
            cells = []
            for k, (tree, part) in enumerate([("dev", 0), ("head", 0),
                                              ("dev", 1), ("head", 1)]):
                txt = []
                for name in names:
                    lo, hi = results[(name, tree, part, crit)]
                    ok = matches(pr[k][0], lo, axis) and \
                        matches(pr[k][1], hi, axis)
                    txt.append(fmt_range(lo, hi, axis) + ("" if ok else "*"))
                cells.append(" / ".join(dict.fromkeys(txt)))
            print(f"| {label} | " + " | ".join(cells) + " |")
    candidates(probes, results, a.crit.split(","))
    bot_points(probes)


def round_inward(z, digits=2):
    """Round z toward 0 to `digits` significant digits."""
    if z == 0 or math.isinf(z):
        return z
    e = math.floor(math.log10(abs(z))) - digits + 1
    q = math.trunc(z / 10 ** e)
    r = q * 10 ** e
    return float(f"{r:.{digits}g}") if abs(r) <= DBL_MAX else z


def show_point(probes, name, axis, x):
    rows = []
    for tree in ("dev", "head"):
        v, g = probes[tree](name, x)
        okv = check(name, 0, v, x)
        okg = check(name, 1, g, x)
        rows.append(f"    {tree:4s} value {v!r:24} ({okv[2]}, rel "
                    f"{okv[1] if okv[1] is None else f'{okv[1]:.1e}'})  "
                    f"grad {g!r:24} ({okg[2]}, rel "
                    f"{okg[1] if okg[1] is None else f'{okg[1]:.1e}'})")
    ra = REF_A.value_grad(name, x)
    rb = REF_B.value_grad(name, x)
    if ra is not None:
        rows.append(f"    ref A value {mp.nstr(ra[0], 17)}  grad "
                    f"{mp.nstr(ra[1], 17)}")
        rows.append(f"    ref B value {mp.nstr(rb[0], 17)}  grad "
                    f"{mp.nstr(rb[1], 17)}")
        rows.append(f"    route agreement (rel): value "
                    f"{mp.nstr(rel_diff(ra[0], rb[0]), 2)}, grad "
                    f"{mp.nstr(rel_diff(ra[1], rb[1]), 2)}")
    return rows


def candidates(probes, results, crits):
    """Points just inside a head end where the head range is wider."""
    print("\n## Test candidates: head valid, develop not (T, then A)")
    done = set()
    for label, names, axis, pr in ROWS:
        for name in names:
            for part in (0, 1):
                for side, idx in ((-1, 0), (1, 1)):
                    for crit in [c for c in ("T", "A") if c in crits]:
                        eh = results[(name, "head", part, crit)][idx]
                        ed = results[(name, "dev", part, crit)][idx]
                        if eh["end"] is None:
                            continue
                        zh = zcoord(axis, eh["end"])
                        zd = (zcoord(axis, ed["end"]) if ed["end"] is not None
                              else 0.0)
                        if abs(zh) <= abs(zd) * (1 + 1e-6):
                            continue
                        x = None
                        for digits in (2, 3, 4, 6):
                            zc = round_inward(zh, digits)
                            xc = zc if axis == "z" else float(mp.exp(zc))
                            okh = check(name, part,
                                        probes["head"](name, xc)[part], xc)
                            okd = check(name, part,
                                        probes["dev"](name, xc)[part], xc)
                            if okh[0][crit] and not okd[0][crit]:
                                x = xc
                                break
                        pn = ["value", "grad"][part]
                        if x is None:
                            print(f"{name} {pn} side {side:+d} {crit}: no "
                                  f"rounded point near {zh:.3g} separates")
                            continue
                        key = (name, x)
                        tag = (f"{name} {pn} {'lower' if side < 0 else 'upper'}"
                               f" end ({crit}): dev {zd:.3g}, head {zh:.3g}")
                        if key in done:
                            print(tag + f"; same point y = {x!r}")
                            continue
                        done.add(key)
                        print(tag)
                        print(f"  y = {x!r} ({x.hex()}), z = "
                              f"{zcoord(axis, x):.6g}")
                        for line in show_point(probes, name, axis, x):
                            print(line)


def bot_points(probes):
    """Points that the bot names in the review of the range table."""
    print("\n## Points from the bot review")
    pts = [("std_normal_lcdf", -1.5e154), ("std_normal_lcdf", -1e154),
           ("std_normal_lcdf", -1e150), ("std_normal_lcdf", -1.5e308),
           ("normal_lcdf", -1.5e154), ("normal_lcdf", -1.5e308),
           ("std_normal_lccdf", 1.5e154), ("std_normal_lccdf", 1.5e308),
           ("exp_mod_normal_lcdf", -2e5), ("exp_mod_normal_lcdf", -1e154),
           ("exp_mod_normal_lcdf", 745.0), ("exp_mod_normal_lcdf", 746.0)]
    for name, x in pts:
        print(f"{name} y = {x!r}")
        for line in show_point(probes, name, "z", x):
            print(line)


if __name__ == "__main__":
    main()
