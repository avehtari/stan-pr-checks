#!/usr/bin/env python3
"""Summarize the RESULT and CHECK lines of run_bench_roihu.sh.

Usage: python3 summarize.py <job log> [range] [develop commit] [head commit]
Prints, per argument type, the table of the PR text: develop / head ns per
element (speedup = develop / head), minimum over all rounds and passes; then
the min and max speedup per function, the spread between passes, and the
largest relative difference of the CHECK values between the trees.
"""
import collections
import sys

FUNCS = ["normal_lcdf", "normal_lccdf", "normal_cdf", "std_normal_lcdf",
         "std_normal_lccdf", "std_normal_cdf", "Phi", "lognormal_lcdf",
         "lognormal_lccdf", "lognormal_cdf", "exp_mod_normal_lpdf",
         "exp_mod_normal_lcdf", "exp_mod_normal_lccdf", "exp_mod_normal_cdf",
         "skew_normal_lpdf"]
KINDS = ["double", "var", "varmat", "fvar"]
SIZES = [100, 1000, 10000, 100000]


def main():
    log = sys.argv[1]
    rng = sys.argv[2] if len(sys.argv) > 2 else "bulk"
    dev = sys.argv[3] if len(sys.argv) > 3 else "dc66941969"
    head = sys.argv[4] if len(sys.argv) > 4 else "16c1ca18bf"
    # (tree, fn, range, n, kind) -> pass -> list of times
    times = collections.defaultdict(lambda: collections.defaultdict(list))
    checks = collections.defaultdict(dict)
    with open(log) as f:
        for line in f:
            p = line.split()
            if not p:
                continue
            if p[0] == "RESULT":
                tree, pas = p[1].rsplit("_p", 1)
                fn, r, n = p[2], p[3], int(p[4])
                for i in range(6, len(p), 2):
                    times[(tree, fn, r, n, p[i])][pas].append(float(p[i + 1]))
            elif p[0] == "CHECK":
                tree = p[1].rsplit("_p", 1)[0]
                fn, r, n = p[2], p[3], int(p[4])
                for i in range(5, len(p), 3):
                    checks[(fn, r, n, p[i])][tree] = float(p[i + 1])

    def best(tree, fn, n, kind):
        d = times.get((tree, fn, rng, n, kind))
        if not d:
            return None
        return min(min(v) for v in d.values())

    def spread(tree, fn, n, kind):
        d = times.get((tree, fn, rng, n, kind))
        mins = [min(v) for v in d.values()]
        return max(mins) / min(mins) - 1

    print(f"range {rng}: develop {dev} / head {head}, ns per element "
          "(speedup = develop / head), min over rounds and passes")
    extremes = collections.defaultdict(list)
    for kind in KINDS:
        print(f"\n### {kind}\n")
        print("| function | " + " | ".join(f"N = {n}" for n in SIZES) + " |")
        print("|---" * (len(SIZES) + 1) + "|")
        for fn in FUNCS:
            cells = []
            for n in SIZES:
                a, b = best(dev, fn, n, kind), best(head, fn, n, kind)
                if a is None or b is None:
                    cells.append("n/a")
                    continue
                cells.append(f"{a:.1f} / {b:.1f} ({a / b:.2f}x)")
                extremes[fn].append((a / b, kind, n))
            print(f"| {fn} | " + " | ".join(cells) + " |")

    print("\n### min and max speedup per function\n")
    print("| function | min | case | max | case |")
    print("|---|---|---|---|---|")
    for fn in FUNCS:
        e = sorted(extremes[fn])
        if not e:
            continue
        lo, hi = e[0], e[-1]
        print(f"| {fn} | {lo[0]:.2f}x | {lo[1]}, N = {lo[2]} | "
              f"{hi[0]:.2f}x | {hi[1]}, N = {hi[2]} |")

    sp = [(spread(t, fn, n, k), t, fn, n, k)
          for (t, fn, r, n, k) in times if r == rng]
    sp.sort()
    vals = [s[0] for s in sp]
    vals_sorted = sorted(vals)
    med = vals_sorted[len(vals_sorted) // 2]
    print(f"\nspread between passes (max/min of the per-pass minimum - 1), "
          f"{len(vals)} cells: median {100 * med:.1f} %, "
          f"max {100 * sp[-1][0]:.1f} % ({sp[-1][1]} {sp[-1][2]} "
          f"N = {sp[-1][3]} {sp[-1][4]})")
    over5 = [s for s in sp if s[0] > 0.05]
    print(f"cells above 5 %: {len(over5)}")
    for s in over5[-10:]:
        print(f"  {100 * s[0]:.1f} % {s[1]} {s[2]} N = {s[3]} {s[4]}")

    worst = collections.defaultdict(lambda: (0.0, None))
    for (fn, r, n, kind), d in checks.items():
        if r != rng or dev not in d or head not in d:
            continue
        a, b = d[dev], d[head]
        rel = abs(a - b) / max(abs(a), abs(b), 1e-300)
        if rel > worst[fn][0]:
            worst[fn] = (rel, (kind, n, a, b))
    print("\nlargest relative difference of the CHECK value, develop vs head:")
    for fn in FUNCS:
        rel, case = worst[fn]
        print(f"  {fn}: {rel:.2e} {case if case else ''}")


if __name__ == "__main__":
    main()
