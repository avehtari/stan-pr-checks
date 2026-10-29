#!/usr/bin/env python3
"""
table_bench.py -- the run-time table of a run_bench_tree.sh output.  Claims
of stan-dev/math#3442: the run-time ratios "this PR / develop".

Usage: table_bench.py BENCH_OUT [BASE_TREE NEW_TREE]
  BENCH_OUT: <probe>_bench.out, rows <tree> <round> BENCH <set> <mode> <ns>
  <calls> <checksum>.  Default trees base3442 head3442.

Per set and mode: the median ns per call of each tree over the rounds, the
spread of each tree (max / min - 1 over the rounds), and the ratio
new / base of the medians (< 1: the new tree is faster).
"""
import statistics
import sys


def main():
    base, new = sys.argv[2:4] if len(sys.argv) > 3 else ("base3442", "head3442")
    ns = {}
    order = []
    for line in open(sys.argv[1]):
        f = line.rstrip("\n").split("\t")
        if len(f) < 8 or f[2] != "BENCH":
            continue
        tree, key = f[0], (f[3], f[4])
        if key not in order:
            order.append(key)
        ns.setdefault((tree, key), []).append(float(f[5]))
    print("%-12s %-4s %-6s %12s %7s %12s %7s %7s" % (
        "set", "mode", "rounds", base + " ns", "spread", new + " ns", "spread",
        "ratio"))
    for key in order:
        cells = []
        for t in (base, new):
            x = ns.get((t, key), [])
            cells.append((statistics.median(x) if x else float("nan"),
                          max(x) / min(x) - 1 if x else float("nan"), len(x)))
        ratio = cells[1][0] / cells[0][0]
        print("%-12s %-4s %-6d %12.1f %6.1f%% %12.1f %6.1f%% %7.3f" % (
            key[0], key[1], min(cells[0][2], cells[1][2]), cells[0][0],
            100 * cells[0][1], cells[1][0], 100 * cells[1][1], ratio))


if __name__ == "__main__":
    main()
