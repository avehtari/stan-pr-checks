#!/usr/bin/env python3
"""Median time per kernel and tree from a glm_bench log, and ratios.

Usage: summarize.py glm_bench_<job>.log

Lines read: "round<r> <tree> TIME <kernel> <ns>". The trees are taken in
their order in the log; the first is the baseline (develop). Prints, per
kernel, the median ns of each tree, the ratio of each other tree to the
baseline, the ratio of the last tree to the one before it, and the spread
(max/min - 1) of the rounds. Also checks that the VALUE lines (value and
two adjoints in hex) are equal in every round, and lists the trees whose
VALUE lines differ.
"""
import collections
import statistics
import sys

times = collections.defaultdict(list)
values = collections.defaultdict(set)
order = []
trees = []
for line in open(sys.argv[1]):
    p = line.split()
    if len(p) < 4 or not p[0].startswith("round"):
        continue
    if p[1] not in trees:
        trees.append(p[1])
    if p[2] == "TIME":
        if p[3] not in order:
            order.append(p[3])
        times[(p[3], p[1])].append(float(p[4]))
    elif p[2] == "VALUE":
        values[(p[3], p[1])].add(" ".join(p[4:]))

base = trees[0]
cols = ["%s ns" % t for t in trees]
cols += ["%s/%s" % (t, base) for t in trees[1:]]
cols += ["%s/%s" % (trees[-1], trees[-2])] if len(trees) > 2 else []
print("%-28s " % "kernel" + " ".join("%10s" % c for c in cols)
      + "  spread")
for k in order:
    med = {t: statistics.median(times[(k, t)]) if times[(k, t)] else None
           for t in trees}

    def r(a, b):
        return "%.2f" % (med[a] / med[b]) if med[a] and med[b] else "-"

    cells = ["%.1f" % med[t] if med[t] else "-" for t in trees]
    cells += [r(t, base) for t in trees[1:]]
    if len(trees) > 2:
        cells.append(r(trees[-1], trees[-2]))
    spread = []
    for t in trees:
        v = times[(k, t)]
        spread.append("%.1f%%" % (100 * (max(v) / min(v) - 1)) if v else "-")
    print("%-28s " % k + " ".join("%10s" % c for c in cells) + "  "
          + "/".join(spread))

bad = [k for k, v in values.items() if len(v) != 1]
print("VALUE lines constant over rounds: %s" % ("yes" if not bad else bad))
for k in sorted(set(x for x, _ in values)):
    vs = {t: next(iter(values[(k, t)])) for t in trees if values[(k, t)]}
    if len(set(vs.values())) == 1:
        print("  %-22s all trees equal" % k)
    else:
        print("  %-22s differs:" % k)
        for t, v in vs.items():
            print("    %-5s %s" % (t, v))
