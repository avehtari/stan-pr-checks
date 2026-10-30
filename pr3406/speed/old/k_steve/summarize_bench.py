"""Summarize s3406_bench RESULT lines: minimum ns per element per
(tree, flags, function, range, type), and speed-ups (old time / new time,
bigger is better) for each of SteveBronder's commits and for the PR head
against develop.

Usage: python3 summarize_bench.py s3406_bench_<job>.log
"""
import re
import sys
from collections import defaultdict

DEV, PRE, MID, HEAD = "e4745c87ff", "e8bb9a0b38", "e1fea0651c", "cf177ef96b"


def main():
    best = defaultdict(lambda: float("inf"))
    spread = defaultdict(list)
    funcs = []
    for line in open(sys.argv[1]):
        m = re.match(r"RESULT (\S+)_(default|native) (\S+) (bulk|wide) \d+ \d+ "
                     r"double (\S+) var (\S+)", line)
        if not m:
            continue
        tree, flags, fn, rng, td, tv = m.groups()
        if fn not in funcs:
            funcs.append(fn)
        for typ, t in (("d", float(td)), ("v", float(tv))):
            key = (tree, flags, fn, rng, typ)
            best[key] = min(best[key], t)
            spread[key].append(t)
    worst_spread = max((max(v) / min(v) - 1) for v in spread.values())
    print("max spread (max/min - 1) over passes and rounds: %.1f %%" % (100 * worst_spread))
    for flags in ("default", "native"):
        print("\n## flags: %s; speed-up = old / new (d = double, v = var)" % flags)
        print("| function | range | head ns d / v | e1fea / pre d / v | cf177 / e1fea d / v "
              "| head / pre d / v | head / dev d / v |")
        print("|---|---|---|---|---|---|---|")
        for fn in funcs:
            for rng in ("bulk", "wide"):
                def t(tree, typ):
                    return best[(tree, flags, fn, rng, typ)]

                def su(old, new):
                    return "%.2f / %.2f" % (t(old, "d") / t(new, "d"), t(old, "v") / t(new, "v"))

                print("| %s | %s | %.1f / %.1f | %s | %s | %s | %s |" % (
                    fn, rng, t(HEAD, "d"), t(HEAD, "v"), su(PRE, MID), su(MID, HEAD),
                    su(PRE, HEAD), su(DEV, HEAD)))


if __name__ == "__main__":
    main()
