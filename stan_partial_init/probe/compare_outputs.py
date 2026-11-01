#!/usr/bin/env python3
"""Compare probe_init outputs of develop and the PR.

Usage: compare_outputs.py <develop output> <pr output>
Prints one summary row per case and build, then the bitwise checks:
the parameters of the runs without partial inits, and the next RNG draw
after initialize for every case and seed.
"""
import re
import sys


def parse(path):
    rows, params = {}, {}
    for line in open(path):
        m = re.match(r"(\S+) seed\s+(\d+) params (.*)", line)
        if m:
            params[(m.group(1), int(m.group(2)))] = m.group(3).split()
            continue
        m = re.match(r"(\S+) seed\s+(\d+) status (\S+)(.*)", line)
        if m:
            fields = m.group(4).split()
            vals = dict(zip(fields[0::2], fields[1::2]))
            vals["status"] = m.group(3)
            rows[(m.group(1), int(m.group(2)))] = vals
    return rows, params


def summary(rows, label):
    cases = []
    for case, _ in rows:
        if case not in cases:
            cases.append(case)
    print(f"{'case':<12} {'build':<8} {'ok':>3} {'throw':>5} "
          f"{'abs_missing':>11} {'rel_missing':>11} {'rel_supplied':>12} {'max_abs_u':>9}")
    for case in cases:
        rs = [v for (c, _), v in rows.items() if c == case]
        ok = [r for r in rs if r["status"] == "ok"]
        n_throw = len(rs) - len(ok)

        def mx(key):
            return max((float(r[key]) for r in ok), default=float("nan"))

        print(f"{case:<12} {label:<8} {len(ok):>3} {n_throw:>5} "
              f"{mx('max_abs_missing'):>11.3g} {mx('max_rel_missing'):>11.3g} "
              f"{mx('max_rel_supplied'):>12.3g} {mx('max_abs_unconstrained'):>9.3g}")


def main():
    dev_rows, dev_params = parse(sys.argv[1])
    pr_rows, pr_params = parse(sys.argv[2])
    summary(dev_rows, "develop")
    summary(pr_rows, "PR")
    same = sum(dev_params[k] == pr_params.get(k) for k in dev_params)
    print(f"\nruns without partial inits: {same} of {len(dev_params)} "
          "parameter vectors bitwise equal")
    keys = [k for k in dev_rows if "rng_next" in dev_rows[k] and "rng_next" in pr_rows.get(k, {})]
    same_rng = sum(dev_rows[k]["rng_next"] == pr_rows[k]["rng_next"] for k in keys)
    print(f"next RNG draw after initialize: {same_rng} of {len(keys)} equal")


if __name__ == "__main__":
    main()
