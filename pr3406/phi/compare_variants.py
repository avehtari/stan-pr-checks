#!/usr/bin/env python3
"""Claim N3: check that the copies in phi_variants.hpp equal Phi of the trees.
Compares each column of probe_phi_variants output ("x dev pr d1 d2 e e2")
with "x Phi" files of probe_phi_tree, bitwise, and counts the differences.

Usage: compare_variants.py variants.txt NAME=tree_file.txt [NAME=...]
"""
import sys

COLS = ["dev", "pr", "d1", "d2", "e", "e2"]


def load(path):
    with open(path) as fh:
        return [line.split() for line in fh if line.strip()]


def main():
    var = load(sys.argv[1])
    print("tree           " + "".join("%9s" % c for c in COLS) + "   (points that differ)")
    for arg in sys.argv[2:]:
        name, path = arg.split("=", 1)
        tree = load(path)
        assert len(tree) == len(var)
        assert all(a[0] == b[0] for a, b in zip(var, tree)), "x grids differ"
        counts = [sum(float.fromhex(v[j + 1]) != float.fromhex(t[1])
                      for v, t in zip(var, tree)) for j in range(len(COLS))]
        print("%-15s" % name + "".join("%9d" % c for c in counts))
    print("points: %d" % len(var))


if __name__ == "__main__":
    main()
