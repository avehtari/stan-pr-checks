#!/usr/bin/env python3
"""Claim M5 of README_emg.md: compare the rows that the two generators print
with the constants of exp_mod_normal_lcdf_reference.hpp.

Reads every "{...}," row of 14 numbers from each file, in order.  The first
13 rows of the header are exp_mod_normal_lcdf_reference, the last 6 are
exp_mod_normal_lcdf_deep_tail_reference.  Prints the number of numbers that
differ and the largest relative difference.

Usage: compare_header.py <header> <rows of 13> <rows of 6>
"""
import re
import sys

ROW = re.compile(r"\{\s*([-+0-9.eE,\s]+?)\s*\}", re.S)


def rows(path):
    text = re.sub(r"//[^\n]*", "", open(path).read())
    out = []
    for m in ROW.finditer(text):
        nums = [float(x) for x in m.group(1).replace("\n", " ").split(",")
                if x.strip()]
        if len(nums) == 14:
            out.append(nums)
    return out


def compare(name, got, ref):
    n_diff, worst = 0, 0.0
    for r, g in zip(ref, got):
        for a, b in zip(r, g):
            if a != b:
                n_diff += 1
                worst = max(worst, abs(a - b) / max(abs(a), abs(b)))
    print("%s: %d rows in the header, %d generated; %d of %d numbers differ;"
          " largest relative difference %.1e"
          % (name, len(ref), len(got), n_diff, 14 * len(ref), worst))


def main(header, rows13, rows6):
    h = rows(header)
    compare("reference (13 rows)", rows(rows13), h[:13])
    compare("deep tail (6 rows)", rows(rows6), h[13:])


if __name__ == "__main__":
    main(*sys.argv[1:4])
