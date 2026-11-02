#!/usr/bin/env python3
"""Score the output of probe_ocl.cpp: the OpenCL device value of beta(a, b)
against the CPU value (ulp distance) and both against mpmath (50 digits).

Groups: R (the reference rows of the unit tests) and the three grid regions
(0 small, 1 mixed, 2 large).  Points whose reference is below 2.2e-308 are
counted separately and do not enter the ulp and relative columns.

Usage: python3 score_ocl.py < probe_ocl_output.txt
"""
import struct
import sys

import mpmath as mp

DBL_MIN = 2.2250738585072014e-308
NAMES = {"R": "reference rows", "0": "small", "1": "mixed", "2": "large"}


def ordered(x):
    i = struct.unpack("<q", struct.pack("<d", x))[0]
    return i if i >= 0 else -(i & 0x7FFFFFFFFFFFFFFF)


def main():
    stats = {}
    for line in sys.stdin:
        tag, reg, a, b, cpu, dev = line.split()
        a, b = float.fromhex(a), float.fromhex(b)
        cpu, dev = float.fromhex(cpu), float.fromhex(dev)
        key = "R" if tag == "R" else reg
        s = stats.setdefault(key, dict(n=0, under=0, ulp=0, over4=0,
                                       rc=0.0, rd=0.0, at=None))
        s["n"] += 1
        with mp.workdps(50):
            r = mp.beta(mp.mpf(a), mp.mpf(b))
            if r < DBL_MIN:
                s["under"] += 1
                continue
            ec = float(abs(mp.mpf(cpu) - r) / r)
            ed = float(abs(mp.mpf(dev) - r) / r)
        u = abs(ordered(cpu) - ordered(dev))
        if u > s["ulp"]:
            s["ulp"], s["at"] = u, (a, b)
        s["over4"] += u > 4
        s["rc"] = max(s["rc"], ec)
        s["rd"] = max(s["rd"], ed)
    print("# group            points  underflow  max_ulp_cpu_device  at (a, b)"
          "                 points_above_4_ulp  cpu_worst_rel  device_worst_rel")
    for key in ["R", "0", "1", "2"]:
        if key not in stats:
            continue
        s = stats[key]
        at = "(%.4g, %.4g)" % s["at"] if s["at"] else "-"
        print("%-16s %7d  %9d  %18d  %-24s  %18d  %13.2e  %16.2e"
              % (NAMES[key], s["n"], s["under"], s["ulp"], at, s["over4"],
                 s["rc"], s["rd"]))


if __name__ == "__main__":
    main()
