#!/usr/bin/env python3
"""Tabulate probe_speed3 output from one or more job logs.

Usage: python3 table_speed.py <log> [<log> ...]

Every log holds rounds "######## speed <tree>, round <r> ########". For each
kernel (K lines: net ns; D and F lines: ns per observation) prints the
median over all rounds of every tree, the spread (max/min - 1) per tree, and
the ratios ours/steve and dev/steve (above 1: steve is faster).
"""

import os
import re
import statistics
import sys
from collections import defaultdict

TREES = ["dev", "ours", "steve"]


def parse(paths):
  vals = defaultdict(lambda: defaultdict(list))
  order = []
  for path in paths:
    tree = None
    for line in open(path):
      m = re.match(r"######## speed (\w+), round (\d+)", line)
      if m:
        tree = m.group(1)
        continue
      if tree is None:
        continue
      if line.startswith("K "):
        name = line[2:48].strip()
        net = float(line[48:].split()[1])
        key = "K " + name
      elif line.startswith("D ") or line.startswith("F "):
        m = re.match(r"([DF]) (.{30}) a=(\S+)\s+b=(\S+)\s+([\d.]+) ns/obs", line)
        if not m:
          continue
        key = "%s %s a=%s b=%s" % (m.group(1), m.group(2).strip(), m.group(3),
                                   m.group(4))
        net = float(m.group(5))
      else:
        continue
      if key not in order:
        order.append(key)
      vals[key][tree].append(net)
  return order, vals


def main():
  # TREES=dev,old,new python3 table_speed.py ... selects other tree names;
  # the ratios are second/third and first/third
  trees = os.environ.get("TREES", ",".join(TREES)).split(",")
  t0, t1, t2 = trees
  order, vals = parse(sys.argv[1:])
  print("%-52s %10s %10s %10s %6s %6s %7s %7s" %
        ("kernel", t0, t1, t2, "sp_" + t1[0], "sp_" + t2[0],
         t1[0] + "/" + t2[0], t0[0] + "/" + t2[0]))
  for key in order:
    med = {}
    spread = {}
    for t in trees:
      v = vals[key].get(t, [])
      med[t] = statistics.median(v) if v else float("nan")
      spread[t] = (max(v) / min(v) - 1) if v and min(v) > 0 else float("nan")
    print("%-52s %10.1f %10.1f %10.1f %5.0f%% %5.0f%% %7.2f %7.2f" %
          (key[:52], med[t0], med[t1], med[t2],
           100 * spread[t1], 100 * spread[t2],
           med[t1] / med[t2], med[t0] / med[t2]))


if __name__ == "__main__":
  main()
