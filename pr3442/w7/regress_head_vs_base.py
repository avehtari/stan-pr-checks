#!/usr/bin/env python3
# regress_head_vs_base.py -- per point and quantity, flag head error > 4 x base
# error and > 1e-8 over the rerun outputs in results/ (README, "Head against base").
# Run from pr3442/: python3 w7/regress_head_vs_base.py
import math, os, re, sys

D = "results"
N8 = ["y", "a", "t0", "w", "v", "sv", "sw", "st0"]
flags = []


def rel(x, r):
  if not math.isfinite(x):
    return math.inf
  return abs(x - r) / max(abs(r), 1.0)


def absd(x, r):
  return abs(x - r) if math.isfinite(x) else math.inf


def check(tag, pid, qname, eb, eh, extra=""):
  if math.isnan(eh) or math.isnan(eb):
    return
  if eh > 1e-8 and eh > 4 * eb:
    flags.append((tag, pid, qname, eb, eh, extra))


def rows(path):
  for line in open(path):
    yield line.rstrip("\n").split("\t")


# p4/p5 (mpmath refs; p5 lcdf/lccdf values + FD gradients from ref7_p5)
ref4 = {(f[0], f[1]): [float(x) for x in f[2:9]] for f in rows("data/ref_w4_p4p5.out")}
refp5 = {}
for f in rows("data/ref7_p5_cdf_fd.out"):
  if f[0] == "RFD":
    g = [float(x) for x in f[4:20]]
    refp5[("LCDF", f[1])] = [float(f[2])] + g[:6]
    refp5[("LCCDF", f[1])] = [float(f[3])] + g[8:14]
w4 = {}
for t in ("base3442", "head3442"):
  for f in rows(os.path.join(D, "probe_w4_%s.out" % t)):
    if f[9] == "THROW":
      continue
    w4[(t, f[0], f[1], f[2])] = [float(x) for x in f[9:16]]
for (t, fn, st, pid), z in w4.items():
  if t != "head3442":
    continue
  zb = w4.get(("base3442", fn, st, pid))
  if st == "p5" and fn in ("LCDF", "LCCDF"):
    r = refp5.get((fn, pid))
  else:
    r = ref4.get(("LPDF" if fn == "OLD" else fn, pid))
  if r is None or zb is None:
    continue
  tag = "%s %s" % (fn, st)
  check(tag, pid, "value", absd(zb[0], r[0]), absd(z[0], r[0]))
  for k in range(6):
    if math.isnan(r[k + 1]) or (fn == "OLD" and k == 5):
      continue
    check(tag, pid, N8[k], rel(zb[k + 1], r[k + 1]), rel(z[k + 1], r[k + 1]))

# D lpdf
refD = {f[1]: [float(x) for x in f[10:19]] for f in rows("data/ref7_D_lpdf.out")}
dd = {}
for t in ("base3442", "head3442"):
  for f in rows(os.path.join(D, "probe_dacc_%s.out" % t)):
    dd[(t, f[1], f[10])] = [float(x) for x in f[11:20]]
for (t, pid, prec), z in dd.items():
  if t != "head3442":
    continue
  zb, r = dd[("base3442", pid, prec)], refD[pid]
  check("LPDF7 D", pid, "value", absd(zb[0], r[0]), absd(z[0], r[0]))
  for k in range(8):
    check("LPDF7 D", pid, N8[k], rel(zb[k + 1], r[k + 1]), rel(z[k + 1], r[k + 1]))

# D lcdf / lccdf
refC = {}
for i in range(4):
  for f in rows("data/ref7_D_cdf_fd_r%d.out" % i):
    if f[0] == "RFD":
      g = [float(x) for x in f[4:20]]
      refC[("LCDF", f[1])] = [float(f[2])] + g[:8]
      refC[("LCCDF", f[1])] = [float(f[3])] + g[8:]
cc = {}
dims = {}
for t in ("base3442", "head3442"):
  for f in rows(os.path.join(D, "probe_cdfacc_%s.out" % t)):
    if f[11] == "THROW":
      continue
    p = [float(x) for x in f[2:10]]
    dims[f[1]] = (p[5] != 0) + (p[6] != 0) + (p[7] != 0)
    cc[(t, f[0], f[1], f[10])] = [float(x) for x in f[11:20]]
cnt = {}
for (t, fn, pid, prec), z in cc.items():
  r = refC.get((fn, pid))
  if r is None:
    continue
  if absd(z[0], r[0]) > 1e-2:
    cnt[(t, fn, dims[pid])] = cnt.get((t, fn, dims[pid]), 0) + 1
  if t != "head3442":
    continue
  zb = cc.get(("base3442", fn, pid, prec))
  tag = "%s7 D dim%d" % (fn, dims[pid])
  check(tag, pid, "value", absd(zb[0], r[0]), absd(z[0], r[0]))
  for k in range(8):
    if math.isnan(r[k + 1]):
      continue
    check(tag, pid, N8[k], rel(zb[k + 1], r[k + 1]), rel(z[k + 1], r[k + 1]))
print("D cdf points with value error > 1e-2 or -inf:", sorted(cnt.items()))

# swst0
refS = {}
for i in range(4):
  for f in rows("data/ref7_swst0_s%d.out" % i):
    if f[0] == "REF" and len(f) >= 24 and f[10] != "FAIL":
      refS[f[1]] = [float(x) for x in f[10:19]]
ss = {}
for t in ("base3442", "head3442"):
  for f in rows(os.path.join(D, "probe_swst0_%s.out" % t)):
    if f[11] != "THROW":
      ss[(t, f[1], f[10])] = [float(x) for x in f[11:20]]
for (t, pid, prec), z in ss.items():
  if t != "head3442" or pid not in refS:
    continue
  zb, r = ss[("base3442", pid, prec)], refS[pid]
  tag = "LPDF7 swst0 prec %s" % prec
  check(tag, pid, "value", absd(zb[0], r[0]), absd(z[0], r[0]))
  for k in range(8):
    check(tag, pid, N8[k], rel(zb[k + 1], r[k + 1]), rel(z[k + 1], r[k + 1]))

# extreme grid and D rows of probe_shift vs ref7_E_D
refE = {f[1]: [float(x) for x in f[10:19]] for f in rows("data/ref7_E_D.out") if f[0] == "REF"}
ee = {}
for t in ("base3442", "head3442"):
  for f in rows(os.path.join(D, "probe_shift_%s.out" % t)):
    if f[0] in ("E", "D") and len(f) > 19 and f[11] != "THROW":
      try:
        ee[(t, f[1])] = [float(x) for x in f[11:20]]
      except ValueError:
        pass
for (t, pid), z in ee.items():
  if t != "head3442" or pid not in refE:
    continue
  zb, r = ee.get(("base3442", pid)), refE[pid]
  if zb is None or not math.isfinite(r[0]):
    continue
  check("LPDF7 E", pid, "value", absd(zb[0], r[0]), absd(z[0], r[0]))
  for k in range(8):
    check("LPDF7 E", pid, N8[k], rel(zb[k + 1], r[k + 1]), rel(z[k + 1], r[k + 1]))

print("flags: %d" % len(flags))
for fl in sorted(flags):
  print("%-22s %-6s %-6s base %.2e head %.2e  ratio %.1f" % (
      fl[0], fl[1], fl[2], fl[3], fl[4], fl[4] / fl[3] if fl[3] > 0 else math.inf))
