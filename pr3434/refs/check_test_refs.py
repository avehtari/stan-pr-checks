#!/usr/bin/env python3
"""Recompute the mpmath reference constants of the PR tests by two routes.

Claim S2: "All references are mpmath at 80 to 140 digits, and two
independent routes agree: loggamma for the values, closed-form digamma
differences and mp.diff for the gradients, and sums of the pmf for the
cdfs."

Usage: check_test_refs.py <tree with test/>
       (for example a `git archive <commit> test` of the PR head)

Files and routes:
  rev/prob/large_shapes_test.cpp      the {"TAG", {args}, value, {grads}}
      rows. Families of ref_sib.py: route A analytic at 90 digits (cdfs:
      finite sums of the pmf), route B value at 140 digits and mp.diff.
      ST, MST, MSTF, YS, LBETA: ref_sweep.py, route A analytic at 50
      digits, route B mp.diff at 110 digits.
  prim/prob/beta_binomial_test.cpp    values; ref_bb.py route A (80
      digits) and B (130 digits)
  rev/prob/beta_binomial_lpmf_test.cpp  values and alpha d/dalpha,
      beta d/dbeta; route B is mp.diff at 130 digits; the lc gradients
      s (d/dalpha + d/dbeta), s = exp(lc) / 2 in double
  prim/fun/digamma_diff_test.cpp      route A digamma difference, route B
      mp.diff of the loggamma difference
  opencl/device_functions/digamma_diff_test.cpp  {x, d, ref} rows, as above
  opencl/rev/*_test.cpp               every other literal with at least 15
      significant digits must equal a constant checked above
Per file and tag: the worst disagreement of the routes, and the worst
difference between the constant in the test and route A, both relative to
the larger magnitude, and absolute where both are below 1e-40.
"""

import glob
import os
import re
import sys

from mpmath import mp, mpf

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import ref_bb as RB  # noqa: E402
import ref_sib as RS  # noqa: E402
import ref_sweep as RW  # noqa: E402

NUM = r"[-+]?(?:0x[0-9a-fA-F.]+p[-+]?\d+|[0-9.]+(?:e[-+]?\d+)?)"
SWEEP_TAGS = {"ST", "MST", "MSTF", "YS", "LBETA"}
OPENCL_FILES = ["beta_binomial_lpmf", "neg_binomial_lpmf",
                "neg_binomial_2_lpmf", "neg_binomial_2_log_lpmf",
                "neg_binomial_2_log_glm_lpmf", "student_t_lpdf"]
OPENCL_NB = {"neg_binomial_lpmf": "NB", "neg_binomial_2_lpmf": "NB2",
             "neg_binomial_2_log_lpmf": "NB2L"}
checked = set()


def num(tok):
  tok = tok.strip()
  neg = tok.startswith("-")
  body = tok.lstrip("+-")
  v = float.fromhex(body) if body.startswith("0x") else float(body)
  return -v if neg else v


def rel(x, y):
  # absolute below 1e-40: an exact zero (NB2 d/dmu at n = mu, ST d/dsigma
  # at |y - mu| = sigma) is 0 in one route and about 1e-50 in the other
  x, y = mpf(x), mpf(y)
  scale = max(abs(x), abs(y))
  if scale < mpf("1e-40"):
    return abs(x - y)
  return abs(x - y) / scale


class Worst:
  def __init__(self):
    self.rows = {}

  def add(self, key, routes, stored):
    w = self.rows.setdefault(key, [0, mpf(0), mpf(0)])
    w[0] += 1
    w[1] = max(w[1], routes)
    w[2] = max(w[2], stored)

  def show(self, title):
    print(f"\n## {title}\n")
    print(f"{'tag':10s} {'rows':>4s} {'routes A-B':>11s} {'test - A':>9s}")
    for k, (n, r, s) in self.rows.items():
      print(f"{k:10s} {n:4d} {float(r):11.1e} {float(s):9.1e}")


def sweep_routes(tag, args):
  t = "MST" if tag == "MSTF" else tag
  mp.dps = RW.DPS
  va, ga = RW.LP[t](*args), RW.GRAD[t](*args)
  with mp.workdps(RW.DPS_DIFF):
    vb = RW.LP[t](*args)
  gb = RW.grad_route_b(t, args)
  mp.dps = RW.DPS
  return [va] + ga, [vb] + gb


def sib_routes(tag, args):
  va, ga = RS.route_a(tag, args)
  vb, gb = RS.route_b(tag, args)
  mp.dps = RS.DPS_A
  return [va] + list(ga), [vb] + list(gb)


def large_shapes(path):
  w = Worst()
  text = open(path).read()
  row = re.compile(r'\{"([A-Z0-9]+)",\s*\{([^}]*)\},\s*(' + NUM
                   + r'),\s*\{([^}]*)\}\}')
  for m in row.finditer(text):
    tag = m.group(1)
    args = [num(t) for t in m.group(2).split(",")]
    stored = [num(m.group(3))] + [num(t) for t in m.group(4).split(",")]
    a, b = (sweep_routes if tag in SWEEP_TAGS else sib_routes)(tag, args)
    mp.dps = 60
    r = max(rel(x, y) for x, y in zip(a, b))
    s = max(rel(x, y) for x, y in zip(stored, a))
    w.add(tag, r, s)
    checked.update(stored)
  w.show("rev/prob/large_shapes_test.cpp")


def bb_rows(path, title, with_grads):
  w = Worst()
  text = open(path).read()
  n_fields = 7 if with_grads else 5
  pat = (r"\{(\d+),\s*(\d+),\s*(" + NUM + r"),\s*(" + NUM + r")"
         + (r",\s*(" + NUM + r")") * (n_fields - 4) + r"\}")
  for m in re.finditer(pat, text):
    n, N = int(m.group(1)), int(m.group(2))
    a, b = num(m.group(3)), num(m.group(4))
    stored = [num(m.group(i)) for i in range(5, n_fields + 1)]
    mp.dps = RB.DPS
    va = RB.lp_bb(n, N, a, b)
    with mp.workdps(RB.DPS_CHECK):
      vb = RB.lp_bb(n, N, a, b)
    ra, rb = [va], [vb]
    if with_grads:
      ga = RB.grad_bb(n, N, a, b)
      gb = RB.grad_bb_diff(n, N, a, b)
      ra += [mpf(a) * ga[0], mpf(b) * ga[1]]
      rb += [mpf(a) * gb[0], mpf(b) * gb[1]]
    mp.dps = 60
    w.add("BB", max(rel(x, y) for x, y in zip(ra, rb)),
          max(rel(x, y) for x, y in zip(stored, ra)))
    checked.update(stored)
  if with_grads:
    m = re.search(r"lcs = \{([^}]*)\};.*?grads = \{([^}]*)\};", text, re.S)
    if m:
      import math
      for lc, g in zip(m.group(1).split(","), m.group(2).split(",")):
        s = math.exp(float(lc)) / 2
        mp.dps = RB.DPS
        ga = RB.grad_bb(57, 117, s, s)
        gb = RB.grad_bb_diff(57, 117, s, s)
        mp.dps = 60
        ra = mpf(s) * (ga[0] + ga[1])
        rb = mpf(s) * (gb[0] + gb[1])
        w.add("BB d/dlc", rel(ra, rb), rel(num(g), ra))
        checked.add(num(g))
  w.show(title)


def dd_ref(x, d):
  mx, md = mpf(x), mpf(d)
  with mp.workdps(120):
    a = mp.digamma(mx + md) - mp.digamma(mx)
    if x < 1e100 and x + d < 1e300:
      with mp.workdps(160):
        # step relative to x, so that x - h stays positive for tiny x
        b = mp.diff(lambda t: mp.loggamma(t + md) - mp.loggamma(t), mx,
                    h=mx * mpf(10) ** -50)
    else:
      b = a  # mp.diff is not reliable there; route A only
  return a, b


def digamma_diff_tests(path, title, pattern):
  w = Worst()
  text = open(path).read()
  for m in re.finditer(pattern, text):
    x, d, ref = num(m.group(1)), num(m.group(2)), num(m.group(3))
    a, b = dd_ref(x, d)
    mp.dps = 60
    w.add("DD", rel(a, b), rel(ref, a))
    checked.add(ref)
  w.show(title)


def opencl_nb(t):
  """{n, a, b, value, d_a, d_b} rows of the three NB OpenCL tests."""
  w = Worst()
  pat = re.compile(r"\{(\d+),\s*(" + NUM + r"),\s*(" + NUM + r"),\s*(" + NUM
                   + r"),\s*(" + NUM + r"),\s*(" + NUM + r")\}")
  for f, tag in OPENCL_NB.items():
    text = open(os.path.join(t, "opencl", "rev", f + "_test.cpp")).read()
    for m in pat.finditer(text):
      args = [int(m.group(1)), num(m.group(2)), num(m.group(3))]
      stored = [num(m.group(i)) for i in (4, 5, 6)]
      a, b = sib_routes(tag, args)
      mp.dps = 60
      w.add(tag, max(rel(x, y) for x, y in zip(a, b)),
            max(rel(x, y) for x, y in zip(stored, a)))
      checked.update(stored)
  w.show("opencl/rev/neg_binomial{,_2,_2_log}_lpmf_test.cpp rows")


def literals(paths):
  print("\n## OpenCL reference tests: literals with at least 15 significant"
        " digits\n")
  lit = re.compile(r"(?<![\w.])-?\d\.\d{14,}(?:e[-+]?\d+)?")
  for path in paths:
    vals = [float(t) for t in lit.findall(open(path).read())]
    missing = [v for v in vals if v not in checked and -v not in checked]
    print(f"{os.path.basename(path):45s} {len(vals):4d} literals, "
          f"{len(missing):3d} not among the checked constants"
          + (f": {missing[:4]}" if missing else ""))


def main():
  if len(sys.argv) != 2:
    print(__doc__)
    return 2
  t = os.path.join(sys.argv[1], "test", "unit", "math")
  large_shapes(os.path.join(t, "rev", "prob", "large_shapes_test.cpp"))
  bb_rows(os.path.join(t, "prim", "prob", "beta_binomial_test.cpp"),
          "prim/prob/beta_binomial_test.cpp", False)
  bb_rows(os.path.join(t, "rev", "prob", "beta_binomial_lpmf_test.cpp"),
          "rev/prob/beta_binomial_lpmf_test.cpp", True)
  bb_rows(os.path.join(t, "opencl", "rev", "beta_binomial_lpmf_test.cpp"),
          "opencl/rev/beta_binomial_lpmf_test.cpp", True)
  digamma_diff_tests(
      os.path.join(t, "prim", "fun", "digamma_diff_test.cpp"),
      "prim/fun/digamma_diff_test.cpp",
      r"digamma_diff\((" + NUM + r"),\s*(" + NUM + r")\),\s*(" + NUM + r")")
  digamma_diff_tests(
      os.path.join(t, "opencl", "device_functions", "digamma_diff_test.cpp"),
      "opencl/device_functions/digamma_diff_test.cpp",
      r"\{(" + NUM + r"),\s*(" + NUM + r"),\s*(" + NUM + r")\}")
  opencl_nb(t)
  literals([os.path.join(t, "opencl", "rev", f + "_test.cpp")
            for f in OPENCL_FILES])
  return 0


if __name__ == "__main__":
  sys.exit(main())
