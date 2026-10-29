#!/usr/bin/env python3
"""References for the test rows added after jachymb's comment on #3434.

Rows: the negative binomial family at small phi (or alpha), where
binomial_coefficient_log(n + phi - 1, n) rounded phi into n + phi - 1, and
the GLM at theta = 720 and 709, where exp(theta) overflowed in the
partials. Values and partials by the route-A formulas of ref_sib.py at 90
digits; every partial is checked against mp.diff of the value at 140
digits (route B).

Usage: make_small_phi_rows.py   (run in beta_binomial_work/sib)
Prints the rows for test/unit/math/rev/prob/large_shapes_test.cpp and for
the opencl_large_shapes_reference tests.
"""

import math

from mpmath import mp, mpf

import ref_sib as R

LOG1000 = math.log(1000.0)

ROWS = [
    # tag, args, why
    ("NB", (1000, 1e-8, 1e-11), "alpha = 1e-8, rounded in n + alpha - 1"),
    ("NB", (1, 1e-17, 1e-17), "alpha = 1e-17: -inf before"),
    ("NB2", (1000, 1000.0, 1e-8), "phi = 1e-8, rounded in n + phi - 1"),
    ("NB2", (1, 1.0, 1e-17), "phi = 1e-17: -inf before"),
    ("NB2L", (1000, LOG1000, 1e-8), "phi = 1e-8, rounded in n + phi - 1"),
    ("NB2L", (1, 0.0, 1e-17), "phi = 1e-17: -inf before"),
    ("GLM1", (3, 1.0, 0.0, 720.0, 1.5), "theta = 720: exp(theta) overflowed"),
    ("GLM1", (5, 1.0, 0.0, 709.0, 1.5), "theta = 709: mu (y + phi) overflowed"),
    ("GLM1", (1, 1.0, 0.0, 0.0, 1e-17), "phi = 1e-17: -inf before"),
    ("GLM1", (1000, 1.0, 0.0, LOG1000, 1e-8), "phi = 1e-8, rounded"),
]

# indices of the arguments that the partials are taken with respect to
WRT = {"NB": (1, 2), "NB2": (1, 2), "NB2L": (1, 2), "GLM1": (2, 3, 4)}


def route_a(tag, args):
  mp.dps = R.DPS_A
  if tag.startswith("GLM"):
    return R.lp_GLM(*args), R.grad_GLM(*args)
  return getattr(R, f"lp_{tag}")(*args), getattr(R, f"grad_{tag}")(*args)


def route_b(tag, args, i):
  mp.dps = R.DPS_B
  f = R.lp_GLM if tag.startswith("GLM") else getattr(R, f"lp_{tag}")

  def g(v):
    a = list(args)
    a[i] = v
    return f(*a)

  return mp.diff(g, mpf(args[i]))


def cxx(x):
  return mp.nstr(mpf(x), 17, min_fixed=-1, max_fixed=-1, strip_zeros=False)


def main():
  out = []
  for tag, args, why in ROWS:
    value, grads = route_a(tag, args)
    worst = mpf(0)
    for g, i in zip(grads, WRT[tag]):
      d = route_b(tag, args, i)
      worst = max(worst, abs(d - g) / max(1, abs(g)))
    mp.dps = R.DPS_A
    out.append((tag, args, value, grads, why, worst))
    print(f"// {tag} {args}: route A vs B, worst relative {mp.nstr(worst, 3)}")
  print("\n// large_shapes_test.cpp rows")
  for tag, args, value, grads, why, _ in out:
    a = ", ".join(float(v).hex() for v in args)
    g = ", ".join(cxx(v) for v in grads)
    print(f'    {{"{tag}", {{{a}}}, {cxx(value)}, {{{g}}}}},  // {why}')
  print("\n// OpenCL rows {n, a, b, value, d_a, d_b} (GLM: see the struct)")
  for tag, args, value, grads, why, _ in out:
    g = ", ".join(cxx(v) for v in grads)
    if tag.startswith("GLM"):
      y, x, alpha, beta, phi = args
      print(f"      // {tag} {why}\n      {{{{{y}}}, {{{x}}}, {float(alpha).hex()}, "
            f"{float(beta).hex()}, {float(phi).hex()}, {cxx(value)}, {g}}},")
    else:
      n, a, b = args
      print(f"      // {tag} {why}\n      {{{n}, {float(a).hex()}, "
            f"{float(b).hex()}, {cxx(value)}, {g}}},")


if __name__ == "__main__":
  main()
