#!/usr/bin/env python3
"""Reference values for the unit tests of stan-dev/math#3409.

Computes d/da and d/db of the regularized incomplete beta I_z(a, b) at
every point of the prim tests by two independent routes and prints their
agreement in decimal digits:

  route A  mp.diff of I at 110 digits (ref_mpmath.dIda_A, dIdb_A)
  route B  density-form quadrature at 80 digits (ref_mpmath.dIda, dIdb)

Both routes work on the small side: where I > 1/2 they differentiate
1 - I_{1-z}(b, a), so that a derivative of 1e-307 next to I = 1 is not
lost.  It also prints the lcdf gradients of the two distribution tests.

Usage: python3 gen_test_refs.py [workers]   (default 8 workers)
"""

import sys
from multiprocessing import Pool

from mpmath import mp, mpf

import ref_mpmath as R

# (a, b, z) of the prim tests: grad_reg_inc_beta_test, inc_beta_dda_test,
# inc_beta_ddb_test.  At the last two points develop throws: in
# grad_reg_inc_beta (through hypergeometric_2F1) and in inc_beta_dda and
# inc_beta_ddb (iteration limit).
POINTS = [
    (1.0, 1.0, 0.4), (2.0, 3.0, 0.25), (0.5, 0.5, 0.5), (20.0, 20.0, 0.5),
    (1.5, 1.25, 0.6), (600.0, 600.0, 0.5), (352.0, 590.0, 0.758),
    (10.0, 1.0, 0.5), (1.0, 100.0, 0.01), (2500.0, 2.0, 0.999),
    (1.145, 6786.0, 1.8e-4), (446.0, 5.0, 0.76), (15000.0, 1.25, 0.999),
    (15000.0, 12500.0, 0.5), (15000.0, 12500.0, 0.6),
    (1.5, 1.25, 0.001), (1.5, 1.25, 0.5), (1.5, 1.25, 0.999),
    (1.5, 12500.0, 0.001),
    (4.0845138072907865, 150.97164240761546, 0.991606836273885),
    (936.5962745738543, 1294.2043377235907, 0.7877755442980134),
]


def route_a(a, b, z, wrt):
  """mp.diff at 110 digits on the small side."""
  with mp.workdps(R.DPS_DIFF):
    a, b, z = mpf(a), mpf(b), mpf(z)
    upper = R.betainc_reg(a, b, z) > mpf(1) / 2
    if not upper:
      if wrt == "a":
        return mp.diff(lambda t: R.betainc_reg(t, b, z), a)
      return mp.diff(lambda t: R.betainc_reg(a, t, z), b)
    zc = 1 - z
    # I_z(a, b) = 1 - I_{1-z}(b, a): d/da is minus the derivative in the
    # second shape of the reflected problem, and the reverse for d/db
    if wrt == "a":
      return -mp.diff(lambda t: R.betainc_reg(b, t, zc), a)
    return -mp.diff(lambda t: R.betainc_reg(t, a, zc), b)


def one(p):
  a, b, z = p
  with mp.workdps(R.DPS_INTEGRAL):
    i = R.betainc_reg(a, b, z)
  out = {"p": p, "I": i}
  for wrt, fb in (("a", R.dIda), ("b", R.dIdb)):
    vb = fb(a, b, z)
    va = route_a(a, b, z, wrt)
    out[wrt] = (vb, R.rel_digits(va, vb))
  return out


def lcdf_grads():
  """Gradients of the two lcdf tests, as d/dtheta of log I."""
  rows = []
  # the divisions run at the integral precision, not at 15 digits
  with mp.workdps(R.DPS_INTEGRAL):
    # beta_lcdf(0.5 | 600, 600): I_{0.5}(600, 600)
    i = R.betainc_reg(600, 600, mpf("0.5"))
    rows.append(("beta_lcdf(0.5 | 600, 600) d/dalpha",
                 R.dIda(600, 600, 0.5) / i, route_a(600, 600, 0.5, "a") / i))
    rows.append(("beta_lcdf(0.5 | 600, 600) d/dbeta",
                 R.dIdb(600, 600, 0.5) / i, route_a(600, 600, 0.5, "b") / i))
  # neg_binomial_lcdf(132 | 313, 0.068) = log I_p(313, 133),
  # p = beta / (1 + beta) with beta the double nearest 0.068
  with mp.workdps(R.DPS_INTEGRAL + 10):
    beta = mpf(0.068)
    p = beta / (1 + beta)
    i = R.betainc_reg(313, 133, p)
    rows.append(("neg_binomial_lcdf(132 | 313, 0.068) d/dalpha",
                 R.dIda(313, 133, p) / i, route_a(313, 133, p, "a") / i))
  return rows


def main():
  workers = int(sys.argv[1]) if len(sys.argv) > 1 else 8
  with Pool(workers) as pool:
    res = pool.map(one, POINTS)
  print("# d/da and d/db of I_z(a, b); digits = agreement of routes A and B")
  print(f"{'a':>20} {'b':>20} {'z':>18} {'I':>10}  "
        f"{'d/da':>24} {'dig':>5}  {'d/db':>24} {'dig':>5}")
  worst = mp.inf
  for r in res:
    a, b, z = r["p"]
    da, dga = r["a"]
    db, dgb = r["b"]
    worst = min(worst, dga, dgb)
    print(f"{a!r:>20} {b!r:>20} {z!r:>18} {mp.nstr(r['I'], 4):>10}  "
          f"{float(da)!r:>24} {float(dga):>5.1f}  "
          f"{float(db)!r:>24} {float(dgb):>5.1f}")
  print(f"# worst agreement: {float(worst):.1f} digits\n")
  print("# lcdf gradients")
  for name, vb, va in lcdf_grads():
    print(f"{name:<46} {float(vb)!r:>24} "
          f"{float(R.rel_digits(va, vb)):>5.1f}")
  return 0


if __name__ == "__main__":
  sys.exit(main())
