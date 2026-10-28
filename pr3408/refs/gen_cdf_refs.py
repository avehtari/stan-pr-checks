#!/usr/bin/env python3
# Two-route references for incomplete_gamma_cdfs_reference.hpp, the
# 18 cdfs at 18 points (324 rows) and the 10 root_table rows of #3408.
"""Two-route references for the 18 incomplete-gamma cdf functions.

Route A: the closed form of each function in mpmath.gammainc (log P and
log Q from the side without cancellation), gradient by mp.diff with
respect to each argument (survey_cases.py).
Route B: log P or log Q and its a derivative by a quadrature of the scaled
density (upper side for z > a - 1, lower side otherwise), the z derivative
as the density divided by the probability, and the chain rule of each
function.  A row is kept only where value and gradient agree between the
routes to 1e-25 relative.

Usage: gen_cdf_refs.py OUT_HEADER
"""
import sys
import mpmath as mp

import survey_cases as sc
from probe_lccdf_tail import ref_quad
from probe_lower_series import ref_quad_lower

mp.mp.dps = 60
AGREE = mp.mpf(10)**-25

# (a, z) as a function of the arguments, and d(a, z)/d(argument) per argument
MAP = {
    "chi_square": (lambda y, nu, _: (nu / 2, y / 2),
                   lambda y, nu, _: [(0, mp.mpf(1) / 2), (mp.mpf(1) / 2, 0)]),
    "gamma": (lambda y, al, be: (al, be * y),
              lambda y, al, be: [(0, be), (1, 0), (0, y)]),
    "inv_gamma": (lambda y, al, be: (al, be / y),
                  lambda y, al, be: [(0, -be / y**2), (1, 0), (0, 1 / y)]),
    "inv_chi_square": (lambda y, nu, _: (nu / 2, 1 / (2 * y)),
                       lambda y, nu, _: [(0, -1 / (2 * y**2)), (mp.mpf(1) / 2, 0)]),
    "scaled_inv_chi_square": (
        lambda y, nu, s: (nu / 2, nu * s * s / (2 * y)),
        lambda y, nu, s: [(0, -nu * s * s / (2 * y**2)),
                          (mp.mpf(1) / 2, s * s / (2 * y)), (0, nu * s / y)]),
    "poisson": (lambda n, lam, _: (n + 1, lam),
                lambda n, lam, _: [(0, 1)]),
    "poisson_log": (lambda n, al, _: (n + 1, mp.exp(al)),
                    lambda n, al, _: [(0, mp.exp(al))]),
}
INVERSE = {"inv_gamma", "inv_chi_square", "scaled_inv_chi_square", "poisson",
           "poisson_log"}


def split(name):
  fam, kind = name.rsplit("_", 1)
  upper = (kind in ("cdf", "lcdf")) == (fam in INVERSE)
  return fam, kind, upper


def route_b(name, y, p1, p2):
  fam, kind, upper = split(name)
  a_of, jac = MAP[fam]
  a, z = a_of(y, p1, p2)
  if z > a - 1:
    lq, dq = ref_quad(a, z)
    q = mp.exp(lq)
    lp = mp.log1p(-q)
    dp = -q / mp.exp(lp) * dq
  else:
    lp, dp = ref_quad_lower(a, z)
    p = mp.exp(lp)
    lq = mp.log1p(-p)
    dq = -p / mp.exp(lq) * dp
  log_density = (a - 1) * mp.log(z) - z - mp.loggamma(a)
  if upper:
    lv, da, dz = lq, dq, -mp.exp(log_density - lq)
  else:
    lv, da, dz = lp, dp, mp.exp(log_density - lp)
  grads = [ja * da + jz * dz for ja, jz in jac(y, p1, p2)]
  if kind == "cdf":
    v = mp.exp(lv)
    return v, [v * g for g in grads]
  return lv, grads


def route_a(name, y, p1, p2):
  fam, kind, upper = split(name)
  if name == "poisson_cdf":
    ref = lambda n, lam, _: sc.Q(n + 1, lam)
    compl = lambda n, lam, _: sc.P(n + 1, lam)
  elif name in sc.POISSON_LOG:
    ref = sc.POISSON_LOG[name][1]
    compl = ((lambda n, al, _: sc.P(n + 1, mp.exp(al)))
             if name == "poisson_log_cdf" else None)
  else:
    ref = sc.FUNCS[name][1]
    compl = sc.COMPL.get(name)
  args = [y, p1, p2]
  val = ref(*args)
  use_compl = compl is not None and val > 0.5
  g_ref = compl if use_compl else ref
  sign = -1 if use_compl else 1
  grads = []
  for k in range(3):
    if args[k] is None or (k == 0 and fam.startswith("poisson")):
      continue
    def f(t, k=k):
      b = list(args)
      b[k] = t
      return g_ref(*b)
    grads.append(sign * mp.diff(f, args[k]))
  return val, grads


def rel(x, r):
  return abs(x - r) / max(abs(r), mp.mpf(10)**-300)


def args_of_name(name):
  if name == "poisson_cdf":
    return sc.FUNCS["poisson_lcdf"][0]
  if name in sc.POISSON_LOG:
    return sc.POISSON_LOG[name][0]
  return sc.FUNCS[name][0]


# Explicit (function, y, p1, p2) cases: the distribution rows of the
# gradient-root table in the PR text, z = a with the shape at 100 (the
# old root wrong by 0.58), 171 (0.72) and 180 (NaN, tgamma overflows).
ROOT_TABLE_CASES = [
    ("chi_square_lccdf", 200.0, 200.0, None),
    ("chi_square_lccdf", 342.0, 342.0, None),
    ("chi_square_lccdf", 360.0, 360.0, None),
    ("chi_square_lcdf", 200.0, 200.0, None),
    ("chi_square_lcdf", 360.0, 360.0, None),
    ("gamma_cdf", 100.0, 100.0, 1.0),
    ("gamma_cdf", 180.0, 180.0, 1.0),
    ("inv_gamma_lccdf", 1.0 / 100.0, 100.0, 1.0),
    ("inv_gamma_lccdf", 1.0 / 180.0, 180.0, 1.0),
    ("scaled_inv_chi_square_lcdf", 1.0, 360.0, 1.0),
]


def cases(names):
  """(name, y, p1, p2, class) in double, from the points and the table."""
  out = []
  for name in names:
    args_of = args_of_name(name)
    for a, z, cls in sc.POINTS:
      y, p1, p2 = args_of(mp.mpf(a), mp.mpf(z))
      out.append((name, float(y), float(p1),
                  None if p2 is None else float(p2), cls))
  for name, y, p1, p2 in ROOT_TABLE_CASES:
    if name in names:
      out.append((name, y, p1, p2, "root_table"))
  return out


def main():
  out = sys.argv[1]
  # optional: comma-separated function names, header tag and namespace
  names = (sys.argv[2].split(",") if len(sys.argv) > 2
           else list(sc.FUNCS.keys()) + ["poisson_cdf"])
  tag = sys.argv[3] if len(sys.argv) > 3 else "INCOMPLETE_GAMMA_CDFS"
  ns = sys.argv[4] if len(sys.argv) > 4 else "incomplete_gamma_cdfs"
  rows = []
  for name, y, p1, p2, cls in cases(names):
    args = [mp.mpf(y), mp.mpf(p1), None if p2 is None else mp.mpf(p2)]
    try:
      vb, gb = route_b(name, *args)
      va, ga = route_a(name, *args)
    except Exception as e:
      print(f"skip {name} ({y}, {p1}, {p2}): {e}")
      continue
    worst = max([rel(va, vb)] + [rel(x, r) for x, r in zip(ga, gb)])
    if worst > AGREE:
      print(f"skip {name} ({y}, {p1}, {p2}): routes differ by "
            f"{mp.nstr(worst, 3)}")
      continue
    rows.append((name, y, p1, p2, vb, gb, cls))
  with open(out, "w") as f:
    f.write(HEADER.replace("INCOMPLETE_GAMMA_CDFS", tag)
            .replace("incomplete_gamma_cdfs", ns))
    for name, y, p1, p2, v, g, cls in rows:
      gs = list(g) + [0] * (3 - len(g))
      nums = [y, p1, p2 if p2 is not None else 0.0]
      f.write(f'      {{"{name}", {nums[0]!r}, {nums[1]!r}, {nums[2]!r}, '
              f'{mp.nstr(v, 17)}, {{{mp.nstr(gs[0], 17)}, {mp.nstr(gs[1], 17)}, '
              f'{mp.nstr(gs[2], 17)}}}, {len(g)}}},  // {cls}\n')
    f.write(FOOTER.replace("incomplete_gamma_cdfs", ns))
  print(f"{len(rows)} rows")


HEADER = """#ifndef TEST_UNIT_MATH_REV_PROB_INCOMPLETE_GAMMA_CDFS_REFERENCE_HPP
#define TEST_UNIT_MATH_REV_PROB_INCOMPLETE_GAMMA_CDFS_REFERENCE_HPP

#include <string>
#include <vector>

namespace stan {
namespace test {
namespace incomplete_gamma_cdfs {

/**
 * Values and gradients of the cdf, lcdf and lccdf of chi_square, gamma,
 * inv_chi_square, inv_gamma, scaled_inv_chi_square and poisson, computed
 * with mpmath at 60 digits by two routes: the closed form in the
 * regularized incomplete gamma function with numerical differentiation,
 * and a quadrature of the scaled density with the chain rule. Only rows
 * where the routes agree to 1e-25 are kept. Each row is one point (a, z)
 * of the incomplete gamma function, mapped to the arguments of the
 * function; the comment names the region. The rows marked root_table put
 * z = a with the shape at 100, 171 and 180 (beta = 1 and s = 1 where the
 * function has them). The gradient is with respect to (y, p1, p2); for
 * poisson with respect to lambda only.
 */
struct Row {
  std::string function;
  double y;
  double p1;
  double p2;
  double value;
  double grad[3];
  int n_grad;
};

inline const std::vector<Row>& rows() {
  // clang-format off
  // NOLINTBEGIN(whitespace/line_length)
  static const std::vector<Row> r = {
"""

FOOTER = """  };
  // NOLINTEND
  // clang-format on
  return r;
}

}  // namespace incomplete_gamma_cdfs
}  // namespace test
}  // namespace stan

#endif
"""

if __name__ == "__main__":
  main()
