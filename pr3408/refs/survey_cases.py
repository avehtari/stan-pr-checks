#!/usr/bin/env python3
# Cases and mpmath references of the cdf survey (17 functions x 18 points).
# Supports the 2e-13 / 8e-13 / 2.7e-10 claim; run from cdfs/.
"""Survey of the incomplete-gamma distribution functions.

Writes survey_cases.txt (one case per line: func y p1 p2 class) and
survey_refs.txt (reference value and gradient per case, mpmath 40 digits,
mp.diff of the closed form in the incomplete gamma function).

Each case is a point (a, z) of the regularized incomplete gamma function,
mapped back to the arguments of the distribution function.
"""
import mpmath as mp

mp.mp.dps = 60

# (a, z, class): class names the region of the incomplete gamma function
POINTS = [
    (10, 1000, "upper_deep"), (150, 1200, "upper_deep"), (200, 1600, "upper_deep"),
    (10, 300, "upper_tail"), (150, 600, "upper_tail"),
    (5, 7.5, "bulk"), (50, 60, "bulk"), (200, 200, "bulk_big"),
    (200, 300, "bulk_big"), (400, 400, "bulk_big"), (1000, 1100, "bulk_big"),
    (1000, 900, "bulk_big"), (400, 300, "bulk_big"),
    (200, 50, "lower_tail"), (50, 5, "lower_tail"),
    (1000, 100, "lower_deep"), (500, 30, "lower_deep"), (100, 0.001, "lower_deep"),
]


def _p(a, z):
  return mp.gammainc(a, 0, z, regularized=True)


def _q(a, z):
  return mp.gammainc(a, z, mp.inf, regularized=True)


def P(a, z):
  return _p(a, z)


def Q(a, z):
  return _q(a, z)


class LogP:
  """log P from the side that does not cancel."""

  def __call__(self, a, z):
    q = _q(a, z)
    return mp.log1p(-q) if q < 0.5 else mp.log(_p(a, z))


class LogQ:
  def __call__(self, a, z):
    p = _p(a, z)
    return mp.log1p(-p) if p < 0.5 else mp.log(_q(a, z))


log_P = LogP()
log_Q = LogQ()


# name: (args from (a, z), reference function of the arguments)
# args are (y, p1, p2); p2 = None for two-argument functions
FUNCS = {
    "chi_square_lcdf": (lambda a, z: (2 * z, 2 * a, None),
                        lambda y, nu, _: log_P(nu / 2, y / 2)),
    "chi_square_lccdf": (lambda a, z: (2 * z, 2 * a, None),
                         lambda y, nu, _: log_Q(nu / 2, y / 2)),
    "chi_square_cdf": (lambda a, z: (2 * z, 2 * a, None),
                       lambda y, nu, _: P(nu / 2, y / 2)),
    "gamma_lcdf": (lambda a, z: (z / 2, a, 2),
                   lambda y, al, be: log_P(al, be * y)),
    "gamma_lccdf": (lambda a, z: (z / 2, a, 2),
                    lambda y, al, be: log_Q(al, be * y)),
    "gamma_cdf": (lambda a, z: (z / 2, a, 2),
                  lambda y, al, be: P(al, be * y)),
    "inv_gamma_lcdf": (lambda a, z: (2 / z, a, 2),
                       lambda y, al, be: log_Q(al, be / y)),
    "inv_gamma_lccdf": (lambda a, z: (2 / z, a, 2),
                        lambda y, al, be: log_P(al, be / y)),
    "inv_gamma_cdf": (lambda a, z: (2 / z, a, 2),
                      lambda y, al, be: Q(al, be / y)),
    "inv_chi_square_lcdf": (lambda a, z: (1 / (2 * z), 2 * a, None),
                            lambda y, nu, _: log_Q(nu / 2, 1 / (2 * y))),
    "inv_chi_square_lccdf": (lambda a, z: (1 / (2 * z), 2 * a, None),
                             lambda y, nu, _: log_P(nu / 2, 1 / (2 * y))),
    "inv_chi_square_cdf": (lambda a, z: (1 / (2 * z), 2 * a, None),
                           lambda y, nu, _: Q(nu / 2, 1 / (2 * y))),
    "scaled_inv_chi_square_cdf": (
        lambda a, z: (a / z, 2 * a, 1),
        lambda y, nu, s: Q(nu / 2, nu * s * s / (2 * y))),
    "scaled_inv_chi_square_lcdf": (
        lambda a, z: (a / z, 2 * a, 1),
        lambda y, nu, s: log_Q(nu / 2, nu * s * s / (2 * y))),
    "scaled_inv_chi_square_lccdf": (
        lambda a, z: (a / z, 2 * a, 1),
        lambda y, nu, s: log_P(nu / 2, nu * s * s / (2 * y))),
    "poisson_lcdf": (lambda a, z: (a - 1, z, None),
                     lambda n, lam, _: log_Q(n + 1, lam)),
    "poisson_lccdf": (lambda a, z: (a - 1, z, None),
                      lambda n, lam, _: log_P(n + 1, lam)),
}
DISCRETE = {"poisson_lcdf", "poisson_lccdf"}
# poisson_log: (n, alpha) with a = n + 1 and z = exp(alpha); not part of
# the survey output (main() writes FUNCS only), used by gen_cdf_refs.py
POISSON_LOG = {
    "poisson_log_lcdf": (lambda a, z: (a - 1, mp.log(z), None),
                         lambda n, al, _: log_Q(n + 1, mp.exp(al))),
    "poisson_log_lccdf": (lambda a, z: (a - 1, mp.log(z), None),
                          lambda n, al, _: log_P(n + 1, mp.exp(al))),
    "poisson_log_cdf": (lambda a, z: (a - 1, mp.log(z), None),
                        lambda n, al, _: Q(n + 1, mp.exp(al))),
}
# complement of each probability-scale function: differentiate it (and
# negate) where the value is above 1/2, so that the small gradient is
# not lost in 1 - (small)
COMPL = {
    "chi_square_cdf": lambda y, nu, _: Q(nu / 2, y / 2),
    "gamma_cdf": lambda y, al, be: Q(al, be * y),
    "inv_gamma_cdf": lambda y, al, be: P(al, be / y),
    "inv_chi_square_cdf": lambda y, nu, _: P(nu / 2, 1 / (2 * y)),
    "scaled_inv_chi_square_cdf": lambda y, nu, s: P(nu / 2, nu * s * s / (2 * y)),
}


def main():
  with open("survey_cases.txt", "w") as fc, open("survey_refs.txt", "w") as fr:
    for name, (args_of, ref) in FUNCS.items():
      for a, z, cls in POINTS:
        y, p1, p2 = args_of(mp.mpf(a), mp.mpf(z))
        y, p1 = float(y), float(p1)
        p2 = None if p2 is None else float(p2)
        args = [mp.mpf(y), mp.mpf(p1), None if p2 is None else mp.mpf(p2)]
        try:
          val = ref(*args)
          use_compl = name in COMPL and val > 0.5
          g_ref = COMPL[name] if use_compl else ref
          sign = -1 if use_compl else 1
          grads = []
          for k in range(3):
            if args[k] is None or (k == 0 and name in DISCRETE):
              continue
            def f(t, k=k):
              b = list(args)
              b[k] = t
              return g_ref(*b)
            grads.append(sign * mp.diff(f, args[k]))
        except Exception as e:  # mpmath non-convergence
          print(f"skip {name} {a} {z}: {e}")
          continue
        fc.write(f"{name} {y!r} {p1!r} {p2 if p2 is not None else 'nan'} {cls} {a} {z}\n")
        fr.write(" ".join(mp.nstr(v, 20) for v in [val] + grads) + "\n")


if __name__ == "__main__":
  main()
