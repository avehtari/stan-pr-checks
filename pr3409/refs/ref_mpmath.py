#!/usr/bin/env python3
"""High-precision references for the incomplete gamma and beta gradients.

Two independent routes for every quantity, as STAN_MATH_DISTRIBUTIONS.md
section 3 requires.

  route A  numerical differentiation of the regularized function (mp.diff)
  route B  an explicit integral identity, integrated with mp.quad

The two routes share no code path inside mpmath apart from the elementary
functions, so agreement is evidence that both are right.

Quantities
  dQda(a, z)        d/da of the regularized upper incomplete gamma Q(a, z)
  dPda(a, z)        d/da of the regularized lower incomplete gamma P(a, z)
  dIda(a, b, z)     d/da of the regularized incomplete beta I_z(a, b)
  dIdb(a, b, z)     d/db of the same
  dIdz(a, b, z)     d/dz of the same

Run this file to print the route A against route B comparison table.
"""

import math
import sys

from mpmath import mp, mpf

# ── Reference precision ──────────────────────────────────────────────────
# 80 digits for the integrals, 110 for the numerical derivative, which
# loses roughly half of its digits. 50 digits is not sufficient for the
# integrals: at (a, z) = (50, 200) route B then agrees with route A to only
# 18 digits. At 80 digits the two routes are bit-identical there.
DPS_INTEGRAL = 80
DPS_DIFF = 110


# ── Route A: numerical differentiation ───────────────────────────────────
def dQda_A(a, z):
  with mp.workdps(DPS_DIFF):
    a, z = mpf(a), mpf(z)
    return mp.diff(lambda t: mp.gammainc(t, z, mp.inf, regularized=True), a)


def dIda_A(a, b, z):
  with mp.workdps(DPS_DIFF):
    a, b, z = mpf(a), mpf(b), mpf(z)
    return mp.diff(lambda t: betainc_reg(t, b, z), a)


def dIdb_A(a, b, z):
  with mp.workdps(DPS_DIFF):
    a, b, z = mpf(a), mpf(b), mpf(z)
    return mp.diff(lambda t: betainc_reg(a, t, z), b)


def dIdz_A(a, b, z):
  with mp.workdps(DPS_DIFF):
    a, b, z = mpf(a), mpf(b), mpf(z)
    return mp.diff(lambda t: mp.betainc(a, b, 0, t, regularized=True), z)


# ── Route B: explicit integral identities ────────────────────────────────
#
#   Q(a, z)   = Gamma(a, z) / Gamma(a)
#   dQ/da     = (1/Gamma(a)) * int_z^inf t^(a-1) e^-t ln t dt - Q * psi(a)
#
#   I_z(a, b) = B_z(a, b) / B(a, b)
#   dI/da     = (1/B(a,b)) * int_0^z t^(a-1)(1-t)^(b-1) ln t dt
#               - I * (psi(a) - psi(a+b))
#   dI/db     = (1/B(a,b)) * int_0^z t^(a-1)(1-t)^(b-1) ln(1-t) dt
#               - I * (psi(b) - psi(a+b))
#   dI/dz     = z^(a-1) (1-z)^(b-1) / B(a, b)
#
# Both integrals are evaluated in **density form**. The naive form
# integrates t^(a-1)(1-t)^(b-1), which is 1e-300 at its peak when a = b =
# 500 and 1e-3000 a short distance away. mp.quad judges convergence by the
# relative size of the correction, so it returns a stable wrong answer.
# Dividing by the normalizing constant inside the integrand makes the
# integrand O(1) and removes the failure. Splitting the interval at the
# mode and at multiples of the standard deviation gives the peak enough
# quadrature nodes.
def _gamma_nodes(a, z):
  """Interval split points for int_z^inf of the gamma density.

  Two node families are needed. Near the mode the width is sqrt(a). In the
  deep tail, where z is far above the mode, the decay length is 1, so the
  mode nodes all fall below z and leave one huge interval.
  """
  mode = a - 1
  sd = mp.sqrt(a)
  pts = [z]
  for j in (-4, -2, -1, 0, 1, 2, 4, 8, 16):
    t = mode + j * sd
    if t > z:
      pts.append(t)
  # Above the mode the density decays like exp(-(t - z) (z - a) / z), so
  # its decay length is z / (z - a). Nodes at fixed offsets 1, 2, ..., 32
  # left the interval [z + 32, inf) with about 1 % of the mass when a = 1e4
  # and z = 1.2e4, and the result was off by 5e-6. Scale the offsets by the
  # decay length instead, and go out far enough that the remainder is
  # exp(-4096).
  decay = max(mpf(1), z / max(z - a, mpf(1)))
  for j in (1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096):
    pts.append(z + j * decay)
  pts.append(mp.inf)
  return sorted(set(pts), key=lambda v: (v == mp.inf, v))


def _beta_nodes(a, b, z):
  """Interval split points for int_0^z of the beta density.

  Two node families. Near the mode the width is sd. When z is below the
  mode the mass of [0, z] sits at the right endpoint z, where the density
  rises like exp((t - z) * slope) with slope = (a-1)/z - (b-1)/(1-z). At
  (a, b, z) = (4933, 2112, 0.569) the decay length 1/slope is 2.6e-4, all
  mode nodes fall outside (0, z), and mp.quad on the single interval agreed
  with mp.diff to 1.1 digits (found 2026-09-20). Nodes at z - j/slope,
  j up to 4096, resolve the endpoint layer; beyond the last node the
  remainder is exp(-4096).
  """
  denom = a + b - 2
  mode = (a - 1) / denom if denom > 0 else mpf(1) / 2
  sd = mp.sqrt(a * b / ((a + b) ** 2 * (a + b + 1)))
  pts = [mpf(0), z]
  for j in (-8, -4, -2, -1, 0, 1, 2, 4, 8):
    t = mode + j * sd
    if 0 < t < z:
      pts.append(t)
  slope = (a - 1) / z - (b - 1) / (1 - z)
  if slope > 0:
    decay = 1 / slope
    for j in (1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096):
      t = z - j * decay
      if 0 < t < z:
        pts.append(t)
  return sorted(set(pts))


def dQda_B(a, z):
  """d/da of Q(a, z), by integrating the side that does not cancel.

  The identity int_0^inf density * (log t - psi(a)) dt = 0 holds exactly,
  because P + Q = 1 for every a. So the integral over [z, inf) and the
  integral over [0, z] are equal and opposite. One of them is a sum of
  same-sign-dominated contributions and the other is a near-total
  cancellation.

  For z far below the mode, Q is 1 to hundreds of digits and dQ/da is tiny.
  Evaluating it as the integral over [z, inf) then subtracts two O(1)
  quantities to get 1e-207, which needs more than 200 digits. At 80 digits
  the result is noise. Found on 2026-09-18: at (a, z) = (169.6, 3.9) this
  route returned -2.47e-84 where the true value is 1.97e-207, and both Stan
  and Eigen had it right.

  The repair is to integrate [0, z] when z is below the mode, which gives
  dP/da directly, and to negate.
  """
  with mp.workdps(DPS_INTEGRAL):
    a, z = mpf(a), mpf(z)
    lg = mp.loggamma(a)
    dig = mp.digamma(a)
    integrand = lambda t: mp.e ** ((a - 1) * mp.log(t) - t - lg) * (
        mp.log(t) - dig)
    if z < a - 1:
      # Below the mode: integrate [0, z], which gives dP/da, and negate.
      # The integrand there grows like t^(a-1), so all the mass sits at the
      # right endpoint. The substitution t = z * s^(1/a) gives
      # t^(a-1) dt = (z^a / a) ds and removes the concentration, leaving
      # only a log singularity at s = 0, which tanh-sinh handles.
      log_z = mp.log(z)
      scale = mp.e ** (a * log_z - lg) / a

      def sub_integrand(s):
        t = z * s ** (1 / a)
        return scale * mp.e ** (-t) * (log_z + mp.log(s) / a - dig)

      return -mp.quad(sub_integrand, [0, mpf(1) / 2, 1])
    return mp.quad(integrand, _gamma_nodes(a, z))


def betainc_reg(a, b, z):
  """Regularized incomplete beta I_z(a, b) at the working precision.

  mp.betainc evaluates a 2F1 series and fails to converge ("hypsum()
  failed") for a + b above about 14 000 when z is near the mode; found on
  2026-09-20 at (6579, 7823, 0.406). The fallback integrates the density
  over the smaller side, with nodes at the mode, and complements if it
  integrated the upper side. The density form has no cancellation, so this
  is accurate at the working precision.
  """
  a, b, z = mpf(a), mpf(b), mpf(z)
  try:
    return mp.betainc(a, b, 0, z, regularized=True)
  except ValueError:
    pass
  lbeta = mp.log(mp.beta(a, b))
  log_dens = lambda t: (a - 1) * mp.log(t) + (b - 1) * mp.log1p(-t) - lbeta
  denom = a + b - 2
  mode = (a - 1) / denom if denom > 0 else mpf(1) / 2

  def scaled_quad(nodes):
    # peak density factored out; see the note in _beta_grad_B on the
    # absolute error tolerance of mp.quad
    log_scale = max(log_dens(t) for t in nodes[:-1] + nodes[1:]
                    if 0 < t < 1)
    return mp.e ** log_scale * mp.quad(
        lambda t: mp.e ** (log_dens(t) - log_scale), nodes)

  if z <= mode:
    return scaled_quad(_beta_nodes(a, b, z))
  # upper side: the reflected problem, I_z(a, b) = 1 - I_{1-z}(b, a)
  return 1 - scaled_quad([1 - t for t in reversed(_beta_nodes(b, a, 1 - z))])


def _beta_grad_B(a, b, z, wrt):
  """int_0^z density * (log-factor - c) dt, for wrt in {'a', 'b'}.

  For a < 1 the factor t^(a-1) is singular at t = 0 and mp.quad
  under-resolves it. The substitution t = z * s^(1/a) gives
  t^(a-1) dt = (z^a / a) ds, which removes the singularity exactly. Only a
  log(s) endpoint singularity is left, and tanh-sinh handles that.
  """
  lbeta = mp.log(mp.beta(a, b))
  i = betainc_reg(a, b, z)
  if wrt == "a":
    c = mp.digamma(a) - mp.digamma(a + b)
  else:
    c = mp.digamma(b) - mp.digamma(a + b)

  # mp.quad stops on an absolute error estimate near eps of the working
  # precision, so an integrand of size 1e-114 is integrated to only 11
  # relative digits at 80 digits (found 2026-09-20 at (4933, 2112, 0.569);
  # a toy exp(linear) integrand scaled by 1e-114 reproduces it). The
  # constant factor is therefore pulled out of the integrand: the
  # quadrature sees an O(1) integrand and the scale is applied afterwards.
  if a < 1:
    log_z = mp.log(z)
    scale = mp.e ** (a * log_z - lbeta) / a

    def integrand(s):
      t = z * s ** (1 / a)
      log_factor = log_z + mp.log(s) / a if wrt == "a" else mp.log1p(-t)
      return (1 - t) ** (b - 1) * log_factor

    num = scale * mp.quad(integrand, [0, mpf(1) / 2, 1])
  else:
    nodes = _beta_nodes(a, b, z)
    log_dens = lambda t: (a - 1) * mp.log(t) + (b - 1) * mp.log1p(-t) - lbeta
    # the density peaks at the mode, or at z when z is below the mode
    log_scale = max(log_dens(t) for t in nodes[1:])

    def integrand(t):
      log_factor = mp.log(t) if wrt == "a" else mp.log1p(-t)
      return mp.e ** (log_dens(t) - log_scale) * log_factor

    num = mp.e ** log_scale * mp.quad(integrand, nodes)
  return num - i * c


def dIda_B(a, b, z):
  with mp.workdps(DPS_INTEGRAL):
    return _beta_grad_B(mpf(a), mpf(b), mpf(z), "a")


def dIdb_B(a, b, z):
  with mp.workdps(DPS_INTEGRAL):
    return _beta_grad_B(mpf(a), mpf(b), mpf(z), "b")


def dIdz_B(a, b, z):
  with mp.workdps(DPS_INTEGRAL):
    a, b, z = mpf(a), mpf(b), mpf(z)
    return z ** (a - 1) * (1 - z) ** (b - 1) / mp.beta(a, b)


# ── Derivative with respect to the second argument ───────────────────────
#
# This one is closed form, so it needs no quadrature and no second route:
#
#   dP/dz =  z^(a-1) e^(-z) / Gamma(a)     the gamma density
#   dQ/dz = -z^(a-1) e^(-z) / Gamma(a)
#
# It is the `der2` term in the fwd and rev gamma_p and gamma_q code. The
# 2026-09-19 change evaluates it in log space, so its overflow boundary
# moves; that is why it is measured.
def dPdz(a, z):
  with mp.workdps(DPS_INTEGRAL):
    a, z = mpf(a), mpf(z)
    return mp.e ** ((a - 1) * mp.log(z) - z - mp.loggamma(a))


def dQdz(a, z):
  return -dPdz(a, z)


def dPdz_A(a, z):
  """Numerical differentiation, as a check on the closed form."""
  with mp.workdps(DPS_DIFF):
    a, z = mpf(a), mpf(z)
    return mp.diff(lambda t: mp.gammainc(a, 0, t, regularized=True), z)


# ── Route C: the Gautschi series, term by term in high precision ─────────
#
#   dP/da = sum_n exp(-z + (a+n) log z - lgamma(a+n+1)) * (log z - psi(a+n+1))
#
# Every term is negative while psi(a+n+1) > log z, which holds for all n
# when z is at or below the mode, so there is no cancellation there and the
# sum is stable at modest precision. Above the mode the factor changes sign
# at n ~ z - a and the sum cancels; the working precision grows with z - a
# to absorb that, and beyond a few standard deviations the upper-side
# quadrature takes over instead.
def dPda_C(a, z, extra_dps=0):
  with mp.workdps(DPS_INTEGRAL + extra_dps):
    a, z = mpf(a), mpf(z)
    lz = mp.log(z)
    s = mpf(0)
    tiny = mpf(10) ** (-(DPS_INTEGRAL + extra_dps) - 5)
    for n in range(0, 2000000):
      term = mp.e ** (-z + (a + n) * lz - mp.loggamma(a + n + 1)) * (
          lz - mp.digamma(a + n + 1))
      s += term
      if n > 50 and abs(term) < abs(s) * tiny:
        return s
    return s


# ── Public references ────────────────────────────────────────────────────
#
# Three failures of the quadrature route were found while this file was in
# use, each caught by a second route:
#   1. the unnormalized integrand spans 300 decades and mp.quad returned a
#      stable wrong answer; fixed by integrating in density form;
#   2. integrating the cancelling side for z far below the mode; fixed by
#      integrating the lower side there;
#   3. the substitution t = z * s^(1/a), added for the a < 1 endpoint
#      singularity, compresses the density peak into a spike of width
#      exp(-sqrt(a)) at s -> 1 when a is large, and mp.quad cannot see it.
#      At (a, z) = (10104.8, 10103.3) it returned 0.003367 where three
#      independent routes give 0.0039686.
#
# The public function therefore uses the series wherever the series has no
# cancellation, which is z at or below the mode and a few standard
# deviations above it, and the upper-side quadrature only beyond that.
def dQda(a, z):
  af, zf = float(a), float(z)
  sd = math.sqrt(max(af, 1.0))
  if zf <= af + 3.0 * sd:
    # series; cancellation grows with z - a, so grow the precision with it
    extra = int(max(0.0, zf - af) * 0.5) + 10
    return -dPda_C(a, z, extra_dps=extra)
  return dQda_B(a, z)


def dPda(a, z):
  return -dQda_B(a, z)


def _beta_small_side(a, b, z):
  """True when I_z(a, b) > 1/2, so the upper side [z, 1] is the small one.

  The quadrature route computes int_0^z density * log-factor - I * c. When
  I is 1 to many digits, the two terms cancel to the working precision and
  the result is the 80-digit noise floor: on 2026-09-20 the value
  -1.687e-80 was returned at three unrelated points where mp.diff gave 0.
  Same rule as dQda: integrate the side that does not cancel. Because
  I_z(a, b) = 1 - I_{1-z}(b, a), the derivative on the upper side is the
  negated derivative of the reflected problem with the parameters swapped.
  """
  with mp.workdps(DPS_INTEGRAL):
    return betainc_reg(a, b, z) > mpf(1) / 2


def _one_minus(z):
  """1 - z at the integral precision. At the global default of 15 digits
  the complement of a double is rounded to 53 bits, and the reflected
  point is then not 1 - z: the two routes differed at the 14th digit at
  (0.4469, 253.8, 0.0405) until this was fixed on 2026-09-20."""
  with mp.workdps(DPS_INTEGRAL + 10):
    return 1 - mpf(z)


def dIda(a, b, z):
  # negate at the integral precision; at the default 15 digits the
  # negation rounds the reference to a double
  with mp.workdps(DPS_INTEGRAL):
    if _beta_small_side(a, b, z):
      return -dIdb_B(b, a, _one_minus(z))
    return dIda_B(a, b, z)


def dIdb(a, b, z):
  with mp.workdps(DPS_INTEGRAL):
    if _beta_small_side(a, b, z):
      return -dIda_B(b, a, _one_minus(z))
    return dIdb_B(a, b, z)


def dIdz(a, b, z):
  return dIdz_B(a, b, z)


# ── Agreement check ──────────────────────────────────────────────────────
def rel_digits(x, y):
  """Number of agreeing decimal digits between two mpf values."""
  with mp.workdps(60):
    x, y = mpf(x), mpf(y)
    if x == y:
      return mp.inf
    if x == 0 or y == 0:
      return -mp.log10(abs(x - y))
    return -mp.log10(abs((x - y) / x))


GAMMA_POINTS = [
    # (a, z): small a, integer a, z near a, z >> a, deep tail
    (0.1, 0.05), (0.5, 0.25), (1.0, 1.0), (2.0, 2.0), (3.0, 3.0),
    (5.0, 5.0), (8.0, 8.0), (10.0, 30.0), (0.5, 20.0), (2.0, 40.0),
    (20.0, 10.0), (100.0, 100.0), (100.0, 150.0), (0.05, 0.001),
    (1.0, 0.001), (50.0, 200.0), (4.0, 7.9), (4.0, 8.1), (8.0, 8.0),
    (12500.0, 12000.0),
]

BETA_POINTS = [
    (0.5, 0.5, 0.5), (1.0, 1.0, 0.5), (2.0, 3.0, 0.25), (2.0, 3.0, 0.75),
    (0.1, 0.1, 0.5), (5.0, 2.0, 0.9), (5.0, 2.0, 0.999), (500.0, 500.0, 0.5),
    (0.5, 5000.0, 0.001), (2500.0, 2.0, 0.999), (1.0, 100.0, 0.01),
    (100.0, 1.0, 0.99), (12500.0, 12500.0, 0.5), (3.0, 4.0, 0.001),
    (0.25, 0.75, 0.999),
]


def main():
  print("# Route A (mp.diff) against route B (integral identity)")
  print("# value: agreeing decimal digits. inf means bit-identical.\n")

  print("## d/da of Q(a, z), the regularized upper incomplete gamma\n")
  print(f"{'a':>10} {'z':>10} {'digits':>10}   {'reference':>28}")
  worst = mp.inf
  for a, z in GAMMA_POINTS:
    ra, rb = dQda_A(a, z), dQda_B(a, z)
    d = rel_digits(ra, rb)
    worst = min(worst, d)
    print(f"{a:>10g} {z:>10g} {float(d):>10.1f}   {mp.nstr(rb, 20):>28}")
  print(f"\nworst agreement: {float(worst):.1f} digits\n")

  for name, fa, fb in (("d/da", dIda_A, dIda_B), ("d/db", dIdb_A, dIdb_B),
                       ("d/dz", dIdz_A, dIdz_B)):
    print(f"## {name} of I_z(a, b), the regularized incomplete beta\n")
    print(f"{'a':>10} {'b':>10} {'z':>10} {'digits':>10}   "
          f"{'reference':>28}")
    worst = mp.inf
    for a, b, z in BETA_POINTS:
      ra, rb = fa(a, b, z), fb(a, b, z)
      d = rel_digits(ra, rb)
      worst = min(worst, d)
      print(f"{a:>10g} {b:>10g} {z:>10g} {float(d):>10.1f}   "
            f"{mp.nstr(rb, 20):>28}")
    print(f"\nworst agreement: {float(worst):.1f} digits\n")

  return 0


if __name__ == "__main__":
  sys.exit(main())
