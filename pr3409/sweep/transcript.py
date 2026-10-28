#!/usr/bin/env python3
"""Binary64 transcription of the Stan Math incomplete gamma and beta
gradient roots, from develop at 5252d51d47.

The method is the one used for #3406: transcribe the algorithm, measure the
math, and leave the C++ plumbing to a later build. Python floats are IEEE
binary64 and Python arithmetic is the same as C++ double arithmetic, so the
transcription reproduces the rounding of the original.

Two deliberate differences from the C++, both stated in the note:

1. `gamma_p`, `gamma_q` and `inc_beta` come from mpmath, correctly rounded
   to binary64. The C++ calls Boost. A correctly rounded value isolates the
   error of the gradient algorithm from the error of the value function.
   Any error measured here is therefore a lower limit for the C++.
2. `digamma`, `tgamma` and `lgamma` also come from mpmath, correctly
   rounded, for the same reason.

Each function returns (value, iterations, branch_name) so the iteration
count is available for the cost question.
"""

import functools
import math

from mpmath import mp, mpf

mp.dps = 40


# ── Correctly rounded scalar helpers ─────────────────────────────────────
@functools.lru_cache(maxsize=1 << 18)
def gamma_p(a, z):
  """Regularized lower incomplete gamma, correctly rounded."""
  return float(mp.gammainc(mpf(a), 0, mpf(z), regularized=True))


@functools.lru_cache(maxsize=1 << 18)
def gamma_q(a, z):
  """Regularized upper incomplete gamma, correctly rounded."""
  return float(mp.gammainc(mpf(a), mpf(z), mp.inf, regularized=True))


@functools.lru_cache(maxsize=1 << 18)
def inc_beta(a, b, z):
  """Regularized incomplete beta, correctly rounded.

  Goes through ref_mpmath.betainc_reg, which falls back to density-form
  quadrature where mp.betainc fails (a + b above about 14 000 near the
  mode).
  """
  import ref_mpmath
  return float(ref_mpmath.betainc_reg(a, b, z))


@functools.lru_cache(maxsize=1 << 18)
def digamma(x):
  return float(mp.digamma(mpf(x)))


@functools.lru_cache(maxsize=1 << 18)
def tgamma(x):
  """std::tgamma semantics: +inf above the binary64 range (x > ~171.62)."""
  v = mp.gamma(mpf(x))
  if v > mpf('1.7976931348623157e308'):
    return math.inf
  return float(v)


def multiply_log(a, b):
  """stan::math::multiply_log(a, b) = a * log(b), with 0 * log(0) = 0."""
  if b == 0.0 and a == 0.0:
    return 0.0
  return a * math.log(b)


def cexp(x):
  """std::exp semantics: overflow gives +inf, it does not raise.

  math.exp raises OverflowError, so every exp in the transcription goes
  through this wrapper. Without it the transcription cannot reproduce the
  inf and NaN results that the C++ produces for large a.
  """
  if x > 709.782712893384:
    return math.inf
  if x < -745.2:
    return 0.0
  return math.exp(x)


INFTY = math.inf


class NotConverged(Exception):
  pass


# ── grad_reg_inc_gamma: d/da of Q(a, z) ──────────────────────────────────
def grad_reg_inc_gamma(a, z, g=None, dig=None, precision=1e-6,
                       max_steps=100000):
  """prim/fun/grad_reg_inc_gamma.hpp, develop 5252d51d47.

  Two branches. `z >= a and z >= 8` takes a fixed 10-term asymptotic
  expansion, which ignores `precision` completely. Otherwise an alternating
  series runs until the **absolute** term size falls below `precision`.
  """
  if g is None:
    g = tgamma(a)
  if dig is None:
    dig = digamma(a)

  l = math.log(z)

  if z >= a and z >= 8:
    # asymptotic expansion, http://dlmf.nist.gov/8.11#E2 ; 10 fixed terms
    S = 0.0
    a_minus_one_minus_k = a - 1.0
    fac = a_minus_one_minus_k
    dfac = 1.0
    zpow = z
    delta = dfac / zpow
    for k in range(1, 10):
      a_minus_one_minus_k -= 1.0
      S += delta
      zpow *= z
      dfac = a_minus_one_minus_k * dfac + fac
      fac *= a_minus_one_minus_k
      delta = dfac / zpow
      if math.isinf(delta):
        raise NotConverged("grad_reg_inc_gamma is not converging")
    value = gamma_q(a, z) * (l - dig) + cexp(-z + (a - 1.0) * l) * S / g
    return value, 9, "asymptotic"

  # gradient of the series expansion, http://dlmf.nist.gov/8.7#E3
  S = 0.0
  log_s = 0.0
  s_sign = 1.0
  log_z = math.log(z)
  log_delta = log_s - multiply_log(2, a)
  log_precision = math.log(precision)
  for k in range(1, max_steps + 1):
    S += cexp(log_delta) if s_sign >= 0.0 else -cexp(log_delta)
    log_s += log_z - math.log(k)
    s_sign = -s_sign
    log_delta = log_s - multiply_log(2, k + a)
    if math.isinf(log_delta):
      raise NotConverged("grad_reg_inc_gamma is not converging")
    if log_delta <= log_precision:
      value = gamma_p(a, z) * (dig - l) + cexp(a * l) * S / g
      return value, k, "series"
  raise NotConverged("grad_reg_inc_gamma exceeded max_steps")


# ── grad_reg_lower_inc_gamma: d/da of P(a, z) ────────────────────────────
def grad_reg_lower_inc_gamma(a, z, precision=1e-10, max_steps=100000):
  """prim/fun/grad_reg_lower_inc_gamma.hpp, develop 5252d51d47.

  For large z relative to a it forwards to grad_reg_inc_gamma, with its own
  `precision`. Otherwise two Gautschi series, each with an **absolute**
  stopping rule.
  """
  if z == 0.0:
    return 0.0, 0, "zero"

  # The third test is `a < sqrt(-756 - z*z + 60*z)`. In C++ a negative
  # argument gives NaN and the comparison is false, so guard it here.
  disc = -756.0 - z * z + 60.0 * z
  third = a < math.sqrt(disc) if disc >= 0.0 else False
  if (a < 0.8 and z > 15.0) or (a < 12.0 and z > 30.0) or third:
    val, iters, _ = grad_reg_inc_gamma(a, z, tgamma(a), digamma(a), precision,
                                       max_steps)
    return -val, iters, "delegated"

  log_z = math.log(z)
  emz = cexp(-z)

  n = 0
  a_plus_n = a
  sum_a = 0.0
  lgamma_a_plus_1 = math.lgamma(a + 1.0)
  lgamma_a_plus_n_plus_1 = lgamma_a_plus_1
  iters_a = 0
  while True:
    term = cexp(a_plus_n * log_z - lgamma_a_plus_n_plus_1)
    sum_a += term
    if term <= precision:
      break
    if n >= max_steps:
      raise NotConverged("grad_reg_lower_inc_gamma sum_a")
    n += 1
    lgamma_a_plus_n_plus_1 += math.log1p(a_plus_n)
    a_plus_n += 1.0
  iters_a = n

  n = 1
  a_plus_n = a + 1.0
  digammap1 = digamma(a_plus_n)
  sum_b = digammap1 * cexp(a * log_z - lgamma_a_plus_1)
  lgamma_a_plus_n_plus_1 = lgamma_a_plus_1 + math.log(a_plus_n)
  while True:
    digammap1 += 1.0 / a_plus_n
    term = cexp(a_plus_n * log_z - lgamma_a_plus_n_plus_1) * digammap1
    sum_b += term
    if term <= precision:
      return emz * (log_z * sum_a - sum_b), iters_a + n, "gautschi"
    if n >= max_steps:
      raise NotConverged("grad_reg_lower_inc_gamma sum_b")
    n += 1
    lgamma_a_plus_n_plus_1 += math.log1p(a_plus_n)
    a_plus_n += 1.0


# ── inc_beta_dda and inc_beta_ddb ────────────────────────────────────────
def inc_beta_dda(a, b, z, digamma_a=None, digamma_ab=None, threshold=1e-10):
  """prim/fun/inc_beta_dda.hpp, develop 5252d51d47.

  The stopping rule is `fabs(summand) > threshold`, an **absolute** test on
  a series whose terms carry the prefactor (a+1)/(a+b) cubed.
  """
  if digamma_a is None:
    digamma_a = digamma(a)
  if digamma_ab is None:
    digamma_ab = digamma(a + b)

  if b > a:
    if ((0.1 < z <= 0.75 and b > 500) or (0.01 < z <= 0.1 and b > 2500)
        or (0.001 < z <= 0.01 and b > 1e5)):
      v, it, _ = inc_beta_ddb(b, a, 1.0 - z, digamma_a, digamma_ab, threshold)
      return -v, it, "reflected"
  if ((z > 0.75 and a < 500) or (z > 0.9 and a < 2500)
      or (z > 0.99 and a < 1e5) or z > 0.999):
    v, it, _ = inc_beta_ddb(b, a, 1.0 - z, digamma_a, digamma_ab, threshold)
    return -v, it, "reflected"

  a_plus_b = a + b
  a_plus_1 = a + 1.0

  digamma_a += 1.0 / a

  prefactor = (a_plus_1 / a_plus_b) ** 3

  sum_numer = (digamma_ab - digamma_a) * prefactor
  sum_denom = prefactor
  summand = prefactor * z * a_plus_b / a_plus_1

  k = 1.0
  digamma_ab += 1.0 / a_plus_b
  digamma_a += 1.0 / a_plus_1

  while abs(summand) > threshold:
    sum_numer += (digamma_ab - digamma_a) * summand
    sum_denom += summand
    summand *= (1.0 + a_plus_b / k) * (1.0 + k) / (1.0 + a_plus_1 / k)
    digamma_ab += 1.0 / (a_plus_b + k)
    digamma_a += 1.0 / (a_plus_1 + k)
    k += 1.0
    summand *= z / k
    if k > 1e5:
      raise NotConverged("inc_beta_dda")
  return (inc_beta(a, b, z) * (math.log(z) + sum_numer / sum_denom), int(k),
          "series")


def inc_beta_ddb(a, b, z, digamma_b=None, digamma_ab=None, threshold=1e-10):
  """prim/fun/inc_beta_ddb.hpp, develop 5252d51d47."""
  if digamma_b is None:
    digamma_b = digamma(b)
  if digamma_ab is None:
    digamma_ab = digamma(a + b)

  if b > a:
    if ((0.1 < z <= 0.75 and b > 500) or (0.01 < z <= 0.1 and b > 2500)
        or (0.001 < z <= 0.01 and b > 1e5)):
      v, it, _ = inc_beta_dda(b, a, 1.0 - z, digamma_b, digamma_ab, threshold)
      return -v, it, "reflected"
  if ((z > 0.75 and a < 500) or (z > 0.9 and a < 2500)
      or (z > 0.99 and a < 1e5) or z > 0.999):
    v, it, _ = inc_beta_dda(b, a, 1.0 - z, digamma_b, digamma_ab, threshold)
    return -v, it, "reflected"

  a_plus_b = a + b
  a_plus_1 = a + 1.0

  prefactor = (a_plus_1 / a_plus_b) ** 3

  sum_numer = digamma_ab * prefactor
  sum_denom = prefactor
  summand = prefactor * z * a_plus_b / a_plus_1

  k = 1.0
  digamma_ab += 1.0 / a_plus_b

  while abs(summand) > threshold:
    sum_numer += digamma_ab * summand
    sum_denom += summand
    summand *= (1.0 + a_plus_b / k) * (1.0 + k) / (1.0 + a_plus_1 / k)
    digamma_ab += 1.0 / (a_plus_b + k)
    k += 1.0
    summand *= z / k
    if k > 1e5:
      raise NotConverged("inc_beta_ddb")
  return (inc_beta(a, b, z) * (math.log1p(-z) - digamma_b
                               + sum_numer / sum_denom), int(k), "series")


# ── Path B: grad_reg_inc_beta through grad_2F1 ───────────────────────────
def clog(x):
  """std::log semantics: log(0) = -inf, it does not raise."""
  if x == 0.0:
    return -INFTY
  return math.log(x)


def csign(x):
  """stan::math::sign: -1, 0 or 1."""
  return (x > 0) - (x < 0)


def grad_2F1_impl_ab(a1, a2, b1, z, precision=1e-14, max_steps=1000000):
  """prim/fun/grad_2F1.hpp internal::grad_2F1_impl_ab<true, true, true>,
  develop 5252d51d47.

  The stopping rule is `inner_diff > precision || k < min_steps`, where
  `inner_diff` is the largest **absolute** increment of the three gradient
  sums at the current step.

  Returns (dF_a1, dF_a2, dF_b1, iterations).
  """
  grad = [0.0, 0.0, 0.0]
  if z == 0.0:
    return grad[0], grad[1], grad[2], 0

  log_g_old = [-INFTY, -INFTY, -INFTY]
  log_t_old = 0.0
  log_t_new = 0.0
  sign_z = csign(z)
  log_z = math.log(abs(z))

  log_t_new_sign = 1
  log_t_old_sign = 1
  log_g_old_sign = [1, 1, 1]

  sign_zk = sign_z
  k = 0
  min_steps = 5
  inner_diff = 1.0
  g_current = [0.0, 0.0, 0.0]
  params = (a1, a2, b1)

  while (inner_diff > precision or k < min_steps) and k < max_steps:
    p = (a1 + k) * (a2 + k) / ((b1 + k) * (1.0 + k))
    if p == 0.0:
      return grad[0], grad[1], grad[2], k
    log_t_new += clog(abs(p)) + log_z
    log_t_new_sign = csign(p) * log_t_new_sign

    for i in range(3):
      # inv(a1 + k), inv(a2 + k), inv(-(b1 + k))
      recip = 1.0 / (params[i] + k) if i < 2 else 1.0 / (-(params[i] + k))
      term = (log_g_old_sign[i] * log_t_old_sign
              * cexp(log_g_old[i] - log_t_old) + recip)
      log_g_old[i] = log_t_new + clog(abs(term))
      log_g_old_sign[i] = csign(term) * log_t_new_sign
      g_current[i] = log_g_old_sign[i] * cexp(log_g_old[i]) * sign_zk
      grad[i] += g_current[i]

    inner_diff = max(abs(g) for g in g_current)

    log_t_old = log_t_new
    log_t_old_sign = log_t_new_sign
    sign_zk *= sign_z
    k += 1

  if k > max_steps:
    raise NotConverged("grad_2F1")
  return grad[0], grad[1], grad[2], k


def grad_reg_inc_beta(a, b, z, digamma_a=None, digamma_b=None,
                      digamma_sum=None, beta_ab=None):
  """prim/fun/grad_reg_inc_beta.hpp, develop 5252d51d47.

  Returns (g1, g2, iterations, branch). For 0 < z < 1 the 2F1 series
  converges, so `grad_2F1_impl` never applies the Euler transform here.
  `dFz` is computed by the C++ through `hypergeometric_2F1` and then not
  used; the transcription omits it.

  When `C = z^a (1-z)^b / a` underflows to 0 the C++ skips `grad_2F1` and
  the `C * dF` term is dropped. The branch name records this.
  """
  if digamma_a is None:
    digamma_a = digamma(a)
  if digamma_b is None:
    digamma_b = digamma(b)
  if digamma_sum is None:
    digamma_sum = digamma(a + b)
  if beta_ab is None:
    beta_ab = beta_fn(a, b)

  c1 = math.log(z)
  c2 = math.log1p(-z)
  c3 = beta_ab * inc_beta(a, b, z)
  C = cexp(a * c1 + b * c2) / a
  dF1 = dF2 = dF3 = 0.0
  iters = 0
  branch = "2F1"
  if C != 0.0:
    dF1, dF2, dF3, iters = grad_2F1_impl_ab(a + b, 1.0, a + 1.0, z)
  else:
    branch = "C_underflow"

  dBda = (c1 - 1.0 / a) * c3 + C * (dF1 + dF3)
  dBdb = c2 * c3 + C * dF1

  # beta_ab underflows to 0 for a + b above about 1100 (a near b). The C++
  # then divides by 0: 0/0 is NaN, x/0 is inf. Python raises instead, so
  # the division goes through cdiv.
  if beta_ab == 0.0:
    branch = "betaAB_underflow"
  g1 = cdiv(dBda - c3 * (digamma_a - digamma_sum), beta_ab)
  g2 = cdiv(dBdb - c3 * (digamma_b - digamma_sum), beta_ab)
  return g1, g2, iters, branch


def cdiv(x, y):
  """IEEE division: x/0 is +-inf, 0/0 and nan/0 are nan."""
  if y != 0.0:
    return x / y
  if x == 0.0 or math.isnan(x):
    return math.nan
  return math.copysign(math.inf, x) * math.copysign(1.0, y)


@functools.lru_cache(maxsize=1 << 18)
def beta_fn(a, b):
  """beta(a, b), correctly rounded; +inf above the binary64 range, 0 below."""
  v = mp.beta(mpf(a), mpf(b))
  if v > mpf('1.7976931348623157e308'):
    return math.inf
  return float(v)


# ── inc_beta_ddz, generic template form ──────────────────────────────────
def inc_beta_ddz_generic(a, b, z):
  """The templated form in prim/fun/inc_beta_ddz.hpp.

  The `double` specialization calls boost::math::ibeta_derivative instead.
  This form is what fvar and var arguments use.
  """
  return cexp((b - 1.0) * math.log1p(-z) + (a - 1.0) * math.log(z)
                  + math.lgamma(a + b) - math.lgamma(a) - math.lgamma(b))
