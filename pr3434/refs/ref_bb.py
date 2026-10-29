#!/usr/bin/env python3
"""High-precision references for beta_binomial_lpmf (stan-dev/math#3153).

Plan: ../2026-10-06_PLAN_fix_issue3153_beta_binomial_lpmf.md, step 5.1.
Standard: ../STAN_MATH_DISTRIBUTIONS.md section 3. Three routes:

  route A  value by mp.loggamma, gradients by closed-form psi differences,
           at DPS digits; this is the reference that score_bb.py uses
  route B  the value again at DPS_CHECK digits, and mp.diff of the value at
           DPS_CHECK digits for the gradients
  route C  the exact finite sums psi(x + k) - psi(x) = sum_{j<k} 1/(x + j)
           for integer k (only for N <= NMAX_SUM, it is O(N))

Route A loses digits to cancellation: at shape 1e18 each loggamma is about
4e19 and the value is O(1), so 20 of the DPS digits go; the tied-shape
log-concentration gradient is about 1e-35 at shape 1e18 against psi values
of about 41, so 37 digits go. DPS = 80 leaves more than 40 digits in both.
Run this file to print the route comparison table.

Derived quantities, all exact functions of the partials:

  ga = alpha * d/dalpha      gradient in log(alpha)
  gb = beta * d/dbeta        gradient in log(beta)
  gc = ga + gb               gradient in log(concentration) when the mean
                             alpha / (alpha + beta) is held fixed; for tied
                             shapes alpha = beta = exp(lc) / 2 this is the
                             d/dlc of the issue
"""

import sys

from mpmath import mp, mpf

DPS = 80
DPS_CHECK = 130
NMAX_SUM = 2000
mp.dps = DPS


def _f(x):
  return mpf(x)


def lchoose(N, n):
  N, n = _f(N), _f(n)
  return mp.loggamma(N + 1) - mp.loggamma(n + 1) - mp.loggamma(N - n + 1)


def _lbeta(a, b):
  return mp.loggamma(a) + mp.loggamma(b) - mp.loggamma(a + b)


def lp_bb(n, N, alpha, beta):
  """Normalized log pmf (propto = false), route A."""
  n, N, alpha, beta = _f(n), _f(N), _f(alpha), _f(beta)
  return (lchoose(N, n) + _lbeta(n + alpha, N - n + beta)
          - _lbeta(alpha, beta))


def grad_bb(n, N, alpha, beta):
  """d/dalpha and d/dbeta, route A (closed-form psi differences)."""
  n, N, alpha, beta = _f(n), _f(N), _f(alpha), _f(beta)
  ab = alpha + beta
  tail = mp.digamma(ab) - mp.digamma(ab + N)
  d_alpha = mp.digamma(alpha + n) - mp.digamma(alpha) + tail
  d_beta = mp.digamma(beta + N - n) - mp.digamma(beta) + tail
  return d_alpha, d_beta


def _psi_diff_sum(x, k):
  """psi(x + k) - psi(x) for integer k >= 0, as an exact finite sum."""
  return mp.fsum(1 / (x + j) for j in range(int(k)))


def grad_bb_sum(n, N, alpha, beta):
  """Route C: the same gradients from exact finite sums."""
  n, N = int(n), int(N)
  alpha, beta = _f(alpha), _f(beta)
  ab = alpha + beta
  tail = -_psi_diff_sum(ab, N)
  d_alpha = _psi_diff_sum(alpha, n) + tail
  d_beta = _psi_diff_sum(beta, N - n) + tail
  return d_alpha, d_beta


def grad_bb_diff(n, N, alpha, beta):
  """Route B: mp.diff of the value at DPS_CHECK digits."""
  with mp.workdps(DPS_CHECK):
    a0, b0 = mpf(alpha), mpf(beta)
    d_alpha = mp.diff(lambda t: lp_bb(n, N, t, b0), a0)
    d_beta = mp.diff(lambda t: lp_bb(n, N, a0, t), b0)
  return d_alpha, d_beta


def ref_point(n, N, alpha, beta):
  """Every reference quantity of one point, route A, as mpf."""
  alpha, beta = _f(alpha), _f(beta)
  value = lp_bb(n, N, alpha, beta)
  d_alpha, d_beta = grad_bb(n, N, alpha, beta)
  ga = alpha * d_alpha
  gb = beta * d_beta
  return {
      "lchoose": lchoose(N, n),
      "value": value,
      "d_alpha": d_alpha,
      "d_beta": d_beta,
      "ga": ga,
      "gb": gb,
      "gc": ga + gb,
  }


# ── Route comparison ─────────────────────────────────────────────────────

def _tied(lc):
  """alpha = beta = exp(lc) / 2 rounded to double, as Stan computes it."""
  import math
  s = math.exp(lc) / 2
  return s


CHECK_POINTS = [
    # the reporter's setup, tied shapes
    (57, 117, _tied(20.0), _tied(20.0)),
    (57, 117, _tied(32.0), _tied(32.0)),
    (57, 117, _tied(42.0), _tied(42.0)),
    (57, 117, 1e18, 1e18),
    # unequal shapes, ordinary and large
    (3, 10, 0.5, 2.0),
    (0, 10, 0.01, 7.0),
    (10, 10, 30.0, 1e15),
    (1, 1000, 2.5, 1e12),
    (400, 1000, 3e13, 7e12),
    (0, 117, 1e17, 3e17),
    (117, 117, 1e17, 3e17),
    # large N: route C is skipped
    (500000, 1000000, 1e12, 1e12),
    (3, 1000000, 1e6, 1e18),
]


def _rel(a, b):
  """Mixed measure: absolute below 1, relative above."""
  return abs(a - b) / (1 + max(abs(a), abs(b)))


def _rel_pure(a, b):
  scale = max(abs(a), abs(b))
  return abs(a - b) / scale if scale > 0 else abs(a - b)


def main():
  worst_value = mpf(0)
  worst_grad = mpf(0)
  print("# route A (DPS=%d) against route B (DPS=%d, mp.diff) and C (sums)"
        % (DPS, DPS_CHECK))
  print(f"{'point':40s} {'value A-B':>10s} {'grad A-B':>10s} "
        f"{'grad A-C':>10s} {'gc rel A-B':>10s}")
  for n, N, alpha, beta in CHECK_POINTS:
    va = lp_bb(n, N, alpha, beta)
    with mp.workdps(DPS_CHECK):
      vb = lp_bb(n, N, alpha, beta)
    dv = _rel(va, vb)
    worst_value = max(worst_value, dv)

    ga = grad_bb(n, N, alpha, beta)
    gb = grad_bb_diff(n, N, alpha, beta)
    dg = max(_rel_pure(x, y) for x, y in zip(ga, gb))
    worst_grad = max(worst_grad, dg)
    gc_a = _f(alpha) * ga[0] + _f(beta) * ga[1]
    gc_b = _f(alpha) * gb[0] + _f(beta) * gb[1]
    dgc = _rel_pure(gc_a, gc_b)
    worst_grad = max(worst_grad, dgc)

    if N <= NMAX_SUM:
      gcs = grad_bb_sum(n, N, alpha, beta)
      dc = max(_rel_pure(x, y) for x, y in zip(ga, gcs))
      worst_grad = max(worst_grad, dc)
      dc_s = f"{float(dc):10.2e}"
    else:
      dc_s = f"{'-':>10s}"
    point = f"({n}, {N}, {float(alpha):.4g}, {float(beta):.4g})"
    print(f"{point:40s} {float(dv):10.2e} {float(dg):10.2e} {dc_s} "
          f"{float(dgc):10.2e}")
    print(f"{'':40s} value {mp.nstr(va, 20)}  gc {mp.nstr(gc_a, 20)}")
  print(f"\nworst value A-B: {float(worst_value):.2e}")
  print(f"worst gradient (relative, A-B, A-C, gc A-B): {float(worst_grad):.2e}")
  if worst_value > mpf("1e-30") or worst_grad > mpf("1e-30"):
    print("FAIL: the routes disagree by more than 1e-30")
    return 1
  print("PASS")
  return 0


if __name__ == "__main__":
  sys.exit(main())
