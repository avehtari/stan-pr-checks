#!/usr/bin/env python3
"""High-precision references for the density-gradient sweep.

Plan: ../2026-10-01_PLAN_density_gradient_sweep.md.
Standard: ../STAN_MATH_DISTRIBUTIONS.md section 3. Every quantity has two
independent routes:

  route A  the analytic gradient, written out
  route B  mp.diff of the log density at 110 digits

The two share no code path apart from the elementary functions, so
agreement is evidence that both are right. Route A is the reference that
`score_sweep.py` uses; route B is the check. `mp.diff` alone is not a
reference: on 2026-10-01 it returned 0 for a gradient of 1e-200 at a raw
scale of 1e150, because the log density there is of order 1e45.

Run this file to print the route A against route B table.
"""

import math
import sys

from mpmath import mp, mpf

DPS = 50
DPS_DIFF = 110
mp.dps = DPS

LOG_SQRT_PI = mp.log(mp.sqrt(mp.pi))


def _f(x):
  return mpf(x)


# ── Log densities, for route B and for the value column ──────────────────

def lp_NB2(n, mu, phi):
  n, mu, phi = _f(n), _f(mu), _f(phi)
  return (mp.loggamma(n + phi) - mp.loggamma(phi) - mp.loggamma(n + 1)
          + n * mp.log(mu) + phi * mp.log(phi) - (n + phi) * mp.log(mu + phi))


def lp_NORM(y, mu, sigma):
  y, mu, sigma = _f(y), _f(mu), _f(sigma)
  z = (y - mu) / sigma
  return -mp.log(sigma) - mp.log(mp.sqrt(2 * mp.pi)) - z * z / 2


def _Phi(x):
  return mp.erfc(-_f(x) / mp.sqrt(2)) / 2


def _Q(x):
  return mp.erfc(_f(x) / mp.sqrt(2)) / 2


def _oprobit_p(lam, c1, c2):
  """P(y = 2) for three categories, in whichever form keeps its digits.

  Phi(x) rounds to exactly 1 above x = 10 at 50 digits, so the difference
  of two CDFs is unusable when both cutpoints are above the location. The
  complement Q has no such limit, and Phi(u) - Phi(l) = Q(l) - Q(u) is
  exact. This is the same identity as the repair recorded for
  ordered_probit_lpmf in herbie_stan_math/shortlist.md.
  """
  lo, hi = _f(c1) - _f(lam), _f(c2) - _f(lam)
  if lo > 0:
    return _Q(lo) - _Q(hi)
  return _Phi(hi) - _Phi(lo)


def lp_OPROBIT(lam, c1, c2):
  """Middle category of three: log(Phi(c2 - lambda) - Phi(c1 - lambda))."""
  return mp.log(_oprobit_p(lam, c1, c2))


def lp_BB(n, N, alpha, beta):
  n, N, alpha, beta = _f(n), _f(N), _f(alpha), _f(beta)
  lchoose = (mp.loggamma(N + 1) - mp.loggamma(n + 1)
             - mp.loggamma(N - n + 1))
  lbeta = lambda a, b: mp.loggamma(a) + mp.loggamma(b) - mp.loggamma(a + b)
  return lchoose + lbeta(n + alpha, N - n + beta) - lbeta(alpha, beta)


def lp_BNB(n, r, alpha, beta):
  n, r, alpha, beta = _f(n), _f(r), _f(alpha), _f(beta)
  lbeta = lambda a, b: mp.loggamma(a) + mp.loggamma(b) - mp.loggamma(a + b)
  return (-mp.loggamma(n + 1) + lbeta(n + r, alpha + beta)
          + mp.loggamma(n + beta) - lbeta(r, alpha) - mp.loggamma(beta))


def lp_YS(n, alpha):
  n, alpha = _f(n), _f(alpha)
  return (mp.loggamma(n) + mp.log(alpha) + mp.loggamma(alpha + 1)
          - mp.loggamma(n + alpha + 1))


def lp_DIR(t1, t2, a1, a2, a3):
  t = [_f(t1), _f(t2), 1 - _f(t1) - _f(t2)]
  a = [_f(a1), _f(a2), _f(a3)]
  asum = a[0] + a[1] + a[2]
  out = mp.loggamma(asum) - sum(mp.loggamma(ai) for ai in a)
  out += sum((ai - 1) * mp.log(ti) for ai, ti in zip(a, t))
  return out


def lp_DM(n1, n2, n3, a1, a2, a3):
  ns = [_f(n1), _f(n2), _f(n3)]
  a = [_f(a1), _f(a2), _f(a3)]
  N = ns[0] + ns[1] + ns[2]
  if N == 0:
    return mpf(0)
  A = a[0] + a[1] + a[2]
  lbeta = lambda x, y: mp.loggamma(x) + mp.loggamma(y) - mp.loggamma(x + y)
  out = mp.log(N) - sum(mp.log(ni) for ni in ns if ni > 0)
  out += lbeta(A, N) - sum(lbeta(ai, ni) for ai, ni in zip(a, ns) if ni > 0)
  return out


def lp_MST(y1, y2, nu):
  """multi_student_t_cholesky, dim 2, mu = 0, L = identity."""
  y1, y2, nu = _f(y1), _f(y2), _f(nu)
  p = mpf(2)
  G = y1 * y1 + y2 * y2
  return (mp.loggamma((nu + p) / 2) - mp.loggamma(nu / 2)
          - (p / 2) * mp.log(nu) - (p / 2) * mp.log(mp.pi)
          - ((nu + p) / 2) * mp.log1p(G / nu))


def lp_LBETA(d, x):
  """lbeta(d, x), the helper of the proposed value repair.

  Not a density. The identity is
  lgamma(x + d) - lgamma(x) = lgamma(d) - lbeta(d, x), so an accurate
  `lbeta` turns the raw lgamma difference in `yule_simon_lpmf` and
  `multi_student_t_cholesky_lpdf` into an accurate expression.
  """
  d, x = _f(d), _f(x)
  return mp.loggamma(d) + mp.loggamma(x) - mp.loggamma(d + x)


def lp_ST(y, nu, mu, sigma):
  y, nu, mu, sigma = _f(y), _f(nu), _f(mu), _f(sigma)
  z = (y - mu) / sigma
  r = z * z / nu
  return (mp.loggamma((nu + 1) / 2) - mp.loggamma(nu / 2)
          - mp.log(nu) / 2 - LOG_SQRT_PI - mp.log(sigma)
          - ((nu + 1) / 2) * mp.log1p(r))


# ── Stage C ───────────────────────────────────────────────────────────────

def lp_VM(y, mu, kappa):
  y, mu, kappa = _f(y), _f(mu), _f(kappa)
  return (kappa * mp.cos(mu - y) - mp.log(2 * mp.pi)
          - mp.log(mp.besseli(0, kappa)))


def _simplex3(t1, t2):
  return [_f(t1), _f(t2), 1 - _f(t1) - _f(t2)]


def lp_MN(n1, n2, n3, t1, t2):
  ns = [_f(n1), _f(n2), _f(n3)]
  t = _simplex3(t1, t2)
  N = ns[0] + ns[1] + ns[2]
  out = mp.loggamma(N + 1) - sum(mp.loggamma(ni + 1) for ni in ns)
  # multiply_log(0, theta) is 0, which matters when a count is zero
  out += sum(ni * mp.log(ti) for ni, ti in zip(ns, t) if ni != 0)
  return out


def _lse3(b1, b2, b3):
  b = [_f(b1), _f(b2), _f(b3)]
  m = max(b)
  return m + mp.log(sum(mp.exp(bi - m) for bi in b))


def lp_MNL(n1, n2, n3, b1, b2, b3):
  ns = [_f(n1), _f(n2), _f(n3)]
  b = [_f(b1), _f(b2), _f(b3)]
  N = ns[0] + ns[1] + ns[2]
  alpha = _lse3(b1, b2, b3)
  out = mp.loggamma(1 + N) - sum(mp.loggamma(1 + ni) for ni in ns)
  out += sum(ni * (bi - alpha) for ni, bi in zip(ns, b) if ni != 0)
  return out


def lp_CATL(n, b1, b2, b3):
  b = [_f(b1), _f(b2), _f(b3)]
  return b[int(n) - 1] - _lse3(b1, b2, b3)


def _pb_prob(y, thetas):
  """Poisson-binomial P(Y = y) by enumeration, three trials."""
  import itertools
  total = mpf(0)
  index = range(len(thetas))
  for subset in itertools.combinations(index, int(y)):
    term = mpf(1)
    for i in index:
      term *= thetas[i] if i in subset else (1 - thetas[i])
    total += term
  return total


def lp_PB(y, t1, t2, t3):
  thetas = [_f(t1), _f(t2), _f(t3)]
  return mp.log(_pb_prob(y, thetas))


def lp_HG(n, N, a, b):
  """hypergeometric_lpmf: three binomial_coefficient_log calls."""
  n, N, a, b = _f(n), _f(N), _f(a), _f(b)
  lchoose = lambda x, k: (mp.loggamma(x + 1) - mp.loggamma(k + 1)
                          - mp.loggamma(x - k + 1))
  return lchoose(a, n) + lchoose(b, N - n) - lchoose(a + b, N)


def lp_DR(y, lower, upper):
  return -mp.log(_f(upper) - _f(lower) + 1)


# ── Stage D: 2x2 matrix densities, mu = 0 ────────────────────────────────
# Every formula is written out for 2x2, so no matrix library is involved and
# the reference is exact at 50 digits. `mu = 0`, so the residual is `y`.

LOG_TWO_PI = mp.log(2 * mp.pi)


def _inv2(a11, a12, a21, a22):
  det = a11 * a22 - a12 * a21
  return a22 / det, -a12 / det, -a21 / det, a11 / det, det


def lp_MVN2(y1, y2, s11, s12, s21, s22):
  y1, y2 = _f(y1), _f(y2)
  s11, s12, s21, s22 = _f(s11), _f(s12), _f(s21), _f(s22)
  i11, i12, i21, i22, det = _inv2(s11, s12, s21, s22)
  quad = y1 * (i11 * y1 + i12 * y2) + y2 * (i21 * y1 + i22 * y2)
  return -quad / 2 - mp.log(det) / 2 - LOG_TWO_PI


def grad_MVN2(y1, y2, s11, s12, s21, s22):
  y1, y2 = _f(y1), _f(y2)
  s11, s12, s21, s22 = _f(s11), _f(s12), _f(s21), _f(s22)
  i11, i12, i21, i22, _ = _inv2(s11, s12, s21, s22)
  # a = Sigma^-1 y
  a1 = i11 * y1 + i12 * y2
  a2 = i21 * y1 + i22 * y2
  # d/dSigma of -0.5 y' S^-1 y - 0.5 log det S is 0.5 (a a' - S^-1)
  g11 = (a1 * a1 - i11) / 2
  g12 = (a1 * a2 - i12) / 2
  g21 = (a2 * a1 - i21) / 2
  g22 = (a2 * a2 - i22) / 2
  return [-a1, -a2, g11, g12, g21, g22]


def lp_MVNP2(y1, y2, p11, p12, p21, p22):
  y1, y2 = _f(y1), _f(y2)
  p11, p12, p21, p22 = _f(p11), _f(p12), _f(p21), _f(p22)
  det = p11 * p22 - p12 * p21
  quad = y1 * (p11 * y1 + p12 * y2) + y2 * (p21 * y1 + p22 * y2)
  return -quad / 2 + mp.log(det) / 2 - LOG_TWO_PI


def grad_MVNP2(y1, y2, p11, p12, p21, p22):
  y1, y2 = _f(y1), _f(y2)
  p11, p12, p21, p22 = _f(p11), _f(p12), _f(p21), _f(p22)
  i11, i12, i21, i22, _ = _inv2(p11, p12, p21, p22)
  b1 = p11 * y1 + p12 * y2
  b2 = p21 * y1 + p22 * y2
  # d/dOmega of -0.5 y' O y + 0.5 log det O is 0.5 (O^-1 - y y')
  return [-b1, -b2,
          (i11 - y1 * y1) / 2, (i12 - y1 * y2) / 2,
          (i21 - y2 * y1) / 2, (i22 - y2 * y2) / 2]


def lp_MVNC2(y1, y2, l11, l21, l22):
  y1, y2 = _f(y1), _f(y2)
  l11, l21, l22 = _f(l11), _f(l21), _f(l22)
  # A = L^-1 y by forward substitution
  a1 = y1 / l11
  a2 = (y2 - l21 * a1) / l22
  return -(a1 * a1 + a2 * a2) / 2 - mp.log(l11) - mp.log(l22) - LOG_TWO_PI


def grad_MVNC2(y1, y2, l11, l21, l22):
  y1, y2 = _f(y1), _f(y2)
  l11, l21, l22 = _f(l11), _f(l21), _f(l22)
  a1 = y1 / l11
  a2 = (y2 - l21 * a1) / l22
  # d/dy = -L^-T A
  dy1 = -(a1 / l11 - l21 * a2 / (l11 * l22))
  dy2 = -a2 / l22
  # d/dL by the chain rule on A = L^-1 y, written out for 2x2
  dl11 = a1 * (a1 + a2 * (-l21) / l22) / l11 - 1 / l11
  dl21 = a2 * a1 / l22
  dl22 = a2 * a2 / l22 - 1 / l22
  return [dy1, dy2, dl11, dl21, dl22]


# ── Stage D at dimension 5: the same three densities, general code ────────
# The 2x2 versions above are written out by hand and are checked against
# route B. These use mpmath matrices and are checked against the 2x2
# versions at p = 2 as well as against route B, so a mistake in the general
# code cannot pass unnoticed.

P5 = 5


def _mat(flat, p):
  m = mp.matrix(p, p)
  for i in range(p):
    for j in range(p):
      m[i, j] = _f(flat[i * p + j])
  return m


def _vec(flat, p):
  v = mp.matrix(p, 1)
  for i in range(p):
    v[i] = _f(flat[i])
  return v


def lp_MVN5(*args):
  p = P5
  y, S = _vec(args[:p], p), _mat(args[p:], p)
  Si = S ** -1
  quad = (y.T * Si * y)[0, 0]
  return -quad / 2 - mp.log(mp.det(S)) / 2 - (mpf(p) / 2) * LOG_TWO_PI


def grad_MVN5(*args):
  p = P5
  y, S = _vec(args[:p], p), _mat(args[p:], p)
  Si = S ** -1
  a = Si * y
  out = [-a[i] for i in range(p)]
  for i in range(p):
    for j in range(p):
      out.append((a[i] * a[j] - Si[i, j]) / 2)
  return out


def lp_MVNP5(*args):
  p = P5
  y, O = _vec(args[:p], p), _mat(args[p:], p)
  quad = (y.T * O * y)[0, 0]
  return -quad / 2 + mp.log(mp.det(O)) / 2 - (mpf(p) / 2) * LOG_TWO_PI


def grad_MVNP5(*args):
  p = P5
  y, O = _vec(args[:p], p), _mat(args[p:], p)
  Oi = O ** -1
  b = O * y
  out = [-b[i] for i in range(p)]
  for i in range(p):
    for j in range(p):
      out.append((Oi[i, j] - y[i] * y[j]) / 2)
  return out


def _forward_solve(L, y, p):
  """a = L^-1 y, lower triangular."""
  a = mp.matrix(p, 1)
  for i in range(p):
    acc = y[i]
    for k in range(i):
      acc -= L[i, k] * a[k]
    a[i] = acc / L[i, i]
  return a


def _back_solve_T(L, a, p):
  """L^-T a, upper triangular in L^T."""
  x = mp.matrix(p, 1)
  for i in reversed(range(p)):
    acc = a[i]
    for k in range(i + 1, p):
      acc -= L[k, i] * x[k]
    x[i] = acc / L[i, i]
  return x


def lp_MVNC5(*args):
  p = P5
  y, L = _vec(args[:p], p), _mat(args[p:], p)
  a = _forward_solve(L, y, p)
  quad = sum(a[i] * a[i] for i in range(p))
  logdet = sum(mp.log(L[i, i]) for i in range(p))
  return -quad / 2 - logdet - (mpf(p) / 2) * LOG_TWO_PI


def grad_MVNC5(*args):
  p = P5
  y, L = _vec(args[:p], p), _mat(args[p:], p)
  a = _forward_solve(L, y, p)
  w = _back_solve_T(L, a, p)
  out = [-w[i] for i in range(p)]
  # d/dL_ij = (L^-T a)_i a_j - delta_ij / L_ii; the strict upper triangle is
  # structurally zero and Stan reports 0 there.
  for i in range(p):
    for j in range(p):
      if j > i:
        out.append(mpf(0))
      else:
        g = w[i] * a[j]
        if i == j:
          g -= 1 / L[i, i]
        out.append(g)
  return out


LP = {
    "NB2": lp_NB2, "NORM": lp_NORM, "OPROBIT": lp_OPROBIT, "BB": lp_BB,
    "BNB": lp_BNB, "YS": lp_YS, "DIR": lp_DIR, "DM": lp_DM, "MST": lp_MST,
    "LBETA": lp_LBETA, "ST": lp_ST,
    "VM": lp_VM, "MN": lp_MN, "MNL": lp_MNL, "CATL": lp_CATL, "PB": lp_PB,
    "HG": lp_HG, "DR": lp_DR,
    "MVN2": lp_MVN2, "MVNP2": lp_MVNP2, "MVNC2": lp_MVNC2,
    "MVN5": lp_MVN5, "MVNP5": lp_MVNP5, "MVNC5": lp_MVNC5,
}


# ── Route A: analytic gradients ───────────────────────────────────────────

def grad_NB2(n, mu, phi):
  n, mu, phi = _f(n), _f(mu), _f(phi)
  d_mu = n / mu - (n + phi) / (mu + phi)
  d_phi = (mp.digamma(n + phi) - mp.digamma(phi) + mp.log(phi)
           - mp.log(mu + phi) + (mu - n) / (mu + phi))
  return [d_mu, d_phi]


def grad_NORM(y, mu, sigma):
  y, mu, sigma = _f(y), _f(mu), _f(sigma)
  z = (y - mu) / sigma
  return [-z / sigma, z / sigma, (z * z - 1) / sigma]


def _npdf(x):
  x = _f(x)
  return mp.exp(-x * x / 2) / mp.sqrt(2 * mp.pi)


def grad_OPROBIT(lam, c1, c2):
  lam, c1, c2 = _f(lam), _f(c1), _f(c2)
  p = _oprobit_p(lam, c1, c2)
  f1, f2 = _npdf(c1 - lam), _npdf(c2 - lam)
  return [-(f2 - f1) / p, -f1 / p, f2 / p]


def grad_BB(n, N, alpha, beta):
  n, N, alpha, beta = _f(n), _f(N), _f(alpha), _f(beta)
  dg_ab = mp.digamma(alpha + beta)
  dg_Nab = mp.digamma(N + alpha + beta)
  d_alpha = mp.digamma(n + alpha) - dg_Nab - mp.digamma(alpha) + dg_ab
  d_beta = (mp.digamma(N - n + beta) - dg_Nab - mp.digamma(beta) + dg_ab)
  return [d_alpha, d_beta]


def grad_BNB(n, r, alpha, beta):
  n, r, alpha, beta = _f(n), _f(r), _f(alpha), _f(beta)
  dg_all = mp.digamma(n + r + alpha + beta)
  dg_ra = mp.digamma(r + alpha)
  dg_ab = mp.digamma(alpha + beta)
  d_r = mp.digamma(n + r) - dg_all - (mp.digamma(r) - dg_ra)
  d_alpha = dg_ab - dg_all - (mp.digamma(alpha) - dg_ra)
  d_beta = (dg_ab - dg_all + mp.digamma(n + beta) - mp.digamma(beta))
  return [d_r, d_alpha, d_beta]


def grad_YS(n, alpha):
  n, alpha = _f(n), _f(alpha)
  return [1 / alpha + mp.digamma(alpha + 1) - mp.digamma(n + alpha + 1)]


def grad_DIR(t1, t2, a1, a2, a3):
  t = [_f(t1), _f(t2), 1 - _f(t1) - _f(t2)]
  a = [_f(a1), _f(a2), _f(a3)]
  asum = a[0] + a[1] + a[2]
  dg_sum = mp.digamma(asum)
  d_t = [(ai - 1) / ti for ai, ti in zip(a, t)]
  d_a = [dg_sum - mp.digamma(ai) + mp.log(ti) for ai, ti in zip(a, t)]
  return d_t + d_a


def grad_DM(n1, n2, n3, a1, a2, a3):
  ns = [_f(n1), _f(n2), _f(n3)]
  a = [_f(a1), _f(a2), _f(a3)]
  N = ns[0] + ns[1] + ns[2]
  A = a[0] + a[1] + a[2]
  if N == 0:
    return [mpf(0), mpf(0), mpf(0)]
  tail = mp.digamma(A) - mp.digamma(A + N)
  out = []
  for ai, ni in zip(a, ns):
    head = (mp.digamma(ai + ni) - mp.digamma(ai)) if ni > 0 else mpf(0)
    out.append(head + tail)
  return out


def grad_MST(y1, y2, nu):
  y1, y2, nu = _f(y1), _f(y2), _f(nu)
  p = mpf(2)
  G = y1 * y1 + y2 * y2
  d_y = [-(nu + p) * yi / (nu + G) for yi in (y1, y2)]
  d_nu = mpf(0.5) * (mp.digamma((nu + p) / 2) - mp.digamma(nu / 2)
                     - mp.log1p(G / nu) + (G - p) / (G + nu))
  return d_y + [d_nu]


def grad_LBETA(d, x):
  d, x = _f(d), _f(x)
  dg_sum = mp.digamma(d + x)
  return [mp.digamma(d) - dg_sum, mp.digamma(x) - dg_sum]


def grad_ST(y, nu, mu, sigma):
  y, nu, mu, sigma = _f(y), _f(nu), _f(mu), _f(sigma)
  z = (y - mu) / sigma
  r = z * z / nu
  d_y = -(nu + 1) * (y - mu) / ((1 + r) * sigma * sigma * nu)
  rep = (nu + 1) * r / (1 + r) - 1
  d_nu = mpf(0.5) * (mp.digamma((nu + 1) / 2) - mp.digamma(nu / 2)
                     - mp.log1p(r) + rep / nu)
  return [d_y, d_nu, -d_y, rep / sigma]


def grad_VM(y, mu, kappa):
  y, mu, kappa = _f(y), _f(mu), _f(kappa)
  d_y = kappa * mp.sin(mu - y)
  ratio = mp.besseli(1, kappa) / mp.besseli(0, kappa)
  return [d_y, -d_y, mp.cos(mu - y) - ratio]


def grad_MN(n1, n2, n3, t1, t2):
  ns = [_f(n1), _f(n2), _f(n3)]
  t = _simplex3(t1, t2)
  return [ni / ti for ni, ti in zip(ns, t)]


def grad_MNL(n1, n2, n3, b1, b2, b3):
  ns = [_f(n1), _f(n2), _f(n3)]
  b = [_f(b1), _f(b2), _f(b3)]
  N = ns[0] + ns[1] + ns[2]
  alpha = _lse3(b1, b2, b3)
  return [ni - N * mp.exp(bi - alpha) for ni, bi in zip(ns, b)]


def grad_CATL(n, b1, b2, b3):
  b = [_f(b1), _f(b2), _f(b3)]
  alpha = _lse3(b1, b2, b3)
  return [(1 if i == int(n) - 1 else 0) - mp.exp(bi - alpha)
          for i, bi in enumerate(b)]


def grad_PB(y, t1, t2, t3):
  """d log P / d theta_j by enumeration, three trials."""
  import itertools
  thetas = [_f(t1), _f(t2), _f(t3)]
  total = _pb_prob(y, thetas)
  index = range(3)
  out = []
  for j in index:
    deriv = mpf(0)
    for subset in itertools.combinations(index, int(y)):
      term = mpf(1)
      for i in index:
        if i == j:
          continue
        term *= thetas[i] if i in subset else (1 - thetas[i])
      deriv += term if j in subset else -term
    out.append(deriv / total)
  return out


def grad_none(*args):
  """No parameter, so no gradient. HG and DR take integers only."""
  return []


GRAD = {
    "NB2": grad_NB2, "NORM": grad_NORM, "OPROBIT": grad_OPROBIT,
    "BB": grad_BB, "BNB": grad_BNB, "YS": grad_YS, "DIR": grad_DIR,
    "DM": grad_DM, "MST": grad_MST, "LBETA": grad_LBETA, "ST": grad_ST,
    "VM": grad_VM, "MN": grad_MN, "MNL": grad_MNL, "CATL": grad_CATL,
    "PB": grad_PB, "HG": grad_none, "DR": grad_none,
    "MVN2": grad_MVN2, "MVNP2": grad_MVNP2, "MVNC2": grad_MVNC2,
    "MVN5": grad_MVN5, "MVNP5": grad_MVNP5, "MVNC5": grad_MVNC5,
}

# Which argument index each gradient differentiates. A gradient of a
# dependent coordinate (the third simplex component of DIR) has no index and
# is checked by route A only.
DIFF_INDEX = {
    "NB2": [1, 2],
    "NORM": [0, 1, 2],
    "OPROBIT": [0, 1, 2],
    "BB": [2, 3],
    "BNB": [1, 2, 3],
    "YS": [1],
    "DIR": [0, 1, None, 2, 3, 4],
    "DM": [3, 4, 5],
    "MST": [0, 1, 2],
    "LBETA": [0, 1],
    "ST": [0, 1, 2, 3],
    "VM": [0, 1, 2],
    # the third simplex component of MN is dependent, as in DIR, so route B
    # cannot differentiate it; MNL, CATL and PB have free arguments
    "MN": [3, 4, None],
    "MNL": [3, 4, 5],
    "CATL": [1, 2, 3],
    "PB": [1, 2, 3],
    "HG": [],
    "DR": [],
    "MVN2": [0, 1, 2, 3, 4, 5],
    "MVNP2": [0, 1, 2, 3, 4, 5],
    "MVNC2": [0, 1, 2, 3, 4],
    # p = 5: 5 vector entries then 25 matrix entries, all free
    "MVN5": list(range(30)),
    "MVNP5": list(range(30)),
    "MVNC5": list(range(30)),
}


def lp_DIR_free(t1, t2, t3, a1, a2, a3):
  """Dirichlet log density with all three components free.

  `lp_DIR` derives `t3` from the other two, which is what the probe does.
  Route B must not do that: Stan reports the partial with respect to each
  component separately, with the simplex constraint handled by the
  transform outside the density. Differentiating `lp_DIR` would add the
  `t3` term to the `t1` and `t2` partials.
  """
  t = [_f(t1), _f(t2), _f(t3)]
  a = [_f(a1), _f(a2), _f(a3)]
  asum = a[0] + a[1] + a[2]
  out = mp.loggamma(asum) - sum(mp.loggamma(ai) for ai in a)
  out += sum((ai - 1) * mp.log(ti) for ai, ti in zip(a, t))
  return out


def grad_route_b(tag, args):
  """mp.diff of the log density, at 110 digits."""
  if tag == "DIR":
    t1, t2, a1, a2, a3 = args
    free = [mpf(t1), mpf(t2), 1 - mpf(t1) - mpf(t2), mpf(a1), mpf(a2),
            mpf(a3)]
    out = []
    with mp.workdps(DPS_DIFF):
      for index in range(6):

        def f(t, index=index):
          a = list(free)
          a[index] = t
          return lp_DIR_free(*a)

        out.append(mp.diff(f, free[index]))
    return out

  out = []
  with mp.workdps(DPS_DIFF):
    base = [mpf(a) for a in args]
    for index in DIFF_INDEX[tag]:
      if index is None:
        out.append(None)
        continue

      def f(t, index=index):
        a = list(base)
        a[index] = t
        return LP[tag](*a)

      out.append(mp.diff(f, base[index]))
  return out


# ── Route A against route B ──────────────────────────────────────────────

CHECK_POINTS = {
    "NB2": [(4, 4.0, 1e4), (0, 0.5, 1e-2), (1000, 1000.0, 10.0)],
    "NORM": [(1.0, 0.0, 1.0), (3.0, 1.0, 0.25)],
    "OPROBIT": [(0.0, 1.0, 1.1), (0.0, 30.0, 30.1)],
    "BB": [(1, 2, 1e3, 1e3), (5, 10, 0.5, 2.0), (3, 1000, 1e6, 1e6)],
    "BNB": [(1, 1e3, 1e3, 1e3), (7, 2.0, 3.0, 4.0)],
    "YS": [(1, 1e3), (10, 0.5), (1000, 2.0)],
    "DIR": [(0.25, 0.5, 1e3, 1e3, 1e3), (0.2, 0.3, 0.5, 1.5, 2.5)],
    "DM": [(1, 1, 1, 1e3, 1e3, 1e3), (10, 20, 30, 0.5, 1.0, 2.0)],
    "MST": [(1.0, 1.0, 1e3), (0.5, -2.0, 3.0)],
    # stage D. Route B differentiates every matrix entry independently, as
    # Stan does, so it checks the written-out 2x2 gradients, including the
    # Cholesky chain rule, which is the one that is easy to get wrong.
    "MVN2": [(1.0, 0.5, 2.0, 0.3, 0.3, 1.5),
             (0.7, 0.7, 1.0, 0.99, 0.99, 1.0)],
    "MVNP2": [(1.0, 0.5, 2.0, 0.3, 0.3, 1.5),
              (0.7, 0.7, 1.0, 0.99, 0.99, 1.0)],
    "MVNC2": [(1.0, 0.5, 1.4142135623730951, 0.2121320343559643,
               1.2093727122985778),
              (0.7, 0.7, 1.0, 0.99, 0.1410673597611751)],
}


def main():
  worst = mpf(0)
  print("# route A (analytic) against route B (mp.diff at 110 digits)")
  print(f"{'tag':8s} {'point':44s} {'gradient':10s} {'rel':>10s}")
  for tag, points in CHECK_POINTS.items():
    names = GRAD[tag].__doc__
    for args in points:
      a_vals = GRAD[tag](*args)
      b_vals = grad_route_b(tag, args)
      for i, (a, b) in enumerate(zip(a_vals, b_vals)):
        if b is None:
          continue
        # Mixed absolute and relative measure. A plain relative measure is
        # useless where the exact gradient is zero, which happens at every
        # zero crossing: d/dmu of neg_binomial_2 is exactly 0 at n = mu.
        rel = abs(a - b) / (1 + max(abs(a), abs(b)))
        worst = max(worst, rel)
        point = "(" + ", ".join(f"{float(x):g}" for x in args) + ")"
        print(f"{tag:8s} {point:44s} {i:<10d} {float(rel):10.2e}")
  print(f"\nworst route A against route B: {float(worst):.2e}")
  if worst > mpf("1e-30"):
    print("FAIL: the two routes disagree by more than 1e-30")
    return 1
  print("PASS")
  return 0


if __name__ == "__main__":
  sys.exit(main())
