#!/usr/bin/env python3
"""mpmath references for the large-shape spot check (probe_sib.cpp).

Routes per point:
  A  value and analytic gradient (psi differences, or the analytic sum
     sum_k pmf_k * d log pmf_k for the cdfs) at DPS_A digits
  B  value at DPS_B digits and mp.diff of the value at DPS_B digits
  C  an independent second formula where one exists:
       LKJ, LKJC  normalizing constant from the partial-correlation
                  integral, -(2 lbeta(1/2, eta + 1/2) + lbeta(1/2, eta))
       YS*        1 - sum_{k=1}^{n} pmf(k) for the cdf
       BNB* a/ab  the direct tail sum sum_{k>n} pmf(k) for the ccdf

Precision: at shape 1e15 each loggamma is about 3e16 and the result can be
O(1), so 17 digits go to cancellation. DPS_A = 90 and DPS_B = 140. The
beta_neg_binomial ccdf is 1 - cdf with the ccdf as small as 1e-900; there
the working precision is raised by -log10(pmf(n + 1)) + 10 digits.

Usage: ref_sib.py <points file> <refs file> [TAG ...]
Appends a line per point that is not yet in <refs file>.
"""

import math
import sys

from mpmath import mp, mpf

DPS_A = 90
DPS_B = 140
NSTR = 50
mp.dps = DPS_A

LOG_PI = None


def _f(x):
  return mpf(x)


def lgam(x):
  return mp.loggamma(x)


def psi(x):
  return mp.digamma(x)


def lbeta(a, b):
  return lgam(a) + lgam(b) - lgam(a + b)


# ── neg_binomial family ───────────────────────────────────────────────────

def lp_NB(n, a, b):
  n, a, b = _f(n), _f(a), _f(b)
  return (lgam(n + a) - lgam(a) - lgam(n + 1) - a * mp.log1p(1 / b)
          - n * mp.log1p(b))


def grad_NB(n, a, b):
  n, a, b = _f(n), _f(a), _f(b)
  return [psi(n + a) - psi(a) - mp.log1p(1 / b),
          a / (b * (1 + b)) - n / (1 + b)]


def lp_NB2(n, mu, phi):
  n, mu, phi = _f(n), _f(mu), _f(phi)
  return (lgam(n + phi) - lgam(phi) - lgam(n + 1) + n * mp.log(mu)
          + phi * mp.log(phi) - (n + phi) * mp.log(mu + phi))


def grad_NB2(n, mu, phi):
  n, mu, phi = _f(n), _f(mu), _f(phi)
  d_mu = n / mu - (n + phi) / (mu + phi)
  d_phi = (psi(n + phi) - psi(phi) + mp.log(phi) - mp.log(mu + phi)
           + (mu - n) / (mu + phi))
  return [d_mu, d_phi]


def lp_NB2L(n, eta, phi):
  return lp_NB2(n, mp.exp(_f(eta)), phi)


def grad_NB2L(n, eta, phi):
  mu = mp.exp(_f(eta))
  d_mu, d_phi = grad_NB2(n, mu, phi)
  return [mu * d_mu, d_phi]


def _glm_split(args):
  k = (len(args) - 3) // 2
  ys = [_f(v) for v in args[:k]]
  xs = [_f(v) for v in args[k:2 * k]]
  alpha, beta, phi = (_f(v) for v in args[2 * k:])
  return ys, xs, alpha, beta, phi


def lp_GLM(*args):
  ys, xs, alpha, beta, phi = _glm_split(args)
  return mp.fsum(lp_NB2L(y, x * beta + alpha, phi) for y, x in zip(ys, xs))


def grad_GLM(*args):
  ys, xs, alpha, beta, phi = _glm_split(args)
  da = db = dp = mpf(0)
  for y, x in zip(ys, xs):
    d_eta, d_phi = grad_NB2L(y, x * beta + alpha, phi)
    da += d_eta
    db += x * d_eta
    dp += d_phi
  return [da, db, dp]


# ── beta_neg_binomial, dirichlet_multinomial ─────────────────────────────

def lp_BNB(n, r, alpha, beta):
  n, r, alpha, beta = _f(n), _f(r), _f(alpha), _f(beta)
  return (-lgam(n + 1) + lbeta(n + r, alpha + beta) + lgam(n + beta)
          - lbeta(r, alpha) - lgam(beta))


def grad_BNB(n, r, alpha, beta):
  n, r, alpha, beta = _f(n), _f(r), _f(alpha), _f(beta)
  dg_all = psi(n + r + alpha + beta)
  dg_ra = psi(r + alpha)
  dg_ab = psi(alpha + beta)
  d_r = psi(n + r) - dg_all - (psi(r) - dg_ra)
  d_alpha = dg_ab - dg_all - (psi(alpha) - dg_ra)
  d_beta = dg_ab - dg_all + psi(n + beta) - psi(beta)
  return [d_r, d_alpha, d_beta]


def lp_DM(n1, n2, n3, a1, a2, a3):
  ns = [_f(n1), _f(n2), _f(n3)]
  a = [_f(a1), _f(a2), _f(a3)]
  N = mp.fsum(ns)
  A = mp.fsum(a)
  out = mp.log(N) - mp.fsum(mp.log(ni) for ni in ns if ni > 0)
  out += lbeta(A, N) - mp.fsum(lbeta(ai, ni) for ai, ni in zip(a, ns)
                               if ni > 0)
  return out


def grad_DM(n1, n2, n3, a1, a2, a3):
  ns = [_f(n1), _f(n2), _f(n3)]
  a = [_f(a1), _f(a2), _f(a3)]
  N = mp.fsum(ns)
  A = mp.fsum(a)
  tail = psi(A) - psi(A + N)
  return [((psi(ai + ni) - psi(ai)) if ni > 0 else mpf(0)) + tail
          for ai, ni in zip(a, ns)]


# ── LKJ, K = 3 ────────────────────────────────────────────────────────────
# do_lkj_constant for K = 3, eta != 1:
#   2 lgamma(eta + 1) - lgamma(eta + 1/2) - lgamma(eta) - 3/2 log(pi)

def lkj_const(eta):
  eta = _f(eta)
  return (2 * lgam(eta + 1) - lgam(eta + mpf(1) / 2) - lgam(eta)
          - mpf(3) / 2 * mp.log(mp.pi))


def lkj_const_c(eta):
  """Route C: 1 / integral of det^(eta - 1) over 3x3 correlation matrices.

  With r12 = u, r13 = v, r23 = u v + w sqrt((1 - u^2)(1 - v^2)) the
  determinant is (1 - u^2)(1 - v^2)(1 - w^2), and the integral factors into
  B(1/2, eta + 1/2)^2 B(1/2, eta).
  """
  eta = _f(eta)
  half = mpf(1) / 2
  return -(2 * lbeta(half, eta + half) + lbeta(half, eta))


def dlkj_const(eta):
  eta = _f(eta)
  return 2 * psi(eta + 1) - psi(eta + mpf(1) / 2) - psi(eta)


def _logdet3(r12, r13, r23):
  r12, r13, r23 = _f(r12), _f(r13), _f(r23)
  return mp.log(1 - r12**2 - r13**2 - r23**2 + 2 * r12 * r13 * r23)


def lp_LKJ(r12, r13, r23, eta):
  return lkj_const(eta) + (_f(eta) - 1) * _logdet3(r12, r13, r23)


def grad_LKJ(r12, r13, r23, eta):
  return [dlkj_const(eta) + _logdet3(r12, r13, r23)]


def _lkjc_logs(l11, l22):
  return mp.log(_f(l11)), mp.log(_f(l22))


def lp_LKJC(l10, l11, l20, l21, l22, eta):
  g1, g2 = _lkjc_logs(l11, l22)
  eta = _f(eta)
  return lkj_const(eta) + (1 + 2 * eta - 2) * g1 + (2 * eta - 2) * g2


def grad_LKJC(l10, l11, l20, l21, l22, eta):
  g1, g2 = _lkjc_logs(l11, l22)
  return [dlkj_const(eta) + 2 * g1 + 2 * g2]


# ── beta_binomial cdfs: finite sums of the pmf ───────────────────────────

def _bb_terms(N, a, b):
  """pmf_k and d log pmf_k / d(a, b) for k = 0..N, by recurrence."""
  N = int(N)
  a, b = _f(a), _f(b)
  lp0 = lbeta(a, N + b) - lbeta(a, b)
  g_a = psi(a + b) - psi(a + b + N)
  g_b = psi(N + b) - psi(N + a + b) - psi(b) + psi(a + b)
  p = mp.exp(lp0)
  out = []
  for k in range(N + 1):
    out.append((p, g_a, g_b))
    if k < N:
      p = p * (N - k) * (k + a) / ((k + 1) * (N - k - 1 + b))
      g_a += 1 / (k + a)
      g_b -= 1 / (N - k - 1 + b)
  return out


def _bb_parts(n, N, a, b):
  terms = _bb_terms(N, a, b)
  n = int(n)
  lo = terms[:n + 1]
  hi = terms[n + 1:]
  cdf = mp.fsum(t[0] for t in lo)
  ccdf = mp.fsum(t[0] for t in hi)
  dcdf = [mp.fsum(t[0] * t[1] for t in lo), mp.fsum(t[0] * t[2] for t in lo)]
  dccdf = [mp.fsum(t[0] * t[1] for t in hi), mp.fsum(t[0] * t[2] for t in hi)]
  return cdf, ccdf, dcdf, dccdf


def lp_BBC(n, N, a, b):
  return _bb_parts(n, N, a, b)[0]


def lp_BBLC(n, N, a, b):
  return mp.log(_bb_parts(n, N, a, b)[0])


def lp_BBLCC(n, N, a, b):
  return mp.log(_bb_parts(n, N, a, b)[1])


def grad_BBC(n, N, a, b):
  return _bb_parts(n, N, a, b)[2]


def grad_BBLC(n, N, a, b):
  cdf, _, dcdf, _ = _bb_parts(n, N, a, b)
  return [d / cdf for d in dcdf]


def grad_BBLCC(n, N, a, b):
  _, ccdf, _, dccdf = _bb_parts(n, N, a, b)
  return [d / ccdf for d in dccdf]


# ── yule_simon cdfs: closed form ──────────────────────────────────────────
# P(Y > n) = Gamma(n + 1) Gamma(a + 1) / Gamma(n + a + 1) = n B(n, a + 1)

def _ys_ccdf(n, a):
  n, a = _f(n), _f(a)
  return mp.exp(lgam(n + 1) + lgam(a + 1) - lgam(n + a + 1))


def _ys_q(n, a):
  n, a = _f(n), _f(a)
  return psi(a + 1) - psi(n + a + 1)


def lp_YSC(n, a):
  return 1 - _ys_ccdf(n, a)


def lp_YSLC(n, a):
  return mp.log1p(-_ys_ccdf(n, a))


def lp_YSLCC(n, a):
  n, a = _f(n), _f(a)
  return lgam(n + 1) + lgam(a + 1) - lgam(n + a + 1)


def grad_YSC(n, a):
  return [-_ys_ccdf(n, a) * _ys_q(n, a)]


def grad_YSLC(n, a):
  c = _ys_ccdf(n, a)
  return [-c * _ys_q(n, a) / (1 - c)]


def grad_YSLCC(n, a):
  return [_ys_q(n, a)]


def ys_cdf_sum(n, a):
  """Route C: sum_{k=1}^{n} a B(k, a + 1)."""
  a = _f(a)
  return mp.fsum(a * mp.beta(k, a + 1) for k in range(1, int(n) + 1))


# ── beta_neg_binomial cdfs: finite sum, 1 - cdf for the ccdf ─────────────

def _bnb_log_pmf_next(n, r, a, b):
  """(log pmf(n + 1), cdf(n)), low precision, for the precision estimate."""
  with mp.workdps(30):
    r, a, b = _f(r), _f(a), _f(b)
    lp = lbeta(r, a + b) - lbeta(r, a)
    cdf = mpf(0)
    for k in range(int(n) + 1):
      cdf += mp.exp(lp)
      lp += (mp.log(k + r) + mp.log(k + b) - mp.log(k + 1)
             - mp.log(k + r + a + b))
    return lp, cdf


def _bnb_extra_dps(n, r, a, b):
  """Digits lost in 1 - cdf: about -log10(ccdf), and ccdf >= pmf(n + 1)."""
  lp, cdf = _bnb_log_pmf_next(n, r, a, b)
  if cdf < mpf("0.5"):
    return 10
  return max(0, int(math.ceil(float(-lp / mp.log(10))))) + 10


def _bnb_parts_at(n, r, a, b):
  r, a, b = _f(r), _f(a), _f(b)
  s = r + a + b
  lp0 = lbeta(r, a + b) - lbeta(r, a)
  g = [psi(r + a) - psi(s), psi(a + b) - psi(s) - psi(a) + psi(r + a),
       psi(a + b) - psi(s)]
  p = mp.exp(lp0)
  cdf = mpf(0)
  dcdf = [mpf(0)] * 3
  for k in range(int(n) + 1):
    cdf += p
    dcdf = [d + p * gi for d, gi in zip(dcdf, g)]
    p = p * (k + r) * (k + b) / ((k + 1) * (k + s))
    g = [g[0] + 1 / (k + r) - 1 / (k + s), g[1] - 1 / (k + s),
         g[2] + 1 / (k + b) - 1 / (k + s)]
  return cdf, 1 - cdf, dcdf


def _bnb_parts(n, r, a, b):
  extra = _bnb_extra_dps(n, r, a, b)
  with mp.workdps(mp.dps + extra):
    cdf, ccdf, dcdf = _bnb_parts_at(n, r, a, b)
    # keep every digit until the log or the ratio is formed
    return cdf, ccdf, dcdf, extra


def lp_BNBC(n, r, a, b):
  return +_bnb_parts(n, r, a, b)[0]


def lp_BNBLC(n, r, a, b):
  cdf, _, _, extra = _bnb_parts(n, r, a, b)
  with mp.workdps(mp.dps + extra):
    return mp.log(cdf)


def lp_BNBLCC(n, r, a, b):
  _, ccdf, _, extra = _bnb_parts(n, r, a, b)
  with mp.workdps(mp.dps + extra):
    return mp.log(ccdf)


def grad_BNBC(n, r, a, b):
  return [+d for d in _bnb_parts(n, r, a, b)[2]]


def grad_BNBLC(n, r, a, b):
  cdf, _, dcdf, extra = _bnb_parts(n, r, a, b)
  with mp.workdps(mp.dps + extra):
    return [d / cdf for d in dcdf]


def grad_BNBLCC(n, r, a, b):
  _, ccdf, dcdf, extra = _bnb_parts(n, r, a, b)
  with mp.workdps(mp.dps + extra):
    return [-d / ccdf for d in dcdf]


def bnb_ccdf_tail(n, r, a, b, kmax=20000):
  """Route C: sum_{k > n} pmf(k) directly; None if it does not converge."""
  r, a, b = _f(r), _f(a), _f(b)
  s = r + a + b
  lp = lbeta(n + r, a + b) + lgam(n + b) - lgam(n + 1) - lbeta(r, a) - lgam(b)
  p = mp.exp(lp)
  total = mpf(0)
  eps = mpf(10) ** (-mp.dps)
  for k in range(int(n), int(n) + kmax):
    p = p * (k + r) * (k + b) / ((k + 1) * (k + s))
    total += p
    ratio = (k + 1 + r) * (k + 1 + b) / ((k + 2) * (k + 1 + s))
    if ratio < mpf("0.9") and p < eps * total:
      return total
  return None


# ── Tag table ─────────────────────────────────────────────────────────────
# tag: (lp, grad, arg index differentiated per gradient)

SPEC = {
    "NB": (lp_NB, grad_NB, [1, 2]),
    "NB2": (lp_NB2, grad_NB2, [1, 2]),
    "NB2L": (lp_NB2L, grad_NB2L, [1, 2]),
    "GLM1": (lp_GLM, grad_GLM, [2, 3, 4]),
    "GLM5": (lp_GLM, grad_GLM, [10, 11, 12]),
    "BNB": (lp_BNB, grad_BNB, [1, 2, 3]),
    "DM": (lp_DM, grad_DM, [3, 4, 5]),
    "LKJ": (lp_LKJ, grad_LKJ, [3]),
    "LKJC": (lp_LKJC, grad_LKJC, [5]),
    "BBC": (lp_BBC, grad_BBC, [2, 3]),
    "BBLC": (lp_BBLC, grad_BBLC, [2, 3]),
    "BBLCC": (lp_BBLCC, grad_BBLCC, [2, 3]),
    "YSC": (lp_YSC, grad_YSC, [1]),
    "YSLC": (lp_YSLC, grad_YSLC, [1]),
    "YSLCC": (lp_YSLCC, grad_YSLCC, [1]),
    "BNBC": (lp_BNBC, grad_BNBC, [1, 2, 3]),
    "BNBLC": (lp_BNBLC, grad_BNBLC, [1, 2, 3]),
    "BNBLCC": (lp_BNBLCC, grad_BNBLCC, [1, 2, 3]),
}


def route_c(tag, args):
  """Independent value by a second formula, or None."""
  if tag == "LKJ":
    return lkj_const_c(args[3]) + (_f(args[3]) - 1) * _logdet3(*args[:3])
  if tag == "LKJC":
    g1, g2 = _lkjc_logs(args[1], args[4])
    eta = _f(args[5])
    return lkj_const_c(eta) + (2 * eta - 1) * g1 + (2 * eta - 2) * g2
  if tag in ("YSC", "YSLC", "YSLCC"):
    cdf = ys_cdf_sum(*args)
    if tag == "YSC":
      return cdf
    if tag == "YSLC":
      return mp.log(cdf)
    return None  # 1 - cdf would cancel; the closed form is route A
  if tag in ("BNBC", "BNBLC", "BNBLCC") and float(args[2]) >= 100:
    extra = _bnb_extra_dps(*args)
    with mp.workdps(mp.dps + extra):
      tail = bnb_ccdf_tail(*args)
      if tail is None:
        return None
      if tag == "BNBC":
        return 1 - tail
      if tag == "BNBLC":
        return mp.log1p(-tail)
      return mp.log(tail)
  return None


def route_a(tag, args):
  lp, grad, _ = SPEC[tag]
  with mp.workdps(DPS_A):
    v = lp(*args)
    g = grad(*args)
    return +v, [+x for x in g]


def route_b(tag, args):
  lp, _, index = SPEC[tag]
  with mp.workdps(DPS_B):
    base = [mpf(x) for x in args]
    v = lp(*base)
    gs = []
    for i in index:

      def f(t, i=i):
        a = list(base)
        a[i] = t
        return lp(*a)

      gs.append(mp.diff(f, base[i]))
    return v, gs


def parse_points(path):
  out = []
  with open(path) as fh:
    for line in fh:
      parts = line.split()
      if not parts or parts[0].startswith("#"):
        continue
      out.append((parts[0], parts[1], parts[2:]))
  return out


def key_of(tag, label, hexes):
  return " ".join([tag, label] + hexes)


def load_keys(path):
  keys = set()
  try:
    with open(path) as fh:
      for line in fh:
        keys.add(line.split(" => ")[0])
  except FileNotFoundError:
    pass
  return keys


def s(x):
  return "none" if x is None else mp.nstr(x, NSTR)


def main():
  if len(sys.argv) < 3:
    print(__doc__)
    return 2
  points = parse_points(sys.argv[1])
  tags = set(sys.argv[3:])
  done = load_keys(sys.argv[2])
  n_new = 0
  with open(sys.argv[2], "a") as out:
    for tag, label, hexes in points:
      if tags and tag not in tags:
        continue
      key = key_of(tag, label, hexes)
      if key in done:
        continue
      args = [float.fromhex(h) for h in hexes]
      va, ga = route_a(tag, args)
      vb, gb = route_b(tag, args)
      vc = route_c(tag, args)
      out.write(key + " => " + " ".join(s(x) for x in [va] + ga) + " | "
                + " ".join(s(x) for x in [vb] + gb) + " | " + s(vc) + "\n")
      out.flush()
      done.add(key)
      n_new += 1
  print(f"{n_new} new reference lines")
  return 0


if __name__ == "__main__":
  sys.exit(main())
