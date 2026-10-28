#!/usr/bin/env python3
# gamma_lccdf d/dalpha references (ref, ref_quad) and the Python probe of
# the deep upper tail.  Supports the 266-cell grid and the all-modes claim.
"""gamma_lccdf d/dalpha in the deep upper tail (beta = 1, y = r * alpha).

Copies the double algorithm of internal::log_q_gamma_cf (develop and
stable-inc-gamma are identical there) and the var-path gradient
dQ_da / exp(log Q).  Classifies exp(log Q) as normal, subnormal or zero.
Scores a candidate repair: the forward-mode derivative of the log-space
continued fraction, d/da log Q = log z - digamma(a) - f'/f.
Reference: mpmath at 50 digits, d/da of log Q by mp.diff.
"""
import math
import sys
import mpmath as mp

mp.mp.dps = 50
EPS = sys.float_info.epsilon
PREC = 1.49012e-12
LOG_DBL_MIN = math.log(sys.float_info.min)  # -708.4
LOG_TRUE_MIN = math.log(5e-324)             # -744.4

ALPHAS = [5, 10, 20, 30, 50, 70, 100, 150, 200, 300, 500, 700, 850, 1000]
RATIOS = [1.05, 1.1, 1.2, 1.3, 1.5, 2, 3, 4, 5, 6, 8, 10, 15, 20, 30, 50, 75,
          100]


def digamma(a):
  return float(mp.digamma(a))


def log_q_cf_dual(a, z, max_steps=250, stop_on_deriv=False, info=None):
  """Returns (log Q, d/da log Q) of the CF in double, a as a dual number.

  stop_on_deriv: also require that the tangent of log(delta) is below
  PREC relative to the tangent of log f.  info: dict, receives 'steps'.
  """
  # values (v) and tangents (d) with respect to a
  lpre_v = a * math.log(z) - z - math.lgamma(a)
  lpre_d = math.log(z) - digamma(a)
  b0_v, b0_d = z + 1.0 - a, -1.0
  if abs(b0_v) >= EPS:
    c_v, c_d = b0_v, b0_d
  else:
    c_v, c_d = EPS, 0.0
  d_v, d_d = 0.0, 0.0
  f_v, f_d = c_v, c_d
  for i in range(1, max_steps + 1):
    an_v, an_d = -i * (i - a), float(i)
    b_v, b_d = b0_v + 2.0 * i, b0_d
    # D = b + an * D
    d_v, d_d = b_v + an_v * d_v, b_d + an_d * d_v + an_v * d_d
    if abs(d_v) < EPS:
      d_v, d_d = EPS, 0.0
    # C = b + an / C
    c_v, c_d = (b_v + an_v / c_v,
                b_d + (an_d * c_v - an_v * c_d) / (c_v * c_v))
    if abs(c_v) < EPS:
      c_v, c_d = EPS, 0.0
    # D = 1 / D
    d_v, d_d = 1.0 / d_v, -d_d / (d_v * d_v)
    dl_v, dl_d = c_v * d_v, c_d * d_v + c_v * d_d
    f_v, f_d = f_v * dl_v, f_d * dl_v + f_v * dl_d
    if abs(dl_v - 1.0) < PREC:
      if not stop_on_deriv:
        break
      if abs(dl_d / dl_v) <= PREC * max(1.0, abs(f_d / f_v)):
        break
  if info is not None:
    info["steps"] = i
  return lpre_v - math.log(f_v), lpre_d - f_d / f_v


def ref_gammainc(a, z):
  a, z = mp.mpf(a), mp.mpf(z)
  lq = lambda s: mp.log(mp.gammainc(s, z, mp.inf, regularized=True))
  return lq(a), mp.diff(lq, a)


def ref_quad(a, z):
  """Second route: integral over t = z + s, scaled to 1 at s = 0 (z > a - 1).

  log Q = (a-1) log z - z - lgamma(a) + log I0,
  d/da log Q = log z - digamma(a) + I1 / I0,
  I0 = int w, I1 = int w log(1 + s/z), w = exp((a-1) log(1 + s/z) - s).
  """
  a, z = mp.mpf(a), mp.mpf(z)
  lw = lambda s: (a - 1) * mp.log1p(s / z) - s
  dec = 1 / (1 - (a - 1) / z)
  nodes = [mp.mpf(0)] + [dec * 2**k for k in range(-4, 12)] + [mp.inf]
  i0 = mp.quad(lambda s: mp.exp(lw(s)), nodes)
  i1 = mp.quad(lambda s: mp.exp(lw(s)) * mp.log1p(s / z), nodes)
  lq = (a - 1) * mp.log(z) - z - mp.loggamma(a) + mp.log(i0)
  return lq, mp.log(z) - mp.digamma(a) + i1 / i0


def ref(a, z):
  """gammainc route where it converges; the routes are compared there."""
  qa, qd = ref_quad(a, z)
  try:
    ga, gd = ref_gammainc(a, z)
  except mp.libmp.libhyper.NoConvergence:
    return qa, qd, float("nan")
  agree = max(abs((qa - ga) / ga), abs((qd - gd) / gd))
  return ga, gd, float(agree)


BIG_ALPHAS = [1e4, 1e5, 1e6]
BIG_RATIOS = [1.01, 1.03, 1.1, 1.13, 1.2, 1.5, 2]


def row(a, z, r):
  lq_ref, dlq_ref, agree = ref(a, z)
  i1, i2 = {}, {}
  lq, dlq = log_q_cf_dual(float(a), float(z), info=i1)
  lq2, dlq2 = log_q_cf_dual(float(a), float(z), max_steps=100000,
                            stop_on_deriv=True, info=i2)
  qv = math.exp(lq)
  cls = ("zero" if qv == 0.0 else
         "subnorm" if qv < sys.float_info.min else "normal")
  e_v = float(abs((lq - lq_ref) / lq_ref))
  e_c = float(abs((dlq - dlq_ref) / dlq_ref))
  e_v2 = float(abs((lq2 - lq_ref) / lq_ref))
  e_c2 = float(abs((dlq2 - dlq_ref) / dlq_ref))
  print(f"{a:7g} {r:5g} {float(lq_ref):10.2f} {cls:>7} {float(dlq_ref):9.4g} "
        f"{e_v:8.1e} {e_c:8.1e} {i1['steps']:4d} {e_v2:8.1e} {e_c2:8.1e} "
        f"{i2['steps']:5d} {agree:8.1e}")
  return e_c, e_c2, e_v, e_v2, (0.0 if math.isnan(agree) else agree)


def main():
  print("cur = log_q_gamma_cf as is (value stop, 250 steps); "
        "new = value and derivative stop, no step cap")
  print(f"{'alpha':>7} {'r':>5} {'logQ_ref':>10} {'Q dbl':>7} {'dlogQ':>9} "
        f"{'cur_val':>8} {'cur_der':>8} {'stp':>4} {'new_val':>8} "
        f"{'new_der':>8} {'stp':>5} {'routes':>8}")
  w = [0.0] * 5
  for alphas, ratios in ((ALPHAS, RATIOS), (BIG_ALPHAS, BIG_RATIOS)):
    for a in alphas:
      for r in ratios:
        z = r * a
        if not z > a + 1.0:
          continue
        w = [max(x, y) for x, y in zip(w, row(a, z, r))]
  print(f"worst rel err: cur_der {w[0]:.2e}  new_der {w[1]:.2e}  "
        f"cur_val {w[2]:.2e}  new_val {w[3]:.2e}  route_disagree {w[4]:.2e}")


if __name__ == "__main__":
  main()
