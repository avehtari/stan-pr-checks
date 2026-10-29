#!/usr/bin/env python3
"""Generate the hex point files for the density-gradient sweep.

Plan: ../2026-10-01_PLAN_density_gradient_sweep.md.
Convention copied from ../inc_gamma_beta_work/gen_points.py: one point per
line, every number in hex, `# <region>` headers. Hex is exact in both
directions, so the C++ probe reads exactly the points the scorer uses.

Line format:

    <TAG> <hex> <hex> ...

The number of arguments is fixed per tag and is listed in `ARITY` below,
which `score_sweep.py` and `probe_sweep.cpp` both follow. Integer arguments
are written as hex doubles and cast back to int in C++; every count used
here is below 2^53, so the round trip is exact.

Usage:
  python3 gen_points.py --stage a > points_a.txt
  python3 gen_points.py --stage b --n 200 > points_b.txt
"""

import argparse
import math
import random
import sys

SEED = 20261001

# tag -> (argument names, gradient names). The probe prints the value first,
# then the gradients in this order.
ARITY = {
    "NB2": (["n", "mu", "phi"], ["d_mu", "d_phi"]),
    "NORM": (["y", "mu", "sigma"], ["d_y", "d_mu", "d_sigma"]),
    "OPROBIT": (["lambda", "c1", "c2"], ["d_lambda", "d_c1", "d_c2"]),
    "BB": (["n", "N", "alpha", "beta"], ["d_alpha", "d_beta"]),
    "BNB": (["n", "r", "alpha", "beta"], ["d_r", "d_alpha", "d_beta"]),
    "YS": (["n", "alpha"], ["d_alpha"]),
    "DIR": (["t1", "t2", "a1", "a2", "a3"],
            ["d_t1", "d_t2", "d_t3", "d_a1", "d_a2", "d_a3"]),
    "DM": (["n1", "n2", "n3", "a1", "a2", "a3"], ["d_a1", "d_a2", "d_a3"]),
    "MST": (["y1", "y2", "nu"], ["d_y1", "d_y2", "d_nu"]),
    # Not a density. `lbeta(d, x)` is the helper that the proposed repair of
    # the yule_simon and multi_student_t_cholesky values would use, through
    # the identity lgamma(x + d) - lgamma(x) = lgamma(d) - lbeta(d, x).
    # Measuring it decides whether that repair is sound.
    "LBETA": (["d", "x"], ["d_d", "d_x"]),
    # student_t_lpdf has the same `lgamma(x + d) - lgamma(x)` difference in
    # its value as multi_student_t_cholesky_lpdf. Its gradients were
    # measured by the Herbie prefilter (0.73 bits); its value was not.
    "ST": (["y", "nu", "mu", "sigma"], ["d_y", "d_nu", "d_mu", "d_sigma"]),
    # ── stage C: the scalar and simplex composed densities ──
    "VM": (["y", "mu", "kappa"], ["d_y", "d_mu", "d_kappa"]),
    "MN": (["n1", "n2", "n3", "t1", "t2"], ["d_t1", "d_t2", "d_t3"]),
    "MNL": (["n1", "n2", "n3", "b1", "b2", "b3"], ["d_b1", "d_b2", "d_b3"]),
    "CATL": (["n", "b1", "b2", "b3"], ["d_b1", "d_b2", "d_b3"]),
    "PB": (["y", "t1", "t2", "t3"], ["d_t1", "d_t2", "d_t3"]),
    # value only: every argument is an integer, so there is no gradient
    "HG": (["n", "N", "a", "b"], []),
    "DR": (["y", "lower", "upper"], []),
    # ── stage D: matrix densities, dimension 2, mu = 0 ──
    # The matrix entries travel in the point file as hex, so the reference
    # and the probe see the same matrix bit for bit. Stan treats the matrix
    # argument as independent vars, so all four adjoints are read, including
    # both off-diagonal ones.
    "MVN2": (["y1", "y2", "s11", "s12", "s21", "s22"],
             ["d_y1", "d_y2", "d_s11", "d_s12", "d_s21", "d_s22"]),
    "MVNP2": (["y1", "y2", "p11", "p12", "p21", "p22"],
              ["d_y1", "d_y2", "d_p11", "d_p12", "d_p21", "d_p22"]),
    "MVNC2": (["y1", "y2", "l11", "l21", "l22"],
              ["d_y1", "d_y2", "d_l11", "d_l21", "d_l22"]),
    # dimension 5: five vector entries then the full 5x5 matrix, row major.
    # The full matrix travels even for the Cholesky form, so one probe code
    # path serves all three and the zero adjoints of the strict upper
    # triangle are checked too.
    "MVN5": ([f"y{i}" for i in range(5)]
             + [f"s{i}{j}" for i in range(5) for j in range(5)],
             [f"d_y{i}" for i in range(5)]
             + [f"d_s{i}{j}" for i in range(5) for j in range(5)]),
    "MVNP5": ([f"y{i}" for i in range(5)]
              + [f"p{i}{j}" for i in range(5) for j in range(5)],
              [f"d_y{i}" for i in range(5)]
              + [f"d_p{i}{j}" for i in range(5) for j in range(5)]),
    "MVNC5": ([f"y{i}" for i in range(5)]
              + [f"l{i}{j}" for i in range(5) for j in range(5)],
              [f"d_y{i}" for i in range(5)]
              + [f"d_l{i}{j}" for i in range(5) for j in range(5)]),
}


def hx(x):
  return float(x).hex()


def emit(tag, args):
  names = ARITY[tag][0]
  assert len(args) == len(names), f"{tag}: {len(args)} args, want {len(names)}"
  print(tag, " ".join(hx(a) for a in args))


def log_uniform(rng, lo, hi):
  return math.exp(rng.uniform(math.log(lo), math.log(hi)))


def log_uniform_int(rng, lo, hi):
  """Integer, log-uniform in [lo, hi]. lo may be 0, which is then kept."""
  if lo <= 0:
    if rng.random() < 0.2:
      return 0
    lo = 1
  return int(round(log_uniform(rng, lo, hi)))


# ── Stage A: controls with known answers ─────────────────────────────────
# Each of these reproduces a measurement that is already recorded, so a
# harness that disagrees with the recorded value is broken.
# See ../2026-10-01_PLAN_density_gradient_sweep.md step 1.

def stage_a():
  print("# A1 neg_binomial_2 dphi, recorded rel 3.1e-08 at (4, 4, 1e4)")
  for phi in (1e0, 1e2, 1e4, 1e6, 1e8, 1e10):
    emit("NB2", [4, 4.0, phi])

  print("# A2 normal_lpdf dsigma at the zero crossing, recorded rel 9.93e-02")
  emit("NORM", [1.0000000000000002e-16, 0.0, 1e-16])
  emit("NORM", [1.0, 0.0, 1.0])
  emit("NORM", [1.0000000000000002, 0.0, 1.0])
  emit("NORM", [3.0, 0.0, 1.0])

  print("# A3 ordered_probit, recorded -inf when both cutpoints are 38.5 up")
  for lo in (0.0, 10.0, 30.0, 37.0, 38.5, 40.0, 50.0):
    emit("OPROBIT", [0.0, lo, lo + 0.1])


# ── Stage B: the digamma-difference densities ────────────────────────────
# The hypothesis is in the plan, section 3: digamma(x + k) - digamma(x)
# keeps log10(x / k) fewer digits than its inputs. B3 and B5 are the
# regions that test it; B1, B2 and B4 are the controls around them.

B_REGIONS = (
    ("B1 small shape, small counts", 1e-2, 1.0, 0, 10),
    ("B2 moderate shape and counts", 1.0, 1e2, 0, 100),
    ("B3 large shape, small counts", 1e3, 1e8, 0, 10),
    ("B4 moderate shape, large counts", 1.0, 1e3, 1000, 1000000),
    ("B5 large shape, large counts", 1e4, 1e8, 1000, 1000000),
)


def stage_b(n, rng):
  for name, s_lo, s_hi, k_lo, k_hi in B_REGIONS:
    print(f"# {name}")
    for _ in range(n):
      # beta_binomial: n of N trials, shapes alpha and beta
      N = max(1, log_uniform_int(rng, max(1, k_lo), k_hi))
      k = rng.randint(0, N)
      emit("BB", [k, N, log_uniform(rng, s_lo, s_hi),
                  log_uniform(rng, s_lo, s_hi)])
      # beta_neg_binomial: count n, shapes r, alpha, beta
      emit("BNB", [log_uniform_int(rng, k_lo, k_hi),
                   log_uniform(rng, s_lo, s_hi),
                   log_uniform(rng, max(s_lo, 1e-2), s_hi),
                   log_uniform(rng, max(s_lo, 1.01), s_hi)])
      # yule_simon: count n >= 1, shape alpha
      emit("YS", [max(1, log_uniform_int(rng, max(1, k_lo), k_hi)),
                  log_uniform(rng, s_lo, s_hi)])
      # dirichlet: simplex of 3, shapes alpha
      u = sorted((rng.random(), rng.random()))
      t1, t2 = max(u[0], 1e-12), max(u[1] - u[0], 1e-12)
      emit("DIR", [t1, t2,
                   log_uniform(rng, s_lo, s_hi),
                   log_uniform(rng, s_lo, s_hi),
                   log_uniform(rng, s_lo, s_hi)])
      # dirichlet_multinomial: counts of 3, shapes alpha
      emit("DM", [log_uniform_int(rng, k_lo, k_hi),
                  log_uniform_int(rng, k_lo, k_hi),
                  log_uniform_int(rng, k_lo, k_hi),
                  log_uniform(rng, s_lo, s_hi),
                  log_uniform(rng, s_lo, s_hi),
                  log_uniform(rng, s_lo, s_hi)])
      # multi_student_t_cholesky, dim 2, mu = 0, L = I: sweep nu and y
      emit("MST", [rng.uniform(-5, 5), rng.uniform(-5, 5),
                   log_uniform(rng, max(s_lo, 1e-2), s_hi)])

  print("# B8 student_t over nu decades, deterministic")
  for e in range(0, 15):
    emit("ST", [1.0, 10.0 ** e, 0.0, 1.0])
    emit("ST", [3.0, 10.0 ** e, 0.0, 2.0])
  for _ in range(n):
    emit("ST", [rng.uniform(-5, 5), log_uniform(rng, 1e-2, 1e14),
                rng.uniform(-2, 2), log_uniform(rng, 1e-2, 1e2)])

  print("# B7 lbeta(d, x) over x decades, deterministic")
  for d in (0.5, 1.0, 1.5, 2.0, 3.0):
    for e in range(0, 15):
      emit("LBETA", [d, 10.0 ** e])
  for _ in range(n):
    emit("LBETA", [log_uniform(rng, 1e-2, 1e2), log_uniform(rng, 1e-2, 1e14)])

  # The digamma difference is largest when the count is tiny next to the
  # shape. This region is deterministic so the worst case is always in the
  # file.
  print("# B6 shape decades with count 1, deterministic")
  for e in range(0, 13):
    s = 10.0 ** e
    emit("BB", [1, 2, s, s])
    emit("BNB", [1, s, s, max(s, 1.01)])
    emit("YS", [1, s])
    emit("DIR", [0.25, 0.5, s, s, s])
    emit("DM", [1, 1, 1, s, s, s])
    emit("MST", [1.0, 1.0, s])


# ── Stage C: the scalar and simplex composed densities ───────────────────
# `categorical_lpmf` is not here: it is `log(theta[n - 1])`, one log call, and
# `categorical_logit_glm_lpmf` shares the softmax kernel of `CATL`.

def simplex3(rng):
  u = sorted((rng.random(), rng.random()))
  return max(u[0], 1e-12), max(u[1] - u[0], 1e-12)


def stage_c(n, rng):
  print("# C1 ordinary")
  for _ in range(n):
    emit("VM", [rng.uniform(-math.pi, math.pi), rng.uniform(-math.pi, math.pi),
                log_uniform(rng, 1e-2, 1e2)])
    t1, t2 = simplex3(rng)
    emit("MN", [rng.randint(0, 100), rng.randint(0, 100), rng.randint(0, 100),
                t1, t2])
    emit("MNL", [rng.randint(0, 100), rng.randint(0, 100), rng.randint(0, 100),
                 rng.uniform(-5, 5), rng.uniform(-5, 5), rng.uniform(-5, 5)])
    emit("CATL", [rng.randint(1, 3), rng.uniform(-5, 5), rng.uniform(-5, 5),
                  rng.uniform(-5, 5)])
    emit("PB", [rng.randint(0, 3), rng.uniform(0.01, 0.99),
                rng.uniform(0.01, 0.99), rng.uniform(0.01, 0.99)])
    a, b = rng.randint(1, 100), rng.randint(1, 100)
    N = rng.randint(1, a + b)
    emit("HG", [rng.randint(max(0, N - b), min(N, a)), N, a, b])
    lower = rng.randint(-100, 100)
    upper = lower + rng.randint(0, 1000)
    emit("DR", [rng.randint(lower, upper), lower, upper])

  print("# C2 extremes")
  for _ in range(n):
    # kappa large: log_modified_bessel_first_kind is the question here
    emit("VM", [rng.uniform(-math.pi, math.pi), rng.uniform(-math.pi, math.pi),
                log_uniform(rng, 1e2, 1e8)])
    # counts to 1e6, so lgamma(sum + 1) - sum lgamma(n + 1) cancels
    t1, t2 = simplex3(rng)
    emit("MN", [log_uniform_int(rng, 1000, 1000000),
                log_uniform_int(rng, 1000, 1000000),
                log_uniform_int(rng, 1000, 1000000), t1, t2])
    emit("MNL", [log_uniform_int(rng, 1000, 1000000),
                 log_uniform_int(rng, 1000, 1000000),
                 log_uniform_int(rng, 1000, 1000000),
                 rng.uniform(-700, 700), rng.uniform(-700, 700),
                 rng.uniform(-700, 700)])
    # logits far apart, so log_sum_exp must not overflow
    emit("CATL", [rng.randint(1, 3), rng.uniform(-700, 700),
                  rng.uniform(-700, 700), rng.uniform(-700, 700)])
    # probabilities at the ends, where the recursion can multiply by zero
    emit("PB", [rng.randint(0, 3), log_uniform(rng, 1e-12, 1e-2),
                1 - log_uniform(rng, 1e-12, 1e-2), rng.uniform(0.01, 0.99)])
    a = log_uniform_int(rng, 100, 100000)
    b = log_uniform_int(rng, 100, 100000)
    N = log_uniform_int(rng, 1, min(a + b, 100000))
    lo, hi = max(0, N - b), min(N, a)
    emit("HG", [rng.randint(lo, hi) if hi >= lo else lo, N, a, b])
    lower = -10 ** rng.randint(0, 8)
    upper = lower + 10 ** rng.randint(0, 8)
    emit("DR", [rng.randint(lower, upper), lower, upper])

  print("# C3 kappa decades, deterministic")
  for e in range(-2, 9):
    kappa = 10.0 ** e
    emit("VM", [0.0, 0.0, kappa])
    emit("VM", [1.0, 0.0, kappa])


# ── Stage D: matrix densities over the condition number ──────────────────
# The sweep variable is the condition number of the matrix argument, not a
# parameter magnitude. An ill-conditioned covariance is an ordinary
# modelling situation, not an extreme input, so this is the one stage where
# a finding would not be an extreme-value finding.
#
# Sigma = R(theta) diag(1, 1/cond) R(theta)^T, built in binary64 here and
# shipped as hex, so the reference and the probe share the same matrix. The
# Cholesky factor and the precision matrix come from the same Sigma through
# exact 2x2 formulas, so the three parameterizations describe the same
# density and their log densities must agree.

def spd2(cond, theta):
  """A 2x2 SPD matrix with eigenvalues 1 and 1/cond, rotated by theta."""
  c, s = math.cos(theta), math.sin(theta)
  lo = 1.0 / cond
  return (c * c + lo * s * s, c * s - lo * s * c, s * s + lo * c * c)


def chol2(s11, s12, s22):
  l11 = math.sqrt(s11)
  l21 = s12 / l11
  return l11, l21, math.sqrt(s22 - l21 * l21)


def inv2(s11, s12, s22):
  det = s11 * s22 - s12 * s12
  return s22 / det, -s12 / det, s11 / det


def _matmul(a, b, p):
  return [[sum(a[i][k] * b[k][j] for k in range(p)) for j in range(p)]
          for i in range(p)]


def _givens(p, i, j, angle):
  q = [[1.0 if r == c else 0.0 for c in range(p)] for r in range(p)]
  c, s = math.cos(angle), math.sin(angle)
  q[i][i], q[j][j] = c, c
  q[i][j], q[j][i] = -s, s
  return q


def spd5(cond):
  """5x5 SPD, eigenvalues log-spaced from 1 to 1/cond, fixed rotation."""
  p = 5
  lam = [10.0 ** (-math.log10(cond) * i / (p - 1)) for i in range(p)]
  q = [[1.0 if r == c else 0.0 for c in range(p)] for r in range(p)]
  for (i, j, ang) in ((0, 1, 0.4), (1, 2, 0.7), (2, 3, 1.1), (3, 4, 0.9),
                      (0, 4, 0.3)):
    q = _matmul(q, _givens(p, i, j, ang), p)
  d = [[lam[i] if i == j else 0.0 for j in range(p)] for i in range(p)]
  qt = [[q[j][i] for j in range(p)] for i in range(p)]
  s = _matmul(_matmul(q, d, p), qt, p)
  # force exact symmetry, so the probe and the reference see one matrix
  for i in range(p):
    for j in range(i):
      s[i][j] = s[j][i]
  return s


def chol5(s):
  p = 5
  l = [[0.0] * p for _ in range(p)]
  for i in range(p):
    for j in range(i + 1):
      acc = s[i][j] - sum(l[i][k] * l[j][k] for k in range(j))
      l[i][j] = math.sqrt(acc) if i == j else acc / l[j][j]
  return l


def inv5_from_chol(l):
  """S^-1 = L^-T L^-1: solve L X = I, then L^T S^-1 = X."""
  p = 5
  x = [[0.0] * p for _ in range(p)]
  for col in range(p):
    for i in range(p):
      acc = (1.0 if i == col else 0.0) - sum(l[i][k] * x[k][col]
                                            for k in range(i))
      x[i][col] = acc / l[i][i]
  out = [[0.0] * p for _ in range(p)]
  for col in range(p):
    for i in reversed(range(p)):
      acc = x[i][col] - sum(l[k][i] * out[k][col] for k in range(i + 1, p))
      out[i][col] = acc / l[i][i]
  return out


def stage_d5(n, rng):
  def triple(cond, ys):
    s = spd5(cond)
    l = chol5(s)
    pi = inv5_from_chol(l)
    fs = [s[i][j] for i in range(5) for j in range(5)]
    fl = [l[i][j] for i in range(5) for j in range(5)]
    fp = [pi[i][j] for i in range(5) for j in range(5)]
    for y in ys:
      emit("MVN5", y + fs)
      emit("MVNP5", y + fp)
      emit("MVNC5", y + fl)

  print("# D3 dimension 5, condition number half decades")
  for half_decade in range(0, 29):
    cond = 10.0 ** (half_decade / 2.0)
    base = [1.0, -0.5, 0.3, 0.8, -0.2]
    triple(cond, [base, [0.1 * v for v in base]])

  print("# D4 dimension 5, random condition numbers and directions")
  for _ in range(n):
    triple(log_uniform(rng, 1.0, 1e14),
           [[rng.uniform(-2, 2) for _ in range(5)]])


def stage_d(n, rng):
  thetas = (0.0, 0.3, math.pi / 4)
  print("# D1 condition number half decades, fixed directions")
  for half_decade in range(0, 29):
    cond = 10.0 ** (half_decade / 2.0)
    for theta in thetas:
      s11, s12, s22 = spd2(cond, theta)
      l11, l21, l22 = chol2(s11, s12, s22)
      p11, p12, p22 = inv2(s11, s12, s22)
      for y1, y2 in ((1.0, 0.0), (0.0, 1.0), (0.7, 0.7)):
        emit("MVN2", [y1, y2, s11, s12, s12, s22])
        emit("MVNP2", [y1, y2, p11, p12, p12, p22])
        emit("MVNC2", [y1, y2, l11, l21, l22])

  print("# D2 random condition numbers and directions")
  for _ in range(n):
    cond = log_uniform(rng, 1.0, 1e14)
    theta = rng.uniform(0, math.pi)
    s11, s12, s22 = spd2(cond, theta)
    l11, l21, l22 = chol2(s11, s12, s22)
    p11, p12, p22 = inv2(s11, s12, s22)
    y1, y2 = rng.uniform(-3, 3), rng.uniform(-3, 3)
    emit("MVN2", [y1, y2, s11, s12, s12, s22])
    emit("MVNP2", [y1, y2, p11, p12, p12, p22])
    emit("MVNC2", [y1, y2, l11, l21, l22])


def main():
  parser = argparse.ArgumentParser()
  parser.add_argument("--stage", choices=["a", "b", "c", "d"], required=True)
  parser.add_argument("--n", type=int, default=200)
  args = parser.parse_args()
  rng = random.Random(SEED)
  print(f"# density_sweep_work/gen_points.py --stage {args.stage} "
        f"--n {args.n}  seed {SEED}")
  if args.stage == "a":
    stage_a()
  elif args.stage == "b":
    stage_b(args.n, rng)
  elif args.stage == "c":
    stage_c(args.n, rng)
  else:
    stage_d(args.n, rng)
    stage_d5(args.n, rng)


if __name__ == "__main__":
  main()
