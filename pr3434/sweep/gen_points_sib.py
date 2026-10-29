#!/usr/bin/env python3
"""Generate the hex point files for the large-shape spot check.

Output:
  points.txt       every point that runs fast on develop
  points_slow.txt  beta_neg_binomial cdf points whose 3F2 series needs about
                   `shape` terms (configs all, r, b at shape >= 1e8)

Line format (read by probe_sib.cpp and ref_sib.py):
  <TAG> <LABEL> <hex arg> <hex arg> ...

LABEL is `<config>:<decade>`, for example `all:1e15`. The scorer groups by
tag, config and decade. Argument order per tag is the ARITY table in
probe_sib.cpp and in ref_sib.py.
"""

import math

from mpmath import mp, mpf

DECADES = [1e2, 1e4, 1e6, 1e8, 1e10, 1e12, 1e15]
DEC_NAME = {1e2: "1e2", 1e4: "1e4", 1e6: "1e6", 1e8: "1e8", 1e10: "1e10",
            1e12: "1e12", 1e15: "1e15"}
INT_MAX = 2**31 - 1


def h(x):
  return float(x).hex()


def line(tag, label, args):
  return f"{tag} {label} " + " ".join(h(a) for a in args)


def nb_family(out):
  """NB, NB2, NB2L: two means, counts 0, 1, 5, 57 and the mean."""
  cases = [(3.0, [0, 1, 3, 57]), (50.0, [5, 50])]
  for s in DECADES:
    for mu, counts in cases:
      cfg = f"mu{int(mu)}"
      lab = f"{cfg}:{DEC_NAME[s]}"
      for n in counts:
        out.append(line("NB", lab, [n, s, s / mu]))
        out.append(line("NB2", lab, [n, mu, s]))
        out.append(line("NB2L", lab, [n, math.log(mu), s]))


GLM_X5 = [0.5, -0.25, 1.25, 0.0, 2.0]
GLM_Y5 = [[0, 1, 5, 3, 57], [2, 4, 3, 5, 6]]


def glm(out):
  for s in DECADES:
    lab = f"one:{DEC_NAME[s]}"
    for y in [0, 5, 57]:
      out.append(line("GLM1", lab, [y, 0.5, 1.0, 0.75, s]))
    lab = f"five:{DEC_NAME[s]}"
    for ys in GLM_Y5:
      out.append(line("GLM5", lab, ys + GLM_X5 + [1.0, 0.75, s]))


def bnb_configs(s):
  """(config, r, alpha, beta, mean) for beta_neg_binomial."""
  return [
      ("all", s, s, s, s * s / (s - 1)),
      ("r", s, 5.0, 2.0, s * 2.0 / 4.0),
      ("a", 3.0, s, 2.0, 6.0 / (s - 1)),
      ("b", 3.0, 5.0, s, 3.0 * s / 4.0),
      ("ab", 5.0, s, 2.0 * s, 5.0 * 2.0 * s / (s - 1)),
  ]


def bnb_counts(mean):
  counts = [0, 1, 5, 57]
  m = int(round(mean))
  if 1 < m < INT_MAX and m not in counts:
    counts.append(m)
  return counts


def bnb(out):
  for s in DECADES:
    for cfg, r, a, b, mean in bnb_configs(s):
      lab = f"{cfg}:{DEC_NAME[s]}"
      for n in bnb_counts(mean):
        out.append(line("BNB", lab, [n, r, a, b]))


DM_COUNTS = [(0, 1, 0), (1, 2, 2), (5, 0, 0), (57, 3, 10), (10, 15, 25)]
DM_PROP = [0.2, 0.3, 0.5]


def dm(out):
  for s in DECADES:
    lab = f"k:{DEC_NAME[s]}"
    alpha = [s * p for p in DM_PROP]
    for ns in DM_COUNTS:
      out.append(line("DM", lab, list(ns) + alpha))


LKJ_R = (0.3, -0.2, 0.1)  # r12, r13, r23


def lkj_cholesky():
  """Lower Cholesky factor of the LKJ_R matrix, rounded to double."""
  with mp.workdps(60):
    r12, r13, r23 = (mpf(x) for x in LKJ_R)
    m = mp.matrix([[1, r12, r13], [r12, 1, r23], [r13, r23, 1]])
    L = mp.cholesky(m)
    return [float(L[1, 0]), float(L[1, 1]), float(L[2, 0]), float(L[2, 1]),
            float(L[2, 2])]


def lkj(out):
  L = lkj_cholesky()
  for s in DECADES:
    lab = f"K3:{DEC_NAME[s]}"
    for eta in [s, 3.3 * s]:
      out.append(line("LKJ", lab, list(LKJ_R) + [eta]))
      out.append(line("LKJC", lab, L + [eta]))


def bb_cdf(out):
  N = 117
  for s in DECADES:
    for cfg, a, b, counts in [("tied", s, s, [0, 1, 5, 57]),
                              ("ab3", s, 3.0 * s, [0, 1, 5, 29, 57])]:
      lab = f"{cfg}:{DEC_NAME[s]}"
      for n in counts:
        for tag in ["BBC", "BBLC", "BBLCC"]:
          out.append(line(tag, lab, [n, N, a, b]))


def ys_cdf(out):
  for s in DECADES:
    lab = f"a:{DEC_NAME[s]}"
    for n in [1, 2, 5, 57]:
      for tag in ["YSC", "YSLC", "YSLCC"]:
        out.append(line(tag, lab, [n, s]))


def bnb_cdf(out, slow):
  for s in DECADES:
    for cfg, r, a, b, mean in bnb_configs(s):
      lab = f"{cfg}:{DEC_NAME[s]}"
      if cfg in ("a", "ab"):
        counts = [0, 1, 5, 57] + ([10] if cfg == "ab" else [])
        dest = out
      else:
        # the 3F2 series of develop needs about `shape` terms here
        counts = [0, 5]
        if s > 1e6:
          if s > 1e8:
            continue
          counts = [5]
          dest = slow
        else:
          dest = out
      for n in counts:
        for tag in ["BNBC", "BNBLC", "BNBLCC"]:
          if dest is slow and tag != "BNBC":
            continue
          dest.append(line(tag, lab, [n, r, a, b]))


def extra(out):
  """Second pass, written to points_extra.txt.

  - NB2, NB2L, GLM1 at a non-round phi = 1.37 * shape and mu = 0.7: at
    phi = 1e15 the develop phi partial of NB2 rounds to exactly 0 for the
    first points, which hides the one-ulp error of the psi difference.
  - LKJ and LKJC at the identity matrix: the value is then only the
    normalizing constant, without the (eta - 1) log det term.
  - beta_binomial and beta_neg_binomial cdfs at shapes 1 and 10: the
    baseline at ordinary shapes.
  """
  for s in DECADES:
    lab = f"x1.37:{DEC_NAME[s]}"
    phi = 1.37 * s
    for n in [0, 1, 2, 57]:
      out.append(line("NB2", lab, [n, 0.7, phi]))
      out.append(line("NB2L", lab, [n, math.log(0.7), phi]))
      out.append(line("GLM1", lab, [n, 0.5, -0.2, 0.75, phi]))
    lab = f"I:{DEC_NAME[s]}"
    for eta in [s, 3.3 * s]:
      out.append(line("LKJ", lab, [0.0, 0.0, 0.0, eta]))
      out.append(line("LKJC", lab, [0.0, 1.0, 0.0, 0.0, 1.0, eta]))
  for s, name in [(1.0, "1e0"), (10.0, "1e1")]:
    for cfg, a, b, counts in [("tied", s, s, [0, 1, 5, 57]),
                              ("ab3", s, 3.0 * s, [0, 1, 5, 29, 57])]:
      for n in counts:
        for tag in ["BBC", "BBLC", "BBLCC"]:
          out.append(line(tag, f"{cfg}:{name}", [n, 117, a, b]))
  s = 10.0
  for cfg, r, a, b, mean in bnb_configs(s):
    for n in [0, 1, 5, 57]:
      for tag in ["BNBC", "BNBLC", "BNBLCC"]:
        out.append(line(tag, f"{cfg}:1e1", [n, r, a, b]))


def main():
  ext = []
  extra(ext)
  with open("points_extra.txt", "w") as fh:
    fh.write("\n".join(ext) + "\n")
  print(f"points_extra.txt: {len(ext)} lines")
  out = []
  slow = []
  nb_family(out)
  glm(out)
  bnb(out)
  dm(out)
  lkj(out)
  bb_cdf(out)
  ys_cdf(out)
  bnb_cdf(out, slow)
  with open("points.txt", "w") as fh:
    fh.write("\n".join(out) + "\n")
  with open("points_slow.txt", "w") as fh:
    fh.write("\n".join(slow) + "\n")
  print(f"points.txt: {len(out)} lines, points_slow.txt: {len(slow)} lines")


if __name__ == "__main__":
  main()
