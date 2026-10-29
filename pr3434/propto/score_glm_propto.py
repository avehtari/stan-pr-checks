#!/usr/bin/env python3
"""Score glm_propto output (claim P4).

Usage: score_glm_propto.py <glm_propto output> [<second output> ...]

Under propto, with data phi and parameters alpha and beta, the GLM keeps
  sum_i y_i theta_i - (y_i + phi) log(exp(theta_i) + phi)
(develop's rule; large_shapes_propto_dropped_terms tests it). Two routes:
  A  the formula above at 60 digits
  B  the full log pmf (ref_sib.lp_GLM at 90 digits) minus the dropped terms
     sum_i lgamma(y_i + phi) - lgamma(phi) - lgamma(y_i + 1) + phi log(phi)
Prints per phi: route agreement, the absolute error of the propto value and
of the full value, and eps * N * phi * log(phi) for comparison.
"""

import os
import sys

from mpmath import mp, mpf

sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(
    os.path.abspath(__file__))), "refs"))
import ref_sib as RS  # noqa: E402

Y = [0, 1, 3, 10, 57]
X = [1.0, -0.5, 2.0, 0.25, 1.5]
ALPHA, BETA = 0.5, 0.3
EPS = 2.0 ** -52


def refs(phi):
  mp.dps = 90
  ph = mpf(phi)
  thetas = [mpf(ALPHA) + mpf(x) * mpf(BETA) for x in X]
  a = sum(y * t - (y + ph) * mp.log(mp.exp(t) + ph) for y, t in zip(Y, thetas))
  full = mpf(0)
  dropped = mpf(0)
  for y, x in zip(Y, X):
    full += RS.lp_GLM(y, x, ALPHA, BETA, phi)
    dropped += (mp.loggamma(y + ph) - mp.loggamma(ph) - mp.loggamma(y + 1)
                + ph * mp.log(ph))
  return a, full - dropped, full


def main():
  for path in sys.argv[1:]:
    print(f"== {os.path.basename(path)}")
    print(f"{'phi':>7s} {'routes A-B':>10s} {'propto err':>10s} "
          f"{'eps*|propto|':>12s} {'full err':>9s} "
          f"{'eps*N*phi*log(phi)':>18s}")
    for line in open(path):
      left, right = line.split(" => ")
      phi = float.fromhex(left)
      v_t, v_f = (float.fromhex(t) for t in right.split()[:2])
      a, b, full = refs(phi)
      bound = EPS * len(Y) * phi * float(mp.log(phi))
      print(f"{phi:7.0e} {float(abs(a - b) / abs(a)):10.1e} "
            f"{float(abs(mpf(v_t) - a)):10.1e} {EPS * float(abs(a)):12.1e} "
            f"{float(abs(mpf(v_f) - full)):9.1e} {bound:18.1e}")
  return 0


if __name__ == "__main__":
  sys.exit(main())
