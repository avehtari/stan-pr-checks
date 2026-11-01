#!/usr/bin/env python3
# Score the output of probe_beta_lccdf_dx against mpmath. The x derivative
# of log(1 - I_x(a, b)) is -x^(a-1) (1-x)^(b-1) / B(a, b) / (1 - I_x(a, b)).
# Route A: that closed form with the upper incomplete beta. Route B: mp.diff
# of log1p(-I_x(a, b)) with the lower incomplete beta.
# Usage: python3 score_lccdf_dx.py < probe output
import sys
import mpmath as mp

mp.mp.dps = 40
print(f"{'a':>4s} {'b':>11s} {'b*x':>5s} {'d/dx rel err':>12s} {'A vs B':>8s}")
for line in sys.stdin:
    f = line.split()
    if len(f) != 4:
        continue
    a, b, x = (mp.mpf(v) for v in f[:3])
    s = mp.mpf(f[3])
    upper = mp.betainc(b, a, 0, 1 - x, regularized=True)
    ra = -x ** (a - 1) * (1 - x) ** (b - 1) / mp.beta(a, b) / upper
    rb = mp.diff(lambda v: mp.log1p(-mp.betainc(a, b, 0, v, regularized=True)), x)
    print(f"{mp.nstr(a, 3):>4s} {mp.nstr(b, 4):>11s} {mp.nstr(b * x, 2):>5s}"
          f" {mp.nstr(abs(s - ra) / abs(ra), 2):>12s}"
          f" {mp.nstr(abs(ra - rb) / abs(ra), 2):>8s}")
