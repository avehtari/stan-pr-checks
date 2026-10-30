#!/usr/bin/env python3
"""Compare ref_mpmath.py with Table 1 of Boik and Robinson-Cox (1998).

Boik, R. J. and Robinson-Cox, J. F. (1998). Derivatives of the Incomplete
Beta Function. Journal of Statistical Software 3(1).
https://doi.org/10.18637/jss.v003.i01 (Table 1 on page 16 of the PDF).
The table gives I and its first and second derivatives in p and q to 8
significant digits, from Maple numerical integration at 32 digits.

The first derivatives come from both routes of ref_mpmath.py; the script
prints their agreement.  The second derivatives come from mp.diff of
betainc_reg at 60 digits.  Expected: every value equal to 8 digits.

Usage: python3 refs/check_table1.py
"""

import os
import sys

from mpmath import mp, mpf

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ref_mpmath as R  # noqa: E402
from gen_test_refs import route_a  # noqa: E402

# (p, q, x): I, I^p, I^q, I^pp, I^qq, I^pq, as printed in Table 1
TABLE1 = {
    (1.5, 11.0, 0.001): (8.9170111e-04, -4.5720356e-03, 1.1845673e-04,
                         2.3080438e-02, 5.1418717e-06, -5.3324285e-04),
    (1.5, 11.0, 0.5): (9.9861069e-01, -2.5501997e-03, 9.0824388e-04,
                       -3.5047111e-03, -5.8941710e-04, 1.5603497e-03),
    (1000.0, 1000.0, 0.5): (5.0000000e-01, -8.9224793e-03, 8.9224793e-03,
                            4.4630987e-06, -4.4630987e-06, 0.0),
    (1000.0, 1000.0, 0.55): (9.9999632e-01, -3.6713108e-07, 4.0584118e-07,
                             -3.4809144e-08, -4.2964422e-08, 3.8682578e-08),
}
NAMES = ["I", "I^p", "I^q", "I^pp", "I^qq", "I^pq"]


def main():
  n_equal = 0
  n_total = 0
  for (p, q, x), table in TABLE1.items():
    with mp.workdps(R.DPS_INTEGRAL):
      i = R.betainc_reg(p, q, x)
    ip, iq = R.dIda(p, q, x), R.dIdb(p, q, x)
    dig_p = float(R.rel_digits(route_a(p, q, x, "a"), ip))
    dig_q = float(R.rel_digits(route_a(p, q, x, "b"), iq))
    with mp.workdps(60):
      f = lambda s, t: R.betainc_reg(s, t, mpf(x))  # noqa: E731
      pq = (mpf(p), mpf(q))
      ipp = mp.diff(f, pq, (2, 0))
      iqq = mp.diff(f, pq, (0, 2))
      ipq = mp.diff(f, pq, (1, 1))
    print(f"(p, q, x) = ({p:g}, {q:g}, {x:g}); routes A and B agree to "
          f"{dig_p:.0f} (I^p) and {dig_q:.0f} (I^q) digits")
    for name, t, v in zip(NAMES, table, (i, ip, iq, ipp, iqq, ipq)):
      v = float(v)
      same = f"{t:.7e}" == f"{v:.7e}"
      n_equal += same
      n_total += 1
      print(f"   {name:5s} table {t: .7e}   ref_mpmath {v: .7e}   "
            f"{'equal' if same else 'DIFFER'}")
  print(f"# {n_equal} of {n_total} values equal to 8 significant digits")
  return 0 if n_equal == n_total else 1


if __name__ == "__main__":
  sys.exit(main())
