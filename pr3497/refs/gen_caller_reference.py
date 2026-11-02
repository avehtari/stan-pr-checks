#!/usr/bin/env python3
"""Reference rows of the beta_lcdf caller test (rev/prob/beta_lcdf_test.cpp).

Columns: x, a, b, beta_lcdf(x | a, b), its x derivative.  The values come
from refs() of score_caller.py (two routes of I_x(a, b), 50 digits).
"""
import mpmath as mp

from score_caller import refs

POINTS = [(2.0**-29, 2.0, 2.0**30), (2.0**-32, 2.0, 2.0**31)]

for x, a, b in POINTS:
    lcdf, ddx, agree = refs(a, b, x)
    print("{%r, %r, %r, %r, %r},  // routes agree to %s"
          % (x, a, b, float(lcdf), float(ddx), mp.nstr(agree, 2)))
