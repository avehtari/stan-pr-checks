#!/usr/bin/env python3
"""TensorFlow Probability's betainc partials at three integer-shape points.

Claim: TFP's d/db is off by 4.2e-4 at (2, 3, 0.25) and by 6.7 % at
(10, 1, 0.5), and its d/da by 5.6 % at (1, 100, 0.01).  TFP stops on the
value alone; for integer b the value has converged while the derivative
series has not.  At (1, 100, 0.01) the symmetry relation swaps the shapes,
so the integer shape there is a = 1.  The second row of each point moves
both shapes by 1e-7; the error then goes away, which shows that the
integer shape is the cause.

Install TFP 0.25.0 (numpy substrate, no TensorFlow) into a separate
directory with the Python that runs this script:

  python3 -m pip install --target /tmp/tfp_lib --no-deps \
      tensorflow-probability==0.25.0 absl-py decorator cloudpickle gast dm-tree
  PYTHONPATH=/tmp/tfp_lib python3 tfp/tfp_betainc.py

The reference is mpmath (refs/ref_mpmath.py, two routes).
"""

import os
import sys

import numpy as np
from tensorflow_probability.substrates.numpy.math import special as S

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                "..", "refs"))
import ref_mpmath as R  # noqa: E402

POINTS = [(2.0, 3.0, 0.25), (10.0, 1.0, 0.5), (1.0, 100.0, 0.01)]


def tfp(a, b, z):
  da, db, _ = S._betainc_partials(np.float64(a), np.float64(b),
                                  np.float64(z))
  return float(da), float(db)


def main():
  import tensorflow_probability as tfp_mod
  print(f"# tensorflow_probability {tfp_mod.__version__}, numpy substrate")
  print(f"{'a':>12} {'b':>12} {'z':>6} {'wrt':>4} {'TFP':>24} "
        f"{'reference':>24} {'rel err':>9}")
  for a0, b0, z in POINTS:
    for a, b in ((a0, b0), (a0 + 1e-7, b0 + 1e-7)):
      da, db = tfp(a, b, z)
      for wrt, got, ref in (("a", da, R.dIda(a, b, z)),
                            ("b", db, R.dIdb(a, b, z))):
        ref = float(ref)
        print(f"{a:>12.10g} {b:>12.10g} {z:>6g} {wrt:>4} {got!r:>24} "
              f"{ref!r:>24} {abs(got - ref) / abs(ref):>9.1e}")
  return 0


if __name__ == "__main__":
  sys.exit(main())
