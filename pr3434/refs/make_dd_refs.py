#!/usr/bin/env python3
"""mpmath references for the OpenCL digamma_diff device-function test.

digamma_diff(x, d) = psi(x + d) - psi(x), for exact double inputs.
Route A: mp.digamma at 80 digits. Route B: the same at 120 digits.
Prints C++ rows {x, d, ref} with x and d in hex and ref to 17 digits.
"""

from mpmath import mp, mpf

POINTS = [
    (1e-300, 1.0), (1e-3, 1.0), (0.5, 0.0), (0.5, 1e-10), (0.5, 3.0),
    (2.5, 0.25), (7.3, 57.0), (9.999, 1.0), (10.0, 1.0), (12.5, 1e6),
    (1e3, 1.0), (1e8, 3.0), (1e12, 57.0), (1e15, 3.0), (1e15, 1e15),
    (3e15, 117.0), (1e20, 1.0), (1.0, 1e300), (1e300, 1e300),
]


def ref(x, d, dps):
    mp.dps = dps
    return mp.digamma(mpf(x) + mpf(d)) - mp.digamma(mpf(x))


def main():
    for x, d in POINTS:
        a = ref(x, d, 80)
        b = ref(x, d, 120)
        mp.dps = 120
        rel = abs(a - b) / abs(b) if b != 0 else abs(a - b)
        assert rel < mpf(10) ** -60, (x, d, rel)
        mp.dps = 30
        print("    {%s, %s, %s}," % (float.hex(x), float.hex(d),
                                      mp.nstr(b, 17, min_fixed=-1,
                                              max_fixed=1)
                                      if b != 0 else "0.0"))


if __name__ == "__main__":
    main()
