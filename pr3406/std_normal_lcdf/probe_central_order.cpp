// Claim N2 (comment on 5f33a10d31): 9.6 ulp on a grid of 20001 points.
// Calls internal::erfcx_small(double) of the tree in STAN_MATH, not the kernel.
// Central range |z| < 0.46875 sqrt(2): three summation orders of
//   log Phi(z) = LOG_HALF - z^2/2 + log(erfcx(-z / sqrt 2))
//   a: (LOG_HALF - h) + L     b: (LOG_HALF + L) - h     c: LOG_HALF + (L - h)
// Output: "# central lo hi n" then "z a b c slope" in hex.
#include <stan/math/prim/fun/erfcx.hpp>
#include <cmath>
#include <cstdio>

int main() {
  namespace sm = stan::math;
  const double lo = -0.66291, hi = 0.66291;
  const int n = 20001;
  std::printf("# central %.17g %.17g %d\n", lo, hi, n);
  for (int i = 0; i < n; ++i) {
    double t = (i + 0.5 * std::sin(1.0 + 7.0 * i)) / (n - 1);
    t = t < 0 ? 0 : (t > 1 ? 1 : t);
    const double z = lo + (hi - lo) * t;
    const double e = sm::internal::erfcx_small(-z * sm::INV_SQRT_TWO);
    const double h = (0.5 * z) * z;
    const double l = std::log(e);
    std::printf("%a %a %a %a %a\n", z, (sm::LOG_HALF - h) + l,
                (sm::LOG_HALF + l) - h, sm::LOG_HALF + (l - h),
                sm::SQRT_TWO_OVER_SQRT_PI / e);
  }
  return 0;
}
