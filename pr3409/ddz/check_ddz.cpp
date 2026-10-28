#include <stan/math/mix.hpp>
#include <boost/math/special_functions/beta.hpp>
#include <cmath>
#include <cstdio>
#include <cstring>

// inc_beta_ddz(double) with the overflow policy against Boost's default
// ibeta_derivative: bit-identical wherever the default does not throw, and
// +inf where it throws.  Grid: 60 x 60 log-spaced shapes in
// [exp(-5), exp(6)] and z = 0, 0.02, ..., 1 (183 600 points).
// Build against the PR tree.
int main() {
  long n = 0, same = 0, default_threw = 0;
  for (int i = 0; i < 60; ++i) {
    const double a = std::exp(-5.0 + 11.0 * i / 59.0);
    for (int j = 0; j < 60; ++j) {
      const double b = std::exp(-5.0 + 11.0 * j / 59.0);
      for (int k = 0; k <= 50; ++k) {
        const double z = k / 50.0;
        double ref;
        try {
          ref = boost::math::ibeta_derivative(a, b, z);
        } catch (const std::exception&) {
          ++default_threw;
          const double v = stan::math::inc_beta_ddz(a, b, z);
          if (!std::isinf(v)) {
            std::printf("not inf where default throws: %g %g %g -> %g\n", a,
                        b, z, v);
          }
          continue;
        }
        const double v = stan::math::inc_beta_ddz(a, b, z);
        ++n;
        same += std::memcmp(&ref, &v, sizeof(double)) == 0;
      }
    }
  }
  std::printf("points %ld, bit-identical %ld, default threw %ld\n", n, same,
              default_threw);
  std::printf("inc_beta_ddz(1.05, 0.35, 1) = %g, (0.5, 2, 0) = %g\n",
              stan::math::inc_beta_ddz(1.05, 0.35, 1.0),
              stan::math::inc_beta_ddz(0.5, 2.0, 0.0));
  return 0;
}
