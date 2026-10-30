// Claim M4 of README_emg.md: non-finite results of exp_mod_normal cdf, lcdf
// and lccdf (var, value and four partials) on the 105-point deep-tail grid
// (15 z from -1e300 to 10, 7 v from 1e-300 to 1e10).  Sections 4 and 5 of
// orig/probe_steve.cpp, with one line per non-finite point.
// Counting rule as in the original: cdf, every component must be finite;
// lcdf and lccdf, the value must not be NaN (-inf is permitted) and the four
// partials must be finite.
#include <stan/math/rev.hpp>
#include <cmath>
#include <cstdio>
#include <exception>

using stan::math::var;

int main() {
  const double zs[] = {-1e300, -1e200, -1e160, -1e155, -1e154, -1e150, -1e100,
                       -1e20,  -1e5,   -40.0,  -10.0,  -1.0,   0.0,    1.0,
                       10.0};
  const double vs[] = {1e-300, 1e-100, 1e-10, 0.1, 1.0, 10.0, 1e10};
  const char* names[3] = {"cdf", "lcdf", "lccdf"};
  for (int which = 0; which < 3; ++which) {
    int bad = 0, total = 0;
    for (double z : zs) {
      for (double v : vs) {
        var y = z, mu = 0.0, sigma = 1.0, lambda = v;
        double o[5];
        try {
          var c = which == 0
                      ? stan::math::exp_mod_normal_cdf(y, mu, sigma, lambda)
                  : which == 1
                      ? stan::math::exp_mod_normal_lcdf(y, mu, sigma, lambda)
                      : stan::math::exp_mod_normal_lccdf(y, mu, sigma, lambda);
          c.grad();
          o[0] = c.val();
          o[1] = y.adj();
          o[2] = mu.adj();
          o[3] = sigma.adj();
          o[4] = lambda.adj();
        } catch (const std::exception&) {
          for (double& x : o) {
            x = std::nan("");
          }
        }
        stan::math::recover_memory();
        const bool value_ok
            = which == 0 ? std::isfinite(o[0]) : !std::isnan(o[0]);
        const bool ok = value_ok && std::isfinite(o[1]) && std::isfinite(o[2])
                        && std::isfinite(o[3]) && std::isfinite(o[4]);
        bad += !ok;
        ++total;
        if (!ok) {
          std::printf("NONFINITE %s z=%g v=%g val=%g dy=%g dmu=%g dsigma=%g"
                      " dlambda=%g\n",
                      names[which], z, v, o[0], o[1], o[2], o[3], o[4]);
        }
      }
    }
    std::printf("COUNT emg_%s_nonfinite %d of %d\n", names[which], bad, total);
  }
  return 0;
}
