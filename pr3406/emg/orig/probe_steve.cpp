// Probe of SteveBronder's #3406 commits e1fea0651c and cf177ef96b.
// Output lines (%.17g), scored against mpmath by score_steve.py:
//   ERFCX x value
//   LCDF z value slope            (std_normal_lcdf, var)
//   EDGE name value g1 g2 ...     (edge cases, checked by eye)
#include <stan/math/rev.hpp>
#include <cmath>
#include <cstdio>
#include <vector>

using stan::math::var;

int main() {
  // 1. erfcx on [4, 2^28], log-spaced, plus the branch edges
  for (int i = 0; i <= 400; ++i) {
    const double x = 4.0 * std::pow(2.0, 26.0 * i / 400.0);
    std::printf("ERFCX %.17g %.17g\n", x, stan::math::erfcx(x));
  }
  for (double x : {4.0, std::nextafter(4.0, 5.0), 0x1p27, std::nextafter(0x1p27, 1e300),
                   0x1p28, 1e10, 1e300}) {
    std::printf("ERFCX %.17g %.17g\n", x, stan::math::erfcx(x));
  }
  // 2. std_normal_lcdf value and slope, var
  std::vector<double> zs;
  for (int i = 0; i <= 300; ++i) {
    zs.push_back(-std::pow(10.0, 300.0 * i / 300.0));  // -1 .. -1e300
  }
  for (int i = 0; i <= 400; ++i) {
    zs.push_back(-40.0 + 80.0 * i / 400.0);
  }
  for (double z : {-4 * std::sqrt(2.0), std::nextafter(-4 * std::sqrt(2.0), 0.0),
                   std::nextafter(-4 * std::sqrt(2.0), -10.0), -0.46875 * std::sqrt(2.0),
                   -4.0 * std::sqrt(2.0) * 0.999999}) {
    zs.push_back(z);
  }
  for (double z : zs) {
    var zv = z;
    var lp = stan::math::std_normal_lcdf(zv);
    lp.grad();
    std::printf("LCDF %.17g %.17g %.17g\n", z, lp.val(), zv.adj());
    stan::math::recover_memory();
  }
  // 3. subnormal sigma in normal_lcdf and lognormal_lcdf (y = mu: log 0.5)
  for (double s : {1e-300, 1e-308, 1e-310, 4.9e-324}) {
    var mu = 0.0, sigma = s;
    var lp = stan::math::normal_lcdf(0.0, mu, sigma);
    lp.grad();
    std::printf("EDGE normal_lcdf_y_eq_mu sigma=%g %.17g %.17g %.17g\n", s, lp.val(),
                mu.adj(), sigma.adj());
    stan::math::recover_memory();
    var mu2 = 0.0, sigma2 = s;
    var lp2 = stan::math::lognormal_lcdf(1.0, mu2, sigma2);
    lp2.grad();
    std::printf("EDGE lognormal_lcdf_y_eq_e^mu sigma=%g %.17g %.17g %.17g\n", s,
                lp2.val(), mu2.adj(), sigma2.adj());
    stan::math::recover_memory();
  }
  // 4. exp_mod_normal_cdf with var far in the lower tail
  int nonfinite = 0, total = 0;
  for (double z : {-1e300, -1e200, -1e160, -1e155, -1e154, -1e150, -1e100, -1e20,
                   -1e5, -40.0, -10.0, -1.0, 0.0, 1.0, 10.0}) {
    for (double v : {1e-300, 1e-100, 1e-10, 0.1, 1.0, 10.0, 1e10}) {
      var y = z, mu = 0.0, sigma = 1.0, lambda = v;
      var c = stan::math::exp_mod_normal_cdf(y, mu, sigma, lambda);
      c.grad();
      const bool ok = std::isfinite(c.val()) && std::isfinite(y.adj())
                      && std::isfinite(mu.adj()) && std::isfinite(sigma.adj())
                      && std::isfinite(lambda.adj());
      nonfinite += !ok;
      ++total;
      if (!ok) {
        std::printf("EDGE emg_cdf_nonfinite z=%g v=%g %.17g %g %g %g %g\n", z, v,
                    c.val(), y.adj(), mu.adj(), sigma.adj(), lambda.adj());
      }
      stan::math::recover_memory();
    }
  }
  std::printf("EDGE emg_cdf_nonfinite_count %d of %d\n", nonfinite, total);
  // 5. the same grid for the lcdf and lccdf (a value of -inf counts as
  // non-finite only if the gradient is also non-finite)
  for (int which = 0; which < 2; ++which) {
    int bad = 0;
    for (double z : {-1e300, -1e200, -1e160, -1e155, -1e154, -1e150, -1e100,
                     -1e20, -1e5, -40.0, -10.0, -1.0, 0.0, 1.0, 10.0}) {
      for (double v : {1e-300, 1e-100, 1e-10, 0.1, 1.0, 10.0, 1e10}) {
        var y = z, mu = 0.0, sigma = 1.0, lambda = v;
        var c = which == 0 ? stan::math::exp_mod_normal_lcdf(y, mu, sigma, lambda)
                           : stan::math::exp_mod_normal_lccdf(y, mu, sigma, lambda);
        c.grad();
        const bool ok = !std::isnan(c.val()) && std::isfinite(y.adj())
                        && std::isfinite(mu.adj()) && std::isfinite(sigma.adj())
                        && std::isfinite(lambda.adj());
        bad += !ok;
        stan::math::recover_memory();
      }
    }
    std::printf("EDGE emg_%s_nonfinite_count %d of 105\n",
                which == 0 ? "lcdf" : "lccdf", bad);
  }
  return 0;
}
