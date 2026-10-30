// Value and partials of std_normal_cdf, normal_cdf, lognormal_cdf and
// exp_mod_normal_cdf (mu = 0.3, sigma = 1.7, lambda = 1) on a z grid.
// Output: "<name> z value d_y d_mu d_sigma d_lambda" in hex ("nan" on throw).
#include <stan/math/rev.hpp>
#include <cmath>
#include <cstdio>
#include <exception>

using stan::math::var;
constexpr double MU = 0.3, SIGMA = 1.7, LAMBDA = 1.0;

template <typename F>
void row(const char* name, double z, double yv, F f) {
  var y = yv, mu = MU, sigma = SIGMA, lambda = LAMBDA;
  double o[5];
  try {
    var p = f(y, mu, sigma, lambda);
    p.grad();
    o[0] = p.val();
    o[1] = y.adj();
    o[2] = mu.adj();
    o[3] = sigma.adj();
    o[4] = lambda.adj();
  } catch (const std::exception&) {
    for (double& t : o) {
      t = std::nan("");
    }
  }
  stan::math::recover_memory();
  std::printf("%s %a %a %a %a %a %a\n", name, z, o[0], o[1], o[2], o[3], o[4]);
}

int main() {
  for (int i = 0; i <= 184; ++i) {
    const double z = -38.0 + 0.25 * i;
    const double y = MU + SIGMA * z;
    row("std_normal", z, z, [](auto& y, auto&, auto&, auto&) {
      return stan::math::std_normal_cdf(y);
    });
    row("normal", z, y, [](auto& y, auto& m, auto& s, auto&) {
      return stan::math::normal_cdf(y, m, s);
    });
    row("lognormal", z, std::exp(y), [](auto& y, auto& m, auto& s, auto&) {
      return stan::math::lognormal_cdf(y, m, s);
    });
    row("exp_mod_normal", z, y, [](auto& y, auto& m, auto& s, auto& l) {
      return stan::math::exp_mod_normal_cdf(y, m, s, l);
    });
  }
  return 0;
}
