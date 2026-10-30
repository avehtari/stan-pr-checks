// Value and d/dy at the points chosen for the range-end unit tests
// (mu = 0, sigma = 1, lambda = 1, alpha = 1).  Output: "name y value dy" (hex).
#include <stan/math/rev.hpp>
#include <cmath>
#include <cstdio>
#include <string>

template <typename F>
void row(const std::string& name, F f, double y) {
  stan::math::var yv = y;
  stan::math::var lp = f(yv);
  lp.grad();
  std::printf("%s %a %a %a\n", name.c_str(), y, lp.val(), yv.adj());
  stan::math::recover_memory();
}

int main() {
  using stan::math::var;
  auto nl = [](const var& y) { return stan::math::normal_lcdf(y, 0.0, 1.0); };
  auto nc = [](const var& y) { return stan::math::normal_lccdf(y, 0.0, 1.0); };
  auto ll = [](const var& y) { return stan::math::lognormal_lcdf(y, 0.0, 1.0); };
  auto lc = [](const var& y) { return stan::math::lognormal_lccdf(y, 0.0, 1.0); };
  auto ep = [](const var& y) {
    return stan::math::exp_mod_normal_lpdf(y, 0.0, 1.0, 1.0);
  };
  auto el = [](const var& y) {
    return stan::math::exp_mod_normal_lcdf(y, 0.0, 1.0, 1.0);
  };
  auto ec = [](const var& y) {
    return stan::math::exp_mod_normal_lccdf(y, 0.0, 1.0, 1.0);
  };
  auto ef = [](const var& y) {
    return stan::math::exp_mod_normal_cdf(y, 0.0, 1.0, 1.0);
  };
  auto sn = [](const var& y) {
    return stan::math::skew_normal_lpdf(y, 0.0, 1.0, 1.0);
  };
  row("normal_lcdf", nl, -1e100);
  row("normal_lcdf", nl, -1.7e308);
  row("normal_lccdf", nc, 1e100);
  row("normal_lccdf", nc, 1.7e308);
  row("lognormal_lcdf", ll, std::exp(-700.0));
  row("lognormal_lccdf", lc, std::exp(700.0));
  row("exp_mod_normal_lpdf", ep, -1e100);
  row("exp_mod_normal_lcdf", el, -1e100);
  row("exp_mod_normal_lcdf", el, -1e300);
  row("exp_mod_normal_lcdf", el, -30.0);
  row("exp_mod_normal_lccdf", ec, 1e300);
  row("exp_mod_normal_lccdf", ec, -30.0);
  row("exp_mod_normal_cdf", ef, -30.0);
  row("skew_normal_lpdf", sn, -1e100);
  row("skew_normal_lpdf", sn, -1e300);
  return 0;
}
