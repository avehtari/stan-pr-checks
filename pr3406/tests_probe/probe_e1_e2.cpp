// Probe for the new test points: normal_lcdf/lccdf gradients at infinite y
// (issue 2881) and Phi at develop's worst point in [-5, -0.67].
#include <stan/math/mix.hpp>
#include <cstdio>
#include <vector>

using stan::math::var;

template <typename F>
void endpoint(const char* name, F f, double y) {
  var mu = 0.3, sigma = 1.7;
  var lp = f(y, mu, sigma);
  stan::math::grad(lp.vi_);
  std::printf("%-22s y=%-5g value=%-10g mu.adj=%-10g sigma.adj=%g\n", name, y,
              lp.val(), mu.adj(), sigma.adj());
  stan::math::recover_memory();
}

template <typename F>
void endpoint_vec(const char* name, F f, double y) {
  std::vector<var> mu{0.3, -2.0}, sigma{1.7, 0.5};
  std::vector<double> yv{y, 1.0};
  var lp = f(yv, mu, sigma);
  stan::math::grad(lp.vi_);
  std::printf("%-22s y={%g,1} value=%-10g mu0.adj=%-10g sigma0.adj=%g\n",
              name, y, lp.val(), mu[0].adj(), sigma[0].adj());
  stan::math::recover_memory();
}

int main() {
  const double inf = stan::math::INFTY;
  auto lcdf = [](auto y, auto m, auto s) { return stan::math::normal_lcdf(y, m, s); };
  auto lccdf = [](auto y, auto m, auto s) { return stan::math::normal_lccdf(y, m, s); };
  endpoint("normal_lcdf", lcdf, inf);
  endpoint("normal_lccdf", lccdf, -inf);
  endpoint_vec("normal_lcdf vec", lcdf, inf);
  endpoint_vec("normal_lccdf vec", lccdf, -inf);
  const double x = -0x1.3f4395a1d3d18p+2;  // -4.988500030550561
  const double ref = 3.042495096713287367514609e-07;
  const double p = stan::math::Phi(x);
  std::printf("Phi(%.17g) = %.17g rel.err %.3g\n", x, p, std::fabs(p / ref - 1));
  return 0;
}
