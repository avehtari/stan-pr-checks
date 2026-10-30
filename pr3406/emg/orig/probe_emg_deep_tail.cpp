#include <stan/math/rev.hpp>
#include <cstdio>
using stan::math::var;
int main() {
  for (double z : {-1e300, -1e200, -1e160, -1e155, -1e154, -1e150, -1e100, -1e20, -1e5, -40.0, -10.0, -1.0, 0.0, 1.0, 10.0})
    for (double v : {1e-300, 1e-100, 1e-10, 0.1, 1.0, 10.0, 1e10})
      for (int w = 0; w < 2; ++w) {
        var y = z, mu = 0.0, sigma = 1.0, lambda = v;
        var c = w == 0 ? stan::math::exp_mod_normal_lcdf(y, mu, sigma, lambda)
                       : stan::math::exp_mod_normal_lccdf(y, mu, sigma, lambda);
        c.grad();
        bool ok = !std::isnan(c.val()) && std::isfinite(y.adj()) && std::isfinite(mu.adj()) && std::isfinite(sigma.adj()) && std::isfinite(lambda.adj());
        if (!ok || (w == 0 && std::isnan(c.val())))
          std::printf("%s z=%g v=%g val=%g dy=%g dsigma=%g dlambda=%g\n", w ? "lccdf" : "lcdf ", z, v, c.val(), y.adj(), sigma.adj(), lambda.adj());
        stan::math::recover_memory();
      }
}
