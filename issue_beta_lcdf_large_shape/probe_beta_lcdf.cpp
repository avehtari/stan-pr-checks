// beta_lcdf(x | a, b) in the bulk at large b: x = t / b. Output per point:
// a b x value d/dx d/da d/db boost_ibeta, for score_beta_lcdf.py.
#include <stan/math/rev.hpp>
#include <boost/math/special_functions/beta.hpp>
#include <cstdio>
#include <vector>

int main() {
  const std::vector<double> as = {0.5, 2, 10};
  const std::vector<double> bs = {1e3, 1e5, 1e6, 1e7, 1e8, 1e9, 2147483648.0};
  const std::vector<double> ts = {0.5, 2, 8};
  for (double a0 : as) {
    for (double b0 : bs) {
      for (double t : ts) {
        const double x0 = t / b0;
        stan::math::var x = x0, a = a0, b = b0;
        stan::math::var lp = stan::math::beta_lcdf(x, a, b);
        lp.grad();
        const double ib
            = boost::math::ibeta(a0, b0, x0, stan::math::boost_policy_t<>());
        std::printf("%.17g %.17g %.17g %.17g %.17g %.17g %.17g %.17g\n", a0, b0,
                    x0, lp.val(), x.adj(), a.adj(), b.adj(), ib);
        stan::math::recover_memory();
      }
    }
  }
  return 0;
}
