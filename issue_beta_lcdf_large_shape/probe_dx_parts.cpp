// The parts of the beta_lcdf x derivative on develop:
// pow(1 - x, b - 1) * pow(x, a - 1) / beta(a, b) / P. Output per point:
// a b x beta(a, b) pow(1 - x, b - 1) exp(lbeta(a, b)).
#include <stan/math/prim.hpp>
#include <cmath>
#include <cstdio>
#include <vector>

int main() {
  for (double a : {0.5, 2.0, 10.0}) {
    for (double b : {1e3, 1e6, 1e9, 2147483648.0}) {
      const double x = 2.0 / b;
      std::printf("%.17g %.17g %.17g %.17g %.17g %.17g\n", a, b, x,
                  stan::math::beta(a, b), std::pow(1 - x, b - 1),
                  std::exp(stan::math::lbeta(a, b)));
    }
  }
  return 0;
}
