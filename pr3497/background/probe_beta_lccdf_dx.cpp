// The x derivative of beta_lccdf(x | a, b) at large b, x = t / b. On
// develop it divides by beta(a, b); at the head of #3433 it uses
// inc_beta_ddz. Output per point: a b x d/dx, for score_lccdf_dx.py.
#include <stan/math/rev.hpp>
#include <cstdio>
#include <vector>

int main() {
  for (double a0 : {0.5, 2.0, 10.0}) {
    for (double b0 : {1e3, 1e6, 1e9, 2147483648.0}) {
      for (double t : {0.5, 2.0, 8.0}) {
        stan::math::var x = t / b0;
        stan::math::var lp = stan::math::beta_lccdf(x, a0, b0);
        lp.grad();
        std::printf("%.17g %.17g %.17g %.17g\n", a0, b0, t / b0, x.adj());
        stan::math::recover_memory();
      }
    }
  }
  return 0;
}
