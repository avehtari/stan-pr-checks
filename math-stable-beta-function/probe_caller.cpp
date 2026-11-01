// beta_lcdf(x | 2, b) and its x derivative at b = 2^k and x = 2^-m, so that
// x, 1 - x and b x are exact.  Scored by refs/score_caller.py.
#include <stan/math/rev.hpp>
#include <cstdio>

int main() {
  using stan::math::var;
  std::printf("# a b x lcdf d_dx\n");
  const double a = 2.0;
  for (int k : {20, 30, 31}) {
    for (int t : {-1, 0, 1, 3}) {
      const double b = std::ldexp(1.0, k);
      const double x = std::ldexp(1.0, t - k);  // b x = 2^t
      var xv = x;
      var lp = stan::math::beta_lcdf(xv, a, b);
      lp.grad();
      std::printf("%a %a %a %a %a\n", a, b, x, lp.val(), xv.adj());
      stan::math::recover_memory();
    }
  }
  return 0;
}
