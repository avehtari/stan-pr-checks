// Bundle: var d/dalpha of gamma_lccdf on the (alpha, y) pairs of stdin.
// Supports the 266-cell grid claim (3.4e-12) of #3408 Part 3.
// d/dalpha gamma_lccdf(y | alpha, 1) with var, for the (alpha, y) pairs on
// stdin; prints alpha, y and the derivative in hex.
#include <stan/math/rev.hpp>
#include <cstdio>

int main() {
  double a, y;
  while (std::scanf("%lf %lf", &a, &y) == 2) {
    stan::math::var av = a;
    stan::math::var lp = stan::math::gamma_lccdf(y, av, 1.0);
    lp.grad();
    std::printf("%a %a %a\n", a, y, av.adj());
    stan::math::recover_memory();
  }
  return 0;
}
