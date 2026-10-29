// probe_1578.cpp --- the case of stan-dev/math#1578: gradients of
// neg_binomial_lpmf(n | alpha, beta) for large alpha. Build against a tree.
#include <stan/math/rev.hpp>
#include <cstdio>
int main() {
  for (int n : {0, 3}) {
    for (double alpha : {1e10 + 1, 1e12, 1e15}) {
      stan::math::var a = alpha;
      stan::math::var b = 0.8;
      stan::math::var lp = stan::math::neg_binomial_lpmf(n, a, b);
      lp.grad();
      std::printf("n=%d alpha=%.17g value=%.17g d/dalpha=%.17g d/dbeta=%.17g\n",
                  n, alpha, lp.val(), a.adj(), b.adj());
      stan::math::recover_memory();
    }
  }
  return 0;
}
