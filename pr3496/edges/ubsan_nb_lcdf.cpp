// Build with -fsanitize=undefined: on develop, neg_binomial_lcdf computes
// digamma(n + 1) in int for every element, so n = INT_MAX overflows.
#include <stan/math/rev.hpp>
#include <climits>
#include <cstdio>
#include <vector>

int main() {
  stan::math::var alpha = 2.5, beta = 1.5;
  stan::math::var lp
      = stan::math::neg_binomial_lcdf(std::vector<int>{3, INT_MAX}, alpha, beta);
  lp.grad();
  std::printf("neg_binomial_lcdf({3, INT_MAX} | 2.5, 1.5) = %.6g, d/dalpha %.6g\n",
              lp.val(), alpha.adj());
  return 0;
}
