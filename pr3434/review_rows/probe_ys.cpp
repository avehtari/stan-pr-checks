// Candidate rows for E1 (#1578 point) and E2 (yule_simon_cdf, _lcdf at
// large alpha). Prints value and log-scale gradients in hex.
#include <stan/math/rev.hpp>
#include <cstdio>

int main() {
  using stan::math::var;
  const int ns[] = {1, 2, 5, 20};
  const double alphas[] = {1e4, 1e6, 1e8, 1e10, 1e12, 1e15};
  for (int tag = 0; tag < 2; ++tag) {
    for (int n : ns) {
      for (double a : alphas) {
        var alpha = a;
        var lp = tag == 0 ? stan::math::yule_simon_cdf(n, alpha)
                          : stan::math::yule_simon_lcdf(n, alpha);
        lp.grad();
        std::printf("%s %d %a %a %a\n", tag == 0 ? "YSC" : "YSLC", n, a,
                    lp.val(), a * alpha.adj());
        stan::math::recover_memory();
      }
    }
  }
  var alpha = 1e10 + 1, beta = 0.8;
  var lp = stan::math::neg_binomial_lpmf(0, alpha, beta);
  lp.grad();
  std::printf("NB 0 %a %a %a %a %a\n", alpha.val(), beta.val(), lp.val(),
              alpha.adj(), beta.adj());
  return 0;
}
