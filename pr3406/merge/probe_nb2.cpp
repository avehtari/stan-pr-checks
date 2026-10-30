// Probe: neg_binomial_2_lpmf(int, var, std::vector<var>) after the merge of
// develop; run with AddressSanitizer to find a dangling reference.
#include <stan/math/rev.hpp>
#include <cstdio>
#include <vector>

int main() {
  using stan::math::var;
  var mu(1.5);
  std::vector<var> phi{2.0, 4.0, 6.0, 8.0};
  var val = stan::math::neg_binomial_2_lpmf(100, mu, phi);
  val.grad();
  std::printf("val %.17g dmu %.17g dphi0 %.17g\n", val.val(), mu.adj(),
              phi[0].adj());
  std::vector<double> phi_d{2.0, 4.0, 6.0, 8.0};
  std::printf("double %.17g\n", stan::math::neg_binomial_2_lpmf(100, 1.5, phi_d));
  return 0;
}
