// glm_propto.cpp --- neg_binomial_2_log_glm_lpmf with var alpha and beta and
// data phi, with propto = true and false (claim P4: "For the GLM with data
// phi, this brings back a rounding error of about eps * phi log(phi) in the
// propto value").
//
// Five observations y = (0, 1, 3, 10, 57), x = (1, -0.5, 2, 0.25, 1.5),
// alpha = 0.5, beta = 0.3, phi = 1e3, 1e4, ..., 1e15.
// Output line: <phi hex> => <propto value hex> <full value hex>
//              <d/dalpha propto hex> <d/dbeta propto hex>

#include <stan/math/rev.hpp>

#include <cstdio>
#include <vector>

using stan::math::var;

int main() {
  const std::vector<int> y{0, 1, 3, 10, 57};
  Eigen::MatrixXd x(5, 1);
  x << 1.0, -0.5, 2.0, 0.25, 1.5;
  for (int k = 3; k <= 15; ++k) {
    const double phi = std::pow(10.0, k);
    var alpha = 0.5;
    Eigen::Matrix<var, Eigen::Dynamic, 1> beta(1);
    beta(0) = 0.3;
    var lp_t
        = stan::math::neg_binomial_2_log_glm_lpmf<true>(y, x, alpha, beta, phi);
    lp_t.grad();
    const double v_t = lp_t.val();
    const double d_alpha = alpha.adj();
    const double d_beta = beta(0).adj();
    stan::math::recover_memory();
    var alpha2 = 0.5;
    Eigen::Matrix<var, Eigen::Dynamic, 1> beta2(1);
    beta2(0) = 0.3;
    var lp_f = stan::math::neg_binomial_2_log_glm_lpmf<false>(y, x, alpha2,
                                                               beta2, phi);
    std::printf("%a => %a %a %a %a\n", phi, v_t, lp_f.val(), d_alpha,
                d_beta);
    stan::math::recover_memory();
  }
  return 0;
}
