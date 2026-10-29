// E4 check: fvar<var> second derivative of neg_binomial_2_log_glm_lpmf in
// theta at y = 1, phi = 0.07, theta = 35 and 27.
#include <stan/math/mix.hpp>
#include <cstdio>
#include <vector>

int main() {
  using stan::math::fvar;
  using stan::math::var;
  const std::vector<int> y{1};
  Eigen::MatrixXd x(1, 1);
  x << 1.0;
  Eigen::VectorXd beta(1);
  beta << 0.0;
  const double ref[] = {-4.7225324533500949e-17, -1.4077670835874031e-13};
  const double th[] = {35.0, 27.0};
  for (int i = 0; i < 2; ++i) {
    fvar<var> alpha(var(th[i]), 1.0);
    fvar<var> lp = stan::math::neg_binomial_2_log_glm_lpmf(y, x, alpha, beta,
                                                           0.07);
    lp.d_.grad();
    const double h = alpha.val_.adj();
    std::printf("theta %g h %.17g ref %.17g relerr %.3g\n", th[i], h, ref[i],
                std::fabs(h / ref[i] - 1));
    stan::math::recover_memory();
  }
  return 0;
}
