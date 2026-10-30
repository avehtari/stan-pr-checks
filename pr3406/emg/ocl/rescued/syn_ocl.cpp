#include <stan/math.hpp>
#include <vector>

int main() {
  using stan::math::matrix_cl;
  using stan::math::var;
  using stan::math::var_value;
  Eigen::VectorXd y(3);
  y << -1.0, 0.5, 2.0;
  matrix_cl<double> y_cl(y);
  double a = stan::math::exp_mod_normal_lcdf(y_cl, 0.1, 1.3, 0.7);
  a += stan::math::exp_mod_normal_lccdf(y_cl, y_cl, y_cl, y_cl);
  a += stan::math::exp_mod_normal_cdf(y_cl, 0.1, y_cl, 0.7);
  var_value<matrix_cl<double>> yv(y_cl);
  var mu = 0.1;
  var lp = stan::math::exp_mod_normal_lcdf(yv, mu, 1.3, yv)
           + stan::math::exp_mod_normal_lccdf(yv, mu, yv, 0.7)
           + stan::math::exp_mod_normal_cdf(yv, yv, yv, yv);
  return static_cast<int>(a + lp.val());
}
