// Grid over mu and alpha of the math#3469 example: print the value or
// "throws", to find points where develop gives NaN and math#3438 a number.
#include <stan/math/mix.hpp>
#include <cmath>
#include <cstdio>
#include <sstream>
#include <vector>

struct constant_cov {
  template <typename T>
  Eigen::Matrix<T, -1, -1> operator()(const T& alpha, int n,
                                      std::ostream*) const {
    return stan::math::add_diag(
        Eigen::Matrix<T, -1, -1>::Constant(n, n, alpha * alpha), 1e-8);
  }
};

int main() {
  const std::vector<int> y{3, 1, 4, 1, 5}, idx{1, 2, 3, 4, 5};
  for (double mu : {-5.0, 0.0, 1.0, 2.0, 5.0}) {
    for (double la = 7.0; la <= 12.01; la += 0.25) {
      double alpha = std::pow(10.0, la);
      const Eigen::VectorXd mean = Eigen::VectorXd::Constant(5, mu);
      std::stringstream msgs;
      std::printf("%5.1f %6.2f ", mu, la);
      try {
        double v = stan::math::laplace_marginal_poisson_log_lpmf(
            y, idx, mean, 1, constant_cov{}, std::forward_as_tuple(alpha, 5),
            &msgs);
        std::printf("%.10g\n", v);
      } catch (const std::exception& e) {
        std::printf("throws\n");
      }
    }
  }
}
