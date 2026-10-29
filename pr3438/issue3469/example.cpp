#include <stan/math/mix.hpp>
#include <iomanip>
#include <iostream>
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
  std::cout << std::setprecision(10);
  for (double mu : {-30.0, 1.0}) {
    for (double alpha : {1e3, 1e5, 1e7, 1e9}) {
      const Eigen::VectorXd mean = Eigen::VectorXd::Constant(5, mu);
      std::stringstream msgs;
      std::cout << "mu " << mu << " alpha " << alpha << ": ";
      try {
        std::cout << stan::math::laplace_marginal_poisson_log_lpmf(
            y, idx, mean, 1, constant_cov{}, std::forward_as_tuple(alpha, 5),
            &msgs);
      } catch (const std::exception& e) {
        std::cout << "throws " << e.what();
      }
      std::cout << "\n" << msgs.str();
    }
  }
}
