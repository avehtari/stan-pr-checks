// Probe for the saddle-point claims of stan-dev/math#3438.
// C: the final Newton iterate of solver 3 from (0.01, 3) in example C.
// D: a saddle point where the LU factor of B = I + K W has an odd row
//    permutation and positive pivots, so the old log determinant is finite.
#include <stan/math/mix.hpp>
#include <cmath>
#include <cstdio>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>

namespace {
struct cauchy_normal_likelihood {
  template <typename Theta>
  auto operator()(const Theta& theta, const Eigen::VectorXd& y, double s,
                  double z, double sz, std::ostream* pstream) const {
    return stan::math::cauchy_lpdf(y, theta(0), s)
           + stan::math::normal_lpdf(z, theta(1), sz);
  }
};

struct fixed_covariance {
  Eigen::MatrixXd operator()(const Eigen::MatrixXd& K,
                             std::ostream* pstream) const {
    return K;
  }
};

struct problem {
  const char* name;
  Eigen::MatrixXd K;
  Eigen::VectorXd y;
  double s, z, sz;
};

// W = -d^2/dtheta^2 log p(y, z | theta), analytic.
Eigen::VectorXd w_diag(const problem& p, const Eigen::VectorXd& theta) {
  Eigen::VectorXd w(2);
  w(0) = 0;
  for (int i = 0; i < p.y.size(); ++i) {
    const double r = p.y(i) - theta(0);
    const double q = p.s * p.s + r * r;
    w(0) += 2 * (p.s * p.s - r * r) / (q * q);
  }
  w(1) = 1 / (p.sz * p.sz);
  return w;
}

// Gradient of log p(theta | y) = log p(y, z | theta) - theta' K^-1 theta / 2.
Eigen::VectorXd post_grad(const problem& p, const Eigen::VectorXd& theta) {
  Eigen::VectorXd g(2);
  g(0) = 0;
  for (int i = 0; i < p.y.size(); ++i) {
    const double r = p.y(i) - theta(0);
    g(0) += 2 * r / (p.s * p.s + r * r);
  }
  g(1) = (p.z - theta(1)) / (p.sz * p.sz);
  return g - p.K.ldlt().solve(theta);
}

void describe(const problem& p, const Eigen::VectorXd& theta) {
  const Eigen::VectorXd w = w_diag(p, theta);
  const Eigen::MatrixXd B
      = Eigen::MatrixXd::Identity(2, 2) + p.K * w.asDiagonal();
  Eigen::PartialPivLU<Eigen::MatrixXd> lu(B);
  const Eigen::MatrixXd H = p.K.inverse() + Eigen::MatrixXd(w.asDiagonal());
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(H);
  std::printf("    theta (%.6g, %.6g), |grad log p(theta | y)| %.3g\n",
              theta(0), theta(1), post_grad(p, theta).norm());
  std::printf("    W (%.6g, %.6g), eigenvalues of K^-1 + W (%.6g, %.6g)\n",
              w(0), w(1), es.eigenvalues()(0), es.eigenvalues()(1));
  std::printf("    det(B) %.6g, permutation sign %+.0f, LU diagonal (%.6g, "
              "%.6g)\n",
              B.determinant(),
              static_cast<double>(lu.permutationP().determinant()),
              lu.matrixLU()(0, 0), lu.matrixLU()(1, 1));
}

void run(const problem& p, const Eigen::VectorXd& theta0, double tol = 1e-10,
         int steps = 100) {
  std::printf("%s, start (%g, %g), solver 3, tolerance %g, %d steps\n", p.name,
              theta0(0), theta0(1), tol, steps);
  auto ops = std::make_tuple(theta0, tol, steps, 3, steps, false);
  auto opts = stan::math::internal::tuple_to_laplace_options(ops);
  try {
    auto md = stan::math::internal::laplace_marginal_density_est(
        cauchy_normal_likelihood{}, std::forward_as_tuple(p.y, p.s, p.z, p.sz),
        p.K, opts, nullptr);
    std::printf("  laplace_marginal_density_est: lmd %.10g, solver %d\n",
                md.lmd, md.solver_used);
    describe(p, md.theta);
    if (std::string(p.name) == "D" && std::abs(md.theta(0)) > 0.1) {
      // refs/ref_example_d.py; the two maxima are (m0, m1) and (-m0, -m1).
      const double sign = md.theta(0) < 0 ? 1.0 : -1.0;
      const double m0 = sign * -0.35038406125098763;
      const double m1 = sign * 0.0029862277947527355;
      std::printf("    error vs mpmath: value %.3g, mean (%.3g, %.3g)\n",
                  std::abs(md.lmd - -4.8472237608447798),
                  std::abs(md.theta(0) - m0), std::abs(md.theta(1) - m1));
    }
  } catch (const std::exception& e) {
    std::printf("  laplace_marginal_density_est throws: %s\n", e.what());
  }
  try {
    const double lml = stan::math::laplace_marginal_tol<false>(
        cauchy_normal_likelihood{}, std::forward_as_tuple(p.y, p.s, p.z, p.sz),
        2, fixed_covariance{}, std::make_tuple(p.K), ops, nullptr);
    std::printf("  laplace_marginal_tol: %.10g\n", lml);
  } catch (const std::exception& e) {
    std::printf("  laplace_marginal_tol throws: %s\n", e.what());
  }
}
}  // namespace

int main() {
  problem c{"C", Eigen::MatrixXd(2, 2), Eigen::VectorXd(2), 0.3, 0.0, 0.1};
  c.K << 4.0, 1.0, 1.0, 1.0;
  c.y << -3.0, 3.0;
  // r^2 = 3 s^2 gives the most negative Cauchy curvature, -1 / (4 s^2) per
  // observation.  With K(0, 1) < 0 the row swap leaves both pivots positive.
  problem d{"D", Eigen::MatrixXd(2, 2), Eigen::VectorXd(2), 0.3, 0.0, 0.1};
  d.K << 1.0, -1.5, -1.5, 4.0;
  d.y << -0.3 * std::sqrt(3.0), 0.3 * std::sqrt(3.0);

  Eigen::VectorXd t(2);
  t << 0.01, 3.0;
  run(c, t);
  t << 0.3, 3.0;
  run(c, t);
  for (const auto& st : std::vector<std::pair<double, double>>{
           {0.0, 0.0}, {0.01, 0.01}, {0.05, 1.0}}) {
    t << st.first, st.second;
    run(d, t);
  }
  // The settings of the unit test solver_3_odd_permutation_saddle_point_throws.
  for (const auto& st : std::vector<std::pair<double, double>>{
           {0.0, 0.0}, {-0.3, 0.0}, {0.3, 0.0}, {0.01, 0.01}}) {
    t << st.first, st.second;
    run(d, t, 1e-12, 1000);
  }
  return 0;
}
