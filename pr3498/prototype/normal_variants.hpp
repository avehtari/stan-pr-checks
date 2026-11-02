#ifndef NORMAL_LPDF_SPEED_WORK_NORMAL_VARIANTS_HPP
#define NORMAL_LPDF_SPEED_WORK_NORMAL_VARIANTS_HPP

// Variants of normal_lpdf<propto = true> for one type case:
// y data vector, mu var_value<VectorXd>, sigma scalar var.
//   0 dev           stan::math::normal_lpdf<true>
//   1 f1            input checks only when the sum is not finite
//   2 f2            f1, no y_scaled_sq array; sigma partial from the sum
//   3 f3            f2, mu partial assigned from an expression
//   4 f4            no temporary arrays; work in the arena partial of mu
//   5 compute_only  f3 arithmetic, no partials (returns a constant var)
// f1-f4 keep the checks and their order: on any invalid input they run the
// develop checks in the develop order, so the messages are the same.

#include <stan/math/rev.hpp>
#include <cmath>

namespace nv {

using Eigen::ArrayXd;
using Eigen::VectorXd;
using stan::math::var;
using vm = stan::math::var_value<VectorXd>;

constexpr const char* function = "normal_lpdf";

// The develop checks, in the develop order.
inline void dev_checks(const VectorXd& y, const VectorXd& mu_val,
                       double sigma_val) {
  stan::math::check_not_nan(function, "Random variable", y.array());
  stan::math::check_finite(function, "Location parameter", mu_val.array());
  stan::math::check_positive(function, "Scale parameter", sigma_val);
}

inline void front_checks(const VectorXd& y, const vm& mu, const var& sigma) {
  stan::math::check_consistent_sizes(function, "Random variable", y,
                                     "Location parameter", mu,
                                     "Scale parameter", sigma);
  const double s = sigma.val();
  if (!(s > 0)) {
    // throws; the develop order gives the develop message
    dev_checks(y, mu.val(), s);
  }
}

inline var dev(const VectorXd& y, const vm& mu, const var& sigma) {
  return stan::math::normal_lpdf<true>(y, mu, sigma);
}

inline var f1(const VectorXd& y, const vm& mu, const var& sigma) {
  using stan::math::partials;
  front_checks(y, mu, sigma);
  if (y.size() == 0) {
    return 0;
  }
  const double sigma_val = sigma.val();
  auto ops = stan::math::make_partials_propagator(y, mu, sigma);
  const double inv_sigma = 1.0 / sigma_val;
  ArrayXd y_scaled = (y.array() - mu.val().array()) * inv_sigma;
  ArrayXd y_scaled_sq = y_scaled * y_scaled;
  const double sum_sq = y_scaled_sq.sum();
  if (!std::isfinite(sum_sq)) {
    dev_checks(y, mu.val(), sigma_val);
  }
  const double N = static_cast<double>(y.size());
  const double logp = -0.5 * sum_sq - std::log(sigma_val) * N;
  ArrayXd scaled_diff = inv_sigma * y_scaled;
  partials<2>(ops) = inv_sigma * y_scaled_sq - inv_sigma;
  partials<1>(ops) = std::move(scaled_diff);
  return ops.build(logp);
}

inline var f2(const VectorXd& y, const vm& mu, const var& sigma) {
  using stan::math::partials;
  front_checks(y, mu, sigma);
  if (y.size() == 0) {
    return 0;
  }
  const double sigma_val = sigma.val();
  auto ops = stan::math::make_partials_propagator(y, mu, sigma);
  const double inv_sigma = 1.0 / sigma_val;
  ArrayXd y_scaled = (y.array() - mu.val().array()) * inv_sigma;
  const double sum_sq = y_scaled.square().sum();
  if (!std::isfinite(sum_sq)) {
    dev_checks(y, mu.val(), sigma_val);
  }
  const double N = static_cast<double>(y.size());
  const double logp = -0.5 * sum_sq - std::log(sigma_val) * N;
  ArrayXd scaled_diff = inv_sigma * y_scaled;
  partials<2>(ops)[0] = inv_sigma * (sum_sq - N);
  partials<1>(ops) = std::move(scaled_diff);
  return ops.build(logp);
}

inline var f3(const VectorXd& y, const vm& mu, const var& sigma) {
  using stan::math::partials;
  front_checks(y, mu, sigma);
  if (y.size() == 0) {
    return 0;
  }
  const double sigma_val = sigma.val();
  auto ops = stan::math::make_partials_propagator(y, mu, sigma);
  const double inv_sigma = 1.0 / sigma_val;
  ArrayXd y_scaled = (y.array() - mu.val().array()) * inv_sigma;
  const double sum_sq = y_scaled.square().sum();
  if (!std::isfinite(sum_sq)) {
    dev_checks(y, mu.val(), sigma_val);
  }
  const double N = static_cast<double>(y.size());
  const double logp = -0.5 * sum_sq - std::log(sigma_val) * N;
  partials<2>(ops)[0] = inv_sigma * (sum_sq - N);
  partials<1>(ops) = (inv_sigma * y_scaled).matrix();
  return ops.build(logp);
}

inline var f4(const VectorXd& y, const vm& mu, const var& sigma) {
  using stan::math::edge;
  using stan::math::partials;
  front_checks(y, mu, sigma);
  if (y.size() == 0) {
    return 0;
  }
  const double sigma_val = sigma.val();
  auto ops = stan::math::make_partials_propagator(y, mu, sigma);
  const double inv_sigma = 1.0 / sigma_val;
  auto& p = edge<1>(ops).partials_;
  p = (y - mu.val()) * inv_sigma;  // y_scaled
  const double sum_sq = p.squaredNorm();
  if (!std::isfinite(sum_sq)) {
    dev_checks(y, mu.val(), sigma_val);
  }
  p *= inv_sigma;
  const double N = static_cast<double>(y.size());
  const double logp = -0.5 * sum_sq - std::log(sigma_val) * N;
  partials<2>(ops)[0] = inv_sigma * (sum_sq - N);
  return ops.build(logp);
}

inline var compute_only(const VectorXd& y, const vm& mu, const var& sigma) {
  front_checks(y, mu, sigma);
  const double sigma_val = sigma.val();
  const double inv_sigma = 1.0 / sigma_val;
  ArrayXd y_scaled = (y.array() - mu.val().array()) * inv_sigma;
  const double sum_sq = y_scaled.square().sum();
  const double N = static_cast<double>(y.size());
  VectorXd d = (inv_sigma * y_scaled).matrix();
  const double ds = inv_sigma * (sum_sq - N);
  return var(-0.5 * sum_sq - std::log(sigma_val) * N + 0 * d(0) * ds);
}

}  // namespace nv
#endif
