// Reproduction of the timing table in the #3406 PR text: vector arguments of
// size N, "double" (Eigen::VectorXd y, double parameters, value only) and
// "var" (Eigen::Matrix<var> y, var parameters, value and gradient), ns per
// element.  Inputs: z uniform on [-3, 3] ("bulk") or [-30, 30] ("wide"),
// y = mu + sigma z (lognormal: exp of that).  Compile once per tree.
// Usage: bench_pr_table <label> <rounds>
#include <stan/math/rev.hpp>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <random>
#include <string>
#include <vector>

using stan::math::var;
using VecV = Eigen::Matrix<var, Eigen::Dynamic, 1>;
constexpr double MU = 0.3, SIGMA = 1.7, LAMBDA = 1.0;

#define FUN(NAME, CALL, LOG_Y)                                       \
  struct NAME {                                                      \
    static constexpr const char* name = #NAME;                       \
    static constexpr bool log_y = LOG_Y;                             \
    template <typename Y, typename P>                                \
    static auto f(const Y& y, const P& mu, const P& sigma,           \
                  const P& lambda) {                                 \
      (void)mu; (void)sigma; (void)lambda;                           \
      return CALL;                                                   \
    }                                                                \
  };

FUN(Phi, stan::math::sum(stan::math::Phi(y)), false)
FUN(std_normal_cdf, stan::math::std_normal_cdf(y), false)
FUN(normal_cdf, stan::math::normal_cdf(y, mu, sigma), false)
FUN(std_normal_lcdf, stan::math::std_normal_lcdf(y), false)
FUN(normal_lcdf, stan::math::normal_lcdf(y, mu, sigma), false)
FUN(lognormal_lcdf, stan::math::lognormal_lcdf(y, mu, sigma), true)
FUN(lognormal_cdf, stan::math::lognormal_cdf(y, mu, sigma), true)
FUN(exp_mod_normal_lcdf, stan::math::exp_mod_normal_lcdf(y, mu, sigma, lambda),
    false)
FUN(exp_mod_normal_lpdf, stan::math::exp_mod_normal_lpdf(y, mu, sigma, lambda),
    false)
FUN(exp_mod_normal_lccdf,
    stan::math::exp_mod_normal_lccdf(y, mu, sigma, lambda), false)
FUN(exp_mod_normal_cdf, stan::math::exp_mod_normal_cdf(y, mu, sigma, lambda),
    false)
FUN(skew_normal_lpdf, stan::math::skew_normal_lpdf(y, mu, sigma, lambda), false)
FUN(Phi_approx, stan::math::sum(stan::math::Phi_approx(y)), false)
FUN(Phi_approx_log, stan::math::sum(stan::math::log(stan::math::Phi_approx(y))),
    false)
// log Phi_approx(y) without underflow: log_inv_logit of the cubic.
FUN(Phi_approx_lil,
    stan::math::sum(stan::math::log_inv_logit(
        0.07056 * stan::math::elt_multiply(stan::math::elt_multiply(y, y), y)
        + 1.5976 * y)),
    false)
FUN(Phi_log, stan::math::sum(stan::math::log(stan::math::Phi(y))), false)

// Phi_approx with x * x * x instead of std::pow(x, 3.0) in the double path
// (the var path already multiplies); Log: log of it.
template <bool Log, typename Y>
auto phi_approx_cube_sum(const Y& y) {
  if constexpr (std::is_same_v<stan::scalar_type_t<Y>, double>) {
    double s = 0;
    for (Eigen::Index i = 0; i < y.size(); ++i) {
      const double x = y[i];
      const double p = stan::math::inv_logit(0.07056 * x * x * x + 1.5976 * x);
      s += Log ? std::log(p) : p;
    }
    return s;
  } else if constexpr (Log) {
    return stan::math::sum(stan::math::log(stan::math::Phi_approx(y)));
  } else {
    return stan::math::sum(stan::math::Phi_approx(y));
  }
}
FUN(Phi_approx_cube, phi_approx_cube_sum<false>(y), false)
FUN(Phi_approx_cube_log, phi_approx_cube_sum<true>(y), false)

// Phi and std_normal_* take the standardized value, the others y.
template <typename F>
std::vector<double> inputs(int n, double width, std::mt19937_64& rng) {
  std::uniform_real_distribution<double> u(-width, width);
  const bool standard = std::string(F::name).rfind("Phi", 0) == 0
                        || std::string(F::name).rfind("std_normal", 0) == 0;
  std::vector<double> ys(n);
  for (auto& v : ys) {
    const double z = u(rng);
    v = standard ? z : MU + SIGMA * z;
    if (F::log_y) {
      v = std::exp(v);
    }
  }
  return ys;
}

template <typename F>
double time_double(const std::vector<double>& ys, int reps, double& sink) {
  Eigen::VectorXd y(ys.size());
  double acc = sink;
  const auto t0 = std::chrono::steady_clock::now();
  for (int r = 0; r < reps; ++r) {
    const double bump = acc * 1e-300;
    for (size_t i = 0; i < ys.size(); ++i) {
      y[i] = ys[i] + bump;
    }
    // A non-finite result must not reach acc: the bump would then turn all
    // later inputs into inf or NaN.
    const double res = F::f(y, MU, SIGMA, LAMBDA);
    acc += std::isfinite(res) ? res : 1.0;
  }
  const auto t1 = std::chrono::steady_clock::now();
  sink = acc;
  return std::chrono::duration<double, std::nano>(t1 - t0).count()
         / (static_cast<double>(reps) * ys.size());
}

template <typename F>
double time_var(const std::vector<double>& ys, int reps, double& sink) {
  double acc = sink;
  const auto t0 = std::chrono::steady_clock::now();
  for (int r = 0; r < reps; ++r) {
    const double bump = acc * 1e-300;
    VecV y(ys.size());
    for (size_t i = 0; i < ys.size(); ++i) {
      y[i] = ys[i] + bump;
    }
    var mu = MU, sigma = SIGMA, lambda = LAMBDA;
    var lp = F::f(y, mu, sigma, lambda);
    lp.grad();
    const double res = lp.val() + mu.adj() + y[0].adj();
    acc += std::isfinite(res) ? res : 1.0;
    stan::math::recover_memory();
  }
  const auto t1 = std::chrono::steady_clock::now();
  sink = acc;
  return std::chrono::duration<double, std::nano>(t1 - t0).count()
         / (static_cast<double>(reps) * ys.size());
}

// Input ranges (z on [-w, w]) and sizes; set per mode in main.
std::vector<double> widths{3.0, 30.0};
std::vector<int> sizes{100, 1000, 10000, 100000};

template <typename F>
void run(const char* label, int rounds, double& sink) {
  for (double width : widths) {
    for (int n : sizes) {
      std::mt19937_64 rng(20261020 + n);
      const auto ys = inputs<F>(n, width, rng);
      const int reps = std::max(1, 2000000 / n);
      for (int r = 0; r < rounds; ++r) {
        const double td = time_double<F>(ys, reps, sink);
        const double tv = time_var<F>(ys, std::max(1, reps / 2), sink);
        std::printf("RESULT %s %s %s %d %d double %.3f var %.3f\n", label,
                    F::name, width < 10 ? "bulk" : "wide", n, r + 1, td, tv);
      }
    }
  }
}

int main(int argc, char** argv) {
  const char* label = argc > 1 ? argv[1] : "?";
  const int rounds = argc > 2 ? std::atoi(argv[2]) : 2;
  const std::string mode = argc > 3 ? argv[3] : "full";
  double sink = 0;
  if (mode == "emg") {
    widths = {3.0, 30.0};
    sizes = {1000};
    run<exp_mod_normal_lcdf>(label, rounds, sink);
    run<exp_mod_normal_lccdf>(label, rounds, sink);
    run<exp_mod_normal_cdf>(label, rounds, sink);
    std::printf("sink %g\n", sink);
    return 0;
  }
  if (mode == "approx") {
    widths = {3.0, 30.0};
    sizes = {1000};
    run<Phi>(label, rounds, sink);
    run<Phi_approx>(label, rounds, sink);
    run<Phi_log>(label, rounds, sink);
    run<Phi_approx_log>(label, rounds, sink);
    run<Phi_approx_lil>(label, rounds, sink);
    run<Phi_approx_cube>(label, rounds, sink);
    run<Phi_approx_cube_log>(label, rounds, sink);
    run<std_normal_lcdf>(label, rounds, sink);
    std::printf("sink %g\n", sink);
    return 0;
  }
  if (mode == "flags") {
    widths = {3.0, 30.0};
    sizes = {1000};
    run<Phi>(label, rounds, sink);
    run<std_normal_cdf>(label, rounds, sink);
    run<normal_cdf>(label, rounds, sink);
    run<std_normal_lcdf>(label, rounds, sink);
    run<normal_lcdf>(label, rounds, sink);
    run<lognormal_lcdf>(label, rounds, sink);
    run<lognormal_cdf>(label, rounds, sink);
    run<exp_mod_normal_lcdf>(label, rounds, sink);
    run<exp_mod_normal_lccdf>(label, rounds, sink);
    run<exp_mod_normal_cdf>(label, rounds, sink);
    run<exp_mod_normal_lpdf>(label, rounds, sink);
    run<skew_normal_lpdf>(label, rounds, sink);
    std::printf("sink %g\n", sink);
    return 0;
  }
  if (mode == "quick") {
    widths = {3.0};
    sizes = {1000};
    run<std_normal_lcdf>(label, rounds, sink);
    run<normal_lcdf>(label, rounds, sink);
    run<lognormal_lcdf>(label, rounds, sink);
    run<exp_mod_normal_lcdf>(label, rounds, sink);
    run<exp_mod_normal_lpdf>(label, rounds, sink);
    run<skew_normal_lpdf>(label, rounds, sink);
    std::printf("sink %g\n", sink);
    return 0;
  }
  run<Phi>(label, rounds, sink);
  run<std_normal_cdf>(label, rounds, sink);
  run<normal_cdf>(label, rounds, sink);
  run<std_normal_lcdf>(label, rounds, sink);
  run<normal_lcdf>(label, rounds, sink);
  run<lognormal_lcdf>(label, rounds, sink);
  run<lognormal_cdf>(label, rounds, sink);
  run<exp_mod_normal_lcdf>(label, rounds, sink);
  run<exp_mod_normal_lpdf>(label, rounds, sink);
  run<skew_normal_lpdf>(label, rounds, sink);
  std::printf("sink %g\n", sink);
  return 0;
}
