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

// Phi and std_normal_* take the standardized value, the others y.
template <typename F>
std::vector<double> inputs(int n, double width, std::mt19937_64& rng) {
  std::uniform_real_distribution<double> u(-width, width);
  const bool standard = std::string(F::name) == "Phi"
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
    acc += F::f(y, MU, SIGMA, LAMBDA);
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
    acc += lp.val() + mu.adj() + y[0].adj();
    stan::math::recover_memory();
  }
  const auto t1 = std::chrono::steady_clock::now();
  sink = acc;
  return std::chrono::duration<double, std::nano>(t1 - t0).count()
         / (static_cast<double>(reps) * ys.size());
}

template <typename F>
void run(const char* label, int rounds, double& sink) {
  for (double width : {3.0, 30.0}) {
    for (int n : {100, 1000, 10000, 100000}) {
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
  double sink = 0;
  run<Phi>(label, rounds, sink);
  run<std_normal_cdf>(label, rounds, sink);
  run<normal_cdf>(label, rounds, sink);
  run<std_normal_lcdf>(label, rounds, sink);
  run<normal_lcdf>(label, rounds, sink);
  run<lognormal_lcdf>(label, rounds, sink);
  run<lognormal_cdf>(label, rounds, sink);
  run<exp_mod_normal_lcdf>(label, rounds, sink);
  std::printf("sink %g\n", sink);
  return 0;
}
