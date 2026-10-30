// GPU time of exp_mod_normal_lcdf, _lccdf and _cdf through the Stan OpenCL
// API: y on the device as matrix_cl (double, value only) or
// var_value<matrix_cl> (var, with mu, sigma and lambda as var scalars, value
// and gradient).  z = (y - mu) / sigma on [-3, 3] ("bulk") or [-30, 30]
// ("wide").  ns per element; the minimum over the rounds is the result.
// Usage: bench_ocl_emg <label> <rounds> <N>
#include <stan/math/opencl/rev.hpp>
#include <stan/math.hpp>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>

namespace sm = stan::math;
using sm::matrix_cl;
using sm::var;
using sm::var_value;

constexpr double MU = 0.3, SIGMA = 1.7, LAMBDA = 1.0;
constexpr int REPS = 50;

template <int Kind, typename Y, typename P>
auto call(const Y& y, const P& mu, const P& sigma, const P& lambda) {
  if constexpr (Kind == 0) {
    return sm::exp_mod_normal_lcdf(y, mu, sigma, lambda);
  } else if constexpr (Kind == 1) {
    return sm::exp_mod_normal_lccdf(y, mu, sigma, lambda);
  } else {
    return sm::exp_mod_normal_cdf(y, mu, sigma, lambda);
  }
}

template <int Kind>
double time_double(const matrix_cl<double>& y_cl, int n, double& sink) {
  double acc = sink + call<Kind>(y_cl, MU, SIGMA, LAMBDA);  // kernel build
  const auto t0 = std::chrono::steady_clock::now();
  for (int r = 0; r < REPS; ++r) {
    const double res = call<Kind>(y_cl, MU, SIGMA, LAMBDA);
    acc += std::isfinite(res) ? res : 1.0;
  }
  const auto t1 = std::chrono::steady_clock::now();
  sink = acc;
  return std::chrono::duration<double, std::nano>(t1 - t0).count()
         / (static_cast<double>(REPS) * n);
}

template <int Kind>
double time_var(const matrix_cl<double>& y_cl, int n, double& sink) {
  double acc = sink;
  auto once = [&]() {
    var_value<matrix_cl<double>> y(y_cl);
    var mu = MU, sigma = SIGMA, lambda = LAMBDA;
    var lp = call<Kind>(y, mu, sigma, lambda);
    lp.grad();
    const double res = lp.val() + mu.adj() + sigma.adj() + lambda.adj();
    acc += std::isfinite(res) ? res : 1.0;
    sm::recover_memory();
  };
  once();  // kernel build
  const auto t0 = std::chrono::steady_clock::now();
  for (int r = 0; r < REPS; ++r) {
    once();
  }
  const auto t1 = std::chrono::steady_clock::now();
  sink = acc;
  return std::chrono::duration<double, std::nano>(t1 - t0).count()
         / (static_cast<double>(REPS) * n);
}

template <int Kind>
void run(const char* label, int rounds, int n, double& sink) {
  static const char* names[] = {"exp_mod_normal_lcdf", "exp_mod_normal_lccdf",
                                "exp_mod_normal_cdf"};
  for (double width : {3.0, 30.0}) {
    Eigen::VectorXd y(n);
    for (int i = 0; i < n; ++i) {
      y[i] = MU + SIGMA * width * std::sin(1.0 + 7.0 * i);
    }
    matrix_cl<double> y_cl(y);
    for (int r = 0; r < rounds; ++r) {
      const double td = time_double<Kind>(y_cl, n, sink);
      const double tv = time_var<Kind>(y_cl, n, sink);
      std::printf("RESULT %s %s %s %d %d double %.4f var %.4f\n", label,
                  names[Kind], width < 10 ? "bulk" : "wide", n, r + 1, td, tv);
      std::fflush(stdout);
    }
  }
}

int main(int argc, char** argv) {
  const char* label = argc > 1 ? argv[1] : "?";
  const int rounds = argc > 2 ? std::atoi(argv[2]) : 3;
  const int n = argc > 3 ? std::atoi(argv[3]) : 1000000;
  double sink = 0;
  run<0>(label, rounds, n, sink);
  run<1>(label, rounds, n, sink);
  run<2>(label, rounds, n, sink);
  std::printf("sink %g\n", sink);
  return 0;
}
