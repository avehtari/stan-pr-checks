// Time of one evaluation of normal_lpdf (create the autodiff arguments,
// call, grad(), recover_memory()) for one type case of cases.hpp.
// One binary per (tree, case, propto):
//   -DTREE=0 develop copy (normal_lpdf_dev), -DTREE=1 branch
//   -DCASE=<index into nc::type_cases>  -DPROPTO=0|1
// Output: RESULT <tree> <case> <propto> <N> <median_ns> <min_ns> <reps>

#include "cases.hpp"
#include <algorithm>
#include <chrono>
#include <cstdio>
#include <random>
#include <vector>

#if !defined(TREE) || !defined(CASE) || !defined(PROPTO)
#error "define TREE, CASE and PROPTO"
#endif

using Eigen::VectorXd;
using stan::math::var;

constexpr nc::type_case tc = nc::type_cases[CASE];

struct data_t {
  VectorXd y, mu, sigma;
  explicit data_t(int N) : y(N), mu(N), sigma(N) {
    std::mt19937 rng(7919u * N);
    std::normal_distribution<double> z(0, 1);
    for (int i = 0; i < N; ++i) {
      mu[i] = 0.3 + z(rng);
      y[i] = mu[i] + 0.7 * z(rng);
      sigma[i] = 0.7 * std::exp(0.2 * z(rng));
    }
  }
};

volatile double sink = 0;

template <typename T>
double finish(T& lp) {
  if constexpr (std::is_same_v<T, var>) {
    lp.grad();
    return lp.val();
  } else {
    return lp;
  }
}

static double one(const data_t& d) {
  double v = 0;
  {
    auto y = nc::make<tc.y, var>(d.y);
    auto mu = nc::make<tc.mu, var>(d.mu);
    auto s = nc::make<tc.sigma, var>(d.sigma);
#if TREE == 0
    auto lp = stan::math::normal_lpdf_dev<PROPTO>(y, mu, s);
#else
    auto lp = stan::math::normal_lpdf<PROPTO>(y, mu, s);
#endif
    v = finish(lp);
  }
  stan::math::recover_memory();
  return v;
}

static double time_batch(const data_t& d, int reps) {
  const auto t0 = std::chrono::steady_clock::now();
  double acc = 0;
  for (int r = 0; r < reps; ++r) {
    acc += one(d);
  }
  const auto t1 = std::chrono::steady_clock::now();
  sink = sink + acc;
  return std::chrono::duration<double, std::nano>(t1 - t0).count() / reps;
}

int main() {
  for (int N : {100, 1000, 10000, 100000}) {
    data_t d(N);
    double t = time_batch(d, 3);
    int reps = std::max(1, static_cast<int>(2e7 / t));
    t = time_batch(d, reps);
    reps = std::max(1, static_cast<int>(2e7 / t));
    std::vector<double> ts;
    for (int b = 0; b < 7; ++b) {
      ts.push_back(time_batch(d, reps));
    }
    std::sort(ts.begin(), ts.end());
    std::printf("RESULT %s %s %d %d %.1f %.1f %d\n", TREE == 0 ? "dev" : "br",
                tc.name, PROPTO, N, ts[3], ts[0], reps);
    std::fflush(stdout);
  }
  std::printf("SINK %g\n", static_cast<double>(sink));
  return 0;
}
