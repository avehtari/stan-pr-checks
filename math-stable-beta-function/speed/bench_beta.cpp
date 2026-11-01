// Speed of beta(double, double) in three regions (BENCH_FN=1), or of
// beta_lcdf with var arguments at b = 1e9 and its gradient (BENCH_FN=2).
// Build one binary per function and per tree; run_bench.sh interleaves them.
// Output: "TIME <name> <ns per call>", the minimum over REPEATS blocks.
#include <stan/math.hpp>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <random>
#include <utility>
#include <vector>

#ifndef BENCH_FN
#error "define BENCH_FN as 1 (beta) or 2 (beta_lcdf)"
#endif

static constexpr int N_IN = 4096;
static constexpr int REPEATS = 15;

// log-uniform on [lo, hi]
static double log_uniform(std::mt19937_64& g, double lo, double hi) {
  std::uniform_real_distribution<double> u(std::log(lo), std::log(hi));
  return std::exp(u(g));
}

template <typename F>
static double min_ns_per_call(F&& block, int calls_per_block) {
  double best = 1e300;
  double sink = 0;
  for (int r = 0; r < REPEATS; ++r) {
    const auto t0 = std::chrono::steady_clock::now();
    sink += block();
    const auto t1 = std::chrono::steady_clock::now();
    const double ns = std::chrono::duration<double, std::nano>(t1 - t0).count()
                      / calls_per_block;
    best = std::min(best, ns);
  }
  std::fprintf(stderr, "sink %g\n", sink);
  return best;
}

int main(int argc, char** argv) {
  const int inner = argc > 1 ? std::atoi(argv[1]) : 20;
  std::mt19937_64 g(20261101);
#if BENCH_FN == 1
  const double top = std::ldexp(1.0, 31);
  const char* names[3] = {"beta_small", "beta_mixed", "beta_large"};
  std::vector<std::pair<double, double>> in[3];
  for (int i = 0; i < N_IN; ++i) {
    in[0].emplace_back(log_uniform(g, 1e-3, 10), log_uniform(g, 1e-3, 10));
    const double s = log_uniform(g, 1e-3, 10);
    const double l = log_uniform(g, 10, top);
    in[1].emplace_back(i % 2 ? s : l, i % 2 ? l : s);
    in[2].emplace_back(log_uniform(g, 10, top), log_uniform(g, 10, top));
  }
  for (int k = 0; k < 3; ++k) {
    const auto& v = in[k];
    double acc = 0;
    auto block = [&]() {
      for (int rep = 0; rep < inner; ++rep) {
        for (int i = 0; i < N_IN; ++i) {
          const double bump = acc * 1e-300;
          acc += stan::math::beta(v[i].first + bump, v[i].second);
        }
      }
      return acc;
    };
    std::printf("TIME %s %.2f\n", names[k], min_ns_per_call(block, inner * N_IN));
  }
#else
  std::vector<std::array<double, 3>> in;
  for (int i = 0; i < N_IN; ++i) {
    const double b = 1e9;
    const double a = log_uniform(g, 0.5, 10);
    const double y = log_uniform(g, 0.1, 10) / b;
    in.push_back({y, a, b});
  }
  double acc = 0;
  auto block = [&]() {
    for (int rep = 0; rep < inner; ++rep) {
      for (int i = 0; i < N_IN; ++i) {
        stan::math::nested_rev_autodiff nested;
        const double bump = acc * 1e-300;
        stan::math::var y = in[i][0] + bump;
        stan::math::var a = in[i][1];
        stan::math::var b = in[i][2];
        stan::math::var lp = stan::math::beta_lcdf(y, a, b);
        lp.grad();
        acc += lp.val() + y.adj() * 1e-30 + a.adj() + b.adj();
      }
    }
    return acc;
  };
  std::printf("TIME beta_lcdf_var_b1e9 %.2f\n",
              min_ns_per_call(block, inner * N_IN));
#endif
  return 0;
}
