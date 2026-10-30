// PR #3406, claim S4 (commit 2015557859): scalar exp_mod_normal_lpdf with
// z = (y - mu - lambda sigma^2) / sigma near 9.4, the argument of the
// std_normal lcdf part.  Develop computes log(erfc(-z / sqrt(2))).
// "double": all four arguments double, value only.
// "var": all four arguments var, value and gradient (includes grad() and
// recover_memory() per call).
// Throughput loop over a block of 1000 inputs with the accumulator bump
// (STAN_MATH_DISTRIBUTIONS.md section 4).  Build once per tree.
// Usage: bench_emg_scalar <label>
#include <stan/math/rev.hpp>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <random>
#include <vector>

using stan::math::var;
constexpr double MU = 0.0, SIGMA = 1.0, LAMBDA = 1.0;
constexpr int BLOCK = 1000;

// y = MU + LAMBDA SIGMA^2 + SIGMA z, with z uniform on [lo, hi].
std::vector<double> inputs(double lo, double hi, unsigned seed) {
  std::mt19937_64 rng(seed);
  std::uniform_real_distribution<double> u(lo, hi);
  std::vector<double> ys(BLOCK);
  for (auto& y : ys) {
    y = MU + LAMBDA * SIGMA * SIGMA + SIGMA * u(rng);
  }
  return ys;
}

double time_double(const std::vector<double>& ys, int reps, double& sink) {
  double acc = sink;
  const auto t0 = std::chrono::steady_clock::now();
  for (int r = 0; r < reps; ++r) {
    const double bump = acc * 1e-300;
    for (int i = 0; i < BLOCK; ++i) {
      const double res
          = stan::math::exp_mod_normal_lpdf(ys[i] + bump, MU, SIGMA, LAMBDA);
      acc += std::isfinite(res) ? res : 1.0;
    }
  }
  const auto t1 = std::chrono::steady_clock::now();
  sink = acc;
  return std::chrono::duration<double, std::nano>(t1 - t0).count()
         / (static_cast<double>(reps) * BLOCK);
}

double time_var(const std::vector<double>& ys, int reps, double& sink) {
  double acc = sink;
  const auto t0 = std::chrono::steady_clock::now();
  for (int r = 0; r < reps; ++r) {
    const double bump = acc * 1e-300;
    for (int i = 0; i < BLOCK; ++i) {
      var y = ys[i] + bump, mu = MU, sigma = SIGMA, lambda = LAMBDA;
      var lp = stan::math::exp_mod_normal_lpdf(y, mu, sigma, lambda);
      lp.grad();
      const double res
          = lp.val() + y.adj() + mu.adj() + sigma.adj() + lambda.adj();
      acc += std::isfinite(res) ? res : 1.0;
      stan::math::recover_memory();
    }
  }
  const auto t1 = std::chrono::steady_clock::now();
  sink = acc;
  return std::chrono::duration<double, std::nano>(t1 - t0).count()
         / (static_cast<double>(reps) * BLOCK);
}

int main(int argc, char** argv) {
  const char* label = argc > 1 ? argv[1] : "?";
  struct Regime {
    const char* name;
    double lo, hi;
  };
  const Regime regimes[] = {{"z9.4", 9.3, 9.5}, {"z2", 1.5, 2.5},
                            {"z0", -0.1, 0.1}};
  double sink = 0;
  for (const auto& g : regimes) {
    const auto ys = inputs(g.lo, g.hi, 20261030);
    // Value check: the same inputs give the same sum in every tree.
    double check = 0, check_grad = 0;
    for (int i = 0; i < BLOCK; ++i) {
      check += stan::math::exp_mod_normal_lpdf(ys[i], MU, SIGMA, LAMBDA);
      var y = ys[i];
      var lp = stan::math::exp_mod_normal_lpdf(y, MU, SIGMA, LAMBDA);
      lp.grad();
      check_grad += y.adj();
      stan::math::recover_memory();
    }
    time_double(ys, 200, sink);  // warm-up
    const double td = time_double(ys, 4000, sink);
    const double tv = time_var(ys, 1000, sink);
    std::printf("RESULT %s %s double %.3f var %.3f sum_lpdf %.17g sum_dy %.17g\n",
                label, g.name, td, tv, check, check_grad);
  }
  std::printf("sink %g\n", sink);
  return 0;
}
