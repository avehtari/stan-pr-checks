// Speed of Phi(double) and of Phi(var) with its gradient, per range.
// Compile once per tree (develop, #3406 commit) and alternate the binaries.
// Usage: bench_phi <label> <rounds>
#include <stan/math/rev.hpp>
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <random>
#include <vector>

double time_double(const std::vector<double>& in, int reps, double& sink) {
  double acc = sink;
  const auto t0 = std::chrono::steady_clock::now();
  for (int r = 0; r < reps; ++r) {
    const double bump = acc * 1e-300;
    for (double x : in) {
      acc += stan::math::Phi(x + bump);
    }
  }
  const auto t1 = std::chrono::steady_clock::now();
  sink = acc;
  return std::chrono::duration<double, std::nano>(t1 - t0).count()
         / (static_cast<double>(reps) * in.size());
}

// One reverse pass per block: Phi of every element, summed, then grad.
double time_var(const std::vector<double>& in, int reps, double& sink) {
  double acc = sink;
  const auto t0 = std::chrono::steady_clock::now();
  for (int r = 0; r < reps; ++r) {
    const double bump = acc * 1e-300;
    std::vector<stan::math::var> xs(in.begin(), in.end());
    stan::math::var total = 0;
    for (auto& x : xs) {
      x = x + bump;
      total += stan::math::Phi(x);
    }
    total.grad();
    acc += total.val() + xs[0].adj();
    stan::math::recover_memory();
  }
  const auto t1 = std::chrono::steady_clock::now();
  sink = acc;
  return std::chrono::duration<double, std::nano>(t1 - t0).count()
         / (static_cast<double>(reps) * in.size());
}

int main(int argc, char** argv) {
  const char* label = argc > 1 ? argv[1] : "?";
  const int rounds = argc > 2 ? std::atoi(argv[2]) : 3;
  const int block = 4096;
  const struct {
    const char* name;
    double lo, hi;
  } ranges[] = {{"central[-0.66,0.66]", -0.66, 0.66},
                {"mid[-5,5]", -5.0, 5.0},
                {"lower[-37,-5]", -37.0, -5.0},
                {"upper[5,8]", 5.0, 8.0},
                {"probit[-3,3]", -3.0, 3.0}};
  std::mt19937_64 rng(20261019);
  double sink = 0;
  for (const auto& r : ranges) {
    std::uniform_real_distribution<double> u(r.lo, r.hi);
    std::vector<double> in(block);
    for (auto& v : in) {
      v = u(rng);
    }
    for (int k = 0; k < rounds; ++k) {
      const double td = time_double(in, 600, sink);
      const double tv = time_var(in, 100, sink);
      std::printf("RESULT %s %s %d double %.3f var %.3f\n", label, r.name,
                  k + 1, td, tv);
    }
  }
  std::printf("sink %g\n", sink);
  return 0;
}
