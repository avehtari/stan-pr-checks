// Speed of std_normal_cdf, normal_cdf, lognormal_cdf and exp_mod_normal_cdf,
// per element, for double and var (all arguments var, one reverse pass per
// block of 4096 elements).  K = elements per call (1 = scalar calls, 64 =
// vector argument).  Compile once per tree and alternate the binaries.
// Usage: bench_cdfs <label> <rounds>
#include <stan/math/rev.hpp>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <random>
#include <vector>

using stan::math::var;

constexpr double MU = 0.3, SIGMA = 1.7, LAMBDA = 1.0;

// y for the location-scale cdfs; z is the standardized value.
struct StdNormal {
  static constexpr const char* name = "std_normal_cdf";
  static double y(double z) { return z; }
  template <typename Y, typename P>
  static auto f(const Y& y, const P&, const P&, const P&) {
    return stan::math::std_normal_cdf(y);
  }
};
struct Normal {
  static constexpr const char* name = "normal_cdf";
  static double y(double z) { return MU + SIGMA * z; }
  template <typename Y, typename P>
  static auto f(const Y& y, const P& mu, const P& sigma, const P&) {
    return stan::math::normal_cdf(y, mu, sigma);
  }
};
struct Lognormal {
  static constexpr const char* name = "lognormal_cdf";
  static double y(double z) { return std::exp(MU + SIGMA * z); }
  template <typename Y, typename P>
  static auto f(const Y& y, const P& mu, const P& sigma, const P&) {
    return stan::math::lognormal_cdf(y, mu, sigma);
  }
};
struct ExpModNormal {
  static constexpr const char* name = "exp_mod_normal_cdf";
  static double y(double z) { return MU + SIGMA * z; }
  template <typename Y, typename P>
  static auto f(const Y& y, const P& mu, const P& sigma, const P& lambda) {
    return stan::math::exp_mod_normal_cdf(y, mu, sigma, lambda);
  }
};

template <typename F>
double time_double(const std::vector<double>& ys, int k, int reps,
                   double& sink) {
  double acc = sink;
  std::vector<double> chunk(k);
  const auto t0 = std::chrono::steady_clock::now();
  for (int r = 0; r < reps; ++r) {
    const double bump = acc * 1e-300;
    for (size_t i = 0; i + k <= ys.size(); i += k) {
      if (k == 1) {
        acc += F::f(ys[i] + bump, MU, SIGMA, LAMBDA);
      } else {
        for (int j = 0; j < k; ++j) {
          chunk[j] = ys[i + j] + bump;
        }
        acc += F::f(chunk, MU, SIGMA, LAMBDA);
      }
    }
  }
  const auto t1 = std::chrono::steady_clock::now();
  sink = acc;
  return std::chrono::duration<double, std::nano>(t1 - t0).count()
         / (static_cast<double>(reps) * ys.size());
}

template <typename F>
double time_var(const std::vector<double>& ys, int k, int reps, double& sink) {
  double acc = sink;
  const auto t0 = std::chrono::steady_clock::now();
  for (int r = 0; r < reps; ++r) {
    const double bump = acc * 1e-300;
    var mu = MU + bump, sigma = SIGMA, lambda = LAMBDA;
    var total = 0;
    std::vector<var> chunk(k);
    for (size_t i = 0; i + k <= ys.size(); i += k) {
      if (k == 1) {
        total += F::f(var(ys[i] + bump), mu, sigma, lambda);
      } else {
        for (int j = 0; j < k; ++j) {
          chunk[j] = ys[i + j] + bump;
        }
        total += F::f(chunk, mu, sigma, lambda);
      }
    }
    total.grad();
    acc += total.val() + mu.adj() + sigma.adj() + lambda.adj();
    stan::math::recover_memory();
  }
  const auto t1 = std::chrono::steady_clock::now();
  sink = acc;
  return std::chrono::duration<double, std::nano>(t1 - t0).count()
         / (static_cast<double>(reps) * ys.size());
}

template <typename F>
void run(const char* label, int rounds, double& sink) {
  const struct {
    const char* name;
    double lo, hi;
    int k;
  } cases[] = {{"bulk[-3,3]", -3.0, 3.0, 1},
               {"bulk[-3,3]", -3.0, 3.0, 64},
               {"lower[-8,-3]", -8.0, -3.0, 1},
               {"upper[3,8]", 3.0, 8.0, 1}};
  std::mt19937_64 rng(20261019);
  for (const auto& c : cases) {
    std::uniform_real_distribution<double> u(c.lo, c.hi);
    std::vector<double> ys(4096);
    for (auto& v : ys) {
      v = F::y(u(rng));
    }
    for (int r = 0; r < rounds; ++r) {
      const double td = time_double<F>(ys, c.k, 200, sink);
      const double tv = time_var<F>(ys, c.k, 50, sink);
      std::printf("RESULT %s %s %s %d %d double %.3f var %.3f\n", label,
                  F::name, c.name, c.k, r + 1, td, tv);
    }
  }
}

int main(int argc, char** argv) {
  const char* label = argc > 1 ? argv[1] : "?";
  const int rounds = argc > 2 ? std::atoi(argv[2]) : 3;
  double sink = 0;
  run<StdNormal>(label, rounds, sink);
  run<Normal>(label, rounds, sink);
  run<Lognormal>(label, rounds, sink);
  run<ExpModNormal>(label, rounds, sink);
  std::printf("sink %g\n", sink);
  return 0;
}
