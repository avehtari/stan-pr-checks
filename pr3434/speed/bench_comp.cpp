// bench_comp.cpp --- cost of the parts of the student_t_lpdf nu terms, to
// find the cause of the branch slowdown (job 20815014: scalar value 1.42x,
// scalar var 1.30x). stan-dev/math#3153, decision 7.
//
// Develop forms lgamma(h + 1/2) - lgamma(h) and digamma(h + 1/2) -
// digamma(h), h = nu / 2. The branch forms -lbeta(1/2, h) and
// digamma_diff(h, 1/2). This file times each part per call, in three nu
// ranges, because lbeta and digamma_diff take different code paths below
// and above h = 10. It also times two candidates (asymptotic series in 1/h
// for h >= 20): cost only here, accuracy is not checked by this file.
//
// Method: STAN_MATH_DISTRIBUTIONS.md section 4. The loop carries
// bump = acc * 1e-300, which changes no result bit and stops the compiler
// from hoisting the call. Compiled once against develop's student_t_lpdf
// (overlay) and once against the branch; the components are the same code
// in both binaries and show the noise between the two.
//
// Output: TIME <range> <kernel> <ns per call>

#include <stan/math/rev.hpp>

#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <vector>

using stan::math::var;

namespace {

constexpr int BLOCK = 256;

template <typename F>
double time_it(const char* range, const char* name, int reps, F&& f) {
  double acc = 0.0;
  auto t0 = std::chrono::steady_clock::now();
  for (int r = 0; r < reps; ++r) {
    const double bump = acc * 1e-300;
    for (int i = 0; i < BLOCK; ++i) {
      acc += f(i, bump);
    }
  }
  auto t1 = std::chrono::steady_clock::now();
  const double ns
      = 1e9 * std::chrono::duration<double>(t1 - t0).count() / (reps * BLOCK);
  std::printf("TIME %s %s %.3f\n", range, name, ns);
  std::printf("# %s %s acc=%.17g\n", range, name, acc);
  return ns;
}

// lgamma(h + 1/2) - lgamma(h) = log(h) / 2 - 1/(8h) + 1/(192 h^3)
//   - 1/(640 h^5) + 17/(14336 h^7) - 31/(18432 h^9) + ...
// (from B_{k+1}(1/2) - B_{k+1}(0)); next term 3.8e-3 / h^11.
inline double asym_const(double h) {
  if (h < 20.0) {
    return stan::math::lgamma(h + 0.5) - stan::math::lgamma(h);
  }
  const double z = 1.0 / h;
  const double z2 = z * z;
  return 0.5 * std::log(h)
         + z
               * (-0.125
                  + z2
                        * (1.0 / 192.0
                           + z2
                                 * (-1.0 / 640.0
                                    + z2
                                          * (17.0 / 14336.0
                                             + z2 * (-31.0 / 18432.0)))));
}

// The derivative of the series above in h: digamma(h + 1/2) - digamma(h).
inline double asym_dpsi(double h) {
  if (h < 20.0) {
    return stan::math::digamma_diff(h, 0.5);
  }
  const double z = 1.0 / h;
  const double z2 = z * z;
  return 0.5 * z
         + z2
               * (0.125
                  + z2
                        * (-1.0 / 64.0
                           + z2
                                 * (1.0 / 128.0
                                    + z2
                                          * (-119.0 / 14336.0
                                             + z2 * (279.0 / 18432.0)))));
}

void run_range(const char* range, const std::vector<double>& nu, int reps) {
  std::vector<double> h(BLOCK);
  for (int i = 0; i < BLOCK; ++i) {
    h[i] = 0.5 * nu[i];
  }
  const int n = reps * 20;
  time_it(range, "dev_const", n, [&](int i, double b) {
    const double x = h[i] + b;
    return stan::math::lgamma(x + 0.5) - stan::math::lgamma(x);
  });
  time_it(range, "lbeta_const", n, [&](int i, double b) {
    return stan::math::LOG_SQRT_PI - stan::math::lbeta(0.5, h[i] + b);
  });
  time_it(range, "lgamma_half", n, [&](int i, double b) {
    return stan::math::lgamma(0.5 + b);
  });
  time_it(range, "lgamma_one", n, [&](int i, double b) {
    return stan::math::lgamma(h[i] + b);
  });
  time_it(range, "lgsd_pair", n, [&](int i, double b) {
    const double x = h[i] + b;
    return stan::math::lgamma_stirling_diff(x)
           - stan::math::lgamma_stirling_diff(x + 0.5);
  });
  time_it(range, "log_h", n,
          [&](int i, double b) { return std::log(h[i] + b); });
  time_it(range, "asym_const", n,
          [&](int i, double b) { return asym_const(h[i] + b); });
  time_it(range, "dev_dpsi", n, [&](int i, double b) {
    const double x = h[i] + b;
    return stan::math::digamma(x + 0.5) - stan::math::digamma(x);
  });
  time_it(range, "digamma_one", n, [&](int i, double b) {
    return stan::math::digamma(h[i] + b);
  });
  time_it(range, "digamma_diff", n, [&](int i, double b) {
    return stan::math::digamma_diff(h[i] + b, 0.5);
  });
  time_it(range, "asym_dpsi", n,
          [&](int i, double b) { return asym_dpsi(h[i] + b); });
  time_it(range, "st_value", n, [&](int i, double b) {
    return stan::math::student_t_lpdf(1.0, nu[i] + b, 0.0, 1.0);
  });
  time_it(range, "st_grad", n, [&](int i, double b) {
    var v_nu = nu[i] + b;
    var lp = stan::math::student_t_lpdf(1.0, v_nu, 0.0, 1.0);
    lp.grad();
    const double out = lp.val() + v_nu.adj();
    stan::math::recover_memory();
    return out;
  });
}

}  // namespace

int main(int argc, char** argv) {
  const int reps = argc > 1 ? std::atoi(argv[1]) : 200;
  std::vector<double> nu(BLOCK);
  for (int i = 0; i < BLOCK; ++i) {
    nu[i] = 2.0 + 60.0 * (i % 97) / 97.0;
  }
  run_range("ALL_2-62", nu, reps);
  for (int i = 0; i < BLOCK; ++i) {
    nu[i] = 2.0 + 18.0 * (i % 97) / 97.0;
  }
  run_range("A_2-20", nu, reps);
  for (int i = 0; i < BLOCK; ++i) {
    nu[i] = 20.0 + 42.0 * (i % 97) / 97.0;
  }
  run_range("B_20-62", nu, reps);
  for (int i = 0; i < BLOCK; ++i) {
    nu[i] = std::pow(10.0, 2.0 + 6.0 * (i % 97) / 97.0);
  }
  run_range("C_1e2-1e8", nu, reps);
  return 0;
}
