// bench_st.cpp --- time of student_t_lpdf, scalar value and scalar value
// with the nu gradient, per nu range. One binary for student_t_lpdf only
// (STAN_MATH_DISTRIBUTIONS.md section 4: one binary per function), built
// once per stan/ tree by run_bench_pr3434_roihu.sh.
//
// Ranges (97 points each, repeated over a block of 256):
//   A_2-20       nu from 2 to 20 (the range A of ir_st.cpp)
//   L_2-40       nu from 2 to 40 (develop's code path in the PR: nu < 40)
//   C_1e2-1e8    nu log-spaced from 1e2 to 1e8 (the series in the PR)
//
// Kernels: st_value (all double), st_grad (nu var; value plus adjoint),
// ctrl_lgamma (lgamma(nu / 2); the same code in every tree: the noise
// between binaries).
//
// Method: STAN_MATH_DISTRIBUTIONS.md section 4. The loop carries
// bump = acc * 1e-300, which changes no result bit and stops the compiler
// from hoisting the call.
//
// Output: TIME <range>:<kernel> <ns per call>

#include <stan/math/rev.hpp>

#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>

using stan::math::var;

namespace {

constexpr int BLOCK = 256;

template <typename F>
void time_it(const std::string& name, int reps, F&& f) {
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
  std::printf("TIME %s %.3f\n", name.c_str(), ns);
  std::printf("# %s acc=%.17g\n", name.c_str(), acc);
}

void run_range(const std::string& range, const std::vector<double>& nu,
               int n) {
  time_it(range + ":ctrl_lgamma", n, [&](int i, double b) {
    return stan::math::lgamma(0.5 * (nu[i] + b));
  });
  time_it(range + ":st_value", n, [&](int i, double b) {
    return stan::math::student_t_lpdf(1.0, nu[i] + b, 0.0, 1.0);
  });
  time_it(range + ":st_grad", n, [&](int i, double b) {
    var v_nu = nu[i] + b;
    var lp = stan::math::student_t_lpdf(1.0, v_nu, 0.0, 1.0);
    lp.grad();
    const double out = lp.val() + v_nu.adj();
    stan::math::recover_memory();
    return out;
  });
}

std::vector<double> nu_grid(double lo, double hi, bool log_scale) {
  std::vector<double> v(BLOCK);
  for (int i = 0; i < BLOCK; ++i) {
    const double t = (i % 97) / 97.0;
    v[i] = log_scale ? lo * std::pow(hi / lo, t) : lo + (hi - lo) * t;
  }
  return v;
}

}  // namespace

int main(int argc, char** argv) {
  const int reps = argc > 1 ? std::atoi(argv[1]) : 200;
  const int n = reps * 20;
  run_range("A_2-20", nu_grid(2.0, 20.0, false), n);
  run_range("L_2-40", nu_grid(2.0, 40.0, false), n);
  run_range("C_1e2-1e8", nu_grid(1e2, 1e8, true), n);
  return 0;
}
