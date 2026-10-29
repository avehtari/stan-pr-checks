// bench_st_vec.cpp --- speed of student_t_lpdf, scalar and vectorized calls,
// develop against the branch (stan-dev/math#3153 uniformity change:
// -lbeta(1/2, nu/2) for the constant, digamma_diff for the nu partial).
//
// Method: STAN_MATH_DISTRIBUTIONS.md section 4. The loop carries
// bump = acc * 1e-300, which changes no result bit and stops the compiler
// from hoisting the call. Compiled once against develop (overlay) and once
// against the branch; the job runs the two alternately and reports ratios.
//
// Output: TIME <kernel> <ns per call>

#include <stan/math/rev.hpp>

#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <vector>

using stan::math::var;

namespace {

constexpr int BLOCK = 256;
constexpr int N_OBS = 1000;

template <typename F>
double time_it(const char* name, int reps, F&& f) {
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
  std::printf("TIME %s %.3f\n", name, ns);
  std::printf("# %s acc=%.17g\n", name, acc);
  return ns;
}

}  // namespace

int main(int argc, char** argv) {
  const int reps = argc > 1 ? std::atoi(argv[1]) : 200;
  std::vector<double> nu(BLOCK);
  for (int i = 0; i < BLOCK; ++i) {
    nu[i] = 2.0 + 60.0 * (i % 97) / 97.0;
  }
  Eigen::VectorXd y(N_OBS);
  for (int j = 0; j < N_OBS; ++j) {
    y(j) = -3.0 + 6.0 * (j % 101) / 101.0;
  }

  time_it("scalar_value", reps * 20, [&](int i, double bump) {
    return stan::math::student_t_lpdf(1.0, nu[i] + bump, 0.0, 1.0);
  });
  time_it("scalar_grad", reps * 20, [&](int i, double bump) {
    var v_nu = nu[i] + bump;
    var lp = stan::math::student_t_lpdf(1.0, v_nu, 0.0, 1.0);
    lp.grad();
    const double out = lp.val() + v_nu.adj();
    stan::math::recover_memory();
    return out;
  });
  time_it("vector1000_value", reps, [&](int i, double bump) {
    return stan::math::student_t_lpdf(y, nu[i] + bump, 0.0, 1.0);
  });
  time_it("vector1000_grad_nu", reps, [&](int i, double bump) {
    var v_nu = nu[i] + bump;
    var lp = stan::math::student_t_lpdf(y, v_nu, 0.0, 1.0);
    lp.grad();
    const double out = lp.val() + v_nu.adj();
    stan::math::recover_memory();
    return out;
  });
  return 0;
}
