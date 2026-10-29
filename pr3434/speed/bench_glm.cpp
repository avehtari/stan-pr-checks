// bench_glm.cpp --- time of neg_binomial_2_log_glm_lpmf, N = 1000
// instances, K = 5 predictors, propto = false. Built once per tree by
// run_bench_glm.sh (develop, the #3434 head before the review fixes, the
// current #3434 head).
//
// Kernels:
//   glm_value          all arguments double
//   glm_grad_ab        alpha and beta var, phi data
//   glm_grad_abphi     alpha, beta and a scalar phi var
//   glm_grad_vecphi    alpha, beta and a vector phi var
//   glm_grad_x         x var (matrix), alpha, beta, phi data
//
// Method: STAN_MATH_DISTRIBUTIONS.md section 4. Each call adds
// acc * 1e-300 to alpha; this changes no result bit and stops the compiler
// from hoisting the call. The value and two adjoints are printed in hex so
// that the trees can be compared.
//
// Output: TIME <kernel> <ns per instance>

#include <stan/math/rev.hpp>

#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <random>
#include <utility>
#include <vector>

using stan::math::var;

namespace {

constexpr int N = 1000;
constexpr int K = 5;

struct glm_data {
  std::vector<int> y;
  Eigen::MatrixXd x;
  Eigen::VectorXd beta;
  Eigen::VectorXd phi_v;
  double alpha;
  double phi;
};

glm_data make_data() {
  std::mt19937_64 rng(20261022);
  std::normal_distribution<double> nx(0.0, 0.5);
  std::uniform_real_distribution<double> lphi(-1.0, 3.0);
  glm_data d;
  d.x.resize(N, K);
  d.beta.resize(K);
  d.phi_v.resize(N);
  d.alpha = 1.0;
  d.phi = 2.5;
  for (int k = 0; k < K; ++k) {
    d.beta(k) = 0.3 - 0.15 * k;
  }
  for (int i = 0; i < N; ++i) {
    for (int k = 0; k < K; ++k) {
      d.x(i, k) = nx(rng);
    }
    d.phi_v(i) = std::exp(lphi(rng));
  }
  const Eigen::VectorXd theta = (d.x * d.beta).array() + d.alpha;
  for (int i = 0; i < N; ++i) {
    std::poisson_distribution<int> p(std::exp(theta(i)));
    d.y.push_back(p(rng));
  }
  return d;
}

template <typename F>
void time_it(const char* name, int reps, F&& f) {
  double probe[3] = {0, 0, 0};
  f(0.0, probe);
  std::printf("VALUE %s %a %a %a\n", name, probe[0], probe[1], probe[2]);
  double acc = 0.0;
  auto t0 = std::chrono::steady_clock::now();
  for (int r = 0; r < reps; ++r) {
    acc += f(acc * 1e-300, probe);
  }
  auto t1 = std::chrono::steady_clock::now();
  const double ns
      = 1e9 * std::chrono::duration<double>(t1 - t0).count() / (reps * N);
  std::printf("TIME %s %.3f\n", name, ns);
  std::printf("# %s acc=%.17g\n", name, acc);
}

}  // namespace

int main(int argc, char** argv) {
  using stan::math::neg_binomial_2_log_glm_lpmf;
  const int reps = argc > 1 ? std::atoi(argv[1]) : 2000;
  const glm_data d = make_data();

  time_it("glm_value", reps, [&](double b, double* out) {
    const double lp = neg_binomial_2_log_glm_lpmf(d.y, d.x, d.alpha + b,
                                                  d.beta, d.phi);
    out[0] = lp;
    return lp;
  });

  time_it("glm_grad_ab", reps / 4, [&](double b, double* out) {
    var alpha = d.alpha + b;
    Eigen::Matrix<var, Eigen::Dynamic, 1> beta = d.beta.cast<var>();
    var lp = neg_binomial_2_log_glm_lpmf(d.y, d.x, alpha, beta, d.phi);
    lp.grad();
    out[0] = lp.val();
    out[1] = alpha.adj();
    out[2] = beta(0).adj();
    stan::math::recover_memory();
    return out[0];
  });

  time_it("glm_grad_abphi", reps / 4, [&](double b, double* out) {
    var alpha = d.alpha + b;
    var phi = d.phi;
    Eigen::Matrix<var, Eigen::Dynamic, 1> beta = d.beta.cast<var>();
    var lp = neg_binomial_2_log_glm_lpmf(d.y, d.x, alpha, beta, phi);
    lp.grad();
    out[0] = lp.val();
    out[1] = alpha.adj();
    out[2] = phi.adj();
    stan::math::recover_memory();
    return out[0];
  });

  // the sampling statement y ~ neg_binomial_2_log_glm(...) with parameter
  // phi: propto = true
  time_it("glm_grad_abphi_propto", reps / 4, [&](double b, double* out) {
    var alpha = d.alpha + b;
    var phi = d.phi;
    Eigen::Matrix<var, Eigen::Dynamic, 1> beta = d.beta.cast<var>();
    var lp = neg_binomial_2_log_glm_lpmf<true>(d.y, d.x, alpha, beta, phi);
    lp.grad();
    out[0] = lp.val();
    out[1] = alpha.adj();
    out[2] = phi.adj();
    stan::math::recover_memory();
    return out[0];
  });

  // worst cases for a table of counts: all counts distinct (table used,
  // nothing reused), and counts above 4 N + 64 (no table)
  glm_data d_distinct = d;
  glm_data d_large = d;
  for (int i = 0; i < N; ++i) {
    d_distinct.y[i] = i;
    d_large.y[i] = 5000 + i;
  }
  const std::pair<const char*, const glm_data*> worst[]
      = {{"glm_grad_abphi_propto_distinct", &d_distinct},
         {"glm_grad_abphi_propto_large", &d_large}};
  for (const auto& w : worst) {
    const glm_data& dw = *w.second;
    time_it(w.first, reps / 4, [&](double b, double* out) {
      var alpha = dw.alpha + b;
      var phi = dw.phi;
      Eigen::Matrix<var, Eigen::Dynamic, 1> beta = dw.beta.cast<var>();
      var lp = neg_binomial_2_log_glm_lpmf<true>(dw.y, dw.x, alpha, beta, phi);
      lp.grad();
      out[0] = lp.val();
      out[1] = alpha.adj();
      out[2] = phi.adj();
      stan::math::recover_memory();
      return out[0];
    });
  }

  time_it("glm_grad_vecphi", reps / 4, [&](double b, double* out) {
    var alpha = d.alpha + b;
    Eigen::Matrix<var, Eigen::Dynamic, 1> beta = d.beta.cast<var>();
    Eigen::Matrix<var, Eigen::Dynamic, 1> phi = d.phi_v.cast<var>();
    var lp = neg_binomial_2_log_glm_lpmf(d.y, d.x, alpha, beta, phi);
    lp.grad();
    out[0] = lp.val();
    out[1] = alpha.adj();
    out[2] = phi(0).adj() + phi(N - 1).adj();
    stan::math::recover_memory();
    return out[0];
  });

  time_it("glm_grad_x", reps / 8, [&](double b, double* out) {
    Eigen::Matrix<var, Eigen::Dynamic, Eigen::Dynamic> x = d.x.cast<var>();
    var lp = neg_binomial_2_log_glm_lpmf(d.y, x, d.alpha + b, d.beta, d.phi);
    lp.grad();
    out[0] = lp.val();
    out[1] = x(0, 0).adj();
    out[2] = x(N - 1, K - 1).adj();
    stan::math::recover_memory();
    return out[0];
  });
  return 0;
}
