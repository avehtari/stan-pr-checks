// Speed benchmark behind the timing tables of the stan-dev/math#3406 PR text:
// ns per element for vector arguments of size N, one function per binary.
// Build: STAN_MATH=<tree> EXTRA_CXXFLAGS=-DBENCH_ONLY=<function> ../build.sh
// Usage: bench_one <label> <rounds>
//
// Argument types (the four tables of the PR text):
//   double  Eigen::VectorXd y, double parameters, value only
//   var     Eigen::Matrix<var, -1, 1> y, var parameters, value and gradient
//   varmat  var_value<Eigen::VectorXd> y, var parameters, value and gradient
//   fvar    Eigen::Matrix<fvar<double>, -1, 1> y (tangent 1), fvar<double>
//           parameters (tangent 0), value and tangent
// Inputs: z uniform on [-3, 3] ("bulk") or [-30, 30] ("wide"); Phi and
// std_normal_* take z, the others y = mu + sigma z (lognormal: exp of that),
// with mu = 0.3, sigma = 1.7, lambda = 1.
// Output: one RESULT line per (range, N, round) and one CHECK line per
// (range, N) with the value and one derivative, to compare the trees.
#include <stan/math/mix.hpp>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <random>
#include <string>
#include <vector>

using stan::math::fvar;
using stan::math::var;
using VecV = Eigen::Matrix<var, Eigen::Dynamic, 1>;
using VecF = Eigen::Matrix<fvar<double>, Eigen::Dynamic, 1>;
using VarMat = stan::math::var_value<Eigen::VectorXd>;
constexpr double MU = 0.3, SIGMA = 1.7, LAMBDA = 1.0;

#define FUN(NAME, CALL, LOG_Y)                             \
  struct NAME {                                            \
    static constexpr const char* name = #NAME;             \
    static constexpr bool log_y = LOG_Y;                   \
    template <typename Y, typename P>                      \
    static auto f(const Y& y, const P& mu, const P& sigma, \
                  const P& lambda) {                       \
      (void)mu;                                            \
      (void)sigma;                                         \
      (void)lambda;                                        \
      return CALL;                                         \
    }                                                      \
  };

FUN(Phi, stan::math::sum(stan::math::Phi(y)), false)
FUN(std_normal_cdf, stan::math::std_normal_cdf(y), false)
FUN(std_normal_lcdf, stan::math::std_normal_lcdf(y), false)
FUN(std_normal_lccdf, stan::math::std_normal_lccdf(y), false)
FUN(normal_cdf, stan::math::normal_cdf(y, mu, sigma), false)
FUN(normal_lcdf, stan::math::normal_lcdf(y, mu, sigma), false)
FUN(normal_lccdf, stan::math::normal_lccdf(y, mu, sigma), false)
FUN(lognormal_cdf, stan::math::lognormal_cdf(y, mu, sigma), true)
FUN(lognormal_lcdf, stan::math::lognormal_lcdf(y, mu, sigma), true)
FUN(lognormal_lccdf, stan::math::lognormal_lccdf(y, mu, sigma), true)
FUN(exp_mod_normal_lpdf, stan::math::exp_mod_normal_lpdf(y, mu, sigma, lambda),
    false)
FUN(exp_mod_normal_cdf, stan::math::exp_mod_normal_cdf(y, mu, sigma, lambda),
    false)
FUN(exp_mod_normal_lcdf, stan::math::exp_mod_normal_lcdf(y, mu, sigma, lambda),
    false)
FUN(exp_mod_normal_lccdf,
    stan::math::exp_mod_normal_lccdf(y, mu, sigma, lambda), false)
FUN(skew_normal_lpdf, stan::math::skew_normal_lpdf(y, mu, sigma, lambda), false)

template <typename F>
std::vector<double> inputs(int n, double width, std::mt19937_64& rng) {
  std::uniform_real_distribution<double> u(-width, width);
  const bool standard = std::string(F::name) == "Phi"
                        || std::string(F::name).rfind("std_normal", 0) == 0;
  std::vector<double> ys(n);
  for (auto& v : ys) {
    const double z = u(rng);
    v = standard ? z : MU + SIGMA * z;
    if (F::log_y) {
      v = std::exp(v);
    }
  }
  return ys;
}

// A non-finite result must not reach acc: the bump (acc * 1e-300, which
// keeps the loop from being hoisted) would then turn every later input
// into inf or NaN.
inline void accumulate(double& acc, double res) {
  acc += std::isfinite(res) ? res : 1.0;
}

template <typename F>
double eval_double(const std::vector<double>& ys, double bump,
                   Eigen::VectorXd& y, double* d) {
  for (size_t i = 0; i < ys.size(); ++i) {
    y[i] = ys[i] + bump;
  }
  *d = 0;
  return F::f(y, MU, SIGMA, LAMBDA);
}

template <typename F>
double eval_var(const std::vector<double>& ys, double bump, double* d) {
  VecV y(ys.size());
  for (size_t i = 0; i < ys.size(); ++i) {
    y[i] = ys[i] + bump;
  }
  var mu = MU, sigma = SIGMA, lambda = LAMBDA;
  var lp = F::f(y, mu, sigma, lambda);
  lp.grad();
  const double val = lp.val();
  *d = y[0].adj() + mu.adj();
  stan::math::recover_memory();
  return val;
}

template <typename F>
double eval_varmat(const std::vector<double>& ys, double bump,
                   Eigen::VectorXd& yv, double* d) {
  for (size_t i = 0; i < ys.size(); ++i) {
    yv[i] = ys[i] + bump;
  }
  VarMat y(yv);
  var mu = MU, sigma = SIGMA, lambda = LAMBDA;
  var lp = F::f(y, mu, sigma, lambda);
  lp.grad();
  const double val = lp.val();
  *d = y.adj()(0) + mu.adj();
  stan::math::recover_memory();
  return val;
}

template <typename F>
double eval_fvar(const std::vector<double>& ys, double bump, VecF& y,
                 double* d) {
  for (size_t i = 0; i < ys.size(); ++i) {
    y[i] = fvar<double>(ys[i] + bump, 1.0);
  }
  const fvar<double> mu(MU, 0.0), sigma(SIGMA, 0.0), lambda(LAMBDA, 0.0);
  const fvar<double> lp = F::f(y, mu, sigma, lambda);
  *d = lp.d_;
  return lp.val_;
}

enum class Kind { dbl, var, varmat, fvar };

template <typename F>
double eval(Kind k, const std::vector<double>& ys, double bump,
            Eigen::VectorXd& yd, VecF& yf, double* d) {
  switch (k) {
    case Kind::dbl:
      return eval_double<F>(ys, bump, yd, d);
    case Kind::var:
      return eval_var<F>(ys, bump, d);
    case Kind::varmat:
      return eval_varmat<F>(ys, bump, yd, d);
    default:
      return eval_fvar<F>(ys, bump, yf, d);
  }
}

template <typename F>
double time_kind(Kind k, const std::vector<double>& ys, int reps,
                 double& sink) {
  Eigen::VectorXd yd(ys.size());
  VecF yf(ys.size());
  double acc = sink, d = 0;
  const auto t0 = std::chrono::steady_clock::now();
  for (int r = 0; r < reps; ++r) {
    const double bump = acc * 1e-300;
    const double val = eval<F>(k, ys, bump, yd, yf, &d);
    accumulate(acc, val + d);
  }
  const auto t1 = std::chrono::steady_clock::now();
  sink = acc;
  return std::chrono::duration<double, std::nano>(t1 - t0).count()
         / (static_cast<double>(reps) * ys.size());
}

template <typename F>
void run(const char* label, int rounds, double& sink) {
  const Kind kinds[] = {Kind::dbl, Kind::var, Kind::varmat, Kind::fvar};
  const char* kind_names[] = {"double", "var", "varmat", "fvar"};
  for (double width : {3.0, 30.0}) {
    const char* range = width < 10 ? "bulk" : "wide";
    for (int n : {100, 1000, 10000, 100000}) {
      std::mt19937_64 rng(20261020 + n);
      const auto ys = inputs<F>(n, width, rng);
      // About 2e6 elements per double measurement, 1e6 for the others.
      const int reps = std::max(1, 2000000 / n);
      {
        Eigen::VectorXd yd(n);
        VecF yf(n);
        std::printf("CHECK %s %s %s %d", label, F::name, range, n);
        for (int j = 0; j < 4; ++j) {
          double d = 0;
          const double val = eval<F>(kinds[j], ys, 0.0, yd, yf, &d);
          std::printf(" %s %.17g %.17g", kind_names[j], val, d);
        }
        std::printf("\n");
      }
      for (int r = 0; r < rounds; ++r) {
        std::printf("RESULT %s %s %s %d %d", label, F::name, range, n, r + 1);
        for (int j = 0; j < 4; ++j) {
          const int kr = j == 0 ? reps : std::max(1, reps / 2);
          std::printf(" %s %.3f", kind_names[j],
                      time_kind<F>(kinds[j], ys, kr, sink));
        }
        std::printf("\n");
        std::fflush(stdout);
      }
    }
  }
}

#ifndef BENCH_ONLY
#error "build with -DBENCH_ONLY=<function>: one function per binary"
#endif

int main(int argc, char** argv) {
  const char* label = argc > 1 ? argv[1] : "?";
  const int rounds = argc > 2 ? std::atoi(argv[2]) : 2;
  double sink = 0;
  run<BENCH_ONLY>(label, rounds, sink);
  std::printf("sink %g\n", sink);
  return 0;
}
