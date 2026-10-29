// propto_probe.cpp --- values of every changed lpmf/lpdf of branch
// stable-beta-binomial-lpmf with propto = true and false, for every
// combination of data (double) and parameter (var) arguments
// (stan-dev/math#3153, decision 10: propto must drop the same terms as on
// develop). Build once against develop and once against the branch, then
// compare the two outputs with compare_propto.py.
//
// Output lines: <function> <mask> <propto 0/1> <value>; bit i of mask set
// means real argument i (in the order of the call) is a var.

#include <stan/math/rev.hpp>

#include <cstdio>
#include <type_traits>
#include <utility>
#include <vector>

using stan::math::var;

namespace {

template <bool V>
using S = std::conditional_t<V, var, double>;

template <bool V>
S<V> sc(double x) {
  return S<V>(x);
}

template <bool V, typename M>
auto cv(const M& m) {
  return m.template cast<S<V>>().eval();
}

// recover_memory() frees every var of a mask, so it runs only after all
// calls of that mask.
template <typename F, int... Ms>
void masks_impl(F&& f, std::integer_sequence<int, Ms...>) {
  ((f(std::integral_constant<int, Ms>{}), stan::math::recover_memory()), ...);
}

template <int K, typename F>
void masks(F&& f) {
  masks_impl(f, std::make_integer_sequence<int, (1 << K)>{});
}

void out(const char* fn, int m, int propto, double v) {
  std::printf("%s %d %d %.17g\n", fn, m, propto, v);
}

}  // namespace

#define B(i) (((M >> (i)) & 1) != 0)
#define VAL(x) stan::math::value_of_rec(x)
#define BOTH(name, call_true, call_false) \
  out(name, M, 1, VAL(call_true));        \
  out(name, M, 0, VAL(call_false));

int main() {
  using Eigen::MatrixXd;
  using Eigen::VectorXd;
  using stan::math::beta_binomial_lpmf;
  using stan::math::beta_neg_binomial_lpmf;
  using stan::math::dirichlet_multinomial_lpmf;
  using stan::math::lkj_corr_cholesky_lpdf;
  using stan::math::lkj_corr_lpdf;
  using stan::math::lkj_cov_lpdf;
  using stan::math::multi_student_t_cholesky_lpdf;
  using stan::math::multi_student_t_lpdf;
  using stan::math::neg_binomial_2_log_glm_lpmf;
  using stan::math::neg_binomial_2_log_lpmf;
  using stan::math::neg_binomial_2_lpmf;
  using stan::math::neg_binomial_lpmf;
  using stan::math::student_t_lpdf;
  using stan::math::yule_simon_lpmf;

  const std::vector<int> n2{3, 7};
  const std::vector<int> N2{10, 12};
  VectorXd a2(2);
  a2 << 2.5, 7.0;
  VectorXd b2(2);
  b2 << 1.5, 30.0;
  masks<2>([&](auto mc) {
    constexpr int M = decltype(mc)::value;
    auto a = cv<B(0)>(a2);
    auto b = cv<B(1)>(b2);
    BOTH("beta_binomial", beta_binomial_lpmf<true>(n2, N2, a, b),
         beta_binomial_lpmf<false>(n2, N2, a, b));
  });

  const std::vector<int> nb{3, 5};
  VectorXd r2(2);
  r2 << 2.5, 4.0;
  masks<3>([&](auto mc) {
    constexpr int M = decltype(mc)::value;
    auto r = cv<B(0)>(r2);
    auto a = cv<B(1)>(a2);
    auto b = cv<B(2)>(b2);
    BOTH("beta_neg_binomial", beta_neg_binomial_lpmf<true>(nb, r, a, b),
         beta_neg_binomial_lpmf<false>(nb, r, a, b));
  });

  VectorXd inv2(2);
  inv2 << 0.7, 1.5;
  masks<2>([&](auto mc) {
    constexpr int M = decltype(mc)::value;
    auto a = cv<B(0)>(r2);
    auto b = cv<B(1)>(inv2);
    BOTH("neg_binomial", neg_binomial_lpmf<true>(nb, a, b),
         neg_binomial_lpmf<false>(nb, a, b));
  });

  VectorXd eta2(2);
  eta2 << 0.3, 1.2;
  masks<2>([&](auto mc) {
    constexpr int M = decltype(mc)::value;
    auto mu = cv<B(0)>(r2);
    auto phi = cv<B(1)>(b2);
    BOTH("neg_binomial_2", neg_binomial_2_lpmf<true>(nb, mu, phi),
         neg_binomial_2_lpmf<false>(nb, mu, phi));
    auto eta = cv<B(0)>(eta2);
    BOTH("neg_binomial_2_log", neg_binomial_2_log_lpmf<true>(nb, eta, phi),
         neg_binomial_2_log_lpmf<false>(nb, eta, phi));
  });

  const std::vector<int> yg{3, 5, 2};
  MatrixXd xg(3, 2);
  xg << 0.5, -0.2, 1.1, 0.3, -0.4, 0.8;
  VectorXd bg(2);
  bg << 0.4, -0.3;
  VectorXd ag(3);
  ag << 0.2, 0.5, -0.1;
  VectorXd pg(3);
  pg << 1.5, 30.0, 4.0;
  masks<4>([&](auto mc) {
    constexpr int M = decltype(mc)::value;
    auto x = cv<B(0)>(xg);
    auto alpha = sc<B(1)>(0.3);
    auto beta = cv<B(2)>(bg);
    auto phi = sc<B(3)>(3.5);
    BOTH("glm_scalar_phi",
         neg_binomial_2_log_glm_lpmf<true>(yg, x, alpha, beta, phi),
         neg_binomial_2_log_glm_lpmf<false>(yg, x, alpha, beta, phi));
    auto alpha_v = cv<B(1)>(ag);
    auto phi_v = cv<B(3)>(pg);
    BOTH("glm_vector_phi",
         neg_binomial_2_log_glm_lpmf<true>(yg, x, alpha_v, beta, phi_v),
         neg_binomial_2_log_glm_lpmf<false>(yg, x, alpha_v, beta, phi_v));
    BOTH("glm_scalar_y",
         neg_binomial_2_log_glm_lpmf<true>(4, x, alpha, beta, phi),
         neg_binomial_2_log_glm_lpmf<false>(4, x, alpha, beta, phi));
  });

  const std::vector<int> nd{3, 5, 2};
  VectorXd ad(3);
  ad << 2.5, 1.5, 7.0;
  masks<1>([&](auto mc) {
    constexpr int M = decltype(mc)::value;
    auto alpha = cv<B(0)>(ad);
    BOTH("dirichlet_multinomial", dirichlet_multinomial_lpmf<true>(nd, alpha),
         dirichlet_multinomial_lpmf<false>(nd, alpha));
  });

  MatrixXd corr(3, 3);
  corr << 1.0, 0.3, -0.2, 0.3, 1.0, 0.1, -0.2, 0.1, 1.0;
  MatrixXd lcorr = corr.llt().matrixL();
  VectorXd mu3(3);
  mu3 << 0.1, -0.2, 0.3;
  VectorXd sd3(3);
  sd3 << 1.5, 0.7, 2.0;
  MatrixXd cov = sd3.asDiagonal() * corr * sd3.asDiagonal();
  for (double eta_d : {2.5, 1.0}) {
    masks<2>([&](auto mc) {
      constexpr int M = decltype(mc)::value;
      auto y = cv<B(0)>(corr);
      auto l = cv<B(0)>(lcorr);
      auto eta = sc<B(1)>(eta_d);
      BOTH(eta_d == 1.0 ? "lkj_corr_eta1" : "lkj_corr",
           lkj_corr_lpdf<true>(y, eta), lkj_corr_lpdf<false>(y, eta));
      BOTH(eta_d == 1.0 ? "lkj_corr_cholesky_eta1" : "lkj_corr_cholesky",
           lkj_corr_cholesky_lpdf<true>(l, eta),
           lkj_corr_cholesky_lpdf<false>(l, eta));
    });
    masks<4>([&](auto mc) {
      constexpr int M = decltype(mc)::value;
      // develop defect, not in this PR: with a data eta == 1, lkj_cov_lpdf
      // passes the covariance matrix to lkj_corr_lpdf, which throws
      if (eta_d == 1.0 && !B(3)) {
        return;
      }
      auto y = cv<B(0)>(cov);
      auto mu = cv<B(1)>(mu3);
      auto sigma = cv<B(2)>(sd3);
      auto eta = sc<B(3)>(eta_d);
      BOTH(eta_d == 1.0 ? "lkj_cov_eta1" : "lkj_cov",
           lkj_cov_lpdf<true>(y, mu, sigma, eta),
           lkj_cov_lpdf<false>(y, mu, sigma, eta));
    });
  }

  VectorXd yt(3);
  yt << 0.3, -1.2, 2.0;
  VectorXd nut(3);
  // develop loses eps * lgamma(nu / 2) in the value; keep nu small enough
  // that this stays below the comparison tolerance
  nut << 3.5, 55.0, 500.0;
  VectorXd mut(3);
  mut << 0.1, 0.2, -0.3;
  masks<4>([&](auto mc) {
    constexpr int M = decltype(mc)::value;
    auto y = cv<B(0)>(yt);
    auto nu = cv<B(1)>(nut);
    auto mu = cv<B(2)>(mut);
    auto sigma = cv<B(3)>(sd3);
    BOTH("student_t", student_t_lpdf<true>(y, nu, mu, sigma),
         student_t_lpdf<false>(y, nu, mu, sigma));
    for (double nu_d : {3.5, 55.0}) {
      auto nu_s = sc<B(1)>(nu_d);
      BOTH(nu_d < 40 ? "student_t_scalar_nu_small" : "student_t_scalar_nu_large",
           student_t_lpdf<true>(y, nu_s, mu, sigma),
           student_t_lpdf<false>(y, nu_s, mu, sigma));
    }
  });

  VectorXd y2(2);
  y2 << 0.4, -1.1;
  VectorXd m2(2);
  m2 << 0.1, 0.2;
  MatrixXd sig2(2, 2);
  sig2 << 1.5, 0.3, 0.3, 0.8;
  MatrixXd l2 = sig2.llt().matrixL();
  for (double nu_d : {4.5, 60.0}) {
    masks<4>([&](auto mc) {
      constexpr int M = decltype(mc)::value;
      auto y = cv<B(0)>(y2);
      auto nu = sc<B(1)>(nu_d);
      auto mu = cv<B(2)>(m2);
      auto sigma = cv<B(3)>(sig2);
      auto l = cv<B(3)>(l2);
      BOTH("multi_student_t", multi_student_t_lpdf<true>(y, nu, mu, sigma),
           multi_student_t_lpdf<false>(y, nu, mu, sigma));
      BOTH("multi_student_t_cholesky",
           multi_student_t_cholesky_lpdf<true>(y, nu, mu, l),
           multi_student_t_cholesky_lpdf<false>(y, nu, mu, l));
    });
  }

  const std::vector<int> ny{1, 3, 7};
  VectorXd ay(3);
  ay << 2.5, 0.7, 30.0;
  masks<1>([&](auto mc) {
    constexpr int M = decltype(mc)::value;
    auto alpha = cv<B(0)>(ay);
    BOTH("yule_simon", yule_simon_lpmf<true>(ny, alpha),
         yule_simon_lpmf<false>(ny, alpha));
  });
  return 0;
}
