// Claim B1 of README_emg.md: f4cd10ee4a, c30ee6f9ef and 927099efb8 leave all
// results bitwise identical.  Prints in hex the values and all partials of
// the touched functions for a fixed grid and several argument types; build it
// against a commit and its parent and diff the outputs.
#include <stan/math/mix.hpp>
#include <cmath>
#include <cstdio>
#include <exception>
#include <string>
#include <vector>

using stan::math::fvar;
using stan::math::var;
using stan::math::var_value;
using Eigen::VectorXd;
using vec_var = Eigen::Matrix<var, Eigen::Dynamic, 1>;

namespace {

// Standardized grid, z from the far lower to the far upper tail.
const std::vector<double> Z = {-40, -12, -5, -2, -0.7, 0, 0.4, 1.5, 4, 9, 40};
const double MU = 0.3, SIG = 1.7, LAM = 0.8;
const double MU_V[3] = {0.3, -0.2, 1.1};
const double SIG_V[3] = {1.7, 0.4, 2.5};
const double LAM_V[3] = {0.8, 3.0, 0.05};

void emit(const std::string& tag, const std::vector<double>& out) {
  std::printf("%s", tag.c_str());
  for (double x : out) {
    std::printf(" %a", x);
  }
  std::printf("\n");
}

template <typename T>
void push_adj(std::vector<double>& out, const T& x) {
  for (int i = 0; i < static_cast<int>(x.size()); ++i) {
    out.push_back(x[i].adj());
  }
}

// One struct per function: nargs, whether y > 0, and the call.
#define DEF_FN(NAME, NARGS, POS, EXPR)                                  \
  struct NAME {                                                         \
    static constexpr int nargs = NARGS;                                 \
    static constexpr bool positive = POS;                               \
    static const char* name() { return #NAME; }                         \
    template <typename Y, typename M, typename S, typename L>           \
    static auto call(const Y& y, const M& mu, const S& s, const L& l) { \
      (void)l;                                                          \
      return EXPR;                                                      \
    }                                                                   \
  };

DEF_FN(lognormal_cdf, 3, true, stan::math::lognormal_cdf(y, mu, s))
DEF_FN(lognormal_lcdf, 3, true, stan::math::lognormal_lcdf(y, mu, s))
DEF_FN(lognormal_lccdf, 3, true, stan::math::lognormal_lccdf(y, mu, s))
DEF_FN(normal_lcdf, 3, false, stan::math::normal_lcdf(y, mu, s))
DEF_FN(normal_lccdf, 3, false, stan::math::normal_lccdf(y, mu, s))
DEF_FN(exp_mod_normal_lcdf, 4, false,
       stan::math::exp_mod_normal_lcdf(y, mu, s, l))
DEF_FN(exp_mod_normal_lccdf, 4, false,
       stan::math::exp_mod_normal_lccdf(y, mu, s, l))
DEF_FN(exp_mod_normal_cdf, 4, false,
       stan::math::exp_mod_normal_cdf(y, mu, s, l))
DEF_FN(skew_normal_lpdf_full, 4, false,
       stan::math::skew_normal_lpdf<false>(y, mu, s, l))
DEF_FN(skew_normal_lpdf_propto, 4, false,
       stan::math::skew_normal_lpdf<true>(y, mu, s, l))

template <typename Fn>
std::vector<double> y_grid() {
  std::vector<double> ys;
  for (double z : Z) {
    ys.push_back(Fn::positive ? std::exp(MU + SIG * z / 8) : MU + SIG * z);
  }
  return ys;
}

// Runs body(); on a throw, prints the tag with "throw".
template <typename B>
void guarded(const std::string& tag, B body) {
  try {
    emit(tag, body());
  } catch (const std::exception&) {
    std::printf("%s throw\n", tag.c_str());
  }
  stan::math::recover_memory();
}

template <typename Fn>
void run() {
  const std::string f = Fn::name();
  const std::vector<double> ys = y_grid<Fn>();
  const std::size_t n = ys.size();
  // double scalars
  for (std::size_t i = 0; i < n; ++i) {
    guarded(f + " double " + std::to_string(i), [&] {
      return std::vector<double>{Fn::call(ys[i], MU, SIG, LAM)};
    });
  }
  // std::vector<double> y, scalar parameters
  guarded(f + " std_vector_double", [&] {
    return std::vector<double>{Fn::call(ys, MU, SIG, LAM)};
  });
  // Eigen vectors of y, mu and lambda, scalar sigma
  guarded(f + " eigen_double", [&] {
    VectorXd y(n), mu(n), l(n);
    for (std::size_t i = 0; i < n; ++i) {
      y[i] = ys[i];
      mu[i] = MU_V[i % 3];
      l[i] = LAM_V[i % 3];
    }
    return std::vector<double>{Fn::call(y, mu, SIG, l)};
  });
  // var scalars, all partials
  for (std::size_t i = 0; i < n; ++i) {
    guarded(f + " var " + std::to_string(i), [&] {
      var y = ys[i], mu = MU, s = SIG, l = LAM;
      var lp = Fn::call(y, mu, s, l);
      lp.grad();
      return std::vector<double>{lp.val(), y.adj(), mu.adj(), s.adj(),
                                 l.adj()};
    });
  }
  // double y, var parameters
  for (std::size_t i = 0; i < n; ++i) {
    guarded(f + " double_y_var_params " + std::to_string(i), [&] {
      var mu = MU, s = SIG, l = LAM;
      var lp = Fn::call(ys[i], mu, s, l);
      lp.grad();
      return std::vector<double>{lp.val(), mu.adj(), s.adj(), l.adj()};
    });
  }
  // std::vector<var> y, sigma and lambda, var mu
  guarded(f + " std_vector_var", [&] {
    std::vector<var> y, s, l;
    for (std::size_t i = 0; i < n; ++i) {
      y.emplace_back(ys[i]);
      s.emplace_back(SIG_V[i % 3]);
      l.emplace_back(LAM_V[i % 3]);
    }
    var mu = MU;
    var lp = Fn::call(y, mu, s, l);
    lp.grad();
    std::vector<double> out{lp.val(), mu.adj()};
    push_adj(out, y);
    push_adj(out, s);
    if (Fn::nargs == 4) {
      push_adj(out, l);
    }
    return out;
  });
  // Eigen vectors of var for y and mu, var sigma and lambda
  guarded(f + " eigen_var", [&] {
    vec_var y(n), mu(n);
    for (std::size_t i = 0; i < n; ++i) {
      y[i] = ys[i];
      mu[i] = MU_V[i % 3];
    }
    var s = SIG, l = LAM;
    var lp = Fn::call(y, mu, s, l);
    lp.grad();
    std::vector<double> out{lp.val(), s.adj(), l.adj()};
    push_adj(out, y);
    push_adj(out, mu);
    return out;
  });
  // var_value<VectorXd> y and sigma, var mu and lambda
  guarded(f + " var_value_vector", [&] {
    VectorXd yd(n), sd(n);
    for (std::size_t i = 0; i < n; ++i) {
      yd[i] = ys[i];
      sd[i] = SIG_V[i % 3];
    }
    var_value<VectorXd> y(yd), s(sd);
    var mu = MU, l = LAM;
    var lp = Fn::call(y, mu, s, l);
    lp.grad();
    std::vector<double> out{lp.val(), mu.adj(), l.adj()};
    for (std::size_t i = 0; i < n; ++i) {
      out.push_back(y.adj()[i]);
    }
    for (std::size_t i = 0; i < n; ++i) {
      out.push_back(s.adj()[i]);
    }
    return out;
  });
  // fvar<double> scalars, one tangent direction per argument
  for (std::size_t i = 0; i < n; ++i) {
    guarded(f + " fvar " + std::to_string(i), [&] {
      std::vector<double> out;
      for (int k = 0; k < Fn::nargs; ++k) {
        fvar<double> y(ys[i], k == 0), mu(MU, k == 1), s(SIG, k == 2),
            l(LAM, k == 3);
        const fvar<double> r = Fn::call(y, mu, s, l);
        if (k == 0) {
          out.push_back(r.val_);
        }
        out.push_back(r.d_);
      }
      return out;
    });
  }
  // std::vector<fvar<double>> y, tangent 1 in each element
  guarded(f + " std_vector_fvar", [&] {
    std::vector<fvar<double>> y;
    for (std::size_t i = 0; i < n; ++i) {
      y.emplace_back(ys[i], 1.0);
    }
    const fvar<double> r = Fn::call(y, MU, SIG, LAM);
    return std::vector<double>{r.val_, r.d_};
  });
}

}  // namespace

int main() {
  run<lognormal_cdf>();
  run<lognormal_lcdf>();
  run<lognormal_lccdf>();
  run<normal_lcdf>();
  run<normal_lccdf>();
  run<exp_mod_normal_lcdf>();
  run<exp_mod_normal_lccdf>();
  run<exp_mod_normal_cdf>();
  run<skew_normal_lpdf_full>();
  run<skew_normal_lpdf_propto>();
  return 0;
}
