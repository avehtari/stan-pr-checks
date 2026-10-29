// glm_compare.cpp --- values and partials of neg_binomial_2_log_glm_lpmf
// for every argument kind, to compare two trees (compare_glm.py). Built once
// against the #3434 head and once against the per-count/shared-e version.
//
// Cases: regimes of theta (alpha from -30 to 715) and phi (1e-17 to 1e15),
// counts from a small range and with a large count; N = 200 (table path)
// and N = 5 with a count of 1000 (fallback path, 1000 > 4 N + 64).
//
// Output: <label> <value> <partials...>, in %a.

#include <stan/math/mix.hpp>

#include <cstdio>
#include <random>
#include <string>
#include <vector>

using stan::math::fvar;
using stan::math::var;
using Eigen::Dynamic;
using Eigen::Matrix;
using Eigen::MatrixXd;
using Eigen::VectorXd;

namespace {

struct glm_case {
  std::string tag;
  std::vector<int> y;
  MatrixXd x;
  VectorXd beta;
  double alpha;
  VectorXd phi_v;
  double phi;
};

void print(const std::string& label, const std::vector<double>& v) {
  std::printf("%s", label.c_str());
  for (double d : v) {
    std::printf(" %a", d);
  }
  std::printf("\n");
}

std::vector<glm_case> make_cases() {
  std::vector<glm_case> cases;
  std::mt19937_64 rng(20261022);
  std::normal_distribution<double> nx(0.0, 0.5);
  std::uniform_real_distribution<double> lu(-1.0, 1.0);
  const double alphas[] = {-30.0, 0.0, 1.0, 5.0, 30.0, 709.0, 715.0};
  const double phis[] = {1e-17, 1e-8, 0.07, 1.5, 40.0, 1e3, 1e8, 1e15};
  const int sizes[] = {200, 5};
  for (int N : sizes) {
    for (double a : alphas) {
      for (double p : phis) {
        glm_case c;
        c.tag = "N" + std::to_string(N) + "_a" + std::to_string(a) + "_p"
                + std::to_string(p);
        c.x.resize(N, 3);
        c.beta.resize(3);
        c.beta << 0.3, -0.2, 0.1;
        c.phi_v.resize(N);
        for (int i = 0; i < N; ++i) {
          for (int k = 0; k < 3; ++k) {
            c.x(i, k) = nx(rng);
          }
          c.phi_v(i) = p * std::exp(lu(rng));
          int yi = (i * 7) % 21;
          if (i == 3) {
            yi = 1000;
          }
          if (i == 1) {
            yi = 0;
          }
          c.y.push_back(yi);
        }
        c.alpha = a;
        c.phi = p;
        cases.push_back(c);
      }
    }
  }
  return cases;
}

template <bool propto>
void run_case(const glm_case& c) {
  using stan::math::neg_binomial_2_log_glm_lpmf;
  const std::string pt = propto ? "_pt1" : "_pt0";
  const int N = c.y.size();
  // all double
  print(c.tag + pt + "_double",
        {stan::math::value_of(neg_binomial_2_log_glm_lpmf<propto>(
            c.y, c.x, c.alpha, c.beta, c.phi))});
  // var alpha, beta, scalar phi
  {
    var alpha = c.alpha;
    var phi = c.phi;
    Matrix<var, Dynamic, 1> beta = c.beta.cast<var>();
    var lp = neg_binomial_2_log_glm_lpmf<propto>(c.y, c.x, alpha, beta, phi);
    lp.grad();
    print(c.tag + pt + "_var_scalar_phi",
          {lp.val(), alpha.adj(), beta(0).adj(), phi.adj()});
    stan::math::recover_memory();
  }
  // var alpha, data phi (only theta partials)
  {
    var alpha = c.alpha;
    var lp
        = neg_binomial_2_log_glm_lpmf<propto>(c.y, c.x, alpha, c.beta, c.phi);
    lp.grad();
    print(c.tag + pt + "_var_alpha_only", {lp.val(), alpha.adj()});
    stan::math::recover_memory();
  }
  // var vector phi
  {
    var alpha = c.alpha;
    Matrix<var, Dynamic, 1> phi = c.phi_v.cast<var>();
    var lp
        = neg_binomial_2_log_glm_lpmf<propto>(c.y, c.x, alpha, c.beta, phi);
    lp.grad();
    print(c.tag + pt + "_var_vector_phi",
          {lp.val(), alpha.adj(), phi(0).adj(), phi(N - 1).adj()});
    stan::math::recover_memory();
  }
  // var x
  {
    Matrix<var, Dynamic, Dynamic> x = c.x.cast<var>();
    var lp
        = neg_binomial_2_log_glm_lpmf<propto>(c.y, x, c.alpha, c.beta, c.phi);
    lp.grad();
    print(c.tag + pt + "_var_x", {lp.val(), x(0, 0).adj(), x(N - 1, 2).adj()});
    stan::math::recover_memory();
  }
  // row-vector x (one row shared by all y)
  {
    var phi = c.phi;
    Matrix<var, 1, Dynamic> xr = c.x.row(0).cast<var>();
    var lp
        = neg_binomial_2_log_glm_lpmf<propto>(c.y, xr, c.alpha, c.beta, phi);
    lp.grad();
    print(c.tag + pt + "_var_rowx", {lp.val(), xr(0).adj(), phi.adj()});
    stan::math::recover_memory();
  }
  // scalar y with a single-row x
  {
    var phi = c.phi;
    MatrixXd x1 = c.x.topRows(1);
    var lp = neg_binomial_2_log_glm_lpmf<propto>(c.y[0], x1, c.alpha, c.beta,
                                                 phi);
    lp.grad();
    print(c.tag + pt + "_var_scalar_y", {lp.val(), phi.adj()});
    stan::math::recover_memory();
  }
  // fvar<double>: tangent in phi
  {
    fvar<double> phi(c.phi, 1.0);
    fvar<double> lp
        = neg_binomial_2_log_glm_lpmf<propto>(c.y, c.x, c.alpha, c.beta, phi);
    print(c.tag + pt + "_fvar_phi", {lp.val_, lp.d_});
  }
  // fvar<var>: second derivative in alpha
  {
    fvar<var> alpha(c.alpha, 1.0);
    fvar<var> lp
        = neg_binomial_2_log_glm_lpmf<propto>(c.y, c.x, alpha, c.beta, c.phi);
    lp.d_.grad();
    print(c.tag + pt + "_fvar_var_alpha",
          {lp.val_.val(), lp.d_.val(), alpha.val_.adj()});
    stan::math::recover_memory();
  }
  // fvar<fvar<double>>: second derivative in phi
  {
    fvar<fvar<double>> phi(fvar<double>(c.phi, 1.0), fvar<double>(1.0, 0.0));
    fvar<fvar<double>> lp
        = neg_binomial_2_log_glm_lpmf<propto>(c.y, c.x, c.alpha, c.beta, phi);
    print(c.tag + pt + "_ffd_phi", {lp.val_.val_, lp.d_.val_, lp.d_.d_});
  }
}

}  // namespace

int main(int argc, char** argv) {
  // "dump <tag>": print y, the rows of x, beta, alpha and phi of one case
  if (argc > 2 && std::string(argv[1]) == "dump") {
    for (const auto& c : make_cases()) {
      if (c.tag != argv[2]) {
        continue;
      }
      std::printf("alpha %a\nphi %a\n", c.alpha, c.phi);
      std::printf("beta %a %a %a\n", c.beta(0), c.beta(1), c.beta(2));
      for (size_t i = 0; i < c.y.size(); ++i) {
        std::printf("row %d %a %a %a\n", c.y[i], c.x(i, 0), c.x(i, 1),
                    c.x(i, 2));
      }
    }
    return 0;
  }
  for (const auto& c : make_cases()) {
    run_case<false>(c);
    run_case<true>(c);
  }
  return 0;
}
