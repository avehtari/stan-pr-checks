// probe_points.cpp --- the defect examples of stan-dev/math#3409 at one
// tree.  Uses only functions that exist on develop, so the same source
// compiles against the merge base and against the PR head.
//
// Prints one line per quantity: name, value (%.17g) or THROW <message>.

#include <stan/math/mix.hpp>

#include <cmath>
#include <cstdio>
#include <exception>
#include <functional>
#include <string>

namespace {

void show(const std::string& name, const std::function<double()>& f) {
  try {
    const double v = f();
    std::printf("%-52s % .17g\n", name.c_str(), v);
  } catch (const std::exception& e) {
    std::string msg = e.what();
    if (msg.size() > 70) {
      msg = msg.substr(0, 70);
    }
    std::printf("%-52s THROW %s\n", name.c_str(), msg.c_str());
  }
}

// both partials of the 9-argument root; which = 0 for g1, 1 for g2
double root9(double a, double b, double z, int which) {
  using stan::math::digamma;
  double g1 = 0;
  double g2 = 0;
  stan::math::grad_reg_inc_beta(g1, g2, a, b, z, digamma(a), digamma(b),
                                digamma(a + b), stan::math::beta(a, b));
  return which == 0 ? g1 : g2;
}

double dda(double a, double b, double z) {
  using stan::math::digamma;
  return stan::math::inc_beta_dda(a, b, z, digamma(a), digamma(a + b));
}

double ddb(double a, double b, double z) {
  using stan::math::digamma;
  return stan::math::inc_beta_ddb(a, b, z, digamma(b), digamma(a + b));
}

// gradient of a function of (y, alpha, beta) with var arguments;
// which = 0 y, 1 alpha, 2 beta
template <typename F>
double grad3(F f, double y, double al, double be, int which) {
  stan::math::nested_rev_autodiff nested;
  stan::math::var yv = y, av = al, bv = be;
  stan::math::var lp = f(yv, av, bv);
  lp.grad();
  const double g[3] = {yv.adj(), av.adj(), bv.adj()};
  return g[which];
}

}  // namespace

int main() {
  std::setvbuf(stdout, nullptr, _IOLBF, 0);

  std::printf("## roots\n");
  for (int w = 0; w < 2; ++w) {
    const char* g = w == 0 ? "g1" : "g2";
    show(std::string("root9 (600, 600, 0.5) ") + g,
         [=] { return root9(600, 600, 0.5, w); });
    show(std::string("root9 (352, 590, 0.758) ") + g,
         [=] { return root9(352, 590, 0.758, w); });
    show(std::string("root9 (1, 1, 1) ") + g,
         [=] { return root9(1, 1, 1, w); });
  }
  show("dda (1.145, 6786, 1.8e-4)", [] { return dda(1.145, 6786, 1.8e-4); });
  show("ddb (1.145, 6786, 1.8e-4)", [] { return ddb(1.145, 6786, 1.8e-4); });
  show("dda (446, 5, 0.76)", [] { return dda(446, 5, 0.76); });
  show("ddb (446, 5, 0.76)", [] { return ddb(446, 5, 0.76); });
  show("dda (1.5, 12500, 0.001)", [] { return dda(1.5, 12500, 0.001); });
  show("ddb (1.5, 12500, 0.001)", [] { return ddb(1.5, 12500, 0.001); });
  show("dda (600, 600, 0.5)", [] { return dda(600, 600, 0.5); });
  show("ddb (600, 600, 0.5)", [] { return ddb(600, 600, 0.5); });
  show("dda (8326, 11470, 0.246)", [] { return dda(8326, 11470, 0.246); });
  show("ddb (8326, 11470, 0.246)", [] { return ddb(8326, 11470, 0.246); });
  show("ddz (1.05, 0.35, 1.0)",
       [] { return stan::math::inc_beta_ddz(1.05, 0.35, 1.0); });
  show("ddz (0.5, 2.0, 0.0)",
       [] { return stan::math::inc_beta_ddz(0.5, 2.0, 0.0); });

  std::printf("## distributions, var arguments\n");
  auto beta_lcdf = [](const auto& y, const auto& a, const auto& b) {
    return stan::math::beta_lcdf(y, a, b);
  };
  auto beta_cdf = [](const auto& y, const auto& a, const auto& b) {
    return stan::math::beta_cdf(y, a, b);
  };
  auto nb_lcdf = [](const auto& n, const auto& a, const auto& b) {
    return stan::math::neg_binomial_lcdf(132, a, b);
  };
  auto inc_beta = [](const auto& z, const auto& a, const auto& b) {
    return stan::math::inc_beta(a, b, z);
  };
  const char* names[3] = {"d/dy", "d/dalpha", "d/dbeta"};
  for (int w = 0; w < 3; ++w) {
    show(std::string("beta_lcdf(0.5 | 600, 600) ") + names[w],
         [=] { return grad3(beta_lcdf, 0.5, 600, 600, w); });
    show(std::string("beta_lcdf(0.758 | 352, 590) ") + names[w],
         [=] { return grad3(beta_lcdf, 0.758, 352, 590, w); });
    show(std::string("beta_cdf(1.8e-4 | 1.145, 6786) ") + names[w],
         [=] { return grad3(beta_cdf, 1.8e-4, 1.145, 6786, w); });
    show(std::string("beta_cdf(0.76 | 446, 5) ") + names[w],
         [=] { return grad3(beta_cdf, 0.76, 446, 5, w); });
    show(std::string("inc_beta(600, 600, 0.5) ") + names[w],
         [=] { return grad3(inc_beta, 0.5, 600, 600, w); });
  }
  for (int w = 1; w < 3; ++w) {
    show(std::string("neg_binomial_lcdf(132 | 313, 0.068) ") + names[w],
         [=] { return grad3(nb_lcdf, 132, 313, 0.068, w); });
  }
  show("beta_cdf(0.76 | 446, 5) value",
       [] { return stan::math::beta_cdf(0.76, 446, 5); });
  show("neg_binomial_lcdf(132 | 313, 0.068) value",
       [] { return stan::math::neg_binomial_lcdf(132, 313, 0.068); });
  return 0;
}
