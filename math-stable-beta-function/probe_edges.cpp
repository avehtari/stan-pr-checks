// Edge cases of beta(a, b): the value of stan::math::beta (the tree in
// STAN_MATH), exp(lbeta(a, b)), and the rule of the branch (exp(lbeta) for
// finite positive a and b, else the lgamma formula), on a cross product of
// inputs.
#include <stan/math/prim.hpp>
#include <cmath>
#include <cstdio>
#include <limits>
#include <string>
#include <vector>

static std::string fmt(double v) {
  char buf[64];
  std::snprintf(buf, sizeof(buf), "%.17g", v);
  return buf;
}

int main() {
  const double nan = std::numeric_limits<double>::quiet_NaN();
  const double inf = std::numeric_limits<double>::infinity();
  const std::vector<double> vals
      = {nan,  0.0,  -0.5, -1.0, -2.5, inf,  4.9e-324, 1e-300,
         1e-3, 0.5,  1.0,  2.0,  9.99, 10.0, 1e3,      1e9,     1.7e308};
  std::printf("# a b tree_beta exp_lbeta rule_beta lgamma_formula\n");
  for (double a : vals) {
    for (double b : vals) {
      std::string tree, el, rule, old;
      try {
        tree = fmt(stan::math::beta(a, b));
      } catch (const std::exception& e) {
        tree = "THROW";
      }
      try {
        el = fmt(std::exp(stan::math::lbeta(a, b)));
      } catch (const std::exception& e) {
        el = "THROW";
      }
      const double formula = std::exp(stan::math::lgamma(a)
                                      + stan::math::lgamma(b)
                                      - stan::math::lgamma(a + b));
      old = fmt(formula);
      if (a > 0 && b > 0 && a < inf && b < inf) {
        try {
          rule = fmt(std::exp(stan::math::lbeta(a, b)));
        } catch (const std::exception& e) {
          rule = "THROW";
        }
      } else {
        rule = old;
      }
      std::printf("%s %s %s %s %s %s\n", fmt(a).c_str(), fmt(b).c_str(),
                  tree.c_str(), el.c_str(), rule.c_str(), old.c_str());
    }
  }
  return 0;
}
