// Value and d/dy of the 15 functions of the "valid ranges" table of
// stan-dev/math#3406, as a line server for ranges/probe_ranges.py.
// Input lines: "<function> <y in %a>".  Output: "<value %a> <d/dy %a>".
#include <stan/math/rev.hpp>
#include <cstdio>
#include <cstdlib>
#include <functional>
#include <iostream>
#include <map>
#include <string>

using stan::math::var;

// ── Parameters: mu = 0, sigma = 1; exp_mod_normal lambda = 1; ──
// ── skew_normal xi = 0, omega = 1, alpha = 1.                   ──
template <typename T>
using fun_t = std::function<T(const T&)>;

struct entry {
  fun_t<double> fd;
  fun_t<var> fv;
};

#define ENTRY(EXPR)                                         \
  entry {                                                   \
    [](const double& y) -> double { return EXPR; },         \
        [](const var& y) -> var { return EXPR; }            \
  }

int main() {
  using namespace stan::math;
  std::map<std::string, entry> table{
      {"std_normal_lcdf", ENTRY(std_normal_lcdf(y))},
      {"normal_lcdf", ENTRY(normal_lcdf(y, 0.0, 1.0))},
      {"std_normal_lccdf", ENTRY(std_normal_lccdf(y))},
      {"normal_lccdf", ENTRY(normal_lccdf(y, 0.0, 1.0))},
      {"std_normal_cdf", ENTRY(std_normal_cdf(y))},
      {"normal_cdf", ENTRY(normal_cdf(y, 0.0, 1.0))},
      {"Phi", ENTRY(Phi(y))},
      {"lognormal_lcdf", ENTRY(lognormal_lcdf(y, 0.0, 1.0))},
      {"lognormal_lccdf", ENTRY(lognormal_lccdf(y, 0.0, 1.0))},
      {"lognormal_cdf", ENTRY(lognormal_cdf(y, 0.0, 1.0))},
      {"exp_mod_normal_lpdf", ENTRY(exp_mod_normal_lpdf(y, 0.0, 1.0, 1.0))},
      {"exp_mod_normal_lcdf", ENTRY(exp_mod_normal_lcdf(y, 0.0, 1.0, 1.0))},
      {"exp_mod_normal_lccdf", ENTRY(exp_mod_normal_lccdf(y, 0.0, 1.0, 1.0))},
      {"exp_mod_normal_cdf", ENTRY(exp_mod_normal_cdf(y, 0.0, 1.0, 1.0))},
      {"skew_normal_lpdf", ENTRY(skew_normal_lpdf(y, 0.0, 1.0, 1.0))},
  };
  std::string name, xs;
  while (std::cin >> name >> xs) {
    auto it = table.find(name);
    if (it == table.end()) {
      std::printf("unknown unknown\n");
      std::fflush(stdout);
      continue;
    }
    const double y = std::strtod(xs.c_str(), nullptr);
    double v = NOT_A_NUMBER, g = NOT_A_NUMBER;
    // A throw prints as nan; the scorer counts it as invalid.
    try {
      v = it->second.fd(y);
    } catch (const std::exception&) {
    }
    try {
      var yv = y;
      var lp = it->second.fv(yv);
      lp.grad();
      g = yv.adj();
    } catch (const std::exception&) {
    }
    recover_memory();
    std::printf("%a %a\n", v, g);
    std::fflush(stdout);
  }
  return 0;
}
