// Bundle: value and var gradient of 17 cdfs at the 306 survey cases.
// Supports the 2e-13 / 8e-13 / 2.7e-10 claim of #3408 Part 3.
// Value and var gradient of the incomplete-gamma distribution functions at
// the cases of survey_cases.txt. Output per line: value grads... (hex), or
// THROW. Gradient order: y, p1, p2 (poisson: lambda only).
#include <stan/math/mix.hpp>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <functional>
#include <map>
#include <sstream>
#include <string>
#include <vector>

using stan::math::var;
namespace sm = stan::math;

using F3 = std::function<var(const var&, const var&, const var&)>;
using F2 = std::function<var(const var&, const var&)>;

int main(int argc, char** argv) {
  const std::map<std::string, F3> f3 = {
      {"gamma_lcdf", [](auto& y, auto& a, auto& b) { return sm::gamma_lcdf(y, a, b); }},
      {"gamma_lccdf", [](auto& y, auto& a, auto& b) { return sm::gamma_lccdf(y, a, b); }},
      {"gamma_cdf", [](auto& y, auto& a, auto& b) { return sm::gamma_cdf(y, a, b); }},
      {"inv_gamma_lcdf", [](auto& y, auto& a, auto& b) { return sm::inv_gamma_lcdf(y, a, b); }},
      {"inv_gamma_lccdf", [](auto& y, auto& a, auto& b) { return sm::inv_gamma_lccdf(y, a, b); }},
      {"inv_gamma_cdf", [](auto& y, auto& a, auto& b) { return sm::inv_gamma_cdf(y, a, b); }},
      {"scaled_inv_chi_square_cdf",
       [](auto& y, auto& n, auto& s) { return sm::scaled_inv_chi_square_cdf(y, n, s); }},
      {"scaled_inv_chi_square_lcdf",
       [](auto& y, auto& n, auto& s) { return sm::scaled_inv_chi_square_lcdf(y, n, s); }},
      {"scaled_inv_chi_square_lccdf",
       [](auto& y, auto& n, auto& s) { return sm::scaled_inv_chi_square_lccdf(y, n, s); }},
  };
  const std::map<std::string, F2> f2 = {
      {"chi_square_lcdf", [](auto& y, auto& n) { return sm::chi_square_lcdf(y, n); }},
      {"chi_square_lccdf", [](auto& y, auto& n) { return sm::chi_square_lccdf(y, n); }},
      {"chi_square_cdf", [](auto& y, auto& n) { return sm::chi_square_cdf(y, n); }},
      {"inv_chi_square_cdf", [](auto& y, auto& n) { return sm::inv_chi_square_cdf(y, n); }},
      {"inv_chi_square_lcdf", [](auto& y, auto& n) { return sm::inv_chi_square_lcdf(y, n); }},
      {"inv_chi_square_lccdf", [](auto& y, auto& n) { return sm::inv_chi_square_lccdf(y, n); }},
  };
  std::ifstream in(argc > 1 ? argv[1] : "survey_cases.txt");
  std::string line;
  while (std::getline(in, line)) {
    std::istringstream ss(line);
    std::string name;
    double y, p1, p2;
    ss >> name >> y >> p1 >> p2;
    std::vector<double> out;
    try {
      if (name == "poisson_lcdf" || name == "poisson_lccdf") {
        var lam = p1;
        const int n = static_cast<int>(y);
        var lp = name == "poisson_lcdf" ? sm::poisson_lcdf(n, lam)
                                        : sm::poisson_lccdf(n, lam);
        lp.grad();
        out = {lp.val(), lam.adj()};
      } else if (f3.count(name)) {
        var yv = y, a = p1, b = p2;
        var lp = f3.at(name)(yv, a, b);
        lp.grad();
        out = {lp.val(), yv.adj(), a.adj(), b.adj()};
      } else {
        var yv = y, a = p1;
        var lp = f2.at(name)(yv, a);
        lp.grad();
        out = {lp.val(), yv.adj(), a.adj()};
      }
    } catch (const std::exception& e) {
      std::printf("THROW\n");
      sm::recover_memory();
      continue;
    }
    sm::recover_memory();
    for (double v : out) {
      std::printf("%a ", v);
    }
    std::printf("\n");
  }
  return 0;
}
