// Log Mills ratio log(Phi(-x) / phi(x)) as wiener4_lcdf_unnorm (logMill)
// computes it, std_normal_lcdf(-x) - std_normal_lpdf(x), and in the erfcx
// form log(sqrt(pi/2) * erfcx(x / sqrt(2))).  Score with score_mills.py.
#include <stan/math/prim.hpp>
#include <cmath>
#include <cstdio>
#include <vector>

int main() {
  using stan::math::erfcx;
  using stan::math::std_normal_lcdf;
  using stan::math::std_normal_lpdf;
  std::vector<double> xs{10.0, 30.0, 100.0, 300.0};
  // log grid from 1 to 1e4, 20 points per decade
  for (int k = 0; k <= 80; ++k) {
    xs.push_back(std::pow(10.0, k / 20.0));
  }
  const double sqrt_half_pi = std::sqrt(std::acos(-1.0) / 2);
  std::printf("# x lcdf(-x) lpdf(x) difference erfcx_form\n");
  for (double x : xs) {
    const double lc = std_normal_lcdf(-x);
    const double lp = std_normal_lpdf(x);
    const double diff = lc - lp;
    const double ex = std::log(sqrt_half_pi * erfcx(x / std::sqrt(2.0)));
    std::printf("%a %a %a %a %a\n", x, lc, lp, diff, ex);
  }
  return 0;
}
