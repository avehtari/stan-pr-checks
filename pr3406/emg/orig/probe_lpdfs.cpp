// skew_normal_lpdf and exp_mod_normal_lpdf of a vector of 200 var y with var
// mu, sigma, alpha / lambda: value and all partials in hex, for a bitwise
// comparison between two trees.
#include <stan/math/rev.hpp>
#include <cmath>
#include <cstdio>

using stan::math::var;

template <typename F>
void run(const char* name, F f) {
  for (double width : {3.0, 30.0}) {
    Eigen::Matrix<var, -1, 1> y(200);
    for (int i = 0; i < 200; ++i) {
      y[i] = 0.3 + 1.7 * width * std::sin(1.0 + 7.0 * i);
    }
    var mu = 0.3, sigma = 1.7, shape = 1.3;
    var lp = f(y, mu, sigma, shape);
    lp.grad();
    std::printf("%s %g %a %a %a %a", name, width, lp.val(), mu.adj(),
                sigma.adj(), shape.adj());
    for (int i = 0; i < 200; i += 7) {
      std::printf(" %a", y[i].adj());
    }
    std::printf("\n");
    stan::math::recover_memory();
  }
}

int main() {
  run("skew_normal_lpdf", [](const auto& y, const auto& m, const auto& s,
                             const auto& a) {
    return stan::math::skew_normal_lpdf(y, m, s, a);
  });
  run("exp_mod_normal_lpdf", [](const auto& y, const auto& m, const auto& s,
                                const auto& l) {
    return stan::math::exp_mod_normal_lpdf(y, m, s, l);
  });
  return 0;
}
