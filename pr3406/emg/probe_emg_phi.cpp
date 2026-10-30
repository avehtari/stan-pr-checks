// Claims M1a and M2 of README_emg.md: exp_mod_normal_lcdf / lccdf, value and
// gradients, on the 170-point (z, v) grid; also Phi on a sweep.
// Compile against any tree.
// Grid: mu = 0, sigma = 1, so y = z and lambda = v.
// Output lines:
//   E <kind> z v value dy dmu dsigma dlambda    (hex; "nan" on throw)
//   P x Phi(x)                                  (hex)
#include <stan/math/rev.hpp>
#include <cmath>
#include <cstdio>
#include <exception>

using stan::math::var;

template <bool Upper>
void emg(double z, double v) {
  var y = z, mu = 0.0, sigma = 1.0, lambda = v;
  double out[5];
  try {
    var lp = Upper ? stan::math::exp_mod_normal_lccdf(y, mu, sigma, lambda)
                   : stan::math::exp_mod_normal_lcdf(y, mu, sigma, lambda);
    lp.grad();
    out[0] = lp.val();
    out[1] = y.adj();
    out[2] = mu.adj();
    out[3] = sigma.adj();
    out[4] = lambda.adj();
  } catch (const std::exception&) {
    for (double& o : out) {
      o = std::nan("");
    }
  }
  stan::math::recover_memory();
  std::printf("E %s %a %a %a %a %a %a %a\n", Upper ? "lccdf" : "lcdf", z, v,
              out[0], out[1], out[2], out[3], out[4]);
}

int main() {
  const double zs[] = {-1e4, -1e3, -300, -100, -30, -10, -5, -3, -1,
                       -0.3, 0,    0.3,  1,    3,   10,  30,  100};
  const double vs[] = {1e-8, 1e-6, 1e-4, 1e-2, 0.1, 0.5, 1, 3, 10, 100};
  for (double z : zs) {
    for (double v : vs) {
      emg<false>(z, v);
      emg<true>(z, v);
    }
  }
  for (int i = 0; i <= 940; ++i) {
    const double x = -38.5 + 0.05 * i;
    std::printf("P %a %a\n", x, stan::math::Phi(x));
  }
  return 0;
}
