// Claim M2 of README_emg.md: exp_mod_normal_lcdf / lccdf, value and
// gradients, on the dense grid: 50 z with |z| in 10^-4 .. 10^2, 25 v = lambda
// in 10^-10 .. 10^2 (1250 points, 2500 lines), mu = 0, sigma = 1.
// Output as probe_emg_phi.cpp.
#include <stan/math/rev.hpp>
#include <cmath>
#include <cstdio>
#include <exception>
#include <vector>

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
  std::vector<double> zs, vs;
  for (double lz = -4; lz <= 2.01; lz += 0.25) {
    zs.push_back(-std::pow(10.0, lz));
    zs.push_back(std::pow(10.0, lz));
  }
  for (double lv = -10; lv <= 2.01; lv += 0.5) {
    vs.push_back(std::pow(10.0, lv));
  }
  for (double z : zs) {
    for (double v : vs) {
      emg<false>(z, v);
      emg<true>(z, v);
    }
  }
  return 0;
}
