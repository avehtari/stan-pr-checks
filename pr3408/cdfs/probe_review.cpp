// Probe for the stevieb-bot review of the incomplete gamma PR (#3408):
// series-branch roots, gamma_lccdf (150, 1200), the #2006 fvar example
// and gamma_lcdf(20000 | 800, 0.1).  Build with ../build.sh.
// Prints one line per quantity: tag, arguments, value.
#include <stan/math/mix.hpp>
#include <cmath>
#include <cstdio>

using stan::math::fvar;
using stan::math::var;

static void line(const char* tag, double a, double z, double v) {
  std::printf("%-34s %-12.6g %-12.6g % .17g\n", tag, a, z, v);
}

int main() {
  // T1, E5: grad_reg_inc_gamma in the series branch with a log z > 709
  for (auto az : {std::pair<double, double>{170, 100}, {165, 80}, {160, 120},
                  {150, 150}}) {
    double a = az.first, z = az.second;
    double r = stan::math::grad_reg_inc_gamma(a, z, std::tgamma(a),
                                              stan::math::digamma(a));
    line("T1 grad_reg_inc_gamma", a, z, r);
  }
  // E2: series branch with an exact g, points of the D5 claim
  for (auto az : {std::pair<double, double>{0.5, 0.25}, {8, 4}}) {
    double a = az.first, z = az.second;
    double r = stan::math::grad_reg_inc_gamma(a, z, std::tgamma(a),
                                              stan::math::digamma(a));
    line("E2 grad_reg_inc_gamma", a, z, r);
  }
  // T6: gamma_lccdf shape derivative, var
  for (auto ay : {std::pair<double, double>{150, 1200}, {10, 1000},
                  {200, 1600}}) {
    var alpha = ay.first;
    var lp = stan::math::gamma_lccdf(ay.second, alpha, 1.0);
    lp.grad();
    line("T6 gamma_lccdf d/dalpha", ay.first, ay.second, alpha.adj());
    stan::math::recover_memory();
  }
  // E1: fvar<double> derivative of log(gamma_q(alpha, beta * y)),
  // y = 12, alpha = 3, beta = 3, all tangents 1
  {
    fvar<double> y(12, 1), alpha(3, 1), beta(3, 1);
    fvar<double> f = stan::math::log(stan::math::gamma_q(alpha, beta * y));
    line("E1 log(gamma_q) value", 3, 36, f.val_);
    line("E1 log(gamma_q) tangent", 3, 36, f.d_);
  }
  // E4: gamma_lcdf(20000 | 800, 0.1)
  line("E4 gamma_lcdf(20000|800,0.1)", 800, 2000,
       stan::math::gamma_lcdf(20000.0, 800.0, 0.1));
  return 0;
}
