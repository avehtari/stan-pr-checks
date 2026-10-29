// probe_shift.cpp -- the 7-parameter wiener_lpdf (value and 8 gradients) at the
// #3329 point and scan, the 300 D points and the 128-point extreme grid.
// Claims of stan-dev/math#3442: the #3329 row of the accuracy table; d/dsw 2e3.
//
// Output rows, tab separated:
//   <section> <id> <y a t0 w v sv sw st0> <prec> <value> <g1..g8> | THROW msg
// Sections:
//   P  the #3329 point (3 variants)
//   S  scan in y across the former -inf boundary, sw = 1e-8, st0 = 1e-3
//   R  reference rows: R <id> <y> <log(F5(y - t0) / st0)> <p5 / F5>, by
//      Boost Gauss-Kronrod on the log-scaled 5-parameter density (sw -> 0)
//   D  random points in the ordinary support, precision 1e-4 and 1e-10
//   E  extreme grid from the #3329 sweep locations
#include <stan/math/rev.hpp>
#include <boost/math/quadrature/gauss_kronrod.hpp>
#include <cmath>
#include <cstdio>
#include <random>
#include <string>
#include <vector>

using stan::math::var;

namespace {

void row(const char* sec, const std::string& id, const std::vector<double>& p,
         double prec) {
  std::printf("%s\t%s", sec, id.c_str());
  for (double v : p) {
    std::printf("\t%.17g", v);
  }
  std::printf("\t%.3g", prec);
  try {
    std::vector<var> x(p.begin(), p.end());
    var lp = stan::math::wiener_lpdf(x[0], x[1], x[2], x[3], x[4], x[5], x[6],
                                     x[7], prec);
    std::vector<double> g;
    lp.grad(x, g);
    std::printf("\t%.17g", lp.val());
    for (double v : g) {
      std::printf("\t%.17g", v);
    }
    std::printf("\n");
  } catch (const std::exception& e) {
    std::printf("\tTHROW\t%s\n", e.what());
  }
  stan::math::recover_memory();
}

double lp5(double s, double a, double w, double v, double sv) {
  return stan::math::wiener_lpdf(s, a, 0.0, w, v, sv, 1e-12);
}

double log_F5(double T, double a, double w, double v, double sv) {
  double M = -INFINITY;
  for (int i = 1; i <= 500; ++i) {
    M = std::fmax(M, lp5(T * i / 500.0, a, w, v, sv));
  }
  auto f = [&](double s) {
    return s <= 0 ? 0.0 : std::exp(lp5(s, a, w, v, sv) - M);
  };
  double err = 0;
  return M
         + std::log(boost::math::quadrature::gauss_kronrod<double, 61>::integrate(
             f, 0.0, T, 12, 1e-11, &err));
}

}  // namespace

int main() {
  std::setvbuf(stdout, nullptr, _IONBF, 0);
  const double a = 0.1, t0 = 0.0, w = 1e-4, v = -10.0, sv = 0.1;
  row("P", "sw_only", {1e-6, a, t0, w, v, sv, 1e-8, 0.0}, 1e-8);
  row("P", "st0_only", {1e-6, a, t0, w, v, sv, 0.0, 1e-3}, 1e-8);
  row("P", "both", {1e-6, a, t0, w, v, sv, 1e-8, 1e-3}, 1e-8);
  for (int k = 0; k <= 15; ++k) {
    const double y = 3e-6 * std::pow(10.0, k / 15.0);
    const std::string id = "s" + std::to_string(k);
    row("S", id, {y, a, t0, w, v, sv, 1e-8, 1e-3}, 1e-8);
    const double lF = log_F5(y - t0, a, w, v, sv);
    std::printf("R\t%s\t%.17g\t%.17g\t%.17g\n", id.c_str(), y,
                lF - std::log(1e-3), std::exp(lp5(y, a, w, v, sv) - lF));
  }
  const double lF_point = log_F5(1e-6, a, w, v, sv);
  std::printf("R\tpoint\t%.17g\t%.17g\t%.17g\n", 1e-6, lF_point - std::log(1e-3),
              std::exp(lp5(1e-6, a, w, v, sv) - lF_point));
  std::printf("R\tfive_param\t%.17g\t%.17g\t0\n", 1e-6, lp5(1e-6, a, w, v, sv));

  std::mt19937 rng(20261009);
  std::uniform_real_distribution<double> U(0.0, 1.0);
  for (int k = 0; k < 300; ++k) {
    const double aa = 0.5 + 2.5 * U(rng);
    const double vv = -3.0 + 6.0 * U(rng);
    const double ww = 0.2 + 0.6 * U(rng);
    const double tt0 = 0.1 + 0.3 * U(rng);
    const double svv = 1.5 * U(rng);
    const double sww
        = (k % 3 == 2) ? 0.0 : 0.8 * U(rng) * 2.0 * std::fmin(ww, 1 - ww);
    const double stt0 = (k % 3 == 1) ? 0.0 : 0.3 * U(rng);
    const double yy = tt0 + 0.05 + 2.0 * U(rng);
    const std::vector<double> p{yy, aa, tt0, ww, vv, svv, sww, stt0};
    const std::string id = "d" + std::to_string(k);
    row("D", id, p, 1e-4);
    row("D", id, p, 1e-10);
  }

  int k = 0;
  for (double yd : {1e-6, 1e-4, 10.0, 50.0}) {
    for (double aa : {0.1, 8.0}) {
      for (double vv : {-30.0, -10.0, 10.0, 30.0}) {
        for (double ww : {1e-4, 0.5}) {
          for (double stt0 : {0.0, 1e-3}) {
            const double sww = 0.5 * std::fmin(ww, 1 - ww);
            row("E", "e" + std::to_string(k++),
                {0.2 + yd, aa, 0.2, ww, vv, 0.1, sww, stt0}, 1e-6);
          }
        }
      }
    }
  }
  return 0;
}
