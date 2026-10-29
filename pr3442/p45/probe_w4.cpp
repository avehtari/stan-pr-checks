// probe_w4.cpp -- the 4- and 5-parameter Wiener functions at the p4 and p5
// point sets.  Claims of stan-dev/math#3442: 4-parameter wiener_lpdf 3.6e-9 and
// 5.2e-8; 5-parameter errors to 1e-2; 4-parameter lccdf -inf at 8 of 300.
//
// Point sets (t0 = 0.1, y = t0 + t):
//   p4: 300 points, a 0.5..3, v -5..5, w 0.2..0.8, t 0.05..3.05, sv = 0
//   p5: 300 points, the same ranges and sv 0..1.5
// Rows (value, then the gradients in y a t0 w v sv; 0 where the function
// has no sv; then the median time of 3 calls in microseconds):
//   OLD    wiener_lpdf(y, a, t0, w, v) (the older 4-parameter code), p4
//   LPDF   wiener_lpdf(y, a, t0, w, v, sv), p4 and p5
//   LCDF   wiener_lcdf_unnorm(y, a, t0, w, v), p4;
//          wiener_lcdf_unnorm(y, a, t0, w, v, sv, 0, 0), p5
//   LCCDF  the same for wiener_lccdf_unnorm
// Output: <fn> <set> <id> <y a t0 w v sv> <value> <6 gradients> <us>, or
//         <fn> <set> <id> <params> THROW <msg>. Default precision (1e-4).
#include <stan/math/rev.hpp>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <random>
#include <string>
#include <vector>

using stan::math::var;

namespace {

template <typename F>
void row(const char* fn, const char* set, int id, const std::vector<double>& p,
         const F& f) {
  std::printf("%s\t%s\t%s%d", fn, set, set, id);
  for (double q : p) {
    std::printf("\t%.17g", q);
  }
  try {
    std::vector<double> times, g;
    double val = 0;
    for (int r = 0; r < 3; ++r) {
      const auto t0 = std::chrono::steady_clock::now();
      std::vector<var> x(p.begin(), p.end());
      var lp = f(x);
      lp.grad(x, g);
      val = lp.val();
      stan::math::recover_memory();
      times.push_back(
          std::chrono::duration<double>(std::chrono::steady_clock::now() - t0)
              .count());
    }
    std::sort(times.begin(), times.end());
    std::printf("\t%.17g", val);
    for (double q : g) {
      std::printf("\t%.17g", q);
    }
    std::printf("\t%.3f\n", 1e6 * times[1]);
  } catch (const std::exception& e) {
    stan::math::recover_memory();
    std::printf("\tTHROW\t%s\n", e.what());
  }
}

}  // namespace

int main() {
  std::setvbuf(stdout, nullptr, _IONBF, 0);
  std::mt19937 rng(20261010);
  std::uniform_real_distribution<double> U(0.0, 1.0);
  for (const char* set : {"p4", "p5"}) {
    const bool has_sv = set[1] == '5';
    for (int k = 0; k < 300; ++k) {
      const double a = 0.5 + 2.5 * U(rng);
      const double v = -5.0 + 10.0 * U(rng);
      const double w = 0.2 + 0.6 * U(rng);
      const double t = 0.05 + 3.0 * U(rng);
      const double sv = has_sv ? 1.5 * U(rng) : 0.0;
      const std::vector<double> p{0.1 + t, a, 0.1, w, v, sv};
      if (!has_sv) {
        row("OLD", set, k, p, [](std::vector<var>& x) {
          x[5] = 0.0;
          return stan::math::wiener_lpdf(x[0], x[1], x[2], x[3], x[4]);
        });
      }
      row("LPDF", set, k, p, [](std::vector<var>& x) {
        return stan::math::wiener_lpdf(x[0], x[1], x[2], x[3], x[4], x[5]);
      });
      if (!has_sv) {
        row("LCDF", set, k, p, [](std::vector<var>& x) {
          return stan::math::wiener_lcdf_unnorm(x[0], x[1], x[2], x[3], x[4]);
        });
        row("LCCDF", set, k, p, [](std::vector<var>& x) {
          return stan::math::wiener_lccdf_unnorm(x[0], x[1], x[2], x[3], x[4]);
        });
      } else {
        row("LCDF", set, k, p, [](std::vector<var>& x) {
          return stan::math::wiener_lcdf_unnorm(x[0], x[1], x[2], x[3], x[4],
                                                x[5], 0.0, 0.0);
        });
        row("LCCDF", set, k, p, [](std::vector<var>& x) {
          return stan::math::wiener_lccdf_unnorm(x[0], x[1], x[2], x[3], x[4],
                                                 x[5], 0.0, 0.0);
        });
      }
    }
  }
  return 0;
}
