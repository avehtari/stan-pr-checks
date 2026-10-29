// probe_oldw4.cpp -- the older 4-parameter wiener_lpdf(y, a, t0, w, v) at
// small densities.  Its series lengths come from an absolute error bound
// (1e-6) on the standardized density.  Grid g4 (t0 = 0.1, y = t0 + t):
//   t in {1e-3, 0.01, 0.05, 0.2, 1, 3, 10, 30}, a in {0.5, 1, 2, 4, 8},
//   v in {-20, -5, 0, 5, 20}, w in {0.1, 0.5, 0.9}; 600 points, sv = 0.
// Rows as in probe_w4.cpp (ref_w4_mpmath.py reads the LPDF rows):
//   OLD   wiener_lpdf(y, a, t0, w, v)
//   LPDF  wiener_lpdf(y, a, t0, w, v, 0) (5-parameter code)
// Output: <fn> <set> <id> <y a t0 w v sv> <value> <6 gradients> <us>, or
//         <fn> <set> <id> <params> THROW <msg>.
#include <stan/math/rev.hpp>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <vector>

using stan::math::var;

namespace {

template <typename F>
void row(const char* fn, int id, const std::vector<double>& p, const F& f) {
  std::printf("%s\tg4\tg4%d", fn, id);
  for (double q : p) {
    std::printf("\t%.17g", q);
  }
  try {
    const auto start = std::chrono::steady_clock::now();
    std::vector<var> x(p.begin(), p.end());
    var lp = f(x);
    std::vector<double> g;
    lp.grad(x, g);
    const double us = 1e6
                      * std::chrono::duration<double>(
                            std::chrono::steady_clock::now() - start)
                            .count();
    std::printf("\t%.17g", lp.val());
    for (double q : g) {
      std::printf("\t%.17g", q);
    }
    std::printf("\t%.3f\n", us);
  } catch (const std::exception& e) {
    std::printf("\tTHROW\t%s\n", e.what());
  }
  stan::math::recover_memory();
}

}  // namespace

int main() {
  std::setvbuf(stdout, nullptr, _IONBF, 0);
  int id = 0;
  for (double t : {1e-3, 0.01, 0.05, 0.2, 1.0, 3.0, 10.0, 30.0}) {
    for (double a : {0.5, 1.0, 2.0, 4.0, 8.0}) {
      for (double v : {-20.0, -5.0, 0.0, 5.0, 20.0}) {
        for (double w : {0.1, 0.5, 0.9}) {
          const std::vector<double> p{0.1 + t, a, 0.1, w, v, 0.0};
          row("OLD", id, p, [](std::vector<var>& x) {
            x[5] = 0.0;
            return stan::math::wiener_lpdf(x[0], x[1], x[2], x[3], x[4]);
          });
          row("LPDF", id, p, [](std::vector<var>& x) {
            return stan::math::wiener_lpdf(x[0], x[1], x[2], x[3], x[4], x[5]);
          });
          ++id;
        }
      }
    }
  }
  return 0;
}
