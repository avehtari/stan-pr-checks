// beta(a, b) of the tree in STAN_MATH on the sweep grid (grid.hpp) and on
// the edge inputs of probe_edges.cpp.  Output: one line per point, "G" or
// "E", then a, b and beta(a, b) in hex.  Scored by score_sweep.py.
#include <stan/math/prim.hpp>
#include <cstdio>
#include <limits>
#include <vector>
#include "grid.hpp"

int main() {
  const std::vector<double> axis = sweep_axis();
  for (double a : axis) {
    for (double b : axis) {
      std::printf("G %a %a %a\n", a, b, stan::math::beta(a, b));
    }
  }
  const double nan = std::numeric_limits<double>::quiet_NaN();
  const double inf = std::numeric_limits<double>::infinity();
  const std::vector<double> edges
      = {nan,  0.0, -0.5, -1.0, -2.5, inf,  4.9e-324, 1e-300,
         1e-3, 0.5, 1.0,  2.0,  9.99, 10.0, 1e3,      1e9,     1.7e308};
  for (double a : edges) {
    for (double b : edges) {
      double v;
      try {
        v = stan::math::beta(a, b);
      } catch (const std::exception&) {
        std::printf("E %a %a THROW\n", a, b);
        continue;
      }
      std::printf("E %a %a %a\n", a, b, v);
    }
  }
  return 0;
}
