// beta(a, b) of the tree in STAN_MATH at the points of the PR summary
// table.  Output: a, b, beta(a, b) in hex.  Scored by refs/score_table.py.
#include <stan/math/prim.hpp>
#include <cmath>
#include <cstdio>

int main() {
  const double big = std::ldexp(1.0, 31);
  const double pts[][2] = {{0.5, 1e6}, {0.5, 1e9}, {2.0, 1e9},
                           {2.0, big}, {10.0, 1e9}, {0.5, 1.7e308}};
  for (const auto& p : pts) {
    std::printf("%a %a %a\n", p[0], p[1], stan::math::beta(p[0], p[1]));
  }
  return 0;
}
