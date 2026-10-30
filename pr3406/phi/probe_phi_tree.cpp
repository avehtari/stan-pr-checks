// Claim N3 (Phi tables of 2026-10-19 and 2026-10-20): stan::math::Phi(double)
// of the tree in STAN_MATH on the 47001-point sweep of [-38.5, 8.5] of
// probe_phi_variants.cpp.  Output: lines "x Phi(x)" in hex, for score_phi.py.
#include <stan/math/prim/fun/Phi.hpp>
#include <cmath>
#include <cstdio>

int main() {
  const int n = 47001;
  for (int i = 0; i < n; ++i) {
    double t = (i + 0.5 * std::sin(1.0 + 7.0 * i)) / (n - 1);
    t = t < 0 ? 0 : (t > 1 ? 1 : t);
    const double x = -38.5 + 47.0 * t;
    std::printf("%a %a\n", x, stan::math::Phi(x));
  }
  return 0;
}
