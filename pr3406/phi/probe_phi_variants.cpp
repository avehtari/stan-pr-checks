// Claim N3: the Phi variants of the 2026-10-19/20 tables on a sweep of
// [-38.5, 8.5]: "x dev pr d1 d2 e e2" in hex.  Build against 5f33a10d31.
#include "phi_variants.hpp"
#include <cmath>
#include <cstdio>

int main() {
  namespace pv = phi_variants;
  const int n = 47001;
  for (int i = 0; i < n; ++i) {
    double t = (i + 0.5 * std::sin(1.0 + 7.0 * i)) / (n - 1);
    t = t < 0 ? 0 : (t > 1 ? 1 : t);
    const double x = -38.5 + 47.0 * t;
    std::printf("%a %a %a %a %a %a %a\n", x, pv::dev(x), pv::pr(x), pv::d1(x),
                pv::d2(x), pv::e(x), pv::e2(x));
  }
  return 0;
}
