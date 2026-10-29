// probe_series.cpp --- the two nu terms of student_t_lpdf in double:
//   C(nu) = lgamma((nu + 1) / 2) - lgamma(nu / 2) - log(nu) / 2
//   D(nu) = digamma((nu + 1) / 2) - digamma(nu / 2)
// at 2012 values h = nu / 2 from 20 to 1e15 (log-spaced grid), where
// student_t_lpdf uses the asymptotic series. Claim S8: "Against mpmath, from
// h = 20 to 1e15, the errors are at most 0.9 and 2.4 eps."
// Builds only against a tree that has internal::student_t_nu_constant_fun.
//
// Output line: <nu hex> => <C hex> <D hex>

#include <stan/math/prim.hpp>

#include <cmath>
#include <cstdio>

int main() {
  const int n = 2012;
  const double lo = std::log(20.0);
  const double hi = std::log(1e15);
  for (int i = 0; i < n; ++i) {
    // exp(log(20)) rounds to 19.999999999999996, below the series range;
    // clamp so that every point is at h >= 20
    const double h = std::fmax(20.0, std::exp(lo + (hi - lo) * i / (n - 1)));
    const double nu = 2.0 * h;
    const double c
        = stan::math::internal::student_t_nu_constant_fun::fun(nu);
    const double d = stan::math::internal::student_t_nu_digamma_fun::fun(nu);
    std::printf("%a => %a %a\n", nu, c, d);
  }
  return 0;
}
