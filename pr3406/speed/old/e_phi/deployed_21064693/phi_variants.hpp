// Candidate implementations of Phi(double), for the accuracy probe and the
// speed bench.  Compile against the #3406 tree with commit 5f33a10d31.
//   dev  develop: erfc below -5, 0.5 (1 + erf) above (copy)
//   pr   #3406: exp(log Phi)
//   d1   0.5 erfc(-x / sqrt 2) everywhere
//   d2   0.5 erfc(-x / sqrt 2) below -0.67, 0.5 + 0.5 erf(x / sqrt 2) above
#ifndef PHI_VARIANTS_HPP
#define PHI_VARIANTS_HPP
#include <stan/math/prim/fun/std_normal_lcdf_impl.hpp>
#include <cmath>

namespace phi_variants {

inline double dev(double x) {
  const double c = stan::math::INV_SQRT_TWO;
  if (x < -37.5) {
    return 0;
  } else if (x < -5.0) {
    return 0.5 * std::erfc(-c * x);
  } else if (x > 8.25) {
    return 1;
  }
  return 0.5 * (1.0 + std::erf(c * x));
}

inline double pr(double x) {
  return std::exp(stan::math::internal::std_normal_lcdf_value_grad<false>(x));
}

inline double d1(double x) {
  return 0.5 * std::erfc(-stan::math::INV_SQRT_TWO * x);
}

inline double d2(double x) {
  const double u = stan::math::INV_SQRT_TWO * x;
  return x < -0.67 ? 0.5 * std::erfc(-u) : 0.5 + 0.5 * std::erf(u);
}

}  // namespace phi_variants
#endif
