// Claim N3: candidate implementations of Phi(double), for the accuracy probe
// and the speed bench.  Compile against the #3406 tree with commit 5f33a10d31.
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

// No C-library erf/erfc: Stan's erfcx pieces and one exp.
// Phi(x) = 0.5 exp(-u^2) erfcx(-u), u = x / sqrt 2; for x > 0 the
// complement 1 - 0.5 exp(-u^2) erfcx(u) has no cancellation (Phi >= 1/2).
inline double e(double x) {
  namespace smi = stan::math::internal;
  const double u = stan::math::INV_SQRT_TWO * x;
  const double a = u < 0 ? -u : u;
  // -(0.5 x) x rounds once; -u * u would round twice.
  const double g = 0.5 * std::exp(-(0.5 * x) * x);
  if (a < 0.46875) {
    return g * smi::erfcx_small(-u);
  }
  const double c = a >= 4.0 ? smi::erfcx_cody_tail(a) : smi::erfcx_cody_middle(a);
  return x < 0 ? g * c : 1.0 - g * c;
}

// e2: as e, but 0.5 + 0.5 erf(u) in the central range, with the Cody (1969)
// erf rational for |u| < 0.46875 (no exp; Phi is near 1/2, no cancellation).
inline double e2(double x) {
  namespace smi = stan::math::internal;
  const double u = stan::math::INV_SQRT_TWO * x;
  const double a = u < 0 ? -u : u;
  if (a < 0.46875) {
    const double s = u * u;
    const double num
        = (((1.85777706184603153e-1 * s + 3.16112374387056560) * s
            + 1.13864154151050156e2) * s + 3.77485237685302021e2) * s
          + 3.20937758913846947e3;
    const double den
        = (((s + 2.36012909523441209e1) * s + 2.44024637934444173e2) * s
           + 1.28261652607737228e3) * s + 2.84423683343917062e3;
    return 0.5 + 0.5 * (u * num / den);
  }
  const double g = 0.5 * std::exp(-(0.5 * x) * x);
  const double c = a >= 4.0 ? smi::erfcx_cody_tail(a) : smi::erfcx_cody_middle(a);
  return x < 0 ? g * c : 1.0 - g * c;
}

// Value and ratio phi(x) / Phi(x), as internal::std_normal_cdf_value_ratio.
// d1_vr: the helper of commit 0d4ac36e34 (erfc).
inline std::pair<double, double> d1_vr(double x) {
  const double cdf = 0.5 * std::erfc(-stan::math::INV_SQRT_TWO * x);
  const double ratio
      = cdf > 0 ? stan::math::INV_SQRT_TWO_PI * std::exp(-(0.5 * x) * x) / cdf
                : 0.0;
  return {cdf, ratio};
}

// e_vr: erfcx.  For x < 0 the ratio is sqrt(2 / pi) / erfcx(-u), with no
// exp and no 0 / 0 where Phi underflows.
inline std::pair<double, double> e_vr(double x) {
  namespace smi = stan::math::internal;
  const double u = stan::math::INV_SQRT_TWO * x;
  const double a = u < 0 ? -u : u;
  const double g = 0.5 * std::exp(-(0.5 * x) * x);
  if (a < 0.46875 || x < 0) {
    const double c = a < 0.46875 ? smi::erfcx_small(-u)
                     : a >= 4.0  ? smi::erfcx_cody_tail(a)
                                 : smi::erfcx_cody_middle(a);
    return {g * c, stan::math::SQRT_TWO_OVER_SQRT_PI / c};
  }
  const double c = a >= 4.0 ? smi::erfcx_cody_tail(a) : smi::erfcx_cody_middle(a);
  const double cdf = 1.0 - g * c;
  return {cdf, stan::math::SQRT_TWO_OVER_SQRT_PI * g / cdf};
}

}  // namespace phi_variants
#endif
