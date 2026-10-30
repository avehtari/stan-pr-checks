// Claim N2 (comment on 5f33a10d31): the before/after table, 2001 points per
// range.  Value and slope of internal::std_normal_lcdf_value_grad<true> over
// the whole range, in hex.  Compile against the tree in STAN_MATH.
// With -DCENTRAL_ERFCX also print the central-range variant that uses
// erfcx_small instead of the Cody erf rational:
//   log Phi(z) = LOG_HALF - z^2/2 + log(erfcx(-z/sqrt2)),
//   slope      = sqrt(2/pi) / erfcx(-z/sqrt2).
// Output: "# <name> <lo> <hi> <n>" blocks; lines "z value slope [v2 s2]".
#include <stan/math/prim/fun/std_normal_lcdf_impl.hpp>
#include <cmath>
#include <cstdio>
#include <cstdlib>

int main(int argc, char** argv) {
  namespace sm = stan::math;
  const struct {
    const char* name;
    double lo, hi;
  } ranges[] = {{"deep", -300.0, -40.0},        {"lower", -40.0, -5.66},
                {"seam_neg", -5.8, -5.5},       {"neg_mid", -5.66, -0.66},
                {"seam_c_neg", -0.70, -0.62},   {"central", -0.66, 0.66},
                {"seam_c_pos", 0.62, 0.70},     {"pos_mid", 0.66, 5.66},
                {"seam_pos", 5.5, 5.8},         {"upper", 5.66, 37.0}};
  // Optional argument: the points per range (default 2001).
  const int n = argc > 1 ? std::atoi(argv[1]) : 2001;
  for (const auto& r : ranges) {
    std::printf("# %s %.17g %.17g %d\n", r.name, r.lo, r.hi, n);
    for (int i = 0; i < n; ++i) {
      double t = (i + 0.5 * std::sin(1.0 + 7.0 * i)) / (n - 1);
      t = t < 0 ? 0 : (t > 1 ? 1 : t);
      const double z = r.lo + (r.hi - r.lo) * t;
      const auto vg = sm::internal::std_normal_lcdf_value_grad<true>(z);
      std::printf("%a %a %a", z, vg.first, vg.second);
#ifdef CENTRAL_ERFCX
      const double e = sm::internal::erfcx_small(-z * sm::INV_SQRT_TWO);
      std::printf(" %a %a", sm::LOG_HALF - (0.5 * z) * z + std::log(e),
                  sm::SQRT_TWO_OVER_SQRT_PI / e);
#endif
      std::printf("\n");
    }
  }
  return 0;
}
