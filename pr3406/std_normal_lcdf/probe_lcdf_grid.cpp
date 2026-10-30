// Claim N1 (comment of 2026-09-18): worst error of std_normal_lcdf value and
// d/dy per range, develop against the PR.  Calls the public function with
// double (value) and var (value and gradient); prints hex for the scorer.
// Output: "# <name> <lo> <hi> <n>" blocks, lines
// "y value_double value_var grad grad_of_normal_lcdf(y|0,1)".
#include <stan/math/rev.hpp>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <vector>

namespace {

struct range {
  const char* name;
  double lo, hi;
};

void print_point(double y) {
  const double v = stan::math::std_normal_lcdf(y);
  stan::math::var yv = y;
  stan::math::var lp = stan::math::std_normal_lcdf(yv);
  lp.grad();
  const double g_std = yv.adj();
  const double v_var = lp.val();
  stan::math::recover_memory();
  // normal_lcdf(y | 0, 1) scales by 1 / (sigma sqrt 2), not by INV_SQRT_TWO,
  // so near a branch end it can select another branch than std_normal_lcdf.
  stan::math::var yn = y;
  stan::math::var ln = stan::math::normal_lcdf(yn, 0.0, 1.0);
  ln.grad();
  std::printf("%a %a %a %a %a\n", y, v, v_var, g_std, yn.adj());
  stan::math::recover_memory();
}

void print_block(const range& r, int n) {
  std::printf("# %s %.17g %.17g %d\n", r.name, r.lo, r.hi, n);
  for (int i = 0; i < n; ++i) {
    // Uniform grid with both ends included.
    const double y = i == n - 1 ? r.hi : r.lo + (r.hi - r.lo) * i / (n - 1);
    print_point(y);
  }
}

}  // namespace

int main(int argc, char** argv) {
  const int n = argc > 1 ? std::atoi(argv[1]) : 2001;
  const double s2 = std::sqrt(2.0);
  // The nine ranges of the 2026-09-18 comment, in y.
  const range comment[] = {{"c_deep", -300.0, -40.0},   {"c_m40_m10", -40.0, -10.0},
                           {"c_m10_m566", -10.0, -5.66}, {"c_seam_m58_m55", -5.8, -5.5},
                           {"c_m566_m1", -5.66, -1.0},   {"c_seam_m07_m06", -0.7, -0.6},
                           {"c_m06_06", -0.6, 0.6},      {"c_06_4", 0.6, 4.0},
                           {"c_4_10", 4.0, 10.0}};
  // The gradient branches of develop, given in s = y / sqrt(2), as y ranges.
  const struct {
    const char* name;
    double slo, shi;
  } branches[] = {{"d_AS_s_gt_2.9", 2.9, 6.0},       {"d_Taylor_2.7", 2.5, 2.9},
                  {"d_Taylor_2.3", 2.1, 2.5},        {"d_Taylor_1.85", 1.5, 2.1},
                  {"d_Taylor_1.15", 0.8, 1.5},       {"d_Taylor_0.45", 0.1, 0.8},
                  {"d_AS_nofit", -2.1, 0.1},         {"d_AS_fit_m2.1", -3.9, -2.1},
                  {"d_AS_fit_m3.9", -7.0, -3.9},     {"d_AS_fit_m7", -17.0, -7.0},
                  {"d_AS_fit_m17", -29.0, -17.0},    {"d_Mills", -60.0, -29.0}};
  for (const auto& r : comment) {
    print_block(r, n);
  }
  for (const auto& b : branches) {
    print_block({b.name, b.slo * s2, b.shi * s2}, n);
  }
  // The 15 points of expect_branch_accuracy (normal_lcdf_tail_test_helpers.hpp).
  const std::vector<double> test_points
      = {3.5496760415564683, 3.818376618407357,  4.1012193308819755,
         2.9839906166072305, 3.2526911934581184, 3.5355339059327378,
         2.135462479183374,  2.5455844122715714, 2.9698484809835,
         1.145512985522207,  1.6263455967290592, 2.121320343559643,
         0.15556349186104046, 0.6363961030678928, 1.1313708498984762};
  std::printf("# test_points 0 0 %zu\n", test_points.size());
  for (double y : test_points) {
    print_point(y);
  }
  return 0;
}
