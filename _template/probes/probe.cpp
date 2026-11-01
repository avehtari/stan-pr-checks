// <Purpose in one line, and the claim of the PR text that this probe backs.>
// Output: one line per point, "<inputs> value d/d<arg1> ... status", in
// %.17g, for refs/score.py.  Build against each tree of out/trees/ with
// build.sh; the same source must compile against develop and the PR.
//
// Rules learned on earlier PRs (STEVIEB_BOT_REPLIES.md section 5c):
// - call each function in its own try block, so that a throw marks only
//   its own column;
// - print values with %.17g (never %.6g), so that equality claims such as
//   "equals log(2^-53)" can be checked;
// - use exactly representable inputs where a claim depends on rounding
//   (for example x = 2^-29 instead of 2 / 1e9).
#include <stan/math/rev.hpp>
#include <cmath>
#include <cstdio>
#include <vector>

namespace {

// Replace with the function under test.  Arguments are var, so that one
// reverse pass gives every partial.
stan::math::var f(const stan::math::var& x, const stan::math::var& a) {
  return stan::math::lgamma(x + a) - stan::math::lgamma(x);
}

void point(double x0, double a0) {
  double value = NAN, dx = NAN, da = NAN;
  const char* status = "ok";
  try {
    stan::math::var x = x0, a = a0;
    stan::math::var y = f(x, a);
    y.grad();
    value = y.val();
    dx = x.adj();
    da = a.adj();
  } catch (const std::exception& e) {
    status = "throw";
  }
  stan::math::recover_memory();
  std::printf("%.17g %.17g %.17g %.17g %.17g %s\n", x0, a0, value, dx, da,
              status);
}

}  // namespace

int main() {
  // A log grid over the regions of the claim; add the edge points
  // (0, denormals, inf, NaN, negative values) separately.
  for (double x : {1e-3, 0.5, 1.0, 10.0, 1e3, 1e6, 1e9}) {
    for (double a : {0.5, 2.0, 10.0}) {
      point(x, a);
    }
  }
  return 0;
}
