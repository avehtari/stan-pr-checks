// neg_binomial_2_lcdf at large finite counts, where the second beta_lcdf
// shape n + 1 is large. The PR does not change this path; the probe
// measures the limit that the PR text names. Output: n mu phi value d/dmu
// d/dphi, one line per point, for score_large_n.py.
#include <stan/math/rev.hpp>
#include <cstdio>
#include <vector>

int main() {
  const double phi0 = 2.0;
  const std::vector<std::pair<double, double>> points = {
      {10, 4},     {1e3, 4},       {1e6, 4},   {1e9, 4},
      {2147483646, 4}, {1e3, 1e3}, {1e6, 1e6}, {1e9, 1e9},
      {2147483646, 2147483646}};
  for (const auto& pt : points) {
    const int n = static_cast<int>(pt.first);
    stan::math::var mu = pt.second, phi = phi0;
    stan::math::var lp = stan::math::neg_binomial_2_lcdf(n, mu, phi);
    lp.grad();
    std::printf("%d %.17g %.17g %.17g %.17g %.17g\n", n, pt.second, phi0,
                lp.val(), mu.adj(), phi.adj());
    stan::math::recover_memory();
  }
  return 0;
}
