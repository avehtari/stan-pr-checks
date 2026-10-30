// Compare the device source (compiled as C on the host) with the prim
// kernel: worst relative difference per output over a dense (z, v) grid.
#include <stan/math/prim.hpp>
#include <cmath>
#include <cstdio>
#include <vector>

extern "C" double exp_mod_normal_lcdf_part(double z, double v, int upper,
                                           int which);

template <bool Upper>
void cmp(const std::vector<double>& zs, const std::vector<double>& vs) {
  double worst[4] = {0, 0, 0, 0};
  double wz[4], wv[4];
  int bad = 0;
  for (double z : zs) {
    for (double v : vs) {
      const auto r
          = stan::math::internal::exp_mod_normal_lcdf_kernel<Upper, true>(z,
                                                                         v);
      const double p[4] = {r.value, r.d_z, r.d_v, r.d_s};
      for (int k = 0; k < 4; ++k) {
        const double o = exp_mod_normal_lcdf_part(z, v, Upper, k);
        if (std::isfinite(p[k]) != std::isfinite(o)) {
          ++bad;
          continue;
        }
        const double s = std::fabs(p[k]);
        const double e = p[k] == o ? 0 : std::fabs(o - p[k]) / (s > 0 ? s : 1);
        if (e > worst[k]) {
          worst[k] = e;
          wz[k] = z;
          wv[k] = v;
        }
      }
    }
  }
  std::printf("%s finite mismatch %d\n", Upper ? "lccdf" : "lcdf", bad);
  for (int k = 0; k < 4; ++k) {
    std::printf("  out %d worst rel %.2e at z %g v %g\n", k, worst[k],
                worst[k] > 0 ? wz[k] : 0, worst[k] > 0 ? wv[k] : 0);
  }
}

int main() {
  std::vector<double> zs, vs;
  for (double lz = -4; lz <= 2.01; lz += 0.05) {
    zs.push_back(-std::pow(10.0, lz));
    zs.push_back(std::pow(10.0, lz));
  }
  zs.push_back(0.0);
  for (double lv = -10; lv <= 2.01; lv += 0.1) {
    vs.push_back(std::pow(10.0, lv));
  }
  cmp<false>(zs, vs);
  cmp<true>(zs, vs);
  return 0;
}
