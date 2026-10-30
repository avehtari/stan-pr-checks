// Claim M3 of README_emg.md: the OpenCL device code of the exp_mod_normal
// lcdf / lccdf, compiled on the host as C++, against the CPU kernel of the
// same tree.  Worst relative difference per output (value, d_z, d_v, d_s).
//
// Build: the device source comes from gen_dev_src (file dev_src.inc in a
// directory given with EXTRA_CXXFLAGS="-I <dir>").
// Usage: cmp_host [dense|fine]
//   dense: the 1250 (z, v) points of probe_emg_dense.cpp (2500 with both
//          functions); fine: the grid of rescued/cmp.cpp (243 z x 121 v).
#include <stan/math/prim.hpp>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <vector>

// OpenCL qualifiers and the one OpenCL extension pragma, for the host.
#define __kernel
#define __global
namespace dev {
using std::exp;
using std::expm1;
using std::fabs;
using std::fma;
using std::fmax;
using std::fmin;
using std::isfinite;
using std::isinf;
using std::isnan;
using std::log;
using std::log1p;
using std::sqrt;
#include "dev_src.inc"
}  // namespace dev

namespace internal = stan::math::internal;

// The CPU kernel; the template has (upper, calc_grad, log_value) at the PR
// head and (upper, calc_grad) at 39886e0eb8.
template <bool Upper>
auto prim_kernel(double z, double v, int)
    -> decltype(internal::exp_mod_normal_lcdf_kernel<Upper, true, true>(z, v)) {
  return internal::exp_mod_normal_lcdf_kernel<Upper, true, true>(z, v);
}
template <bool Upper>
auto prim_kernel(double z, double v, long)
    -> decltype(internal::exp_mod_normal_lcdf_kernel<Upper, true>(z, v)) {
  return internal::exp_mod_normal_lcdf_kernel<Upper, true>(z, v);
}

template <bool Upper>
void cmp(const std::vector<double>& zs, const std::vector<double>& vs) {
  const char* names[4] = {"value", "d_z", "d_v", "d_s"};
  double worst[4] = {0, 0, 0, 0};
  double wz[4] = {0, 0, 0, 0}, wv[4] = {0, 0, 0, 0};
  int mismatch = 0, n = 0;
  for (double z : zs) {
    for (double v : vs) {
      ++n;
      const auto r = prim_kernel<Upper>(z, v, 0);
      const double p[4] = {r.value, r.d_z, r.d_v, r.d_s};
      for (int k = 0; k < 4; ++k) {
        const double o = dev::exp_mod_normal_lcdf_part(z, v, Upper, k);
        if (std::isfinite(p[k]) != std::isfinite(o)) {
          ++mismatch;
          std::printf("  finite mismatch %s at z %g v %g: cpu %g device %g\n",
                      names[k], z, v, p[k], o);
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
  std::printf("%s: %d points, finite mismatch %d\n", Upper ? "lccdf" : "lcdf",
              n, mismatch);
  for (int k = 0; k < 4; ++k) {
    std::printf("  %-5s worst rel %.2e at z %g v %g\n", names[k], worst[k],
                wz[k], wv[k]);
  }
}

int main(int argc, char** argv) {
  const bool fine = argc > 1 && std::strcmp(argv[1], "fine") == 0;
  std::vector<double> zs, vs;
  if (fine) {
    for (double lz = -4; lz <= 2.01; lz += 0.05) {
      zs.push_back(-std::pow(10.0, lz));
      zs.push_back(std::pow(10.0, lz));
    }
    zs.push_back(0.0);
    for (double lv = -10; lv <= 2.01; lv += 0.1) {
      vs.push_back(std::pow(10.0, lv));
    }
  } else {
    for (double lz = -4; lz <= 2.01; lz += 0.25) {
      zs.push_back(-std::pow(10.0, lz));
      zs.push_back(std::pow(10.0, lz));
    }
    for (double lv = -10; lv <= 2.01; lv += 0.5) {
      vs.push_back(std::pow(10.0, lv));
    }
  }
  std::printf("grid %s\n", fine ? "fine" : "dense");
  cmp<false>(zs, vs);
  cmp<true>(zs, vs);
  return 0;
}
