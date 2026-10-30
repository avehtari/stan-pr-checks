// Claim M3 of README_emg.md: print the OpenCL device source of the
// exp_mod_normal lcdf functions (erfcx first, as the kernel generator
// includes it) from the tree in STAN_MATH, for cmp_host.cpp.
#define STAN_OPENCL
#include <stan/math/opencl/kernels/device_functions/erfcx.hpp>
#include <stan/math/opencl/kernels/device_functions/exp_mod_normal_lcdf.hpp>
#include <cstdio>

int main() {
  std::printf("%s\n%s\n", stan::math::opencl_kernels::erfcx_device_function,
              stan::math::opencl_kernels::exp_mod_normal_lcdf_device_function);
  return 0;
}
