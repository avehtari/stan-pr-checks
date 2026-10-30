// Print the device source of the exp_mod_normal lcdf functions, as the
// kernel generator includes it (erfcx first).
#include <stan/math/opencl/kernels/device_functions/erfcx.hpp>
#include <stan/math/opencl/kernels/device_functions/exp_mod_normal_lcdf.hpp>
#include <cstdio>

int main() {
  std::printf("#pragma OPENCL EXTENSION cl_khr_fp64 : enable\n%s\n%s\n",
              stan::math::opencl_kernels::erfcx_device_function,
              stan::math::opencl_kernels::exp_mod_normal_lcdf_device_function);
  std::printf(
      "__kernel void k(__global double* o, __global const double* z,"
      " __global const double* v) { int i = get_global_id(0);"
      " o[i] = exp_mod_normal_lcdf_value(z[i], v[i])"
      " + exp_mod_normal_lccdf_d_s(z[i], v[i]); }\n");
  return 0;
}
