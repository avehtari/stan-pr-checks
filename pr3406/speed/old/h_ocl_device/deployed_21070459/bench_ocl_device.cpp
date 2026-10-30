// GPU time of the device functions Phi, std_normal_lcdf_impl and
// std_normal_lcdf_derivative, with the device-function strings of the tree
// on the include path.  STAN_MATH_DISTRIBUTIONS.md section 4: each work item
// walks through the range in registers (a pseudo-random sequence), so the
// kernel measures the function and lane divergence, not memory; a baseline
// kernel with the same walk is subtracted.  Device time from OpenCL event
// profiling.  HAVE_LCDF: the tree has the #3406 std_normal_lcdf device
// function (not on develop).
// Usage: bench_ocl_device <label> <rounds>
#include <CL/cl.h>
#define STAN_OPENCL
#ifdef HAVE_ERFCX
#include <stan/math/opencl/kernels/device_functions/erfcx.hpp>
#endif
#ifdef HAVE_LCDF
#include <stan/math/opencl/kernels/device_functions/std_normal_lcdf.hpp>
#endif
#include <stan/math/opencl/kernels/device_functions/Phi.hpp>
#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>

#define CHECK(call)                                                        \
  do {                                                                     \
    cl_int err_ = (call);                                                  \
    if (err_ != CL_SUCCESS) {                                              \
      std::fprintf(stderr, "OpenCL error %d at line %d\n", err_, __LINE__); \
      std::exit(1);                                                        \
    }                                                                      \
  } while (0)

namespace {
const char* WALK = R"CLC(
#define WALK_BEGIN                                                        \
  size_t i = get_global_id(0);                                            \
  ulong s = i * 6364136223846793005UL + 1442695040888963407UL;            \
  double acc = 0.0;                                                       \
  for (int k = 0; k < iters; ++k) {                                       \
    s = s * 6364136223846793005UL + 1442695040888963407UL;                \
    double x = lo + (hi - lo) * (double)(s >> 11) * 0x1.0p-53;
#define WALK_END                                                          \
  }                                                                       \
  out[i] = acc;
__kernel void walk_base(double lo, double hi, int iters, __global double* out) {
  WALK_BEGIN acc += x; WALK_END
}
__kernel void walk_phi(double lo, double hi, int iters, __global double* out) {
  WALK_BEGIN acc += Phi(x); WALK_END
}
#ifdef HAVE_LCDF
__kernel void walk_lcdf(double lo, double hi, int iters, __global double* out) {
  WALK_BEGIN acc += std_normal_lcdf_impl(x); WALK_END
}
__kernel void walk_dlcdf(double lo, double hi, int iters, __global double* out) {
  WALK_BEGIN acc += std_normal_lcdf_derivative(x); WALK_END
}
#endif
)CLC";

double run_ns(cl_command_queue q, cl_kernel k, double lo, double hi, int iters,
              size_t items, cl_mem out) {
  CHECK(clSetKernelArg(k, 0, sizeof(double), &lo));
  CHECK(clSetKernelArg(k, 1, sizeof(double), &hi));
  CHECK(clSetKernelArg(k, 2, sizeof(int), &iters));
  CHECK(clSetKernelArg(k, 3, sizeof(cl_mem), &out));
  double best = 1e300;
  for (int rep = 0; rep < 5; ++rep) {
    cl_event ev;
    CHECK(clEnqueueNDRangeKernel(q, k, 1, nullptr, &items, nullptr, 0,
                                 nullptr, &ev));
    CHECK(clWaitForEvents(1, &ev));
    cl_ulong t0, t1;
    CHECK(clGetEventProfilingInfo(ev, CL_PROFILING_COMMAND_START,
                                  sizeof(t0), &t0, nullptr));
    CHECK(clGetEventProfilingInfo(ev, CL_PROFILING_COMMAND_END, sizeof(t1),
                                  &t1, nullptr));
    clReleaseEvent(ev);
    best = std::min(best, static_cast<double>(t1 - t0));
  }
  return best;
}
}  // namespace

int main(int argc, char** argv) {
  const char* label = argc > 1 ? argv[1] : "?";
  const int rounds = argc > 2 ? std::atoi(argv[2]) : 3;
  cl_platform_id platform;
  CHECK(clGetPlatformIDs(1, &platform, nullptr));
  cl_device_id device;
  CHECK(clGetDeviceIDs(platform, CL_DEVICE_TYPE_GPU, 1, &device, nullptr));
  cl_int err;
  cl_context ctx = clCreateContext(nullptr, 1, &device, nullptr, nullptr, &err);
  CHECK(err);
  cl_command_queue q
      = clCreateCommandQueue(ctx, device, CL_QUEUE_PROFILING_ENABLE, &err);
  CHECK(err);

  std::string src = "#pragma OPENCL EXTENSION cl_khr_fp64 : enable\n";
#ifdef HAVE_ERFCX
  src += stan::math::opencl_kernels::erfcx_device_function;
#endif
#ifdef HAVE_LCDF
  src += stan::math::opencl_kernels::std_normal_lcdf_device_function;
  src += "\n#define HAVE_LCDF\n";
#endif
  src += stan::math::opencl_kernels::phi_device_function;
  src += WALK;
  const char* s = src.c_str();
  cl_program prog = clCreateProgramWithSource(ctx, 1, &s, nullptr, &err);
  CHECK(err);
  if (clBuildProgram(prog, 1, &device, "", nullptr, nullptr) != CL_SUCCESS) {
    size_t len;
    clGetProgramBuildInfo(prog, device, CL_PROGRAM_BUILD_LOG, 0, nullptr, &len);
    std::string log(len, ' ');
    clGetProgramBuildInfo(prog, device, CL_PROGRAM_BUILD_LOG, len, &log[0],
                          nullptr);
    std::fprintf(stderr, "build failed:\n%s\n", log.c_str());
    return 1;
  }
  const size_t items = 1 << 20;
  const int iters = 256;
  cl_mem out = clCreateBuffer(ctx, CL_MEM_WRITE_ONLY, items * sizeof(double),
                              nullptr, &err);
  CHECK(err);
  std::vector<std::string> names{"walk_phi"};
#ifdef HAVE_LCDF
  names.push_back("walk_lcdf");
  names.push_back("walk_dlcdf");
#endif
  cl_kernel base = clCreateKernel(prog, "walk_base", &err);
  CHECK(err);
  const struct {
    const char* name;
    double lo, hi;
  } ranges[] = {{"central[-0.66,0.66]", -0.66, 0.66},
                {"probit[-3,3]", -3.0, 3.0},
                {"mid[-5,5]", -5.0, 5.0},
                {"lower[-37,-5]", -37.0, -5.0},
                {"upper[5,8]", 5.0, 8.0}};
  const double evals = static_cast<double>(items) * iters;
  for (int r = 0; r < rounds; ++r) {
    for (const auto& rg : ranges) {
      const double tb = run_ns(q, base, rg.lo, rg.hi, iters, items, out);
      for (const auto& nm : names) {
        cl_kernel k = clCreateKernel(prog, nm.c_str(), &err);
        CHECK(err);
        const double t = run_ns(q, k, rg.lo, rg.hi, iters, items, out);
        clReleaseKernel(k);
        std::printf("RESULT %s %s %s %d ps %.3f\n", label, nm.c_str(), rg.name,
                    r + 1, 1000.0 * (t - tb) / evals);
      }
    }
  }
  return 0;
}
