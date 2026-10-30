// GPU speed of the normal_lcdf value and gradient: origin/develop against
// the erfcx identities, on one OpenCL device.
//
// Method, from STAN_MATH_DISTRIBUTIONS.md section 4: a simple map kernel over
// global memory is bandwidth-bound and measures the memory system, not the
// function. Each work item instead runs an in-register PRNG walk over the
// range. Neighboring lanes therefore land in different branches, which is the
// divergence effect under test. A baseline kernel runs the same walk and
// accumulates the input alone, so the harness overhead is subtracted.
//
// The erfcx device function comes from the branch under test through
// stan/math/opencl/kernels/device_functions/erfcx.hpp. It is not copied here.
// The develop branches are copied, because the point is to compare against
// them; their faithfulness is checked against the CPU result at startup.

#include <CL/cl.h>

#define STAN_OPENCL
#include <stan/math/opencl/kernels/device_functions/erfcx.hpp>

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>

#define CHECK(call)                                                       \
  do {                                                                    \
    cl_int err_ = (call);                                                 \
    if (err_ != CL_SUCCESS) {                                             \
      std::fprintf(stderr, "OpenCL error %d at line %d\n", err_,          \
                   __LINE__);                                             \
      return 1;                                                           \
    }                                                                     \
  } while (0)

namespace {

const char* KERNELS = R"CLC(
#define INV_SQRT_PI       0.56418958354775628694807945156077
#define SQRT_TWO          1.4142135623730950488016887242097
#define INV_SQRT_TWO      0.70710678118654752440084436210485
#define INV_SQRT_TWO_PI_C 0.39894228040143267793994605993438
#define TWO_OVER_SQRT_PI  1.1283791670955125738961589031215
#define LOG_HALF         -0.69314718055994530941723212145818

// ── the gradient on origin/develop: 9 branches ──
double dnlcdf_develop(double s) {
  double x2 = s * s;
  if (s > 2.9) {
    double t = 1.0 / (1.0 + 0.3275911 * s);
    double t2 = t * t;
    double t4 = t2 * t2;
    double e = exp(-x2);
    return INV_SQRT_PI * e
           / (1.0 - e * (0.254829592 - 0.284496736 * t + 1.421413741 * t2
                         - 1.453152027 * t2 * t + 1.061405429 * t4));
  }
  if (s > 2.5) {
    double t = s - 2.7; double t2 = t * t; double t4 = t2 * t2;
    return 0.0003849882382 - 0.002079084702 * t + 0.005229340880 * t2
           - 0.008029540137 * t2 * t + 0.008232190507 * t4
           - 0.005692364250 * t4 * t + 0.002399496363 * pown(t, 6);
  }
  if (s > 2.1) {
    double t = s - 2.3; double t2 = t * t; double t4 = t2 * t2;
    return 0.002846135439 - 0.01310032351 * t + 0.02732189391 * t2
           - 0.03326906904 * t2 * t + 0.02482478940 * t4
           - 0.009883071924 * t4 * t - 0.0002771362254 * pown(t, 6);
  }
  if (s > 1.5) {
    double t = s - 1.85; double t2 = t * t; double t4 = t2 * t2;
    return 0.01849212058 - 0.06876280470 * t + 0.1099906382 * t2
           - 0.09274533184 * t2 * t + 0.03543327418 * t4
           + 0.005644855518 * t4 * t - 0.01111434424 * pown(t, 6);
  }
  if (s > 0.8) {
    double t = s - 1.15; double t2 = t * t; double t4 = t2 * t2;
    return 0.1585747034 - 0.3898677543 * t + 0.3515963775 * t2
           - 0.09748053605 * t2 * t - 0.04347986191 * t4
           + 0.02182506378 * t4 * t + 0.01074751427 * pown(t, 6);
  }
  if (s > 0.1) {
    double t = s - 0.45; double t2 = t * t; double t4 = t2 * t2;
    return 0.6245634904 - 0.9521866949 * t + 0.3986215682 * t2
           + 0.04700850676 * t2 * t - 0.03478651979 * t4
           - 0.01772675404 * t4 * t + 0.0006577254811 * pown(t, 6);
  }
  if (s < -29.0) {
    double iv = 1.0 / x2;
    return -2.0 * s / (1.0 + iv * (-0.5 + iv * (0.75 + iv * -1.875)));
  }
  double t = 1.0 / (1.0 - 0.3275911 * s);
  double t2 = t * t;
  double t4 = t2 * t2;
  double v = 2.0 * INV_SQRT_PI
             / (0.254829592 * t - 0.284496736 * t2 + 1.421413741 * t2 * t
                - 1.453152027 * t4 + 1.061405429 * t4 * t);
  if (s < -17.0) {
    v += 0.0001263257217272 * x2 * s + 0.0123586859488623 * x2
         - 0.0860505264736028 * s - 1.252783383752970;
  } else if (s < -7.0) {
    v += 0.000471585349920831 * x2 * s + 0.0296839305424034 * x2
         + 0.207402143352332 * s + 0.425316974683324;
  } else if (s < -3.9) {
    v += -0.0006972280656443 * x2 * s + 0.0068218494628567 * x2
         + 0.0585761964460277 * s + 0.1034397670201370;
  } else if (s < -2.1) {
    v += -0.0018742199480885 * x2 * s - 0.0097119598291202 * x2
         - 0.0170137970924080 * s - 0.0100428567412041;
  }
  return v;
}

// ── the value on origin/develop: erfc, then the Cody rational ──
double lcdf_develop(double s) {
  double x2 = s * s;
  if (s > 0.0) return log1p(-0.5 * erfc(s));
  if (s > -4.0) return log(erfc(-s)) + LOG_HALF;
  if (10.0 * log(fabs(s)) < log(DBL_MAX)) {
    double x4 = pown(s, 4);
    double x6 = pown(s, 6);
    double x8 = pown(s, 8);
    double x10 = pown(s, 10);
    double p = 0.000658749161529837803157 + 0.0160837851487422766278 / x2
               + 0.125781726111229246204 / x4 + 0.360344899949804439429 / x6
               + 0.305326634961232344035 / x8 + 0.0163153871373020978498 / x10;
    double q = -0.00233520497626869185443 - 0.0605183413124413191178 / x2
               - 0.527905102951428412248 / x4 - 1.87295284992346047209 / x6
               - 2.56852019228982242072 / x8 - 1.0 / x10;
    return LOG_HALF + log(INV_SQRT_PI + (p / q) / x2) - log(-s) - x2;
  }
  return -INFINITY;
}

// ── the erfcx identities ──
double dnlcdf_erfcx(double s) { return TWO_OVER_SQRT_PI / erfcx(-s); }

double lcdf_erfcx(double s) {
  if (s > 0.0) return log1p(-0.5 * erfc(s));
  return LOG_HALF + log(erfcx(-s)) - s * s;
}

double both_develop(double s) { return lcdf_develop(s) + dnlcdf_develop(s); }

double both_erfcx(double s) {
  double ex = erfcx(-s);
  if (s > 0.0) return log1p(-0.5 * erfc(s)) + TWO_OVER_SQRT_PI / ex;
  return LOG_HALF + log(ex) - s * s + TWO_OVER_SQRT_PI / ex;
}

// ── map kernels, for the agreement check ──
__kernel void map_grad_develop(__global const double* x, __global double* y) {
  int i = get_global_id(0); y[i] = dnlcdf_develop(x[i]);
}
__kernel void map_grad_erfcx(__global const double* x, __global double* y) {
  int i = get_global_id(0); y[i] = dnlcdf_erfcx(x[i]);
}
__kernel void map_value_develop(__global const double* x, __global double* y) {
  int i = get_global_id(0); y[i] = lcdf_develop(x[i]);
}
__kernel void map_value_erfcx(__global const double* x, __global double* y) {
  int i = get_global_id(0); y[i] = lcdf_erfcx(x[i]);
}

// ── timing: in-register PRNG walk, compute-bound, divergent by design ──
// ── PR #3406, transcribed from the normal-lcdf-opt branch ──
#define HALF_LOG_TWO_PI 0.91893853320467274178032973640562
#define LOG_HALF_C     -0.69314718055994530941723212145818

double p3406_erf_small(double x) {
  const double a[] = {3.16112374387056560, 1.13864154151050156e2,
                      3.77485237685302021e2, 3.20937758913846947e3,
                      1.85777706184603153e-1};
  const double b[] = {2.36012909523441209e1, 2.44024637934444173e2,
                      1.28261652607737228e3, 2.84423683343917062e3};
  double x2 = x * x, num = a[4] * x2, den = x2;
  for (int i = 0; i < 3; ++i) { num = (num + a[i]) * x2; den = (den + b[i]) * x2; }
  return x * (num + a[3]) / (den + b[3]);
}
double p3406_tail_corr(double r) {
  const double p[] = {0.000658749161529837803157, 0.0160837851487422766278,
                      0.125781726111229246204, 0.360344899949804439429,
                      0.305326634961232344035, 0.0163153871373020978498};
  const double q[] = {-0.00233520497626869185443, -0.0605183413124413191178,
                      -0.527905102951428412248, -1.87295284992346047209,
                      -2.56852019228982242072, -1.0};
  double num = p[5] * r + p[4], den = q[5] * r + q[4];
  for (int i = 3; i >= 0; --i) { num = num * r + p[i]; den = den * r + q[i]; }
  return (num / den) / INV_SQRT_PI;
}
double p3406_erfcx(double x) {
  if (x > 4.0) {
    double r = 1.0 / (x * x);
    return (1.0 + r * p3406_tail_corr(r)) * INV_SQRT_PI / x;
  }
  const double c[] = {5.64188496988670089e-1, 8.88314979438837594,
      6.61191906371416295e1, 2.98635138197400131e2, 8.81952221241769090e2,
      1.71204761263407058e3, 2.05107837782607147e3, 1.23033935479799725e3,
      2.15311535474403846e-8};
  const double d[] = {1.57449261107098347e1, 1.17693950891312499e2,
      5.37181101862009858e2, 1.62138957456669019e3, 3.29079923573345963e3,
      4.36261909014324716e3, 3.43936767414372164e3, 1.23033935480374942e3};
  double num = c[8] * x, den = x;
  for (int i = 0; i < 7; ++i) { num = (num + c[i]) * x; den = (den + d[i]) * x; }
  return (num + c[7]) / (den + d[7]);
}
double both_pr3406(double z) {
  if (z > 40.0) return 0.0;
  if (z <= -4.0 * SQRT_TWO) {
    double a = -z, r = 2.0 / (a * a), rc = r * p3406_tail_corr(r);
    return -(0.5 * z) * z - HALF_LOG_TWO_PI - log(a) + log1p(rc) + a / (1.0 + rc);
  }
  double x = fabs(z) * INV_SQRT_TWO;
  if (x < 0.46875) {
    double e = p3406_erf_small(x);
    return LOG_HALF_C + log1p(z < 0.0 ? -e : e)
         + TWO_OVER_SQRT_PI * INV_SQRT_TWO * exp(-x * x) / (z < 0.0 ? 1.0 - e : 1.0 + e);
  }
  double ex = p3406_erfcx(x);
  if (z < 0.0) return -(0.5 * z) * z + LOG_HALF_C + log(ex)
                    + TWO_OVER_SQRT_PI * INV_SQRT_TWO / ex;
  double dens = exp(-(0.5 * z) * z), tail = 0.5 * dens * ex;
  return log1p(-tail) + INV_SQRT_TWO_PI_C * dens / (1.0 - tail);
}
double both_pr3406_fma(double z) {
  double x = fabs(z) * INV_SQRT_TWO;
  if (z > 40.0 || z <= -4.0 * SQRT_TWO || x < 0.46875 || z < 0.0)
    return both_pr3406(z);
  double ex = p3406_erfcx(x);
  double h = z * z, err = fma(z, z, -h);
  double dens = exp(-0.5 * h) * (1.0 - 0.5 * err);
  double tail = 0.5 * dens * ex;
  return log1p(-tail) + INV_SQRT_TWO_PI_C * dens / (1.0 - tail);
}

#define WALK(EXPR)                                    \
  int gid = get_global_id(0);                         \
  uint s = (uint)gid * 747796405u + 2891336453u;      \
  double acc = 0.0;                                   \
  for (int k = 0; k < iters; ++k) {                   \
    s = s * 747796405u + 2891336453u;                 \
    double u = (double)(s >> 8) * (1.0 / 16777216.0); \
    double xi = lo + (hi - lo) * u;                   \
    acc += (EXPR);                                    \
  }                                                   \
  out[gid] = acc;

__kernel void bench_base(double lo, double hi, int iters,
                         __global double* out) { WALK(xi) }
__kernel void bench_grad_develop(double lo, double hi, int iters,
                                 __global double* out) { WALK(dnlcdf_develop(xi)) }
__kernel void bench_grad_erfcx(double lo, double hi, int iters,
                               __global double* out) { WALK(dnlcdf_erfcx(xi)) }
__kernel void bench_value_develop(double lo, double hi, int iters,
                                  __global double* out) { WALK(lcdf_develop(xi)) }
__kernel void bench_value_erfcx(double lo, double hi, int iters,
                                __global double* out) { WALK(lcdf_erfcx(xi)) }
__kernel void bench_both_develop(double lo, double hi, int iters,
                                 __global double* out) { WALK(both_develop(xi)) }
__kernel void bench_both_erfcx(double lo, double hi, int iters,
                               __global double* out) { WALK(both_erfcx(xi)) }
__kernel void bench_both_pr3406(double lo, double hi, int iters,
                                __global double* out) { WALK(both_pr3406(xi)) }
__kernel void bench_both_pr3406_fma(double lo, double hi, int iters,
                                    __global double* out) { WALK(both_pr3406_fma(xi)) }
)CLC";

struct Range {
  const char* name;
  double lo;
  double hi;
};

const Range RANGES[] = {{"far tail [-40,-10]", -40.0, -10.0},
                        {"Cody     [-10,-4] ", -10.0, -4.0},
                        {"Cody edge[-4,-2.5]", -4.0, -2.5},
                        {"interior [-2.5,-0.1]", -2.5, -0.1},
                        {"interior [0.1,2.9]", 0.1, 2.9},
                        {"A&S      [2.9,6]  ", 2.9, 6.0}};

}  // namespace

int main(int argc, char** argv) {
  const int iters = (argc > 1) ? std::atoi(argv[1]) : 2000;
  const size_t N = 1u << 20;
  const int REPEATS = 5;

  cl_uint nplat = 0;
  CHECK(clGetPlatformIDs(0, nullptr, &nplat));
  std::vector<cl_platform_id> plats(nplat);
  CHECK(clGetPlatformIDs(nplat, plats.data(), nullptr));

  cl_device_id dev = nullptr;
  for (cl_uint i = 0; i < nplat && !dev; ++i) {
    cl_uint nd = 0;
    if (clGetDeviceIDs(plats[i], CL_DEVICE_TYPE_GPU, 0, nullptr, &nd)
            == CL_SUCCESS
        && nd > 0) {
      std::vector<cl_device_id> ds(nd);
      CHECK(clGetDeviceIDs(plats[i], CL_DEVICE_TYPE_GPU, nd, ds.data(),
                           nullptr));
      dev = ds[0];
    }
  }
  if (!dev) {
    std::fprintf(stderr, "no OpenCL GPU device found\n");
    return 1;
  }
  char name[256] = {0};
  clGetDeviceInfo(dev, CL_DEVICE_NAME, sizeof name, name, nullptr);
  std::printf("device: %s\n", name);

  cl_int err = CL_SUCCESS;
  cl_context ctx = clCreateContext(nullptr, 1, &dev, nullptr, nullptr, &err);
  CHECK(err);
  const cl_queue_properties props[]
      = {CL_QUEUE_PROPERTIES, CL_QUEUE_PROFILING_ENABLE, 0};
  cl_command_queue q = clCreateCommandQueueWithProperties(ctx, dev, props,
                                                          &err);
  CHECK(err);

  std::string src;
  src += stan::math::opencl_kernels::erfcx_device_function;
  src += KERNELS;
  const char* csrc = src.c_str();
  const size_t slen = src.size();
  cl_program prog = clCreateProgramWithSource(ctx, 1, &csrc, &slen, &err);
  CHECK(err);
  if (clBuildProgram(prog, 1, &dev, "", nullptr, nullptr) != CL_SUCCESS) {
    size_t n = 0;
    clGetProgramBuildInfo(prog, dev, CL_PROGRAM_BUILD_LOG, 0, nullptr, &n);
    std::vector<char> log(n + 1, 0);
    clGetProgramBuildInfo(prog, dev, CL_PROGRAM_BUILD_LOG, n, log.data(),
                          nullptr);
    std::fprintf(stderr, "build failed:\n%s\n", log.data());
    return 1;
  }

  cl_mem bout = clCreateBuffer(ctx, CL_MEM_WRITE_ONLY, N * sizeof(double),
                               nullptr, &err);
  CHECK(err);

  auto time_kernel = [&](const char* kname, double lo, double hi) -> double {
    cl_int e = CL_SUCCESS;
    cl_kernel k = clCreateKernel(prog, kname, &e);
    if (e != CL_SUCCESS) {
      std::fprintf(stderr, "no kernel %s\n", kname);
      return -1.0;
    }
    clSetKernelArg(k, 0, sizeof(double), &lo);
    clSetKernelArg(k, 1, sizeof(double), &hi);
    clSetKernelArg(k, 2, sizeof(int), &iters);
    clSetKernelArg(k, 3, sizeof(cl_mem), &bout);
    size_t g = N;
    double best = 1e300;
    for (int rep = 0; rep < REPEATS; ++rep) {
      cl_event ev;
      clEnqueueNDRangeKernel(q, k, 1, nullptr, &g, nullptr, 0, nullptr, &ev);
      clFinish(q);
      cl_ulong t0 = 0;
      cl_ulong t1 = 0;
      clGetEventProfilingInfo(ev, CL_PROFILING_COMMAND_START, sizeof t0, &t0,
                              nullptr);
      clGetEventProfilingInfo(ev, CL_PROFILING_COMMAND_END, sizeof t1, &t1,
                              nullptr);
      const double per = static_cast<double>(t1 - t0)
                         / (static_cast<double>(N) * iters);
      if (per < best) {
        best = per;
      }
      clReleaseEvent(ev);
    }
    clReleaseKernel(k);
    return best;
  };

  std::printf("ns per call, %zu work items x %d iters, baseline subtracted\n",
              N, iters);
  struct pair {
    const char* what;
    const char* old_kernel;
    const char* new_kernel;
  };
  const pair pairs[] = {{"gradient", "bench_grad_develop", "bench_grad_erfcx"},
                        {"value", "bench_value_develop", "bench_value_erfcx"},
                        {"both", "bench_both_develop", "bench_both_erfcx"},
                        {"develop vs PR 3406", "bench_both_develop",
                         "bench_both_pr3406"},
                        {"PR 3406 vs the erfcx identity", "bench_both_pr3406",
                         "bench_both_erfcx"},
                        {"PR 3406 vs PR 3406 + fma on z*z", "bench_both_pr3406",
                         "bench_both_pr3406_fma"}};
  for (const pair& p : pairs) {
    std::printf("\n=== %s ===\n", p.what);
    std::printf("%-22s %10s %10s %8s\n", "range (scaled)", "develop", "erfcx",
                "ratio");
    for (const Range& r : RANGES) {
      const double base = time_kernel("bench_base", r.lo, r.hi);
      const double t_old = time_kernel(p.old_kernel, r.lo, r.hi) - base;
      const double t_new = time_kernel(p.new_kernel, r.lo, r.hi) - base;
      std::printf("%-22s %10.4f %10.4f %8.2f\n", r.name, t_old, t_new,
                  t_old / t_new);
    }
  }

  clReleaseMemObject(bout);
  clReleaseProgram(prog);
  clReleaseCommandQueue(q);
  clReleaseContext(ctx);
  return 0;
}
