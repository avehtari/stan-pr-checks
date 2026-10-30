// Speed of the normal_lcdf value and gradient kernels: origin/develop
// against the erfcx identities.
//
// Both variants are compiled into one binary and run alternately on the same
// node in the same job. That removes any node-to-node or session-to-session
// frequency difference from the comparison.
//
// This measures the scalar math only. It does not include autodiff overhead.
// A second, end-to-end measurement through std_normal_lcdf itself follows
// after the implementation lands, by the header-swap method of
// STAN_MATH_DISTRIBUTIONS.md section 8.
//
// Method notes, from STAN_MATH_DISTRIBUTIONS.md section 4:
//   - a dependency chain runs through the input, so the compiler cannot
//     hoist the call or vectorize the loop;
//   - an empty baseline loop with the same chain is subtracted;
//   - five alternating repetitions expose drift.
//
// Do not run this on the workstation. See AGENTS.md.

#include <stan/math/prim/fun/erfcx.hpp>

#include <chrono>
#include <cmath>
#include <cstdio>
#include <limits>

namespace {

constexpr double INV_SQRT_TWO = 0.70710678118654752440084436210485;
constexpr double INV_SQRT_PI = 0.56418958354775628694807945156077;
constexpr double TWO_OVER_SQRT_PI = 1.1283791670955125738961589031215;
constexpr double LOG_HALF = -0.69314718055994530941723212145818;

// ── The gradient on origin/develop: 9 branches, 5 Taylor series, 4 fits ──

double dnlcdf_develop(double s) {
  const double x2 = s * s;
  if (s > 2.9) {
    const double t = 1.0 / (1.0 + 0.3275911 * s);
    const double t2 = t * t;
    const double t4 = t2 * t2;
    const double exp_m_x2 = std::exp(-x2);
    return INV_SQRT_PI * exp_m_x2
           / (1.0
              - exp_m_x2
                    * (0.254829592 - 0.284496736 * t + 1.421413741 * t2
                       - 1.453152027 * t2 * t + 1.061405429 * t4));
  }
  if (s > 2.5) {
    const double t = s - 2.7;
    const double t2 = t * t;
    const double t4 = t2 * t2;
    return 0.0003849882382 - 0.002079084702 * t + 0.005229340880 * t2
           - 0.008029540137 * t2 * t + 0.008232190507 * t4
           - 0.005692364250 * t4 * t + 0.002399496363 * std::pow(t, 6);
  }
  if (s > 2.1) {
    const double t = s - 2.3;
    const double t2 = t * t;
    const double t4 = t2 * t2;
    return 0.002846135439 - 0.01310032351 * t + 0.02732189391 * t2
           - 0.03326906904 * t2 * t + 0.02482478940 * t4
           - 0.009883071924 * t4 * t - 0.0002771362254 * std::pow(t, 6);
  }
  if (s > 1.5) {
    const double t = s - 1.85;
    const double t2 = t * t;
    const double t4 = t2 * t2;
    return 0.01849212058 - 0.06876280470 * t + 0.1099906382 * t2
           - 0.09274533184 * t2 * t + 0.03543327418 * t4
           + 0.005644855518 * t4 * t - 0.01111434424 * std::pow(t, 6);
  }
  if (s > 0.8) {
    const double t = s - 1.15;
    const double t2 = t * t;
    const double t4 = t2 * t2;
    return 0.1585747034 - 0.3898677543 * t + 0.3515963775 * t2
           - 0.09748053605 * t2 * t - 0.04347986191 * t4
           + 0.02182506378 * t4 * t + 0.01074751427 * std::pow(t, 6);
  }
  if (s > 0.1) {
    const double t = s - 0.45;
    const double t2 = t * t;
    const double t4 = t2 * t2;
    return 0.6245634904 - 0.9521866949 * t + 0.3986215682 * t2
           + 0.04700850676 * t2 * t - 0.03478651979 * t4
           - 0.01772675404 * t4 * t + 0.0006577254811 * std::pow(t, 6);
  }
  if (s < -29.0) {
    const double inv_x2 = 1.0 / x2;
    return -2.0 * s
           / (1.0 + inv_x2 * (-0.5 + inv_x2 * (0.75 + inv_x2 * -1.875)));
  }
  const double t = 1.0 / (1.0 - 0.3275911 * s);
  const double t2 = t * t;
  const double t4 = t2 * t2;
  double value = 2.0 * INV_SQRT_PI
                 / (0.254829592 * t - 0.284496736 * t2 + 1.421413741 * t2 * t
                    - 1.453152027 * t4 + 1.061405429 * t4 * t);
  if (s < -17.0) {
    value += 0.0001263257217272 * x2 * s + 0.0123586859488623 * x2
             - 0.0860505264736028 * s - 1.252783383752970;
  } else if (s < -7.0) {
    value += 0.000471585349920831 * x2 * s + 0.0296839305424034 * x2
             + 0.207402143352332 * s + 0.425316974683324;
  } else if (s < -3.9) {
    value += -0.0006972280656443 * x2 * s + 0.0068218494628567 * x2
             + 0.0585761964460277 * s + 0.1034397670201370;
  } else if (s < -2.1) {
    value += -0.0018742199480885 * x2 * s - 0.0097119598291202 * x2
             - 0.0170137970924080 * s - 0.0100428567412041;
  }
  return value;
}

// ── The value on origin/develop: erfc, then the Cody rational ────────────

double lcdf_develop(double s) {
  const double x2 = s * s;
  if (s > 0.0) {
    return std::log1p(-0.5 * std::erfc(s));
  }
  if (s > -4.0) {
    return std::log(std::erfc(-s)) + LOG_HALF;
  }
  if (10.0 * std::log(std::fabs(s))
      < std::log(std::numeric_limits<double>::max())) {
    const double x4 = std::pow(s, 4);
    const double x6 = std::pow(s, 6);
    const double x8 = std::pow(s, 8);
    const double x10 = std::pow(s, 10);
    const double temp_p = 0.000658749161529837803157
                          + 0.0160837851487422766278 / x2
                          + 0.125781726111229246204 / x4
                          + 0.360344899949804439429 / x6
                          + 0.305326634961232344035 / x8
                          + 0.0163153871373020978498 / x10;
    const double temp_q = -0.00233520497626869185443
                          - 0.0605183413124413191178 / x2
                          - 0.527905102951428412248 / x4
                          - 1.87295284992346047209 / x6
                          - 2.56852019228982242072 / x8 - 1.0 / x10;
    return LOG_HALF + std::log(INV_SQRT_PI + (temp_p / temp_q) / x2)
           - std::log(-s) - x2;
  }
  return -std::numeric_limits<double>::infinity();
}

// ── The erfcx identities ─────────────────────────────────────────────────

double dnlcdf_erfcx(double s) {
  return TWO_OVER_SQRT_PI / stan::math::erfcx(-s);
}

double lcdf_erfcx(double s) {
  if (s > 0.0) {
    return std::log1p(-0.5 * std::erfc(s));
  }
  return LOG_HALF + std::log(stan::math::erfcx(-s)) - s * s;
}

// ── PR #3406: Cody rationals, internal erfcx, no libm erf/erfc ───────────
//
// Transcribed from stan/math/prim/fun/std_normal_lcdf_impl.hpp on branch
// normal-lcdf-opt. Returns value + gradient so it is comparable with
// both_develop and both_erfcx below.

namespace pr3406 {

double erf_small(double x) {
  static constexpr double a[] = {3.16112374387056560, 1.13864154151050156e2,
                                 3.77485237685302021e2, 3.20937758913846947e3,
                                 1.85777706184603153e-1};
  static constexpr double b[] = {2.36012909523441209e1, 2.44024637934444173e2,
                                 1.28261652607737228e3, 2.84423683343917062e3};
  const double x2 = x * x;
  double num = a[4] * x2;
  double den = x2;
  for (int i = 0; i < 3; ++i) {
    num = (num + a[i]) * x2;
    den = (den + b[i]) * x2;
  }
  return x * (num + a[3]) / (den + b[3]);
}

double tail_correction(double r) {
  static constexpr double p[] = {0.000658749161529837803157,
                                 0.0160837851487422766278,
                                 0.125781726111229246204,
                                 0.360344899949804439429,
                                 0.305326634961232344035,
                                 0.0163153871373020978498};
  static constexpr double q[] = {-0.00233520497626869185443,
                                 -0.0605183413124413191178,
                                 -0.527905102951428412248,
                                 -1.87295284992346047209,
                                 -2.56852019228982242072,
                                 -1.0};
  double num = p[5] * r + p[4];
  double den = q[5] * r + q[4];
  for (int i = 3; i >= 0; --i) {
    num = num * r + p[i];
    den = den * r + q[i];
  }
  return (num / den) / INV_SQRT_PI;
}

double erfcx(double x) {
  if (x > 4) {
    const double r = 1.0 / (x * x);
    return (1 + r * tail_correction(r)) * INV_SQRT_PI / x;
  }
  static constexpr double c[]
      = {5.64188496988670089e-1, 8.88314979438837594,   6.61191906371416295e1,
         2.98635138197400131e2,  8.81952221241769090e2, 1.71204761263407058e3,
         2.05107837782607147e3,  1.23033935479799725e3, 2.15311535474403846e-8};
  static constexpr double d[]
      = {1.57449261107098347e1, 1.17693950891312499e2, 5.37181101862009858e2,
         1.62138957456669019e3, 3.29079923573345963e3, 4.36261909014324716e3,
         3.43936767414372164e3, 1.23033935480374942e3};
  double num = c[8] * x;
  double den = x;
  for (int i = 0; i < 7; ++i) {
    num = (num + c[i]) * x;
    den = (den + d[i]) * x;
  }
  return (num + c[7]) / (den + d[7]);
}

constexpr double HALF_LOG_TWO_PI = 0.91893853320467274178032973640562;
constexpr double SQRT_TWO = 1.4142135623730950488016887242097;
constexpr double INV_SQRT_TWO_PI = 0.39894228040143267793994605993438;

// value + gradient, as the header returns them as a pair
double both(double z) {
  if (z > 40) {
    return 0.0;
  }
  if (z <= -4 * SQRT_TWO) {
    const double a = -z;
    const double r = 2 * (1.0 / (a * a));
    const double rc = r * tail_correction(r);
    return -(0.5 * z) * z - HALF_LOG_TWO_PI - std::log(a) + std::log1p(rc)
           + a / (1 + rc);
  }
  const double x = (z < 0 ? -z : z) * INV_SQRT_TWO;
  if (x < 0.46875) {
    const double e = erf_small(x);
    return LOG_HALF + std::log1p(z < 0 ? -e : e)
           + TWO_OVER_SQRT_PI * INV_SQRT_TWO * std::exp(-x * x)
                 / (z < 0 ? 1 - e : 1 + e);
  }
  const double ex = erfcx(x);
  if (z < 0) {
    return -(0.5 * z) * z + LOG_HALF + std::log(ex)
           + TWO_OVER_SQRT_PI * INV_SQRT_TWO / ex;
  }
  const double density = std::exp(-(0.5 * z) * z);
  const double tail = 0.5 * density * ex;
  return std::log1p(-tail) + INV_SQRT_TWO_PI * density / (1 - tail);
}

// the same, with the rounding error of z*z recovered on the z > 0 branch
double both_fma(double z) {
  if (z > 40 || z <= -4 * SQRT_TWO || (z < 0 ? -z : z) * INV_SQRT_TWO < 0.46875
      || z < 0) {
    return both(z);
  }
  const double x = z * INV_SQRT_TWO;
  const double ex = erfcx(x);
  const double h = z * z;
  const double err = std::fma(z, z, -h);
  const double density = std::exp(-0.5 * h) * (1.0 - 0.5 * err);
  const double tail = 0.5 * density * ex;
  return std::log1p(-tail) + INV_SQRT_TWO_PI * density / (1 - tail);
}

}  // namespace pr3406

double both_pr3406(double z) { return pr3406::both(z); }
double both_pr3406_fma(double z) { return pr3406::both_fma(z); }

// ── Combined value and gradient ──────────────────────────────────────────
//
// This is what an autodiff call computes. It is the number that decides the
// question, because the erfcx route calls erfcx once and uses the result
// for both, while develop evaluates two separate piecewise structures.

double both_develop(double s) { return lcdf_develop(s) + dnlcdf_develop(s); }

double both_erfcx(double s) {
  if (s > 0.0) {
    const double ex = stan::math::erfcx(-s);
    return std::log1p(-0.5 * std::erfc(s)) + TWO_OVER_SQRT_PI / ex;
  }
  const double ex = stan::math::erfcx(-s);
  return LOG_HALF + std::log(ex) - s * s + TWO_OVER_SQRT_PI / ex;
}

// ── Timing ───────────────────────────────────────────────────────────────

using clock_type = std::chrono::steady_clock;

// A dependency chain runs through `x`, so the call cannot be hoisted and the
// loop cannot be vectorized. `tiny` is 1e-300, so it changes no result bit
// but the compiler cannot prove that.
template <typename F>
double time_ns(F kernel, double lo, double hi, long n, double* checksum) {
  const double step = (hi - lo) / static_cast<double>(n);
  const double tiny = 1e-300;
  double x = lo;
  double acc = 0.0;
  const auto start = clock_type::now();
  for (long i = 0; i < n; ++i) {
    const double value = kernel(x);
    acc += value;
    x += step + value * tiny;
    if (x > hi) {
      x = lo;
    }
  }
  const auto stop = clock_type::now();
  *checksum = acc;
  const double elapsed
      = std::chrono::duration<double, std::nano>(stop - start).count();
  return elapsed / static_cast<double>(n);
}

double identity(double x) { return x; }

// Throughput: the calls are independent, as they are when normal_lcdf runs
// over a vector of observations. The input block is 4096 doubles, 32 kB, so
// it stays in L1 and the loop does not measure the memory system.
template <typename F>
double time_ns_throughput(F kernel, double lo, double hi, long n,
                          double* checksum) {
  constexpr int BLOCK = 4096;
  static double input[BLOCK];
  for (int i = 0; i < BLOCK; ++i) {
    input[i] = lo + (hi - lo) * i / (BLOCK - 1.0);
  }
  const long passes = n / BLOCK + 1;
  double acc = 0.0;
  const auto start = clock_type::now();
  for (long p = 0; p < passes; ++p) {
    // The input block does not change between passes, so the whole inner
    // loop is loop-invariant and the compiler may evaluate it once and reuse
    // the sum. `bump` is acc * 1e-300, which changes no result bit, but the
    // compiler cannot prove that, so the hoist is blocked.
    //
    // Measured 2026-09-18: without this line the normal_cdf erf branch
    // reported 1.56 ns per call, which is 3.9 cycles on a 2.5 GHz part. With
    // it the same kernel costs 39.5 ns. Any throughput number measured
    // without this line is not trustworthy.
    const double bump = acc * 1e-300;
    for (int i = 0; i < BLOCK; ++i) {
      acc += kernel(input[i] + bump);
    }
  }
  const auto stop = clock_type::now();
  *checksum = acc;
  const double elapsed
      = std::chrono::duration<double, std::nano>(stop - start).count();
  return elapsed / static_cast<double>(passes * BLOCK);
}

struct range {
  const char* name;
  double lo;
  double hi;
};

const range RANGES[] = {{"far tail  [-40,-10]", -40.0, -10.0},
                        {"Cody      [-10,-4] ", -10.0, -4.0},
                        {"Cody edge [-4,-2.5]", -4.0, -2.5},
                        {"interior  [-2.5,-0.1]", -2.5, -0.1},
                        {"interior  [0.1,2.9]", 0.1, 2.9},
                        {"A&S       [2.9,6]  ", 2.9, 6.0}};

void report(const char* what, double (*old_kernel)(double),
            double (*new_kernel)(double), long n) {
  std::printf("\n=== %s ===\n", what);
  std::printf("%-22s %8s %8s %6s   %8s %8s %6s\n", "range (scaled)", "dev lat",
              "erfcx", "ratio", "dev thr", "erfcx", "ratio");
  for (const range& r : RANGES) {
    double sink = 0.0;
    const double base = time_ns(identity, r.lo, r.hi, n, &sink);
    const double lat_old = time_ns(old_kernel, r.lo, r.hi, n, &sink) - base;
    const double lat_new = time_ns(new_kernel, r.lo, r.hi, n, &sink) - base;
    const double base_t
        = time_ns_throughput(identity, r.lo, r.hi, n, &sink);
    const double thr_old
        = time_ns_throughput(old_kernel, r.lo, r.hi, n, &sink) - base_t;
    const double thr_new
        = time_ns_throughput(new_kernel, r.lo, r.hi, n, &sink) - base_t;
    std::printf("%-22s %8.2f %8.2f %6.2f   %8.2f %8.2f %6.2f\n", r.name,
                lat_old, lat_new, lat_old / lat_new, thr_old, thr_new,
                thr_old / thr_new);
  }
}

// Guard: if the two kernels disagree badly, the benchmark is measuring the
// wrong thing. The develop gradient is only accurate to about 6e-5.
void check_agreement() {
  double worst_grad = 0.0;
  double worst_value = 0.0;
  for (int i = 0; i <= 4000; ++i) {
    const double s = -40.0 + 46.0 * i / 4000.0;
    const double g_old = dnlcdf_develop(s);
    const double g_new = dnlcdf_erfcx(s);
    const double v_old = lcdf_develop(s);
    const double v_new = lcdf_erfcx(s);
    if (std::isfinite(g_old) && g_old != 0.0) {
      worst_grad = std::fmax(worst_grad, std::fabs(g_new / g_old - 1.0));
    }
    if (std::isfinite(v_old) && v_old != 0.0) {
      worst_value = std::fmax(worst_value, std::fabs(v_new / v_old - 1.0));
    }
  }
  std::printf("agreement: gradient %.2e, value %.2e\n", worst_grad,
              worst_value);
  if (worst_grad > 1e-3 || worst_value > 1e-6) {
    std::printf("FATAL: kernels disagree, the comparison is not valid\n");
  }
}

}  // namespace

int main(int argc, char** argv) {
  const long n = (argc > 1) ? std::atol(argv[1]) : 2000000L;
  std::printf("points per range: %ld\n", n);
  check_agreement();
  report("gradient  d/ds log Phi(s * sqrt(2))", dnlcdf_develop, dnlcdf_erfcx,
         n);
  report("value     log Phi(s * sqrt(2))", lcdf_develop, lcdf_erfcx, n);
  report("both, one shared erfcx call (the autodiff case)", both_develop,
         both_erfcx, n);
  report("develop vs PR 3406", both_develop, both_pr3406, n);
  report("PR 3406 vs the erfcx identity", both_pr3406, both_erfcx, n);
  report("PR 3406 vs PR 3406 + fma on z*z", both_pr3406, both_pr3406_fma, n);
  return 0;
}
