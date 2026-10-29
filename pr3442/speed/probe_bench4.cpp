// probe_bench4.cpp -- time the 4- and 5-parameter Wiener functions on one
// tree.  Claims of stan-dev/math#3442: run-time ratios of the 4-parameter lcdf
// and lccdf and of the 5-parameter lpdf.
//
// Points: p4 and p5 of probe_w4.cpp; default precision.
// Sets: old_p4 (wiener_lpdf(y, a, t0, w, v), the older code), lpdf_p4,
// lpdf_p5 (wiener_lpdf(y, a, t0, w, v, sv)), lcdf_p4, lccdf_p4 (4-parameter
// lcdf/lccdf), lcdf_p5, lccdf_p5 (7-argument lcdf/lccdf with sw = st0 = 0).
// Modes: dbl (value only) and var (value and gradient). Each (set, mode) is
// timed over whole passes until at least 1 s has elapsed, after one warm-up
// pass. Each timed call adds sum * 1e-300 to y (sum: the running checksum),
// which changes no bit of y but makes the loop depend on the previous
// results, so that no call can be hoisted.
// Output: BENCH <set> <mode> <ns per call> <calls> <checksum>
#include <stan/math/rev.hpp>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <functional>
#include <random>
#include <string>
#include <vector>

using stan::math::var;

namespace {

using Pt = std::vector<double>;

template <typename F>
double call(const F& f, const Pt& p0, bool use_var, double bump) {
  static thread_local Pt p;
  p = p0;
  p[0] += bump;
  try {
    if (!use_var) {
      return stan::math::value_of(f(p));
    }
    std::vector<var> x(p.begin(), p.end());
    var lp = f(x);
    std::vector<double> g;
    lp.grad(x, g);
    const double out = lp.val() + 1e-30 * g[0];
    stan::math::recover_memory();
    return out;
  } catch (const std::exception&) {
    stan::math::recover_memory();
    return NAN;
  }
}

template <typename F>
void bench(const std::string& name, const std::vector<Pt>& pts, const F& f,
           bool use_var) {
  double check = 0;
  for (const auto& p : pts) {
    call(f, p, use_var, 0.0);
  }
  long calls = 0;
  const auto t0 = std::chrono::steady_clock::now();
  double elapsed = 0;
  while (elapsed < 1.0) {
    for (const auto& p : pts) {
      const double v = call(f, p, use_var, check * 1e-300);
      if (std::isfinite(v)) {
        check += v;
      }
      ++calls;
    }
    elapsed = std::chrono::duration<double>(std::chrono::steady_clock::now()
                                            - t0)
                  .count();
  }
  std::printf("BENCH\t%s\t%s\t%.1f\t%ld\t%.10g\n", name.c_str(),
              use_var ? "var" : "dbl", 1e9 * elapsed / calls, calls, check);
}

}  // namespace

int main() {
  std::setvbuf(stdout, nullptr, _IONBF, 0);
  std::vector<Pt> p4, p5;
  std::mt19937 rng(20261010);
  std::uniform_real_distribution<double> U(0.0, 1.0);
  for (int s = 0; s < 2; ++s) {
    for (int k = 0; k < 300; ++k) {
      const double a = 0.5 + 2.5 * U(rng);
      const double v = -5.0 + 10.0 * U(rng);
      const double w = 0.2 + 0.6 * U(rng);
      const double t = 0.05 + 3.0 * U(rng);
      const double sv = s == 1 ? 1.5 * U(rng) : 0.0;
      (s == 0 ? p4 : p5).push_back({0.1 + t, a, 0.1, w, v, sv});
    }
  }
  const auto old_f = [](const auto& x) {
    return stan::math::wiener_lpdf(x[0], x[1], x[2], x[3], x[4]);
  };
  const auto lpdf_f = [](const auto& x) {
    return stan::math::wiener_lpdf(x[0], x[1], x[2], x[3], x[4], x[5]);
  };
  const auto lcdf4_f = [](const auto& x) {
    return stan::math::wiener_lcdf_unnorm(x[0], x[1], x[2], x[3], x[4]);
  };
  const auto lccdf4_f = [](const auto& x) {
    return stan::math::wiener_lccdf_unnorm(x[0], x[1], x[2], x[3], x[4]);
  };
  const auto lcdf5_f = [](const auto& x) {
    return stan::math::wiener_lcdf_unnorm(x[0], x[1], x[2], x[3], x[4], x[5],
                                          0.0, 0.0);
  };
  const auto lccdf5_f = [](const auto& x) {
    return stan::math::wiener_lccdf_unnorm(x[0], x[1], x[2], x[3], x[4], x[5],
                                           0.0, 0.0);
  };
  for (bool use_var : {false, true}) {
    bench("old_p4", p4, old_f, use_var);
    bench("lpdf_p4", p4, lpdf_f, use_var);
    bench("lpdf_p5", p5, lpdf_f, use_var);
    bench("lcdf_p4", p4, lcdf4_f, use_var);
    bench("lccdf_p4", p4, lccdf4_f, use_var);
    bench("lcdf_p5", p5, lcdf5_f, use_var);
    bench("lccdf_p5", p5, lccdf5_f, use_var);
  }
  return 0;
}
