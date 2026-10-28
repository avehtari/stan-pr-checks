// Bundle: cdf speed, var value plus all gradients (Roihu job 2430495).
// Supports the cdf speed table of #3408 Part 3.  Not run again here.
// Time value plus all gradients (var) of incomplete-gamma cdf functions.
// Usage: bench_cdfs LABEL
// Output lines: LABEL FUNCTION SET ns_per_call
// Each point is (a, z) of the incomplete gamma function, mapped to the
// arguments of the function as in survey_cases.py.
#include <stan/math/rev.hpp>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <functional>
#include <string>
#include <vector>

using stan::math::var;
namespace sm = stan::math;
using Pts = std::vector<std::array<double, 2>>;

namespace {

double sink = 0;

struct Fn {
  const char* name;
  std::function<double(double, double)> call;  // (a, z) -> sum of gradients
};

double grads3(var lp, var& y, var& p1, var& p2) {
  lp.grad();
  const double g = y.adj() + p1.adj() + p2.adj();
  sm::recover_memory();
  return g;
}

double grads2(var lp, var& y, var& p1) {
  lp.grad();
  const double g = y.adj() + p1.adj();
  sm::recover_memory();
  return g;
}

const std::vector<Fn>& fns() {
  static const std::vector<Fn> f = {
      {"gamma_cdf", [](double a, double z) {
         var y = z / 2, al = a, be = 2.0;
         return grads3(sm::gamma_cdf(y, al, be), y, al, be);
       }},
      {"gamma_lcdf", [](double a, double z) {
         var y = z / 2, al = a, be = 2.0;
         return grads3(sm::gamma_lcdf(y, al, be), y, al, be);
       }},
      {"gamma_lccdf", [](double a, double z) {
         var y = z / 2, al = a, be = 2.0;
         return grads3(sm::gamma_lccdf(y, al, be), y, al, be);
       }},
      {"chi_square_cdf", [](double a, double z) {
         var y = 2 * z, nu = 2 * a;
         return grads2(sm::chi_square_cdf(y, nu), y, nu);
       }},
      {"chi_square_lccdf", [](double a, double z) {
         var y = 2 * z, nu = 2 * a;
         return grads2(sm::chi_square_lccdf(y, nu), y, nu);
       }},
      {"inv_gamma_lcdf", [](double a, double z) {
         var y = 2 / z, al = a, be = 2.0;
         return grads3(sm::inv_gamma_lcdf(y, al, be), y, al, be);
       }},
      {"inv_gamma_cdf", [](double a, double z) {
         var y = 2 / z, al = a, be = 2.0;
         return grads3(sm::inv_gamma_cdf(y, al, be), y, al, be);
       }},
      {"scaled_inv_chi_square_lcdf", [](double a, double z) {
         var y = a / z, nu = 2 * a, s = 1.0;
         return grads3(sm::scaled_inv_chi_square_lcdf(y, nu, s), y, nu, s);
       }},
      {"poisson_lcdf", [](double a, double z) {
         var lam = z;
         var lp = sm::poisson_lcdf(static_cast<int>(a) - 1, lam);
         lp.grad();
         const double g = lam.adj();
         sm::recover_memory();
         return g;
       }},
  };
  return f;
}

double run(const Fn& f, const Pts& pts, long reps) {
  double acc = 0;
  const auto t0 = std::chrono::steady_clock::now();
  for (long r = 0; r < reps; ++r) {
    for (const auto& p : pts) {
      // depends on acc, changes no bit of the input
      const double g = f.call(p[0], p[1] + acc * 1e-300);
      acc += std::isfinite(g) ? g : 1.0;
    }
  }
  const auto t1 = std::chrono::steady_clock::now();
  sink += acc;
  return std::chrono::duration<double, std::nano>(t1 - t0).count()
         / (static_cast<double>(reps) * pts.size());
}

}  // namespace

int main(int argc, char** argv) {
  const std::string label = argc > 1 ? argv[1] : "x";
  // a integer so that the poisson count is exact
  const Pts small = {{5, 1}, {5, 3}, {5, 4.5}, {5, 7}, {5, 10}, {5, 15}};
  const Pts mid = {{50, 30}, {50, 40}, {50, 48}, {50, 55}, {50, 65}, {50, 80}};
  for (const auto& f : fns()) {
    for (const auto& s : {std::make_pair("a5", &small), std::make_pair("a50", &mid)}) {
      const double ns0 = run(f, *s.second, 200);
      const long reps = std::max(200L, static_cast<long>(2e8 / (ns0 * s.second->size())));
      std::printf("%s %s %s %.1f\n", label.c_str(), f.name, s.first,
                  run(f, *s.second, reps));
    }
  }
  std::fprintf(stderr, "sink %g\n", sink);
  return 0;
}
