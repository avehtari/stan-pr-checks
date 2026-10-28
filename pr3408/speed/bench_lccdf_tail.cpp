// Bundle: gamma_lccdf shape-derivative speed (Roihu job 2430495).
// Supports the gamma_lccdf speed table of #3408 Part 3.  Not run again.
// Time gamma_lccdf(y | alpha, 1) value plus alpha derivative per call.
// Usage: bench_lccdf_tail LABEL ROUNDS
// Output lines: LABEL MODE SET ns_per_call
#include <stan/math/mix.hpp>
#include <array>
#include <chrono>
#include <cstdio>
#include <string>
#include <vector>

using stan::math::fvar;
using stan::math::gamma_lccdf;
using stan::math::var;
using Pts = std::vector<std::array<double, 2>>;

namespace {

double sink = 0;

double call_var(double a, double y) {
  var av = a;
  var lp = gamma_lccdf(y, av, 1.0);
  lp.grad();
  const double g = av.adj();
  stan::math::recover_memory();
  return g;
}

double call_fvar(double a, double y) {
  return gamma_lccdf(y, fvar<double>(a, 1.0), 1.0).d_;
}

double call_fvar_var(double a, double y) {
  var av = a;
  fvar<var> lp = gamma_lccdf(y, fvar<var>(av, 1.0), 1.0);
  lp.d_.grad();
  const double h = av.adj();
  stan::math::recover_memory();
  return h;
}

template <typename F>
double run(F f, const Pts& pts, long reps) {
  double acc = 0;
  const auto t0 = std::chrono::steady_clock::now();
  for (long r = 0; r < reps; ++r) {
    for (const auto& p : pts) {
      // depends on acc, changes no bit of the input; the old code returns
      // non-finite derivatives in the deep tail
      const double g = f(p[0] + acc * 1e-300, p[1]);
      acc += std::isfinite(g) ? g : 1.0;
    }
  }
  const auto t1 = std::chrono::steady_clock::now();
  sink += acc;
  return std::chrono::duration<double, std::nano>(t1 - t0).count()
         / (static_cast<double>(reps) * pts.size());
}

template <typename F>
double time_it(F f, const Pts& pts) {
  const double ns = run(f, pts, 200);
  const long reps = std::max(200L, static_cast<long>(3e8 / (ns * pts.size())));
  return run(f, pts, reps);
}

}  // namespace

int main(int argc, char** argv) {
  const std::string label = argc > 1 ? argv[1] : "x";
  const int rounds = argc > 2 ? std::atoi(argv[2]) : 1;
  const Pts bulk_cf = {{5, 7.5},    {10, 30},    {50, 100},  {100, 150},
                       {300, 400},  {1000, 1200}, {20, 200}, {2, 40}};
  const Pts deep = {{150, 1200}, {10, 1000}, {200, 1600}, {1000, 3000}};
  const Pts log1m = {{5, 5}, {100, 100}, {20, 15}, {1000, 990}};
  const std::vector<std::pair<std::string, const Pts*>> sets
      = {{"bulk_cf", &bulk_cf}, {"deep", &deep}, {"log1m", &log1m}};
  for (int r = 0; r < rounds; ++r) {
    for (const auto& s : sets) {
      std::printf("%s var %s %.1f\n", label.c_str(), s.first.c_str(),
                  time_it(call_var, *s.second));
      std::printf("%s fvar %s %.1f\n", label.c_str(), s.first.c_str(),
                  time_it(call_fvar, *s.second));
      std::printf("%s fvar_var %s %.1f\n", label.c_str(), s.first.c_str(),
                  time_it(call_fvar_var, *s.second));
    }
  }
  std::fprintf(stderr, "sink %g\n", sink);
  return 0;
}
