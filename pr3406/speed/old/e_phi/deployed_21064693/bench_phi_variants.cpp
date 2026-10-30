// Speed of the Phi(double) variants of phi_variants.hpp, per range, with the
// variants interleaved in each round.  Usage: bench_phi_variants <rounds>
#include "phi_variants.hpp"
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <random>
#include <vector>

template <double (*F)(double)>
double time_ns(const std::vector<double>& in, int reps, double& sink) {
  double acc = sink;
  const auto t0 = std::chrono::steady_clock::now();
  for (int r = 0; r < reps; ++r) {
    const double bump = acc * 1e-300;
    for (double x : in) {
      acc += F(x + bump);
    }
  }
  const auto t1 = std::chrono::steady_clock::now();
  sink = acc;
  return std::chrono::duration<double, std::nano>(t1 - t0).count()
         / (static_cast<double>(reps) * in.size());
}

int main(int argc, char** argv) {
  namespace pv = phi_variants;
  const int rounds = argc > 1 ? std::atoi(argv[1]) : 3;
  const struct {
    const char* name;
    double lo, hi;
  } ranges[] = {{"central[-0.66,0.66]", -0.66, 0.66},
                {"probit[-3,3]", -3.0, 3.0},
                {"mid[-5,5]", -5.0, 5.0},
                {"lower[-37,-5]", -37.0, -5.0},
                {"upper[5,8]", 5.0, 8.0}};
  std::mt19937_64 rng(20261019);
  double sink = 0;
  for (const auto& r : ranges) {
    std::uniform_real_distribution<double> u(r.lo, r.hi);
    std::vector<double> in(4096);
    for (auto& v : in) {
      v = u(rng);
    }
    for (int k = 0; k < rounds; ++k) {
      const double a = time_ns<pv::dev>(in, 600, sink);
      const double b = time_ns<pv::pr>(in, 600, sink);
      const double c = time_ns<pv::d1>(in, 600, sink);
      const double d = time_ns<pv::d2>(in, 600, sink);
      std::printf("RESULT %s %d dev %.3f pr %.3f d1 %.3f d2 %.3f\n", r.name,
                  k + 1, a, b, c, d);
    }
  }
  std::printf("sink %g\n", sink);
  return 0;
}
