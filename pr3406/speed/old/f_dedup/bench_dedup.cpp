// Speed of internal::std_normal_lcdf_value_grad (value and slope, and value
// only) for double, per range.  Compile once per tree (PR, prototype, ...)
// and alternate the binaries.
//
// Usage: bench_dedup <label> <rounds>
#include <stan/math/prim/fun/std_normal_lcdf_impl.hpp>
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <random>
#include <vector>

namespace smi = stan::math::internal;

template <bool Grad>
double time_ns(const std::vector<double>& in, int reps, double& sink) {
  double acc = sink;
  const auto t0 = std::chrono::steady_clock::now();
  for (int r = 0; r < reps; ++r) {
    const double bump = acc * 1e-300;
    for (double z : in) {
      if constexpr (Grad) {
        const auto p = smi::std_normal_lcdf_value_grad<true>(z + bump);
        acc += p.first + p.second;
      } else {
        acc += smi::std_normal_lcdf_value_grad<false>(z + bump);
      }
    }
  }
  const auto t1 = std::chrono::steady_clock::now();
  sink = acc;
  return std::chrono::duration<double, std::nano>(t1 - t0).count()
         / (static_cast<double>(reps) * in.size());
}

int main(int argc, char** argv) {
  const char* label = argc > 1 ? argv[1] : "?";
  const int rounds = argc > 2 ? std::atoi(argv[2]) : 3;
  const int block = 4096;
  const int reps = 600;
  const struct {
    const char* name;
    double lo, hi;
  } ranges[] = {{"deep[-40,-5.66]", -40.0, -5.66},
                {"neg[-5.66,-0.66]", -5.66, -0.663},
                {"central[-0.66,0.66]", -0.66, 0.66},
                {"pos[0.66,5.66]", 0.663, 5.66},
                {"upper[5.66,37]", 5.66, 37.0},
                {"mixed[-6,6]", -6.0, 6.0}};
  std::mt19937_64 rng(20261019);
  double sink = 0;
  for (const auto& r : ranges) {
    std::uniform_real_distribution<double> u(r.lo, r.hi);
    std::vector<double> in(block);
    for (auto& v : in) {
      v = u(rng);
    }
    for (int k = 0; k < rounds; ++k) {
      const double tg = time_ns<true>(in, reps, sink);
      const double tv = time_ns<false>(in, reps, sink);
      std::printf("RESULT %s %s %d grad %.3f value %.3f\n", label, r.name,
                  k + 1, tg, tv);
    }
  }
  std::printf("sink %g\n", sink);
  return 0;
}
