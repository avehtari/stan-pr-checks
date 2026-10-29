// ir_st.cpp --- instruction count per student_t_lpdf call, for valgrind
// --tool=cachegrind (stan-dev/math#3153, decision 11). Instruction counts
// do not depend on code layout, so they separate a real cost from a
// layout effect when two timing binaries disagree.
//
// Usage: ir_st <range A|B|C> <kernel value|grad> <calls>
// Run with calls = 0 and calls = K; (I refs(K) - I refs(0)) / K is the
// count per call.

#include <stan/math/rev.hpp>

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <vector>

using stan::math::var;

int main(int argc, char** argv) {
  if (argc != 4) {
    std::fprintf(stderr, "usage: ir_st <A|B|C> <value|grad> <calls>\n");
    return 2;
  }
  const char range = argv[1][0];
  const bool grad = std::strcmp(argv[2], "grad") == 0;
  const long calls = std::atol(argv[3]);
  std::vector<double> nu(97);
  for (int i = 0; i < 97; ++i) {
    const double t = i / 97.0;
    nu[i] = range == 'A'   ? 2.0 + 18.0 * t
            : range == 'B' ? 20.0 + 42.0 * t
                           : std::pow(10.0, 2.0 + 6.0 * t);
  }
  double acc = 0;
  for (long k = 0; k < calls; ++k) {
    const double v = nu[k % 97];
    if (grad) {
      var v_nu = v;
      var lp = stan::math::student_t_lpdf(1.0, v_nu, 0.0, 1.0);
      lp.grad();
      acc += lp.val() + v_nu.adj();
      stan::math::recover_memory();
    } else {
      acc += stan::math::student_t_lpdf(1.0, v, 0.0, 1.0);
    }
  }
  std::printf("acc %.17g\n", acc);
  return 0;
}
