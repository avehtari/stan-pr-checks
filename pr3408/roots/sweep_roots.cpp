// Bundle: evaluates both roots and Eigen igamma_der_a at the 1 800 points.
// Supports the nine-region table of #3408 Part 1 (build with ../build.sh).
// sweep_roots.cpp --- evaluate the develop gamma gradient roots and Eigen at
// every point of a hex point file, and print the results in hex.
//
// The host scores the output against mpmath. Nothing is judged here.
//
// Input line:  <a hex> <z hex>
// Output line: <a hex> <z hex> <grad_reg_inc_gamma> <grad_reg_lower_inc_gamma>
//              <Eigen igamma_der_a> <status>
//
// status is a three-character field, one character per function, in the same
// order: '.' finished, 'T' threw, 'N' returned a non-finite value.

#include <stan/math/prim/fun/digamma.hpp>
#include <stan/math/prim/fun/grad_reg_inc_gamma.hpp>
#include <stan/math/prim/fun/grad_reg_lower_inc_gamma.hpp>
#include <stan/math/prim/fun/tgamma.hpp>

#include <unsupported/Eigen/SpecialFunctions>

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>

namespace {

char run_one(double (*fn)(double, double), double a, double z, double* out) {
  try {
    *out = fn(a, z);
  } catch (const std::exception&) {
    *out = std::nan("");
    return 'T';
  }
  return std::isfinite(*out) ? '.' : 'N';
}

double stan_upper(double a, double z) {
  const double g = stan::math::tgamma(a);
  const double dig = stan::math::digamma(a);
  return stan::math::grad_reg_inc_gamma(a, z, g, dig);
}

double stan_lower(double a, double z) {
  return stan::math::grad_reg_lower_inc_gamma(a, z);
}

double eigen_lower(double a, double z) {
  return Eigen::numext::igamma_der_a(a, z);
}

}  // namespace

int main(int argc, char** argv) {
  std::setvbuf(stdout, nullptr, _IOLBF, 0);
  const char* path = (argc > 1) ? argv[1] : "points.txt";
  FILE* fh = std::fopen(path, "r");
  if (fh == nullptr) {
    std::printf("# ERROR: cannot open %s\n", path);
    return 1;
  }

  char line[512];
  long count = 0;
  while (std::fgets(line, sizeof(line), fh)) {
    if (line[0] == '#' || line[0] == '\n') {
      std::fputs(line, stdout);  // pass region headers through unchanged
      continue;
    }
    char a_tok[128], z_tok[128];
    if (std::sscanf(line, "%127s %127s", a_tok, z_tok) != 2) {
      continue;
    }
    const double a = std::strtod(a_tok, nullptr);
    const double z = std::strtod(z_tok, nullptr);

    double up = 0.0, lo = 0.0, ei = 0.0;
    char s0 = run_one(stan_upper, a, z, &up);
    char s1 = run_one(stan_lower, a, z, &lo);
    char s2 = run_one(eigen_lower, a, z, &ei);

    std::printf("%a %a %a %a %a %c%c%c\n", a, z, up, lo, ei, s0, s1, s2);
    ++count;
  }
  std::fclose(fh);
  std::printf("# evaluated %ld points\n", count);
  return 0;
}
