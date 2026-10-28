// probe_cmp_beta.cpp --- the incomplete beta gradient roots at every point
// of the shared hex point file, in one or four scalar configurations.
//
// Compiled from this one source against the pristine develop tree and
// against the branch, so it must not name anything that exists only on the
// branch. Path A is inc_beta_dda / inc_beta_ddb; path B is
// grad_reg_inc_beta. On the branch both paths run the same code.
//
// Usage: probe_cmp_beta <points file> [double|all]
//
//   double : 4 columns. Build with -DPROBE_PRIM_ONLY to avoid mix.hpp.
//   all    : 16 columns, four scalar configurations per column.
//
// Output line: <a> <b> <z> then the values in hex, then a status field with
// one character per value: '.' finished, 'T' threw.
//
//   columns, mode=double:  A_da  A_db  B_g1  B_g2
//   columns, mode=all:     for each of A_da, A_db, B_g1, B_g2 the four
//                          configurations double, var, fvar<double>,
//                          fvar<fvar<double>> (value only)

#ifdef PROBE_PRIM_ONLY
#include <stan/math/prim.hpp>
#else
#include <stan/math/mix.hpp>
#endif

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>

namespace {

// One root per call, so that a throw marks only its own column.
// which: 0 A_da, 1 A_db, 2 B_g1, 3 B_g2
template <typename T>
T root(int which, const T& a, const T& b, const T& z) {
  using stan::math::digamma;
  switch (which) {
    case 0:
      return stan::math::inc_beta_dda(a, b, z, digamma(a), digamma(a + b));
    case 1:
      return stan::math::inc_beta_ddb(a, b, z, digamma(b), digamma(a + b));
    default: {
      T g1 = 0;
      T g2 = 0;
      stan::math::grad_reg_inc_beta(g1, g2, a, b, z, digamma(a), digamma(b),
                                    digamma(a + b), stan::math::beta(a, b));
      return which == 2 ? g1 : g2;
    }
  }
}

// config: 0 double, 1 var, 2 fvar, 3 fvar<fvar>
double call_one(int which, int config, double a, double b, double z) {
  switch (config) {
    case 0:
      return root(which, a, b, z);
#ifndef PROBE_PRIM_ONLY
    case 1: {
      stan::math::nested_rev_autodiff nested;
      stan::math::var av = a, bv = b, zv = z;
      return root(which, av, bv, zv).val();
    }
    case 2: {
      stan::math::fvar<double> av = a, bv = b, zv = z;
      av.d_ = 1.0;
      return root(which, av, bv, zv).val_;
    }
    default: {
      stan::math::fvar<stan::math::fvar<double>> av = a, bv = b, zv = z;
      av.d_ = 1.0;
      return root(which, av, bv, zv).val_.val_;
    }
#else
    default:
      return std::nan("");
#endif
  }
}

}  // namespace

int main(int argc, char** argv) {
  std::setvbuf(stdout, nullptr, _IOLBF, 0);
  const char* path = (argc > 1) ? argv[1] : "points_beta.txt";
  const bool all_modes = (argc > 2) && std::strcmp(argv[2], "all") == 0;
  const int nconf = all_modes ? 4 : 1;

  FILE* fh = std::fopen(path, "r");
  if (fh == nullptr) {
    std::printf("# ERROR: cannot open %s\n", path);
    return 1;
  }
  std::printf("# probe_cmp_beta mode=%s\n", all_modes ? "all" : "double");

  char line[512];
  long count = 0;
  while (std::fgets(line, sizeof(line), fh)) {
    if (line[0] == '#' || line[0] == '\n') {
      std::fputs(line, stdout);
      continue;
    }
    char ta[128], tb[128], tz[128];
    if (std::sscanf(line, "%127s %127s %127s", ta, tb, tz) != 3) {
      continue;
    }
    const double a = std::strtod(ta, nullptr);
    const double b = std::strtod(tb, nullptr);
    const double z = std::strtod(tz, nullptr);

    double v[16];
    char s[17];
    int k = 0;
    for (int which = 0; which < 4; ++which) {
      for (int c = 0; c < nconf; ++c, ++k) {
        try {
          v[k] = call_one(which, c, a, b, z);
          s[k] = '.';
        } catch (const std::exception&) {
          v[k] = std::nan("");
          s[k] = 'T';
        }
      }
    }
    s[k] = '\0';
    std::printf("%a %a %a", a, b, z);
    for (int i = 0; i < k; ++i) {
      std::printf(" %a", v[i]);
    }
    std::printf(" %s\n", s);
    ++count;
  }
  std::fclose(fh);
  std::printf("# evaluated %ld points\n", count);
  return 0;
}
