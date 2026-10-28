// Bundle: the roots at 49 990 points (double) or 2 000 points (four modes).
// Supports the non-finite counts, the worst errors and the autodiff spread.
// probe_cmp_roots.cpp --- the two gradient roots at every point of a shared
// hex point file, in one or four scalar configurations.
//
// Compiled twice from this one source, once against the pristine develop
// tree and once against the branch, so it must not name anything that
// exists only on the branch.
//
// Usage: probe_cmp_roots <points file> [double|all]
//
//   double : columns 1 and 5 only. Safe over the whole grid.
//   all    : all eight columns. Develop does not terminate for large z in
//            the fvar configurations, because its fwd gamma_q series has no
//            iteration cap, so run this mode on a capped subset only.
//
// Output line: <a> <z> then the values in hex, then a status field with one
// character per value, '.' finished or 'T' threw.
//
//   1 grad_reg_inc_gamma        double
//   2 grad_reg_inc_gamma        var, value
//   3 grad_reg_inc_gamma        fvar<double>, value
//   4 grad_reg_inc_gamma        fvar<fvar<double>>, value
//   5 grad_reg_lower_inc_gamma  double
//   6 grad_reg_lower_inc_gamma  var, value
//   7 grad_reg_lower_inc_gamma  fvar<double>, value
//   8 grad_reg_lower_inc_gamma  fvar<fvar<double>>, value
//
// Columns 1 to 4 are dQ/da, columns 5 to 8 are dP/da = -dQ/da. Any spread
// within a group is an inconsistency between autodiff orders.

#include <stan/math/mix.hpp>

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>

namespace {

using stan::math::fvar;
using stan::math::var;

constexpr int NCOL = 8;

double call_one(int which, double a, double z) {
  switch (which) {
    case 0:
      return stan::math::grad_reg_inc_gamma(a, z, stan::math::tgamma(a),
                                            stan::math::digamma(a));
    case 1: {
      stan::math::nested_rev_autodiff nested;
      var av = a, zv = z;
      return stan::math::grad_reg_inc_gamma(av, zv, stan::math::tgamma(av),
                                            stan::math::digamma(av))
          .val();
    }
    case 2: {
      fvar<double> av = a, zv = z;
      av.d_ = 1.0;
      return stan::math::grad_reg_inc_gamma(av, zv, stan::math::tgamma(av),
                                            stan::math::digamma(av))
          .val_;
    }
    case 3: {
      fvar<fvar<double>> av = a, zv = z;
      av.d_ = 1.0;
      return stan::math::grad_reg_inc_gamma(av, zv, stan::math::tgamma(av),
                                            stan::math::digamma(av))
          .val_.val_;
    }
    case 4:
      return stan::math::grad_reg_lower_inc_gamma(a, z);
    case 5: {
      stan::math::nested_rev_autodiff nested;
      var av = a, zv = z;
      return stan::math::grad_reg_lower_inc_gamma(av, zv).val();
    }
    case 6: {
      fvar<double> av = a, zv = z;
      av.d_ = 1.0;
      return stan::math::grad_reg_lower_inc_gamma(av, zv).val_;
    }
    default: {
      fvar<fvar<double>> av = a, zv = z;
      av.d_ = 1.0;
      return stan::math::grad_reg_lower_inc_gamma(av, zv).val_.val_;
    }
  }
}

}  // namespace

int main(int argc, char** argv) {
  std::setvbuf(stdout, nullptr, _IOLBF, 0);

  const char* path = (argc > 1) ? argv[1] : "points_cmp.txt";
  const bool all_modes = (argc > 2) && std::strcmp(argv[2], "all") == 0;
  const int cols[8] = {0, 1, 2, 3, 4, 5, 6, 7};
  const int ncols = all_modes ? 8 : 2;
  const int double_only[2] = {0, 4};

  FILE* fh = std::fopen(path, "r");
  if (fh == nullptr) {
    std::printf("# ERROR: cannot open %s\n", path);
    return 1;
  }

  std::printf("# probe_cmp_roots mode=%s\n", all_modes ? "all" : "double");

  char line[512];
  long count = 0;
  while (std::fgets(line, sizeof(line), fh)) {
    if (line[0] == '#' || line[0] == '\n') {
      std::fputs(line, stdout);
      continue;
    }
    char a_tok[128], z_tok[128];
    if (std::sscanf(line, "%127s %127s", a_tok, z_tok) != 2) {
      continue;
    }
    const double a = std::strtod(a_tok, nullptr);
    const double z = std::strtod(z_tok, nullptr);

    double v[NCOL];
    char s[NCOL + 1];
    for (int i = 0; i < ncols; ++i) {
      const int w = all_modes ? cols[i] : double_only[i];
      try {
        v[i] = call_one(w, a, z);
        s[i] = '.';
      } catch (const std::exception&) {
        v[i] = std::nan("");
        s[i] = 'T';
      }
    }
    s[ncols] = '\0';

    std::printf("%a %a", a, z);
    for (int i = 0; i < ncols; ++i) {
      std::printf(" %a", v[i]);
    }
    std::printf(" %s\n", s);
    ++count;
  }
  std::fclose(fh);
  std::printf("# evaluated %ld points\n", count);
  return 0;
}
