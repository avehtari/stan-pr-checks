// Bundle: times each root call; on develop the fvar call at a = 0.168,
// z = 4 748 does not return.  Supports the forward-mode gamma_q claim.
// probe_diag_slow.cpp --- which (point, scalar configuration) is slow?
//
// The develop build of probe_cmp_roots wrote 6 rows in 55 minutes and then
// hit the walltime, with 1.4 GB resident. This probe finds the call
// responsible. It writes an unbuffered marker to stderr *before* each call
// and the elapsed time *after* it, so the last marker in the log names the
// call that hung, even if the process is killed.
//
// Usage: probe_diag_slow <points file> [n_points] [report_ms]

#include <stan/math/mix.hpp>

#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>

namespace {

using stan::math::fvar;
using stan::math::var;

const char* kNames[8] = {"upper_double", "upper_var", "upper_fd", "upper_ffd",
                         "lower_double", "lower_var", "lower_fd", "lower_ffd"};

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
  std::setvbuf(stderr, nullptr, _IONBF, 0);

  const char* path = (argc > 1) ? argv[1] : "points_cmp.txt";
  const long nmax = (argc > 2) ? std::atol(argv[2]) : 200;
  const double report_ms = (argc > 3) ? std::atof(argv[3]) : 50.0;

  FILE* fh = std::fopen(path, "r");
  if (fh == nullptr) {
    std::fprintf(stderr, "cannot open %s\n", path);
    return 1;
  }

  char line[512];
  long n = 0;
  while (std::fgets(line, sizeof(line), fh) && n < nmax) {
    if (line[0] == '#' || line[0] == '\n') {
      continue;
    }
    char a_tok[128], z_tok[128];
    if (std::sscanf(line, "%127s %127s", a_tok, z_tok) != 2) {
      continue;
    }
    const double a = std::strtod(a_tok, nullptr);
    const double z = std::strtod(z_tok, nullptr);
    ++n;

    for (int w = 0; w < 8; ++w) {
      // marker first: if the call never returns, this line is the evidence
      std::fprintf(stderr, "START n=%ld %s a=%.17g z=%.17g\n", n, kNames[w], a,
                   z);
      const auto t0 = std::chrono::steady_clock::now();
      double v = std::nan("");
      const char* status = "ok";
      try {
        v = call_one(w, a, z);
      } catch (const std::exception& e) {
        status = "throw";
      }
      const double ms = std::chrono::duration<double, std::milli>(
                            std::chrono::steady_clock::now() - t0)
                            .count();
      std::fprintf(stderr, "DONE  n=%ld %s %.1f ms %s\n", n, kNames[w], ms,
                   status);
      if (ms >= report_ms) {
        std::printf("SLOW %8.1f ms  n=%-5ld %-13s a=%-18.10g z=%-18.10g "
                    "value=%.10g %s\n",
                    ms, n, kNames[w], a, z, v, status);
      }
    }
  }
  std::fclose(fh);
  std::printf("# examined %ld points\n", n);
  return 0;
}
