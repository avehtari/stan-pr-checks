// probe_bb.cpp --- beta_binomial_lpmf value and reverse-mode shape partials
// at every point of a hex point file (stan-dev/math#3153).
//
// Plan: ../2026-10-06_PLAN_fix_issue3153_beta_binomial_lpmf.md, step 5.3.
// Compiled from this one source against develop and against the branch
// (build_probe.sh), so it must not name anything that exists only on the
// branch.
//
// Usage: probe_bb <points file>
//
// Input line:   <TAG> <n> <N> <alpha hex> <beta hex>
// Output line:  <TAG> <n> <N> <alpha hex> <beta hex> =>
//               <value double> <value var> <d_alpha> <d_beta> <status>
//
// status: one character per call, '.' finished, 'T' threw. Both calls use
// propto = false, so the value is the normalized log pmf.

#include <stan/math/rev.hpp>

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <string>

int main(int argc, char** argv) {
  if (argc != 2) {
    std::fprintf(stderr, "usage: probe_bb <points file>\n");
    return 2;
  }
  std::ifstream in(argv[1]);
  if (!in) {
    std::fprintf(stderr, "cannot open %s\n", argv[1]);
    return 2;
  }
  std::string line;
  while (std::getline(in, line)) {
    if (line.empty()) {
      continue;
    }
    if (line[0] == '#') {
      std::printf("%s\n", line.c_str());
      continue;
    }
    std::istringstream ss(line);
    std::string tag, tn, tN, ta, tb;
    if (!(ss >> tag >> tn >> tN >> ta >> tb)) {
      std::fprintf(stderr, "bad line: %s\n", line.c_str());
      return 2;
    }
    const int n = std::atoi(tn.c_str());
    const int N = std::atoi(tN.c_str());
    const double a = std::strtod(ta.c_str(), nullptr);
    const double b = std::strtod(tb.c_str(), nullptr);

    double v_dbl = std::nan("");
    double v_var = std::nan("");
    double d_a = std::nan("");
    double d_b = std::nan("");
    char st[3] = {'.', '.', '\0'};
    try {
      v_dbl = stan::math::beta_binomial_lpmf(n, N, a, b);
    } catch (const std::exception&) {
      st[0] = 'T';
    }
    try {
      stan::math::var av = a, bv = b;
      stan::math::var lp = stan::math::beta_binomial_lpmf(n, N, av, bv);
      lp.grad();
      v_var = lp.val();
      d_a = av.adj();
      d_b = bv.adj();
    } catch (const std::exception&) {
      st[1] = 'T';
    }
    stan::math::recover_memory();
    std::printf("%s %d %d %a %a => %a %a %a %a %s\n", tag.c_str(), n, N, a, b,
                v_dbl, v_var, d_a, d_b, st);
  }
  return 0;
}
