// probe_pr3442.cpp -- the candidate rows of the new tests of stan-dev/math#3442
// on one tree: value, the 8 gradients (reverse mode) and the time per call in
// microseconds.  Claim: the new test rows fail at the merge base.
// Output: <fn> <id> <y a t0 w v sv sw st0> <prec> <value> <g1..g8> <us>
//         or ... THROW <msg>.  prec 0 means the default precision.
// Sets:
//   LPDF e*  extreme-grid points of probe_shift.cpp (issue regions; d/dsw
//            with st0 = 0; d/dst0 with sw and st0)
//   LCDF/LCCDF d65, d178  D points with 2 integration dimensions
//            (coordinates of probe_cdfacc.cpp)
//   LCCDF k*  drift variability only, log ccdf -39 and near -700
//   neginf_*  the two calls below the smallest double (value -inf)
#include <stan/math/rev.hpp>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <string>
#include <utility>
#include <vector>

using stan::math::var;

namespace {

enum class Fn { lpdf, lcdf, lccdf };

void row(Fn fn, const std::string& id, const std::vector<double>& p,
         double prec) {
  const char* name
      = fn == Fn::lpdf ? "LPDF" : fn == Fn::lcdf ? "LCDF" : "LCCDF";
  std::printf("%s\t%s", name, id.c_str());
  for (double q : p) {
    std::printf("\t%.17g", q);
  }
  std::printf("\t%.3g", prec);
  try {
    const auto start = std::chrono::steady_clock::now();
    std::vector<var> x(p.begin(), p.end());
    var lp;
    if (fn == Fn::lpdf) {
      lp = prec > 0 ? stan::math::wiener_lpdf(x[0], x[1], x[2], x[3], x[4],
                                              x[5], x[6], x[7], prec)
                    : stan::math::wiener_lpdf(x[0], x[1], x[2], x[3], x[4],
                                              x[5], x[6], x[7]);
    } else if (fn == Fn::lcdf) {
      lp = prec > 0 ? stan::math::wiener_lcdf_unnorm(x[0], x[1], x[2], x[3],
                                                     x[4], x[5], x[6], x[7],
                                                     prec)
                    : stan::math::wiener_lcdf_unnorm(x[0], x[1], x[2], x[3],
                                                     x[4], x[5], x[6], x[7]);
    } else {
      lp = prec > 0 ? stan::math::wiener_lccdf_unnorm(x[0], x[1], x[2], x[3],
                                                      x[4], x[5], x[6], x[7],
                                                      prec)
                    : stan::math::wiener_lccdf_unnorm(x[0], x[1], x[2], x[3],
                                                      x[4], x[5], x[6], x[7]);
    }
    std::vector<double> g;
    lp.grad(x, g);
    const double us = 1e6
                      * std::chrono::duration<double>(
                            std::chrono::steady_clock::now() - start)
                            .count();
    std::printf("\t%.17g", lp.val());
    for (double q : g) {
      std::printf("\t%.17g", q);
    }
    std::printf("\t%.1f\n", us);
  } catch (const std::exception& e) {
    std::printf("\tTHROW\t%s\n", e.what());
  }
  stan::math::recover_memory();
}

}  // namespace

int main() {
  std::setvbuf(stdout, nullptr, _IONBF, 0);
  // y, a, t0, w, v, sv, sw, st0
  const std::vector<std::pair<std::string, std::vector<double>>> lpdf_pts
      = {{"e65", {10.2, 0.1, 0.2, 1e-4, -30, 0.1, 5e-5, 1e-3}},
         {"e81", {10.2, 8, 0.2, 1e-4, -30, 0.1, 5e-5, 1e-3}},
         {"e100", {50.2, 0.1, 0.2, 1e-4, -10, 0.1, 5e-5, 0}},
         {"e113", {50.2, 8, 0.2, 1e-4, -30, 0.1, 5e-5, 1e-3}},
         {"e84", {10.2, 8, 0.2, 1e-4, -10, 0.1, 5e-5, 0}},
         {"e88", {10.2, 8, 0.2, 1e-4, 10, 0.1, 5e-5, 0}},
         {"e87", {10.2, 8, 0.2, 0.5, -10, 0.1, 0.25, 1e-3}},
         {"e91", {10.2, 8, 0.2, 0.5, 10, 0.1, 0.25, 1e-3}}};
  for (const auto& pt : lpdf_pts) {
    for (double prec : {0.0, 1e-6, 1e-8}) {
      row(Fn::lpdf, pt.first, pt.second, prec);
    }
  }
  const std::vector<std::pair<std::string, std::vector<double>>> cdf_pts
      = {{"d65",
          {0.20753147706690653, 0.96248012425856511, 0.11408699836935123,
           0.20587307426609611, 2.6537214716508419, 0.80835488843424108, 0,
           0.21064846346879404}},
         {"d178",
          {0.54520091427708517, 1.7495584716477295, 0.30505332833140775,
           0.29376778741039961, 2.2034724592953872, 0.56060550137195397,
           0.14611719349304031, 0}}};
  for (const auto& pt : cdf_pts) {
    for (double prec : {0.0, 1e-6}) {
      row(Fn::lcdf, pt.first, pt.second, prec);
      row(Fn::lccdf, pt.first, pt.second, prec);
    }
  }
  for (double y : {3.1, 9.0, 22.1, 45.0, 50.0, 51.0, 52.0, 53.0, 54.0}) {
    row(Fn::lccdf, "k" + std::to_string(y),
        {y, 2.0, 0.1, 0.5, 5.0, 0.01, 0, 0}, 0.0);
  }
  row(Fn::lccdf, "neginf_ccdf", {60.0, 2.0, 0.1, 0.5, 5.0, 0.01, 0, 0}, 0.0);
  row(Fn::lcdf, "neginf_cdf", {1e-6, 0.1, 0.0, 1e-4, -10.0, 0.1, 0, 0}, 0.0);
  return 0;
}
