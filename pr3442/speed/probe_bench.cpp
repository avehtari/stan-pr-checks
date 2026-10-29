// probe_bench.cpp -- time the 7-parameter wiener_lpdf on one tree.  Claims of
// stan-dev/math#3442: 7-parameter lpdf run-time ratios 0.07, 0.09 and 1.05;
// "11 to 14 times faster when st0 > 0".
//
// Sets:
//   ord_both, ord_sw, ord_st0  the 300 D points of probe_shift.cpp, split by
//                              which of sw and st0 are nonzero
//   ext                        the #3329 point (3 variants) and the S scan
// Modes: dbl (double arguments, value only), var (all 8 arguments var, value
// and gradient). Precision 1e-4 (the default) for ord, 1e-8 for ext.
// Each (set, mode) is timed over whole passes until at least MIN_SECONDS
// have elapsed, after one warm-up pass.  Each timed call adds sum * 1e-300 to
// y (sum: the running checksum), which changes no bit of y but makes the loop
// depend on the previous results, so that no call can be hoisted.
// Output: BENCH <set> <mode> <ns per call> <calls> <checksum of finite values>
#include <stan/math/rev.hpp>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <random>
#include <string>
#include <vector>

using stan::math::var;

namespace {

constexpr double MIN_SECONDS = 1.0;

using Pt = std::vector<double>;

double call(const Pt& p, double prec, bool use_var, double bump) {
  const double y = p[0] + bump;
  try {
    if (!use_var) {
      return stan::math::wiener_lpdf(y, p[1], p[2], p[3], p[4], p[5], p[6],
                                     p[7], prec);
    }
    std::vector<var> x(p.begin(), p.end());
    x[0] = y;
    var lp = stan::math::wiener_lpdf(x[0], x[1], x[2], x[3], x[4], x[5], x[6],
                                     x[7], prec);
    std::vector<double> g;
    lp.grad(x, g);
    const double out = lp.val() + 1e-30 * g[0];
    stan::math::recover_memory();
    return out;
  } catch (const std::exception&) {
    stan::math::recover_memory();
    return NAN;
  }
}

void bench(const std::string& name, const std::vector<Pt>& pts, double prec,
           bool use_var) {
  double check = 0;
  for (const auto& p : pts) {
    call(p, prec, use_var, 0.0);
  }
  long calls = 0;
  const auto t0 = std::chrono::steady_clock::now();
  double elapsed = 0;
  while (elapsed < MIN_SECONDS) {
    for (const auto& p : pts) {
      const double v = call(p, prec, use_var, check * 1e-300);
      if (std::isfinite(v)) {
        check += v;
      }
      ++calls;
    }
    elapsed = std::chrono::duration<double>(std::chrono::steady_clock::now()
                                            - t0)
                  .count();
  }
  std::printf("BENCH\t%s\t%s\t%.1f\t%ld\t%.10g\n", name.c_str(),
              use_var ? "var" : "dbl", 1e9 * elapsed / calls, calls, check);
}

}  // namespace

int main() {
  std::setvbuf(stdout, nullptr, _IONBF, 0);
  std::vector<Pt> both, sw_only, st0_only, ext;
  std::mt19937 rng(20261009);
  std::uniform_real_distribution<double> U(0.0, 1.0);
  for (int k = 0; k < 300; ++k) {
    const double aa = 0.5 + 2.5 * U(rng);
    const double vv = -3.0 + 6.0 * U(rng);
    const double ww = 0.2 + 0.6 * U(rng);
    const double tt0 = 0.1 + 0.3 * U(rng);
    const double svv = 1.5 * U(rng);
    const double sww
        = (k % 3 == 2) ? 0.0 : 0.8 * U(rng) * 2.0 * std::fmin(ww, 1 - ww);
    const double stt0 = (k % 3 == 1) ? 0.0 : 0.3 * U(rng);
    const double yy = tt0 + 0.05 + 2.0 * U(rng);
    const Pt p{yy, aa, tt0, ww, vv, svv, sww, stt0};
    (k % 3 == 0 ? both : k % 3 == 1 ? sw_only : st0_only).push_back(p);
  }
  const double a = 0.1, t0 = 0.0, w = 1e-4, v = -10.0, sv = 0.1;
  ext.push_back({1e-6, a, t0, w, v, sv, 1e-8, 0.0});
  ext.push_back({1e-6, a, t0, w, v, sv, 0.0, 1e-3});
  ext.push_back({1e-6, a, t0, w, v, sv, 1e-8, 1e-3});
  for (int k = 0; k <= 15; ++k) {
    const double y = 3e-6 * std::pow(10.0, k / 15.0);
    ext.push_back({y, a, t0, w, v, sv, 1e-8, 1e-3});
  }
  // ORD_PRECS = "1e-4/1e-6": precisions for the ord sets (set names get a
  // suffix @<prec> if given); NO_EXT = 1 skips the ext set
  const char* ord_precs = std::getenv("ORD_PRECS");
  std::vector<double> precs{1e-4};
  if (ord_precs) {
    precs.clear();
    std::string s(ord_precs);
    size_t q = 0;
    while (q < s.size()) {
      size_t r = s.find('/', q);
      if (r == std::string::npos) {
        r = s.size();
      }
      precs.push_back(std::strtod(s.substr(q, r - q).c_str(), nullptr));
      q = r + 1;
    }
  }
  for (bool use_var : {false, true}) {
    for (double prec : precs) {
      char suffix[32] = "";
      if (ord_precs) {
        std::snprintf(suffix, sizeof(suffix), "@%.0e", prec);
      }
      bench(std::string("ord_both") + suffix, both, prec, use_var);
      bench(std::string("ord_sw") + suffix, sw_only, prec, use_var);
      bench(std::string("ord_st0") + suffix, st0_only, prec, use_var);
    }
    if (!std::getenv("NO_EXT")) {
      bench("ext", ext, 1e-8, use_var);
    }
  }
  return 0;
}
