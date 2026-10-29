// probe_cdfacc.cpp -- the 7-parameter lcdf and lccdf at the 300 D points.
// Claims of stan-dev/math#3442: the 7-parameter lcdf and lccdf rows, the
// 3-dimensional limit (4.3e-5) and the lcdf and lccdf run-time ratios.
//
// wiener_lcdf_unnorm and wiener_lccdf_unnorm (value and the 8 gradients,
// reverse mode) at the D points of probe_shift.cpp, at the precisions in
// PRECS ("1e-4/1e-6", default 1e-4), with the time per call in
// microseconds (median of 3).  PTS (see main) replaces the D points.
// Output: <LCDF|LCCDF> <id> <y a t0 w v sv sw st0> <prec> <value> <g1..g8>
//         <us>, or ... THROW <msg>.
#include <stan/math/rev.hpp>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <random>
#include <string>
#include <vector>

using stan::math::var;

namespace {

std::vector<double> split(const std::string& s, char sep) {
  std::vector<double> out;
  size_t q = 0;
  while (q <= s.size()) {
    size_t r = s.find(sep, q);
    if (r == std::string::npos) {
      r = s.size();
    }
    if (r > q) {
      out.push_back(std::strtod(s.substr(q, r - q).c_str(), nullptr));
    }
    q = r + 1;
  }
  return out;
}

template <bool Ccdf>
void row(const std::string& id, const std::vector<double>& p, double prec) {
  std::printf("%s\t%s", Ccdf ? "LCCDF" : "LCDF", id.c_str());
  for (double q : p) {
    std::printf("\t%.17g", q);
  }
  std::printf("\t%.3g", prec);
  try {
    std::vector<double> times, g;
    double val = 0;
    for (int r = 0; r < 3; ++r) {
      const auto t0 = std::chrono::steady_clock::now();
      std::vector<var> x(p.begin(), p.end());
      var lp = Ccdf ? stan::math::wiener_lccdf_unnorm(x[0], x[1], x[2], x[3],
                                                      x[4], x[5], x[6], x[7],
                                                      prec)
                    : stan::math::wiener_lcdf_unnorm(x[0], x[1], x[2], x[3],
                                                     x[4], x[5], x[6], x[7],
                                                     prec);
      lp.grad(x, g);
      val = lp.val();
      stan::math::recover_memory();
      times.push_back(
          std::chrono::duration<double>(std::chrono::steady_clock::now() - t0)
              .count());
    }
    std::sort(times.begin(), times.end());
    std::printf("\t%.17g", val);
    for (double q : g) {
      std::printf("\t%.17g", q);
    }
    std::printf("\t%.2f\n", 1e6 * times[1]);
  } catch (const std::exception& e) {
    stan::math::recover_memory();
    std::printf("\tTHROW\t%s\n", e.what());
  }
}

}  // namespace

int main() {
  std::setvbuf(stdout, nullptr, _IONBF, 0);
  const char* precs = std::getenv("PRECS");
  const std::vector<double> prec_list = split(precs ? precs : "1e-4", '/');
  // PTS = "y/a/t0/w/v/sv/sw/st0:...": these points (ids l0, l1, ...) only
  if (const char* pts = std::getenv("PTS")) {
    std::string s(pts);
    size_t pos = 0;
    int id = 0;
    while (pos < s.size()) {
      size_t end = s.find(':', pos);
      if (end == std::string::npos) {
        end = s.size();
      }
      const std::vector<double> p = split(s.substr(pos, end - pos), '/');
      for (double prec : prec_list) {
        row<false>("l" + std::to_string(id), p, prec);
        row<true>("l" + std::to_string(id), p, prec);
      }
      ++id;
      pos = end + 1;
    }
    return 0;
  }
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
    const std::vector<double> p{yy, aa, tt0, ww, vv, svv, sww, stt0};
    const std::string id = "d" + std::to_string(k);
    for (double prec : prec_list) {
      row<false>(id, p, prec);
      row<true>(id, p, prec);
    }
  }
  return 0;
}
