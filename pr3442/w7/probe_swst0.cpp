// probe_swst0.cpp -- the 7-parameter wiener_lpdf at 240 points with small sw
// or st0.  Claims of stan-dev/math#3442: d/dsw error 1.1 at sw = 1e-8 on
// develop, 3.3e-9 at the PR; end-point error about 3e-10/sw.
//
// Value and the 8 gradients (reverse mode) at the points of swst0_points.hpp,
// at the precisions in PRECS ("1e-4/1e-8" by default), with the time per
// call in microseconds (median of 3).
// Output: LPDF <id> <y a t0 w v sv sw st0> <prec> <value> <g1..g8> <us>, or
//         LPDF <id> <params> <prec> THROW <msg>.
#include <stan/math/rev.hpp>
#include "swst0_points.hpp"
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
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

void row(const std::string& id, const std::vector<double>& p, double prec) {
  std::printf("LPDF\t%s", id.c_str());
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
      var lp = stan::math::wiener_lpdf(x[0], x[1], x[2], x[3], x[4], x[5],
                                       x[6], x[7], prec);
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
    std::printf("\t%.2f\n", times[1] * 1e6);
  } catch (const std::exception& e) {
    stan::math::recover_memory();
    std::printf("\tTHROW\t%s\n", e.what());
  }
}

}  // namespace

int main() {
  std::setvbuf(stdout, nullptr, _IONBF, 0);
  const char* precs = std::getenv("PRECS");
  const std::vector<double> prec_list = split(precs ? precs : "1e-4/1e-8", '/');
  // PTS = "y/a/t0/w/v/sv/sw/st0:...": these points (ids l0, l1, ...) only
  if (const char* pts = std::getenv("PTS")) {
    const std::string s(pts);
    size_t pos = 0;
    int id = 0;
    while (pos < s.size()) {
      size_t end = s.find(':', pos);
      if (end == std::string::npos) {
        end = s.size();
      }
      const std::vector<double> p = split(s.substr(pos, end - pos), '/');
      for (double prec : prec_list) {
        row("l" + std::to_string(id), p, prec);
      }
      ++id;
      pos = end + 1;
    }
    return 0;
  }
  for (const auto& pt : swst0_points()) {
    for (double prec : prec_list) {
      row(pt.first, pt.second, prec);
    }
  }
  return 0;
}
