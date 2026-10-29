// swst0_points.hpp -- the 240 points with small sw or st0 (probe_swst0.cpp,
// probe_ref7.cpp REF_SWST0).  Claim of stan-dev/math#3442: "240 points with sw
// or st0 from 1e-8 to 0.1".
// 12 base points (a 0.5-3, v -3..3, w 0.2-0.8, t0 0.1-0.4, sv 0-1.5,
// d = y - t0 0.05-2.05), each with
//   k<i>w<j>a: sw = SW[j], st0 = 0
//   k<i>w<j>b: sw = SW[j], st0 = 0.2 d
//   k<i>t<j>a: sw = 0,     st0 = ST0[j]   (st0 < d)
//   k<i>t<j>b: sw = 0.1,   st0 = ST0[j]
// SW = ST0 = {1e-8, 1e-6, 1e-4, 1e-2, 0.1} (ST0[4] = 0.5 d).
#ifndef SWST0_POINTS_HPP
#define SWST0_POINTS_HPP

#include <random>
#include <string>
#include <utility>
#include <vector>

inline std::vector<std::pair<std::string, std::vector<double>>>
swst0_points() {
  std::vector<std::pair<std::string, std::vector<double>>> out;
  std::mt19937 rng(20261011);
  std::uniform_real_distribution<double> U(0.0, 1.0);
  const double vals[5] = {1e-8, 1e-6, 1e-4, 1e-2, 0.1};
  for (int i = 0; i < 12; ++i) {
    const double a = 0.5 + 2.5 * U(rng);
    const double v = -3.0 + 6.0 * U(rng);
    const double w = 0.2 + 0.6 * U(rng);
    const double t0 = 0.1 + 0.3 * U(rng);
    const double sv = 1.5 * U(rng);
    const double d = 0.05 + 2.0 * U(rng);
    const double y = t0 + d;
    const std::string k = "k" + std::to_string(i);
    for (int j = 0; j < 5; ++j) {
      const std::string jj = std::to_string(j);
      out.push_back({k + "w" + jj + "a", {y, a, t0, w, v, sv, vals[j], 0.0}});
      out.push_back(
          {k + "w" + jj + "b", {y, a, t0, w, v, sv, vals[j], 0.2 * d}});
    }
    for (int j = 0; j < 5; ++j) {
      const std::string jj = std::to_string(j);
      const double st0 = j == 4 ? 0.5 * d : vals[j];
      out.push_back({k + "t" + jj + "a", {y, a, t0, w, v, sv, 0.0, st0}});
      out.push_back({k + "t" + jj + "b", {y, a, t0, w, v, sv, 0.1, st0}});
    }
  }
  return out;
}

#endif
