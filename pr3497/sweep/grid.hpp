// The sweep grid of beta(a, b): a and b each take the values
// 10^(-3 + k / 40) below 2^31, and 2^31 (495 values, 245 025 pairs).
#ifndef BETA_SWEEP_GRID_HPP
#define BETA_SWEEP_GRID_HPP

#include <cmath>
#include <vector>

inline std::vector<double> sweep_axis() {
  std::vector<double> v;
  const double top = std::ldexp(1.0, 31);
  for (int k = 0;; ++k) {
    const double x = std::pow(10.0, -3.0 + k / 40.0);
    if (x >= top) {
      break;
    }
    v.push_back(x);
  }
  v.push_back(top);
  return v;
}

// Region of a point: 0 small (both below 10), 1 mixed, 2 large (both >= 10).
inline int sweep_region(double a, double b) {
  return (a >= 10) + (b >= 10);
}

#endif
