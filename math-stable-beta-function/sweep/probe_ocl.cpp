// beta(a, b) on the OpenCL device and on the CPU: the reference rows of
// test/unit/math/prim/fun/beta_reference.hpp, and 1000 points of each
// region of the sweep grid (grid.hpp).  Output per point: R or G, the
// region, a, b, the CPU value and the device value in hex.  Scored by
// score_ocl.py.
#include <stan/math.hpp>
#include <test/unit/math/prim/fun/beta_reference.hpp>
#include <cstdio>
#include <vector>
#include "grid.hpp"

int main() {
  std::vector<double> as, bs;
  std::vector<char> tags;
  for (const auto& r : stan::test::beta_large_reference) {
    as.push_back(r[0]);
    bs.push_back(r[1]);
    tags.push_back('R');
  }
  const std::vector<double> axis = sweep_axis();
  std::vector<std::vector<std::pair<double, double>>> by_region(3);
  for (double a : axis) {
    for (double b : axis) {
      by_region[sweep_region(a, b)].emplace_back(a, b);
    }
  }
  for (const auto& pts : by_region) {
    const size_t stride = pts.size() / 1000;
    for (size_t i = 0; i < 1000; ++i) {
      as.push_back(pts[i * stride].first);
      bs.push_back(pts[i * stride].second);
      tags.push_back('G');
    }
  }
  const int n = as.size();
  Eigen::VectorXd a = Eigen::Map<Eigen::VectorXd>(as.data(), n);
  Eigen::VectorXd b = Eigen::Map<Eigen::VectorXd>(bs.data(), n);
  stan::math::matrix_cl<double> a_cl(a);
  stan::math::matrix_cl<double> b_cl(b);
  stan::math::matrix_cl<double> res_cl = stan::math::beta(a_cl, b_cl);
  Eigen::VectorXd dev = stan::math::from_matrix_cl<Eigen::VectorXd>(res_cl);
  for (int i = 0; i < n; ++i) {
    std::printf("%c %d %a %a %a %a\n", tags[i], sweep_region(as[i], bs[i]),
                as[i], bs[i], stan::math::beta(as[i], bs[i]), dev[i]);
  }
  return 0;
}
