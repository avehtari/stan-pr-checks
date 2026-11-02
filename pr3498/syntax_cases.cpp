// Compile-only instantiation of normal_lpdf for many argument types.
// g++ -std=c++17 -fsyntax-only -I <tree> -isystem <eigen> ... syntax_cases.cpp

#include <stan/math/mix.hpp>
#include <vector>

using stan::math::fvar;
using stan::math::var;
using stan::math::var_value;
using Eigen::Matrix;
using Eigen::RowVectorXd;
using Eigen::VectorXd;

template <typename T>
using vec = Matrix<T, -1, 1>;
template <typename T>
using rvec = Matrix<T, 1, -1>;

template <typename Y, typename M, typename S>
void inst(const Y& y, const M& mu, const S& sigma) {
  auto a = stan::math::normal_lpdf<true>(y, mu, sigma);
  auto b = stan::math::normal_lpdf<false>(y, mu, sigma);
  auto c = stan::math::normal_lpdf(y, mu, sigma);
  (void)a;
  (void)b;
  (void)c;
}

template <typename T>
void all_for(int n) {
  vec<T> v(n);
  rvec<T> rv(n);
  std::vector<T> sv(n);
  T s = 1.0;
  VectorXd dv(n);
  RowVectorXd drv(n);
  std::vector<double> dsv(n);
  double ds = 1.0;
  // all 3^3 combinations of {scalar, vector, data} for the common shapes
  inst(dv, v, s);
  inst(dv, v, v);
  inst(v, s, s);
  inst(dv, v, ds);
  inst(dv, dv, s);
  inst(dsv, sv, s);
  inst(v, v, v);
  inst(v, v, s);
  inst(v, dv, s);
  inst(s, s, s);
  inst(ds, s, ds);
  inst(s, ds, ds);
  inst(ds, ds, s);
  inst(dv, dv, v);
  inst(s, s, v);
  inst(ds, s, v);
  inst(dv, s, v);
  inst(v, s, v);
  inst(v, dv, dv);
  inst(dv, v, dv);
  inst(v, s, dv);
  inst(dv, s, dv);
  inst(s, dv, ds);
  inst(rv, rv, rv);
  inst(drv, rv, s);
  inst(rv, drv, s);
  inst(sv, sv, sv);
  inst(dsv, dsv, sv);
  inst(sv, s, dsv);
  inst(v, rv, sv);
  inst(1, s, 2);
  inst(dv, v, 2);
  inst(v * 2.0, v + v, s);
}

void varmat_cases(int n) {
  using vm = var_value<VectorXd>;
  using rvm = var_value<RowVectorXd>;
  vm m(VectorXd::Zero(n));
  rvm rm(RowVectorXd::Zero(n));
  VectorXd dv(n);
  RowVectorXd drv(n);
  var s = 1.0;
  inst(dv, m, s);
  inst(dv, m, m);
  inst(m, s, s);
  inst(m, m, m);
  inst(m, dv, s);
  inst(dv, dv, m);
  inst(s, s, m);
  inst(drv, rm, s);
  inst(rm, rm, rm);
  inst(rm, drv, 1.0);
  inst(dv, m, dv);
  inst(dv, s, m);
}

// rvalues, expressions, maps and integer containers
template <typename T>
void extra_cases(int n) {
  using stan::math::normal_lpdf;
  VectorXd dv = VectorXd::Ones(n);
  vec<T> v(n);
  T s = 1.0;
  std::vector<int> iv(n, 1);
  Eigen::Map<vec<T>> mv(v.data(), n);
  auto a = normal_lpdf(VectorXd::Ones(n), vec<T>(v), T(1.0));
  auto b = normal_lpdf<true>(dv.segment(0, 2), v.head(2), s);
  auto c = normal_lpdf(dv, mv, s);
  auto d = normal_lpdf(iv, v, s);
  auto e = normal_lpdf(dv * 2.0, mv, std::vector<T>(n, s));
  auto f = normal_lpdf<true>(v.head(2), dv.tail(2), v.tail(2));
  auto g = normal_lpdf(iv, 0, s);
  (void)a;
  (void)b;
  (void)c;
  (void)d;
  (void)e;
  (void)f;
  (void)g;
}

int main() {
  extra_cases<var>(3);
  extra_cases<fvar<double>>(3);
  extra_cases<fvar<var>>(3);
  all_for<var>(3);
  all_for<fvar<double>>(3);
  all_for<fvar<var>>(3);
  all_for<fvar<fvar<var>>>(3);
  all_for<fvar<fvar<double>>>(3);
  varmat_cases(3);
  VectorXd dv(3);
  inst(dv, dv, 1.0);
  inst(1.0, 0.0, 1.0);
  inst(std::vector<double>(3), dv, dv);
  return 0;
}
