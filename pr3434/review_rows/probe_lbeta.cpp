// 3b check: x * d/dx of lbeta at large arguments, var and fvar<double>,
// as error / tolerance of the new mix test (fails if > 1).
#include <stan/math/mix.hpp>
#include <algorithm>
#include <cmath>
#include <cstdio>

int main() {
  using stan::math::fvar;
  using stan::math::var;
  struct T { double a, b, v, da, db; };
  const T rows[] = {
      {0x1.0000000000000p-1, 0x1.6bcc41e900000p+46, -15.545730708033618,
       -34.199701327938065, -5.0000000000000126e-15},
      {0x1.8000000000000p+0, 0x1.6bcc41e900000p+46, -48.475069190510212,
       -32.199701327938072, -1.4999999999999961e-14},
      {0x1.0000000000000p-1, 0x1.7d78400000000p+26, -8.6379754278014822,
       -20.384190769973788, -5.0000000125000003e-09},
      {0x1.9000000000000p+3, 0x1.2a05f20000000p+33, -269.08878911950677,
       -20.540655279865543, -1.24999999928125e-09}};
  auto r = [](double x, double got, double e) {
    double s = x * e;
    return std::fabs(x * got - s) / (1e-11 * std::max(1.0, std::fabs(s)));
  };
  for (const T& t : rows) {
    for (int sw = 0; sw < 2; ++sw) {
      double x1 = sw ? t.b : t.a, x2 = sw ? t.a : t.b;
      double d1 = sw ? t.db : t.da, d2 = sw ? t.da : t.db;
      var v1 = x1, v2 = x2;
      var lp = stan::math::lbeta(v1, v2);
      lp.grad();
      double rv = std::max(r(x1, v1.adj(), d1), r(x2, v2.adj(), d2));
      double val = std::fabs(lp.val() - t.v) / (1e-13 * std::fabs(t.v));
      stan::math::recover_memory();
      double f1 = stan::math::lbeta(fvar<double>(x1, 1), fvar<double>(x2, 0)).d_;
      double f2 = stan::math::lbeta(fvar<double>(x1, 0), fvar<double>(x2, 1)).d_;
      double rf = std::max(r(x1, f1, d1), r(x2, f2, d2));
      std::printf("x1 %g x2 %g  value %.2g  var %.2g  fvar %.2g\n", x1, x2,
                  val, rv, rf);
    }
  }
  return 0;
}
