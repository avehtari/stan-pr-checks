// Bundle: gamma_lccdf d/dalpha in var, fvar, fvar<var>, fvar<fvar> at 17
// points.  Supports 5.6e-12 / 1.9e-10 / 7e-9 of #3408 Part 3.
// gamma_lccdf(y | alpha, 1): value and d/dalpha in every autodiff mode.
// Prints hex so that the Python scorer can read the values back exactly.
#include <stan/math/mix.hpp>
#include <cstdio>

using stan::math::fvar;
using stan::math::gamma_lccdf;
using stan::math::var;

int main() {
  const double pts[][2] = {
      {5, 6.5},     {5, 7.5},      {10, 12},      {150, 1200}, {300, 1500},
      {10, 1000},   {200, 1600},   {1000, 3000},  {1e4, 10100}, {1e6, 1.1e6},
      {700, 2100},  {20, 1000},    {50, 150},     {100, 120},  {1e5, 101000},
      {200, 200},   {1000, 1000},
  };
  std::printf("# alpha y value d_var d_fvar d_fvar_var d_fvar_fvar d2_fvar_var\n");
  for (const auto& p : pts) {
    const double a = p[0];
    const double y = p[1];
    var av = a;
    var lv = gamma_lccdf(y, av, 1.0);
    lv.grad();
    const double d_var = av.adj();
    const double val = lv.val();
    stan::math::recover_memory();

    fvar<double> af(a, 1.0);
    const double d_fvar = gamma_lccdf(y, af, 1.0).d_;

    var a_inner = a;
    fvar<var> afv(a_inner, var(1.0));
    fvar<var> lfv = gamma_lccdf(y, afv, 1.0);
    const double d_fvar_var = lfv.d_.val();
    lfv.d_.grad();
    const double d2_fvar_var = a_inner.adj();
    stan::math::recover_memory();

    fvar<fvar<double>> aff(fvar<double>(a, 1.0), fvar<double>(1.0, 0.0));
    const fvar<fvar<double>> lff = gamma_lccdf(y, aff, 1.0);
    const double d_fvar_fvar = lff.d_.val_;
    std::fprintf(stderr, "%a %a %a\n", a, y, lff.d_.d_);

    std::printf("%a %a %a %a %a %a %a %a\n", a, y, val, d_var, d_fvar,
                d_fvar_var, d_fvar_fvar, d2_fvar_var);
  }
  return 0;
}
