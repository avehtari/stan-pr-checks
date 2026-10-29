// probe_w5grad.cpp -- value and gradients of the 5-parameter wiener_lpdf at
// points with a small density, and a 20 000-point scan.  Claim of
// stan-dev/math#3442: the three faults of the 5-parameter lpdf.
//
// Rows:
//   PUB <id> <t a w v sv> <value> <d/dy d/da d/dt0 d/dw d/dv d/dsv>
//       public wiener_lpdf(y = t, a, t0 = 0, w, v, sv), reverse mode,
//       precision 1e-4 (default) and 1e-8 (PUB8)
//   REL <id> ... the internal functions with a series tolerance relative to
//       the factor in front of the series, log(1e-15) + min(0, log factor)
#include <stan/math/rev.hpp>
#include <cmath>
#include <cstdio>
#include <random>
#include <string>
#include <vector>

using stan::math::var;

namespace {

void pub(const char* tag, const std::string& id, const std::vector<double>& p,
         double prec) {
  std::vector<var> x{p[0], p[1], 0.0, p[2], p[3], p[4]};
  var lp = stan::math::wiener_lpdf(x[0], x[1], x[2], x[3], x[4], x[5], prec);
  std::vector<double> g;
  lp.grad(x, g);
  std::printf("%s\t%s", tag, id.c_str());
  for (double q : p) {
    std::printf("\t%.17g", q);
  }
  std::printf("\t%.17g", lp.val());
  for (double q : g) {
    std::printf("\t%.17g", q);
  }
  std::printf("\n");
  stan::math::recover_memory();
}

void rel(const std::string& id, const std::vector<double>& p) {
  using stan::math::internal::GradientCalc;
  namespace in = stan::math::internal;
  const double t = p[0], a = p[1], w = p[2], v = p[3], sv = p[4];
  const double let = in::wiener5_compute_log_error_term(t, a, v, w, sv);
  const double err = std::log(1e-15) + std::fmin(let, 0.0);
  const double gt = in::wiener5_grad_t<GradientCalc::OFF>(t, a, v, w, sv, err);
  std::printf("REL\t%s", id.c_str());
  for (double q : p) {
    std::printf("\t%.17g", q);
  }
  std::printf("\t%.17g\t%.17g\t%.17g\t%.17g\t%.17g\t%.17g\t%.17g\n",
              in::wiener5_density<GradientCalc::OFF>(t, a, v, w, sv, err), gt,
              in::wiener5_grad_a<GradientCalc::OFF>(t, a, v, w, sv, err), -gt,
              in::wiener5_grad_w<GradientCalc::OFF>(t, a, v, w, sv, err),
              in::wiener5_grad_v<GradientCalc::OFF>(t, a, v, w, sv, err),
              in::wiener5_grad_sv<GradientCalc::OFF>(t, a, v, w, sv, err));
}

}  // namespace

int main() {
  std::setvbuf(stdout, nullptr, _IONBF, 0);
  // t, a, w, v, sv
  const std::vector<std::pair<std::string, std::vector<double>>> pts{
      {"s1", {2, 2, 0.3, 5, 0}},       {"s2", {2, 2, 0.3, -5, 0}},
      {"s4", {1, 2, 0.2, -6, 0}},      {"s5", {10, 8, 1e-4, -10, 0}},
      {"s6", {3, 3, 0.1, -4, 0}},      {"s7", {0.8, 1.5, 0.3, -5, 0}},
      {"v4", {1, 2, 0.2, -6, 0.5}},    {"v6", {3, 3, 0.1, -4, 1.0}},
      {"vw", {1.40338, 2.99164, 0.217078, -4.91775, 0.0470162}},
      {"i3322", {5.99, 10, 0.1, -3, 0.2}}};
  for (const auto& pt : pts) {
    pub("PUB", pt.first, pt.second, 1e-4);
    pub("PUB8", pt.first, pt.second, 1e-8);
    rel(pt.first, pt.second);
    const auto& p = pt.second;
    if (p[4] == 0) {
      // 4-parameter lcdf (same estimate_with_err_check path as the density)
      std::vector<var> x{p[0], p[1], 0.0, p[2], p[3]};
      var lc = stan::math::wiener_lcdf_unnorm(x[0], x[1], x[2], x[3], x[4]);
      std::vector<double> g;
      lc.grad(x, g);
      std::printf("LCDF\t%s", pt.first.c_str());
      for (double q : p) {
        std::printf("\t%.17g", q);
      }
      std::printf("\t%.17g", lc.val());
      for (double q : g) {
        std::printf("\t%.17g", q);
      }
      std::printf("\n");
      stan::math::recover_memory();
    }
  }
  // scan of the plausible region: public value and
  // gradients against the internal functions with relative tolerance
  using stan::math::internal::GradientCalc;
  namespace in = stan::math::internal;
  std::mt19937 rng(20261009);
  std::uniform_real_distribution<double> U(0.0, 1.0);
  const int n = 20000;
  double worst_v = 0, worst_g = 0;
  int nv6 = 0, nv3 = 0, ng6 = 0, ng3 = 0;
  for (int k = 0; k < n; ++k) {
    const double a = 0.5 + 2.5 * U(rng);
    const double v = -5.0 + 10.0 * U(rng);
    const double w = 0.2 + 0.6 * U(rng);
    const double sv = 1.5 * U(rng);
    const double t = 0.05 + 3.0 * U(rng);
    std::vector<var> x{t, a, 0.0, w, v, sv};
    var lp = stan::math::wiener_lpdf(x[0], x[1], x[2], x[3], x[4], x[5]);
    std::vector<double> g;
    lp.grad(x, g);
    const double let = in::wiener5_compute_log_error_term(t, a, v, w, sv);
    const double err = std::log(1e-15) + std::fmin(let, 0.0);
    const double r[4]
        = {in::wiener5_density<GradientCalc::OFF>(t, a, v, w, sv, err),
           in::wiener5_grad_t<GradientCalc::OFF>(t, a, v, w, sv, err),
           in::wiener5_grad_a<GradientCalc::OFF>(t, a, v, w, sv, err),
           in::wiener5_grad_w<GradientCalc::OFF>(t, a, v, w, sv, err)};
    const double dv = std::fabs(lp.val() - r[0]);
    double dg = 0;
    const int idx[3] = {0, 1, 3};
    for (int j = 0; j < 3; ++j) {
      dg = std::fmax(dg, std::fabs(g[idx[j]] - r[j + 1])
                             / std::fmax(1.0, std::fabs(r[j + 1])));
    }
    worst_v = std::fmax(worst_v, dv);
    worst_g = std::fmax(worst_g, dg);
    nv6 += dv > 1e-6;
    nv3 += dv > 1e-3;
    ng6 += dg > 1e-6;
    ng3 += dg > 1e-3;
    stan::math::recover_memory();
  }
  std::printf(
      "SCAN\t%d\tvalue worst %.3g >1e-6 %d >1e-3 %d\tgrad(y,a,w) worst %.3g "
      ">1e-6 %d >1e-3 %d\n",
      n, worst_v, nv6, nv3, worst_g, ng6, ng3);
  return 0;
}
