// Bundle: value and var gradient of the 18 cdfs at the 324 rows of the
// #3408 reference test, on any tree.  Supports "199 of 324 non-finite",
// "1.6e+93" and "values change only where they were -inf or lost digits".
//
// The rows come from incomplete_gamma_cdfs_reference_324.hpp, the
// reference header of the PR head 22d781f290 (same directory).
// Output, one line per row:
//   <index> <function> <n_grad> <ref value> <ref grad x3> <value> <grad x3>
// All numbers in hex.  A row whose call throws prints THROW after n_grad.
#include <stan/math/rev.hpp>
#include "incomplete_gamma_cdfs_reference_324.hpp"
#include <cmath>
#include <cstdio>
#include <functional>
#include <map>
#include <string>
#include <vector>

namespace {

using stan::math::var;
namespace sm = stan::math;
using F3 = std::function<var(const var&, const var&, const var&)>;
using F2 = std::function<var(const var&, const var&)>;
using FP = std::function<var(int, const var&)>;

const std::map<std::string, F3> three_arg = {
    {"gamma_cdf", [](auto& y, auto& a, auto& b) { return sm::gamma_cdf(y, a, b); }},
    {"gamma_lcdf", [](auto& y, auto& a, auto& b) { return sm::gamma_lcdf(y, a, b); }},
    {"gamma_lccdf", [](auto& y, auto& a, auto& b) { return sm::gamma_lccdf(y, a, b); }},
    {"inv_gamma_cdf", [](auto& y, auto& a, auto& b) { return sm::inv_gamma_cdf(y, a, b); }},
    {"inv_gamma_lcdf", [](auto& y, auto& a, auto& b) { return sm::inv_gamma_lcdf(y, a, b); }},
    {"inv_gamma_lccdf",
     [](auto& y, auto& a, auto& b) { return sm::inv_gamma_lccdf(y, a, b); }},
    {"scaled_inv_chi_square_cdf",
     [](auto& y, auto& n, auto& s) { return sm::scaled_inv_chi_square_cdf(y, n, s); }},
    {"scaled_inv_chi_square_lcdf",
     [](auto& y, auto& n, auto& s) { return sm::scaled_inv_chi_square_lcdf(y, n, s); }},
    {"scaled_inv_chi_square_lccdf",
     [](auto& y, auto& n, auto& s) { return sm::scaled_inv_chi_square_lccdf(y, n, s); }},
};
const std::map<std::string, F2> two_arg = {
    {"chi_square_cdf", [](auto& y, auto& n) { return sm::chi_square_cdf(y, n); }},
    {"chi_square_lcdf", [](auto& y, auto& n) { return sm::chi_square_lcdf(y, n); }},
    {"chi_square_lccdf", [](auto& y, auto& n) { return sm::chi_square_lccdf(y, n); }},
    {"inv_chi_square_cdf", [](auto& y, auto& n) { return sm::inv_chi_square_cdf(y, n); }},
    {"inv_chi_square_lcdf", [](auto& y, auto& n) { return sm::inv_chi_square_lcdf(y, n); }},
    {"inv_chi_square_lccdf",
     [](auto& y, auto& n) { return sm::inv_chi_square_lccdf(y, n); }},
};
const std::map<std::string, FP> poisson = {
    {"poisson_cdf", [](int n, auto& l) { return sm::poisson_cdf(n, l); }},
    {"poisson_lcdf", [](int n, auto& l) { return sm::poisson_lcdf(n, l); }},
    {"poisson_lccdf", [](int n, auto& l) { return sm::poisson_lccdf(n, l); }},
};

}  // namespace

int main() {
  std::setvbuf(stdout, nullptr, _IOLBF, 0);
  const auto& rows = stan::test::incomplete_gamma_cdfs::rows();
  int i = 0;
  for (const auto& r : rows) {
    std::printf("%d %s %d %a %a %a %a ", i++, r.function.c_str(), r.n_grad,
                r.value, r.grad[0], r.grad[1], r.grad[2]);
    std::vector<var> args;
    double out[4] = {std::nan(""), 0.0, 0.0, 0.0};
    try {
      var lp;
      if (poisson.count(r.function)) {
        args = {var(r.p1)};
        lp = poisson.at(r.function)(static_cast<int>(r.y), args[0]);
      } else if (two_arg.count(r.function)) {
        args = {var(r.y), var(r.p1)};
        lp = two_arg.at(r.function)(args[0], args[1]);
      } else {
        args = {var(r.y), var(r.p1), var(r.p2)};
        lp = three_arg.at(r.function)(args[0], args[1], args[2]);
      }
      lp.grad();
      out[0] = lp.val();
      for (std::size_t k = 0; k < args.size(); ++k) {
        out[k + 1] = args[k].adj();
      }
    } catch (const std::exception& e) {
      std::printf("THROW\n");
      sm::recover_memory();
      continue;
    }
    sm::recover_memory();
    std::printf("%a %a %a %a\n", out[0], out[1], out[2], out[3]);
  }
  return 0;
}
