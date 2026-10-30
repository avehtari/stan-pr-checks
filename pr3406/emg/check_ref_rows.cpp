// Claims M1b and M4 of README_emg.md: the 13 reference rows and the 6
// deep-tail rows of the PR's exp_mod_normal_lcdf test, run on any tree with
// the comparison of the test: |expected - actual| <= 1e-12 |expected|.
//
// Output, one line per (table, function, row):
//   ROW <table> <fun> <n> z=<z> v=<v> <status of 6 checks> [details]
// The checks are: the double call, then the var value, dy, dmu, dsigma,
// dlambda.  Status "." = pass, "F" = finite but outside the tolerance,
// "N" = NaN or infinite where the reference is finite, "T" = throw.
// Then the vectorized-sum checks of the test (13-row table only) and a summary.
#include <stan/math/rev.hpp>
#include "../refs_emg/exp_mod_normal_lcdf_reference.hpp"
#include <array>
#include <cmath>
#include <cstdio>
#include <exception>
#include <string>
#include <vector>

using stan::math::var;
using row_t = std::array<double, 14>;

namespace {

const char* const CHECK_NAMES[6] = {"double", "value", "dy", "dmu", "dsigma",
                                    "dlambda"};
const char* const FUN_NAMES[3] = {"lcdf", "lccdf", "cdf"};

// The comparison of expect_rel in the test (EXPECT_NEAR semantics).
bool near(double expected, double actual, double tol) {
  return std::fabs(expected - actual) <= tol;
}

char status(double expected, double actual, double tol) {
  if (near(expected, actual, tol)) {
    return '.';
  }
  return std::isfinite(actual) || !std::isfinite(expected) ? 'F' : 'N';
}

template <typename Y, typename M, typename S, typename L>
auto call(int kind, const Y& y, const M& mu, const S& sigma, const L& lambda) {
  if (kind == 0) {
    return stan::math::exp_mod_normal_lcdf(y, mu, sigma, lambda);
  }
  if (kind == 1) {
    return stan::math::exp_mod_normal_lccdf(y, mu, sigma, lambda);
  }
  return stan::math::exp_mod_normal_cdf(y, mu, sigma, lambda);
}

// The reference value and the four partials of row r for kind (as the test).
std::array<double, 5> expected(int kind, const row_t& r) {
  std::array<double, 5> e;
  const std::size_t k = kind == 1 ? 9 : 4;
  for (std::size_t i = 0; i < 5; ++i) {
    e[i] = r[k + i];
  }
  if (kind == 2) {
    e[0] = std::exp(r[4]);
    for (std::size_t i = 1; i < 5; ++i) {
      e[i] = e[0] * r[4 + i];
    }
  }
  return e;
}

int n_fail_checks = 0;
int n_nonfinite = 0;

void check_rows(const char* table, const std::vector<row_t>& rows) {
  for (int kind = 0; kind < 3; ++kind) {
    for (std::size_t n = 0; n < rows.size(); ++n) {
      const row_t& r = rows[n];
      const auto e = expected(kind, r);
      double got[6];
      bool thrown[6] = {false, false, false, false, false, false};
      try {
        got[0] = call(kind, r[0], r[1], r[2], r[3]);
      } catch (const std::exception&) {
        got[0] = std::nan("");
        thrown[0] = true;
      }
      try {
        var y = r[0], mu = r[1], sigma = r[2], lambda = r[3];
        var lp = call(kind, y, mu, sigma, lambda);
        lp.grad();
        got[1] = lp.val();
        got[2] = y.adj();
        got[3] = mu.adj();
        got[4] = sigma.adj();
        got[5] = lambda.adj();
      } catch (const std::exception&) {
        for (int i = 1; i < 6; ++i) {
          got[i] = std::nan("");
          thrown[i] = true;
        }
      }
      stan::math::recover_memory();
      const double z = (r[0] - r[1]) / r[2];
      const double v = r[3] * r[2];
      std::string st, detail;
      for (int i = 0; i < 6; ++i) {
        const double ex = e[i == 0 ? 0 : i - 1];
        char c = thrown[i] ? 'T' : status(ex, got[i], 1e-12 * std::fabs(ex));
        st += c;
        if (!std::isfinite(got[i])) {
          ++n_nonfinite;
        }
        if (c != '.') {
          ++n_fail_checks;
          char buf[160];
          const double rel = std::fabs(got[i] - ex) / std::fabs(ex);
          std::snprintf(buf, sizeof(buf), " %s:%.17g(ref %.17g, rel %.1e)",
                        CHECK_NAMES[i], got[i], ex, rel);
          detail += buf;
        }
      }
      std::printf("ROW %s %s %zu z=%g v=%g %s%s\n", table, FUN_NAMES[kind], n,
                  z, v, st.c_str(), detail.c_str());
    }
  }
}

// The vectorized-sum checks of the test: y and lambda vary, mu and sigma
// are shared scalars.
void check_vectorized(const std::vector<row_t>& rows) {
  for (int kind = 0; kind < 2; ++kind) {
    std::vector<var> y, lambda;
    double value = 0, value_scale = 0, d_mu = 0, d_mu_scale = 0;
    double d_sigma = 0, d_sigma_scale = 0;
    for (const auto& r : rows) {
      const auto e = expected(kind, r);
      y.emplace_back(r[0]);
      lambda.emplace_back(r[3]);
      value += e[0];
      value_scale += std::fabs(e[0]);
      d_mu += e[2];
      d_mu_scale += std::fabs(e[2]);
      d_sigma += e[3];
      d_sigma_scale += std::fabs(e[3]);
    }
    var mu = rows[0][1];
    var sigma = rows[0][2];
    std::string st = "";
    try {
      var lp = call(kind, y, mu, sigma, lambda);
      lp.grad();
      st += near(value, lp.val(), 1e-12 * value_scale) ? '.' : 'F';
      st += near(d_mu, mu.adj(), 1e-12 * d_mu_scale) ? '.' : 'F';
      st += near(d_sigma, sigma.adj(), 1e-12 * d_sigma_scale) ? '.' : 'F';
      int bad_dy = 0, bad_dl = 0;
      for (std::size_t n = 0; n < rows.size(); ++n) {
        const auto e = expected(kind, rows[n]);
        bad_dy += !near(e[1], y[n].adj(), 1e-12 * std::fabs(e[1]));
        bad_dl += !near(e[4], lambda[n].adj(), 1e-12 * std::fabs(e[4]));
      }
      std::printf("VEC %s sum: value dmu dsigma %s; failing dy %d, dlambda %d"
                  " of %zu\n",
                  FUN_NAMES[kind], st.c_str(), bad_dy, bad_dl, rows.size());
    } catch (const std::exception& ex) {
      std::printf("VEC %s throws: %s\n", FUN_NAMES[kind], ex.what());
    }
    stan::math::recover_memory();
  }
}

}  // namespace

int main() {
  const std::vector<row_t> ref(stan::test::exp_mod_normal_lcdf_reference.begin(),
                               stan::test::exp_mod_normal_lcdf_reference.end());
  const std::vector<row_t> deep(
      stan::test::exp_mod_normal_lcdf_deep_tail_reference.begin(),
      stan::test::exp_mod_normal_lcdf_deep_tail_reference.end());
  std::printf("# checks per line: double value dy dmu dsigma dlambda\n");
  check_rows("ref", ref);
  check_rows("deep", deep);
  check_vectorized(ref);
  std::printf("SUMMARY failing row checks %d, non-finite results %d\n",
              n_fail_checks, n_nonfinite);
  return 0;
}
