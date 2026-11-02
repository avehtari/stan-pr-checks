// f1-f4 against stan::math::normal_lpdf<true>: value, mu and sigma
// gradients (N = 1000), and the outcome (value or exception text) for
// invalid and edge inputs.
// Output: CHECK lines and CHECK_SUMMARY <n_fail>.

#include "normal_variants.hpp"
#include <cmath>
#include <cstdio>
#include <functional>
#include <limits>
#include <random>
#include <string>
#include <vector>

using fn_t = std::function<stan::math::var(const Eigen::VectorXd&,
                                           const nv::vm&,
                                           const stan::math::var&)>;
static int n_fail = 0;

struct out_t {
  std::string what;
  double val = 0, sig_adj = 0;
  Eigen::VectorXd mu_adj;
};

static out_t run(const fn_t& f, const Eigen::VectorXd& y,
                 const Eigen::VectorXd& mu_v, double s) {
  out_t o;
  try {
    stan::math::var sigma(s);
    nv::vm mu(mu_v);
    stan::math::var lp = f(y, mu, sigma);
    lp.grad();
    o.val = lp.val();
    o.sig_adj = sigma.adj();
    o.mu_adj = mu.adj();
    o.what = "value";
  } catch (const std::exception& e) {
    o.what = e.what();
  }
  stan::math::recover_memory();
  return o;
}

static bool same(double a, double b, double scale) {
  if (std::isnan(a) || std::isnan(b)) {
    return std::isnan(a) && std::isnan(b);
  }
  if (std::isinf(a) || std::isinf(b)) {
    return a == b;
  }
  return std::fabs(a - b) <= 1e-12 * scale;
}

int main() {
  const std::vector<std::pair<const char*, fn_t>> vars
      = {{"f1", nv::f1}, {"f2", nv::f2}, {"f3", nv::f3}, {"f4", nv::f4}};
  std::mt19937 rng(42);
  std::normal_distribution<double> z(0, 1);
  const int N = 1000;
  Eigen::VectorXd mu(N), y(N);
  for (int i = 0; i < N; ++i) {
    mu[i] = 0.3 + 2 * z(rng);
    y[i] = mu[i] + 0.7 * z(rng);
  }
  const double inf = std::numeric_limits<double>::infinity();
  const double nan = std::numeric_limits<double>::quiet_NaN();

  struct case_t {
    const char* name;
    Eigen::VectorXd y, mu;
    double sigma;
  };
  std::vector<case_t> cases;
  cases.push_back({"regular", y, mu, 0.7});
  cases.push_back({"small_sigma", y, mu, 1e-3});
  cases.push_back({"sigma_inf", y, mu, inf});
  {
    case_t c{"y_nan", y, mu, 0.7};
    c.y[3] = nan;
    cases.push_back(c);
  }
  {
    case_t c{"y_inf", y, mu, 0.7};
    c.y[3] = inf;
    cases.push_back(c);
  }
  {
    case_t c{"mu_inf", y, mu, 0.7};
    c.mu[5] = -inf;
    cases.push_back(c);
  }
  {
    case_t c{"mu_nan", y, mu, 0.7};
    c.mu[5] = nan;
    cases.push_back(c);
  }
  cases.push_back({"sigma_zero", y, mu, 0.0});
  cases.push_back({"sigma_neg", y, mu, -1.0});
  cases.push_back({"sigma_nan", y, mu, nan});
  {
    case_t c{"y_nan_sigma_neg", y, mu, -1.0};
    c.y[3] = nan;
    cases.push_back(c);
  }
  {
    case_t c{"y_inf_mu_inf", y, mu, 0.7};
    c.y[3] = inf;
    c.mu[3] = inf;
    cases.push_back(c);
  }
  cases.push_back(
      {"size_zero", Eigen::VectorXd(0), Eigen::VectorXd(0), 0.7});
  {
    case_t c{"size_mismatch", y, mu.head(N - 1), 0.7};
    cases.push_back(c);
  }

  for (const auto& c : cases) {
    const out_t ref = run(nv::dev, c.y, c.mu, c.sigma);
    for (const auto& v : vars) {
      const out_t o = run(v.second, c.y, c.mu, c.sigma);
      bool ok = o.what == ref.what;
      double m = 0;
      if (ok && ref.what == "value") {
        const double sc = std::max(1.0, std::fabs(ref.val));
        ok = same(o.val, ref.val, sc)
             && same(o.sig_adj, ref.sig_adj,
                     std::max(1.0, std::fabs(ref.sig_adj)) * 10);
        if (ok && ref.mu_adj.size() == o.mu_adj.size()) {
          for (Eigen::Index i = 0; i < ref.mu_adj.size(); ++i) {
            ok = ok
                 && same(o.mu_adj[i], ref.mu_adj[i],
                         std::max(1.0, std::fabs(ref.mu_adj[i])));
          }
        } else if (ref.mu_adj.size() != o.mu_adj.size()) {
          ok = false;
        }
        m = std::fabs(o.sig_adj - ref.sig_adj)
            / std::max(1.0, std::fabs(ref.sig_adj));
      }
      n_fail += !ok;
      std::printf("CHECK %s %s %s sig_rel %.2g | %s\n", c.name, v.first,
                  ok ? "OK" : "FAIL", m,
                  ref.what == "value" ? "value" : ref.what.substr(0, 70).c_str());
      if (!ok) {
        std::printf("  dev: %s\n  var: %s\n", ref.what.c_str(), o.what.c_str());
      }
    }
  }
  std::printf("CHECK_SUMMARY %d\n", n_fail);
  return n_fail == 0 ? 0 : 1;
}
